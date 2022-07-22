from prepare_bin_table_upload import *
from prepare_run_table_upload import *
from prepare_annotation_table_upload import *

def upload_table_to_db(temp_panda,table_name):

    engine=sqlalchemy.create_engine(f"sqlite:///{database_address}")
    #connection=engine.connect()

    temp_panda.to_sql(
        table_name,
        con=engine,
        if_exists='append',
        index=False,
        chunksize=5000,
        method='multi'
    )

    #connection.close()
    engine.dispose()


if __name__=="__main__":
    '''
    order of events:
    1) run_panda
    2) bin_panda
    3) annotations_panda (which depends on bin_panda)
    unfortunately, this breaks if the bin_panda is created first. unless you want to recreate the bin panda after 
    making the run panda. but, ultimately, somehow the alignment panda is modified when making the bin panda
    and that is taken advantage of when making the annotation panda
    '''

    to_transient_for_pycutter_pipeline=True
    ion_mode='pos'
    if to_transient_for_pycutter_pipeline==True:
        database_address="../../../data/database/transient_bucketbase.db"
        final_alignment_address='../../../data/BRYU005_pipeline_test/step_1_post_pycutter/py_cutter_step_1_output.tsv'
        mapping_file_address='../../../data/three_studies/alignment_individual_mappings/BRYU005_pos_mapping.txt'
        individual_files_directory='../../../data/three_studies/individual_sample_data_subset/unzipped/BRYU005_Bacterial_Supernatant/pos/'
    elif to_transient_for_pycutter_pipeline==False:
        database_address="../../../data/database/bucketbase.db"
        final_alignment_address='../../../data/BRYU005_pipeline_test/step_2_final_alignment/BRYU005_CombineSubmit_June2022_pos.txt'
        database_address='../../../data/database/bucketbase.db'
        mapping_file_address='../../../data/three_studies/alignment_individual_mappings/BRYU005_pos_mapping.txt'
        individual_files_directory='../../../data/three_studies/individual_sample_data_subset/unzipped/BRYU005_Bacterial_Supernatant/pos/'    


    #run_panda
    if to_transient_for_pycutter_pipeline==True:
        alignment_panda=pd.read_csv(
            final_alignment_address,sep='\t',nrows=5,header=None
        )
    elif to_transient_for_pycutter_pipeline==False:
        alignment_panda=pd.read_csv(
            final_alignment_address,sep='\t',nrows=4,header=None
        )    
    run_panda_for_upload=create_run_table_upload(
        alignment_panda,
        to_transient_for_pycutter_pipeline
    )

    #bin_panda
    if to_transient_for_pycutter_pipeline==True:
        alignment_panda=pd.read_csv(final_alignment_address,sep='\t',skiprows=4)
        bin_panda_for_upload=create_bin_table_upload(alignment_panda,database_address,to_transient_for_pycutter_pipeline,ion_mode)
    elif to_transient_for_pycutter_pipeline==False:
        alignment_panda=pd.read_csv(final_alignment_address,sep='\t',skiprows=3)    
        bin_panda_for_upload=create_bin_table_upload(alignment_panda,database_address,to_transient_for_pycutter_pipeline)
    

    #annotation_panda (depends on bin_panda)
    alignment_id_bin_id_panda=get_alignment_id_bin_id_map(alignment_panda,bin_panda_for_upload)
    mapping_panda=pd.read_csv(mapping_file_address,sep='\t',skiprows=4)
    mapping_panda=clean_mapping_panda(mapping_panda,alignment_id_bin_id_panda)
    annotation_panda_for_upload=create_annotation_table_wrapper(individual_files_directory,mapping_panda,database_address)

    upload_table_to_db(run_panda_for_upload,'runs')
    upload_table_to_db(bin_panda_for_upload,'bins')
    upload_table_to_db(annotation_panda_for_upload,'annotations')