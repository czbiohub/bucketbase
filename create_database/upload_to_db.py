from prepare_bin_table_upload import *
from prepare_run_table_upload import *
from prepare_annotation_table_upload import *
import sys

def upload_table_to_db(temp_panda,table_name):

    engine=sqlalchemy.create_engine(f"sqlite:///{database_address}")
    #connection=engine.connect()

    temp_panda.to_sql(
        table_name,
        con=engine,
        if_exists='append',
        index=False,
        #chunksize=5000,
        #method='multi'
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

    study_name=sys.argv[1]
    to_transient_for_pycutter_pipeline=sys.argv[2]
    ion_mode=sys.argv[3]


    
    mapping_file_address=f'../data/{study_name}/input/{ion_mode}/{study_name}_{ion_mode}_mapping.txt'
    individual_files_directory=f'../data/{study_name}/input/{ion_mode}/individual_files/'

    if to_transient_for_pycutter_pipeline=='transient':
        final_alignment_address=f'../data/{study_name}/intermediates/step_1_post_pycutter/pycutter_step_1_output_{ion_mode}.tsv'
        database_address=f"../data/{study_name}/database/{to_transient_for_pycutter_pipeline}_bucketbase.db"
    elif to_transient_for_pycutter_pipeline=='main':
        final_alignment_address=''

    #run_panda
    if to_transient_for_pycutter_pipeline=='transient':
        alignment_panda=pd.read_csv(
            final_alignment_address,sep='\t',nrows=5,header=None
        )
    elif to_transient_for_pycutter_pipeline=='main':
        alignment_panda=pd.read_csv(
            final_alignment_address,sep='\t',nrows=4,header=None
        )    
    run_panda_for_upload=create_run_table_upload(
        alignment_panda,
        to_transient_for_pycutter_pipeline
    )

    #bin_panda
    if to_transient_for_pycutter_pipeline=='transient':
        alignment_panda=pd.read_csv(final_alignment_address,sep='\t',skiprows=4)
        bin_panda_for_upload=create_bin_table_upload(alignment_panda,database_address,to_transient_for_pycutter_pipeline,ion_mode)
    elif to_transient_for_pycutter_pipeline=='main':
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