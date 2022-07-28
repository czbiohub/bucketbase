import numpy as np
import pandas as pd
import sqlalchemy


def create_run_table_upload(alignment_panda,to_transient_for_pycutter_pipeline):
    '''
    steps: 
    1) we pre-plan the set of bin_id that will be new
    2) we check conformity to standards? (to some extent?)
    3) we coerce the bin_panda into a panda for upload to db
    '''

    #drop the columns where the first three rows are na
    alignment_panda.dropna(
        axis='columns',
        how='all',
        subset=[0,1,2],
        inplace=True
    )

    alignment_panda.set_index(
        alignment_panda.columns[0],
        drop=True,
        append=False,
        inplace=True
    )

    #an artifact of this operation seems to be that the column index is named
    #which is weird (in this case, 9, which was a column name previous)
    #doesnt seem to affect things, so will ignore
    alignment_panda=alignment_panda.transpose()

    if to_transient_for_pycutter_pipeline=='transient':
        column_swap_dict={
            'Class':'irrelevant_1',
            'File type':'run_type', 
            'Injection order':'irrelevant_2', 
            'Batch ID':'irrelevant_3',
            'MS/MS spectrum':'run_id'
        }       

    elif to_transient_for_pycutter_pipeline=='main':
        column_swap_dict={
            'Class':'class',
            'Sample Type':'run_type', 
            'Name from Collaborator':'name_from_collaborator', 
            'Polarity/Filename':'run_id'
        }

    alignment_panda.rename(
        mapper=column_swap_dict,
        inplace=True,
        axis='columns'
    )

    #strip all whitespaces
    [alignment_panda[temp_col].apply(str.strip) for temp_col in alignment_panda.columns if alignment_panda[temp_col].dtype==str]

    run_panda=alignment_panda.loc[
        :,
        [
            'run_id','run_type'
        ]
    ]

    run_panda['method_id']='my_dummy_method_id'
    run_panda['sample_des_id']='my_dummy_sample_des_id'

    return run_panda

if __name__=="__main__":
    
    #out of date
    final_alignment_address='../../../data/BRYU005_pipeline_test/step_2_final_alignment/BRYU005_CombineSubmit_June2022_pos.txt'
    database_address='../../../data/database/bucketbase.db'
    alignment_panda=pd.read_csv(
        final_alignment_address,sep='\t',nrows=4,header=None
    )
    create_run_table_upload(alignment_panda)