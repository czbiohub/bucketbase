import numpy as np
import pandas as pd
import sqlalchemy

def find_lowest_bin(database_address):
    '''
    if the bin number is not present, then we must find the lowest valid bin numebr by querying the db
    if the db just got created, infer -1 as the lowest bin number
    so that the numbering starts at 0
    '''
    engine=sqlalchemy.create_engine(f"sqlite:///{database_address}")
    connection=engine.connect()

    temp_cursor=connection.execute(
        '''
        select bin_id from bins
        order by bin_id desc
        limit 1
        '''
    )

    temp_result=temp_cursor.fetchall()

    connection.close()
    engine.dispose()

    if len(temp_result)==0:
        return -1
    elif len(temp_result)==1:
        return temp_result[0][0]



def create_bin_table_upload(alignment_panda,database_address):
    '''
    a bin panda is basically the output of pycutter 2 except that we skip
    the first several (currently 3) rows as well as skip the columns
    with feature magnitues.

    steps: 
    1) we pre-plan the set of bin_id that will be new
    2) we check conformity to standards? (to some extent?)
    3) we coerce the bin_panda into a panda for upload to db
    '''

    column_swap_dict={
        'Alignment ID':'alignment_id',
        'Metabolite name':'english_name',	
        'INCHIKEY':'inchikey',
        'Polarity/Filename':'polarity'
    }

    alignment_panda.rename(
        mapper=column_swap_dict,
        inplace=True,
        axis='columns'
    )
    
    
    #strip all whitespaces
    [alignment_panda[temp_col].apply(str.strip) for temp_col in alignment_panda.columns if alignment_panda[temp_col].dtype==str]

    bin_panda=alignment_panda.loc[
        alignment_panda.bin_id.isna()==True,
        [
            'bin_id','inchikey','adduct','english_name',
            'comment','polarity'
        ]
    ]

    #still need to make the is_known and is_istd columns
    #the bin_id column is there but is full of nan (yum!)
    bin_id_list=list(range(
        1+find_lowest_bin(database_address),len(bin_panda.index)
    ))
    is_istd_list=[
        1 if (x=='Internal Standard') else 0 for x in bin_panda.inchikey.tolist()
    ]
    is_known_list=[
        1 if (x=='Unknown') else 0 for x in bin_panda.english_name.tolist()
    ]
    group_id_list=[
        np.nan for x in range(len(bin_panda.index))
    ]
    valid_for_autocuration_list=[
        np.nan for x in range(len(bin_panda.index))
    ]
    consensus_rt=[
        np.nan for x in range(len(bin_panda.index))
    ]
    consensus_mz=[
        np.nan for x in range(len(bin_panda.index))
    ]
    consensus_spectrum=[
        np.nan for x in range(len(bin_panda.index))
    ]

    bin_panda['bin_id']=bin_id_list
    bin_panda['is_istd']=is_istd_list
    bin_panda['is_known']=is_known_list
    bin_panda['group_id']=group_id_list
    bin_panda['valid_for_autocuration']=valid_for_autocuration_list
    bin_panda['consensus_rt']=consensus_rt
    bin_panda['consensus_mz']=consensus_mz
    bin_panda['consensus_spectrum']=consensus_spectrum

    return bin_panda



if __name__=="__main__":
    final_alignment_address='../../../data/BRYU005_pipeline_test/step_2_final_alignment/BRYU005_CombineSubmit_June2022_pos.txt'
    database_address='../../../data/database/bucketbase.db'
    alignment_panda=pd.read_csv(final_alignment_address,sep='\t',skiprows=3)
    create_bin_table_upload(alignment_panda,database_address)