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



def create_bin_table_upload(alignment_panda,database_address,to_transient_for_pycutter_pipeline,ion_mode=None):
    '''
    a bin panda is basically the output of pycutter 2 except that we skip
    the first several (currently 3) rows as well as skip the columns
    with feature magnitues.

    steps: 
    1) we pre-plan the set of bin_id that will be new
    2) we check conformity to standards? (to some extent?)
    3) we coerce the bin_panda into a panda for upload to db

    In retrospect, some of the things here are little clunky. 
    '''


    if to_transient_for_pycutter_pipeline==True:
        alignment_panda['bin_id']=np.nan
        alignment_panda['comment']='junk comment'
        alignment_panda['polarity']=ion_mode
        column_swap_dict={
            'INCHIKEY':'inchikey',
            'Adduct type':'adduct',
            'Metabolite name':'english_name',
            'Alignment ID':'alignment_id'
        }
    elif to_transient_for_pycutter_pipeline==False:
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
    
    print(alignment_panda)
    hold=input('alignment panda')

    #strip all whitespaces
    #print(alignment_panda[alignment_panda.columns[1]])
    [alignment_panda[temp_col].apply(str.strip) for temp_col in alignment_panda.columns if alignment_panda[temp_col].dtype==str]
    #[print(alignment_panda[temp_col]) for temp_col in alignment_panda.columns]
    
    


    bin_panda=alignment_panda.loc[
        alignment_panda.bin_id.astype(str).str.isdigit()==False,
        [
            'bin_id','inchikey','adduct','english_name',
            'comment','polarity'
        ]
    ]
    bin_panda.reset_index(inplace=True,drop=True)
    print(bin_panda)
    hold=input('see above bin panda')

    print(find_lowest_bin(database_address))

    #still need to make the is_known and is_istd columns
    #the bin_id column is there but is full of nan (yum!)
    bin_id_list=list(range(
        1+find_lowest_bin(database_address),(1+find_lowest_bin(database_address)+len(bin_panda.index))
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
    mzrt_only=[
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
    bin_panda['mzrt_only']=mzrt_only


    print(bin_panda)
    hold=input('verify bin panda')



    return bin_panda



if __name__=="__main__":
    
    #out of date
    final_alignment_address='../../../data/BRYU005_pipeline_test/step_2_final_alignment/BRYU005_CombineSubmit_June2022_pos.txt'
    database_address='../../../data/database/bucketbase.db'
    alignment_panda=pd.read_csv(final_alignment_address,sep='\t',skiprows=3)
    create_bin_table_upload(alignment_panda,database_address)