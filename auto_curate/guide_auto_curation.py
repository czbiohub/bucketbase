import pandas as pd
import sys
sys.path.insert(0, '../utils/')
from utils import execute_query_connection_established
import sqlalchemy

'''
the purpose of this is to provide some level of automatic curation to the post-pycutter-step-1 file
at this point,we have already taken the pycutter step 1 output and put it into a transient databse
we go through every "bin" in the transient database, and attempt curations by comparing with the main databsae
at current, we select the entire bins table from the transient databse to make this faster

i am deliberating on whether to select the entire main databse bin table. i will go against that strategy for now
but i expect that it would be significantly faster
i compromise by making a persistent connection (as i always should have) and then selecting a subset
'''


def select_all_bins(a_database_connection,ion_mode):
    '''
    get the entire bins table from the specified databse
    '''
    query=f'''PRAGMA table_info(bins)'''
    bins_column_reply=execute_query_connection_established(a_database_connection,query)
    bin_columns=[element[1] for element in bins_column_reply]
    
    query=f'''
    select * from bins
    where
    (polarity = "{ion_mode}") AND
    (valid_for_autocuration=1)'''
    bins_reply=execute_query_connection_established(a_database_connection,query)
    bins_df=pd.DataFrame.from_records(
        data=bins_reply,
        columns=bin_columns
    )

    #print(bins_df)
    return bins_df
    #return [element[0] for element in execute_query(database_address,query)]

def select_candidate_bins(
    a_database_connection,
    #i guess we really dont need ion mode if we specify adduct
    ion_mode,
    transient_adduct,
    transient_mz,
    transient_rt,
    mz_tolerance,
    rt_tolerance
):

    query=f'''PRAGMA table_info(bins)'''
    bins_column_reply=execute_query_connection_established(a_database_connection,query)
    bin_columns=[element[1] for element in bins_column_reply]

    mz_min=transient_mz-mz_tolerance
    mz_max=transient_mz+mz_tolerance
    rt_min=transient_rt-rt_tolerance
    rt_max=transient_rt+rt_tolerance   

    query=f'''
    select * from bins
    where
    (valid_for_autocuration=1) AND
    (polarity="{ion_mode}") AND
    (adduct="{transient_adduct}") AND
    ((consensus_mz>{mz_min}) AND (consensus_mz<{mz_max})) AND
    ((consensus_rt>{rt_min}) AND (consensus_rt<{rt_max}))
    '''
    bins_reply=execute_query_connection_established(a_database_connection,query)
    if bins_reply=='no results found':
       return
    bins_df=pd.DataFrame.from_records(
        data=bins_reply,
        columns=bin_columns
    )

    return bins_df

#def select_candidate_bins(database_address,)


#can refactor as an apply if desired once we see what columns we weant
#expect a 4x speedup compared to iterrows
def autocurate_transient_bins(
    transient_bins_panda,
    database_connection,
    ion_mode,
    weight_of_dot_product,
    weight_of_reverse_dot_product,
    mz_tolerance,
    rt_tolerance
):
    '''
    '''

    #select_database_bins=
    for row,series in transient_bins_panda.iterrows():
        candidate_bins_df=select_candidate_bins(
            database_connection,
            ion_mode,
            series['adduct'],
            series['consensus_mz'],
            series['consensus_rt'],
            mz_tolerance,
            rt_tolerance
        )
        print(candidate_bins_df)
        hold=input('above is candidates df')

def guide_auto_curation_wrapper(
        database_connection,
        transient_database_connection,
        ion_mode,
        weight_of_dot_product,
        weight_of_reverse_dot_product,
        mz_tolerance,
        rt_tolerance
    ):
    transient_bins_panda=select_all_bins(transient_connection,ion_mode)
    #database_bins_panda=select_all_bins(database_connection,ion_mode)
    print(transient_bins_panda)
    print(transient_bins_panda.columns)
    autocurate_transient_bins(
        transient_bins_panda,
        database_connection,
        ion_mode,
        weight_of_dot_product,
        weight_of_reverse_dot_product,
        mz_tolerance,
        rt_tolerance
    )




if __name__=="__main__":

    ion_mode='pos'
    weight_of_dot_product=0.5
    weight_of_reverse_dot_product=0.5
    mz_tolerance=0.01
    rt_tolerance=1


    transient_database_address='../../data/database/transient_bucketbase.db'
    database_address='../../data/database/bucketbase.db'

    transient_engine=sqlalchemy.create_engine(f"sqlite:///{transient_database_address}")
    transient_connection=transient_engine.connect()

    database_engine=sqlalchemy.create_engine(f"sqlite:///{database_address}")
    database_connection=database_engine.connect()







    guide_auto_curation_wrapper(
        database_connection,
        transient_connection,
        ion_mode,
        weight_of_dot_product,
        weight_of_reverse_dot_product,
        mz_tolerance,
        rt_tolerance
        )

    
    transient_connection.close()
    transient_engine.dispose()

    database_connection.close()
    database_engine.dispose()