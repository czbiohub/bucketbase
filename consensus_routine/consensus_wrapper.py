import sqlalchemy
import pandas as pd
import sys
sys.path.insert(0, '../utils/')
from utils import execute_query
from valid_for_autocuration_test import *
from generate_mzrt_consensus import *

def aquire_bin_ids_without_autocuration_status(database_address):
    query='''
    select bin_id 
    from bins
    where valid_for_autocuration is null
    '''
    return [element[0] for element in execute_query(database_address,query)]

def aquire_bin_ids_and_spectrum_count_without_autocuration(database_address):
    #this query gives you the number of spectra that are not null
    #corresponding to each bin id that does not have a valid_for_autocurate assigned
    #and where spectra come from samples
    #there must be at least 1 spectrum for a bin_id to show up (group by bin id)
    query='''
    select b.bin_id, count(*)
    from bins b
    inner join (
        select a.bin_id,a.spectrum from
        annotations a
        inner join
        runs r 
        on a.run_id=r.run_id 
        where r.run_type='Sample'
        ) as iq1
    on b.bin_id=iq1.bin_id
    where (b.valid_for_autocuration is null) and (iq1.spectrum is not null)
    group by b.bin_id
    order by b.bin_id
    '''
    #print(execute_query(database_address,query))
    #print(len(execute_query(database_address,query)))
    return [element for element in execute_query(database_address,query)]
   
def set_mzrt_only(bins,value):
    '''
    receives a set or list of bins and a 0 or 1.
    '''
    #actually not terible approach
    #https://stackoverflow.com/questions/38199080/efficient-sqlite-query-based-on-list-of-primary-keys
    #maybe a faster way would be making a table on the fly?
    #https://stackoverflow.com/questions/51760204/unnest-or-similar-in-sqlite
    bins_as_strings=[str(element) for element in bins]
    total_string='\', \''.join(bins_as_strings)
    total_string='(\''+total_string+'\')'

    query=f'''UPDATE bins
    set mzrt_only = {value}
    where bin_id IN '''+total_string

    execute_query(database_address,query,returns_rows=False)


def set_autocurate_valid(bins,value):
    '''
    receives a set or list of bins and a 0 or 1.
    '''
    #actually not terible approach
    #https://stackoverflow.com/questions/38199080/efficient-sqlite-query-based-on-list-of-primary-keys
    #maybe a faster way would be making a table on the fly?
    #https://stackoverflow.com/questions/51760204/unnest-or-similar-in-sqlite
    bins_as_strings=[str(element) for element in bins]
    total_string='\', \''.join(bins_as_strings)
    total_string='(\''+total_string+'\')'

    query=f'''UPDATE bins
    set valid_for_autocuration = {value}
    where bin_id IN '''+total_string

    execute_query(database_address,query,returns_rows=False)

def update_bins_with_spectra(
    dataframe,
    database_address
):
    values_text_list='('+dataframe.bin_id.astype(str)+', '+dataframe.valid_for_autocuration.astype(str)+', \''+dataframe.consensus_spectrum+'\')'
    values_text=',\n'.join(values_text_list.values)
    print(values_text)

    #for the life of me, i tried to rename the columns in temp table. could not get the syntax to work
    query_string='''
    update bins set
        consensus_spectrum = temp_table.column3,
        valid_for_autocuration = temp_table.column2
    from (values
    '''+values_text+'''
    )as temp_table
    where bins.bin_id=temp_table.column1'''

    engine=sqlalchemy.create_engine(f"sqlite:///{database_address}")
    connection=engine.connect()

    temp_cursor=connection.execute(query_string)

    connection.close()
    engine.dispose()


def update_bins_with_mzrt(
    dataframe,
    database_address
):
    values_text_list='('+dataframe.bin_id.astype(str)+', '+dataframe.consensus_mz.astype(str)+', '+dataframe.consensus_rt.astype(str)+')'
    values_text=',\n'.join(values_text_list.values)
    print(values_text)

    #for the life of me, i tried to rename the columns in temp table. could not get the syntax to work
    query_string='''
    update bins set
        consensus_mz = temp_table.column2,
        consensus_rt = temp_table.column3
    from (values
    '''+values_text+'''
    )as temp_table
    where bins.bin_id=temp_table.column1'''

    engine=sqlalchemy.create_engine(f"sqlite:///{database_address}")
    connection=engine.connect()

    temp_cursor=connection.execute(query_string)

    connection.close()
    engine.dispose()

def guide_consensus_routine(database_address,spectrum_cutoff):
    
    ########get bin numbers that are new (need to be curated/checked for curation)#########
    bins_without_autocuration_status=aquire_bin_ids_without_autocuration_status(database_address)
    #bins that have at least one spectrum (and therefore need the valid-for-autocuration test)
    #we do all the different stuff like cluster on all of them for code-base simplicity
    bins_without_autocuration_non_zero_spectrum_count=aquire_bin_ids_and_spectrum_count_without_autocuration(database_address)

    bins_mzrt_only=set(bins_without_autocuration_status).difference(
         #do the comprehension to extract the bin ids on the fly
         {element[0] for element in bins_without_autocuration_non_zero_spectrum_count}
    )
    ######################################################################################

    ############set mzrt only property###################################################
    set_mzrt_only(
        [element[0] for element in bins_without_autocuration_non_zero_spectrum_count],
        0
    )
    set_mzrt_only(bins_mzrt_only,1)
    set_autocurate_valid(bins_mzrt_only,0)
    ######################################################################################

    
    ####################do auto-curation################################################
    #we receive a dataframe in response.
    #we choose a dataframe, rather than doing any updating in this method, so that we can
    #parallelize it as desired
    bins_panda_spectra=valid_for_autocuration_test_wrapper(
        database_address,
        [element[0] for element in bins_without_autocuration_non_zero_spectrum_count],
        'dot_product',
        0.015,
        0.03,
        0.2,
        True,
        60,
        0.5,
        0.9,
        20,
        0.3,
        3
    )
    ######################################################################################
    
    update_bins_with_spectra(
        bins_panda_spectra,
        database_address
    )

    bins_panda_mzrt=generate_mzrt_wrapper(
        database_address,
        bins_without_autocuration_status
    )

    update_bins_with_mzrt(
        bins_panda_mzrt,
        database_address
    )


if __name__=="__main__":
    
    minimum_count_for_auto_curation_possible=20
    to_transient_for_pycutter_pipeline=False
    if to_transient_for_pycutter_pipeline==True:
        database_address="../../data/database/transient_bucketbase.db"
    elif to_transient_for_pycutter_pipeline==False:
        database_address="../../data/database/bucketbase.db"


    guide_consensus_routine(database_address,minimum_count_for_auto_curation_possible)