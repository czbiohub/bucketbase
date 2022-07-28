import sqlalchemy
import pandas as pd
import sys
sys.path.insert(0, '../utils/')
from utils import execute_query_connection_established
from valid_for_autocuration_test import *
from generate_mzrt_consensus import *
import random

def aquire_bin_ids_without_autocuration_status(database_connection):
    query='''
    select bin_id 
    from bins
    where valid_for_autocuration is null
    '''
    return [element[0] for element in execute_query_connection_established(database_connection,query)]

def aquire_bin_ids_and_spectrum_count_without_autocuration(database_connection):
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
    #print(execute_query_connection_established(database_connection,query))
    #print(len(execute_query_connection_established(database_connection,query)))
    return [element for element in execute_query_connection_established(database_connection,query)]
   
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

    execute_query_connection_established(database_connection,query,returns_rows=False)


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

    execute_query_connection_established(database_connection,query,returns_rows=False)

def update_bins_with_spectra(
    dataframe
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
    where bin_id=temp_table.column1'''

    execute_query_connection_established(database_connection,query_string,returns_rows=False)

    # engine=sqlalchemy.create_engine(f"sqlite:///{database_connection}")
    # connection=engine.connect()

    # temp_cursor=connection.execute(query_string)

    # connection.close()
    # engine.dispose()


def update_bins_with_mzrt(
    dataframe
):
    values_text_list='('+dataframe.bin_id.astype(str)+', '+dataframe.consensus_mz.astype(str)+', '+dataframe.consensus_rt.astype(str)+')'
    values_text=',\n'.join(values_text_list.values)
    print(values_text)
    # print('-----------------------------------------')
    # print(database_connection)
    

    #for the life of me, i tried to rename the columns in temp table. could not get the syntax to work
    query_string='''
    update bins set
        consensus_mz = temp_table.column2,
        consensus_rt = temp_table.column3
    from (values
    '''+values_text+'''
    )as temp_table
    where bin_id=temp_table.column1'''


    execute_query_connection_established(database_connection,query_string,returns_rows=False)
    #hold=input('inspect db connect')
    # engine=sqlalchemy.create_engine(f"sqlite:///{database_connection}")
    # connection=engine.connect()

    # temp_cursor=connection.execute(query_string)

    # connection.close()
    # engine.dispose()

def guide_consensus_routine_generate(database_connection,spectrum_cutoff):
    
    ########get bin numbers that are new (need to be curated/checked for curation)#########
    bins_without_autocuration_status=aquire_bin_ids_without_autocuration_status(database_connection)
    #bins that have at least one spectrum (and therefore need the valid-for-autocuration test)
    #we do all the different stuff like cluster on all of them for code-base simplicity
    bins_without_autocuration_non_zero_spectrum_count=aquire_bin_ids_and_spectrum_count_without_autocuration(database_connection)

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
        database_connection,
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
        bins_panda_spectra
    )

    bins_panda_mzrt=generate_mzrt_wrapper(
        database_connection,
        bins_without_autocuration_status
    )

    update_bins_with_mzrt(
        bins_panda_mzrt
    )

def guide_consensus_routine_update(database_connection,final_alignment_address,max_consensus_contributers,
noise_level,ms2_tolerance,similarity_metric,mutual_distance_for_cluster,minimum_percent_present,bin_space_tolerance):
    '''
    '''
    alignment_panda=pd.read_csv(final_alignment_address,sep='\t',skiprows=3)
    bins_for_updating=alignment_panda['bin_id'].loc[
        #alignment_panda['bin_id'].notna()
            alignment_panda.bin_id.astype(str).str.isdigit()==True
        ].astype(int).to_list()
    print(bins_for_updating)

    result_dict={
        'bin_id':[],
        'valid_for_autocuration':[],
        'consensus_spectrum':[]
    }

    for temp_bin in bins_for_updating:
        #a lot of this code is borrowed form valid_for_autocurations_test
        #the difference in this method is taht we dont check whether the autocuration is valid
        #the reason for this is that we made the automatic annotation based on whether the main DB
        #bin is like the bin that we see in this study
        #it would be unfortunate if our sampling, by chance, crossed the "fail" threshold and this bin was
        #frozen in time
        #so we proceed directly to clustering and consensussing


        annotations_and_spectra=select_spectra_for_bin(database_connection,temp_bin)

        # query='select * from bins limit 10'
        # junk_test=execute_query_connection_established(database_connection,query)
        # print(junk_test)
        # hold=input('junk test')

        temp_annotation_ids=[element[0] for element in annotations_and_spectra]
        temp_spectra_text=[element[1] for element in annotations_and_spectra]

        print(len(temp_annotation_ids))
        print('temp annotations ids length')

        update_member_of_consensus(database_connection,temp_annotation_ids,0)

        
        #choose up to max_consensus_contributers randomly
        #we do this because getting the clusters from a large number of spectra (could be come 10k+)
        #is slow. and the consensus spectra probably asymptote
        if len(annotations_and_spectra)>max_consensus_contributers:
            random.shuffle(annotations_and_spectra)
            annotations_and_spectra=annotations_and_spectra[0:max_consensus_contributers]


        print(f'{len(temp_spectra_text)} spectra') #len(temp_spectra_text))
        temp_spectra_paired=parse_text_spectra_return_pairs(temp_spectra_text)
        temp_spectra_paired_cleaned=get_cleaned_spectra(temp_spectra_paired,noise_level,ms2_tolerance)


        update_member_of_consensus(database_connection,temp_annotation_ids,1)


        cluster_assignments=perform_hierarchical_clustering_routine(temp_spectra_paired_cleaned,similarity_metric,ms2_tolerance,mutual_distance_for_cluster)
        #count membership and get cluster percent
        cluster_assignments_sorted_by_membership,biggest_cluster_percent=get_cluster_membership_ordering(cluster_assignments)   

        consensus_spectra_text=generate_consensus_spectra_text_wrapper(
            temp_spectra_paired_cleaned,
            cluster_assignments,
            cluster_assignments_sorted_by_membership,
            ms2_tolerance,
            minimum_percent_present,
            bin_space_tolerance
        )
        result_dict['bin_id'].append(temp_bin)
        #its already valid for autocuration, this changes nothing
        result_dict['valid_for_autocuration'].append(1)
        result_dict['consensus_spectrum'].append(consensus_spectra_text)

    bins_panda_spectra=pd.DataFrame.from_dict(result_dict)

    update_bins_with_spectra(
        bins_panda_spectra
    )

    bins_panda_mzrt=generate_mzrt_wrapper(
        database_connection,
        bins_for_updating
    )

    update_bins_with_mzrt(
        bins_panda_mzrt
    )



if __name__=="__main__":
    #generate or update
    #i would expect that update would be the choice if and only if
    #the database style the "main" database
    
    consensus_style='update'
    #consensus_style='generate'
    max_consensus_contributers=500

    noise_level=0.03
    ms2_tolerance=0.015
    similarity_metric='dot_product'
    mutual_distance_for_cluster=0.2
    minimum_percent_present=0.3
    bin_space_tolerance=3


    minimum_count_for_auto_curation_possible=20
    to_transient_for_pycutter_pipeline=False
    if to_transient_for_pycutter_pipeline==True:
        database_address="../../data/database/transient_bucketbase.db"
    elif to_transient_for_pycutter_pipeline==False:
        database_address="../../data/database/bucketbase.db"


    database_engine=sqlalchemy.create_engine(f"sqlite:///{database_address}")
    database_connection=database_engine.connect()

    

    if consensus_style=='generate':
        guide_consensus_routine_generate(database_connection,minimum_count_for_auto_curation_possible)
    elif consensus_style=='update':
        final_alignment_address='../../data/BRYU005_pipeline_test/step_2_final_alignment/py_cutter_step_2_output_auto_curated.tsv'
        guide_consensus_routine_update(database_connection,final_alignment_address,max_consensus_contributers,
        noise_level,ms2_tolerance,similarity_metric,mutual_distance_for_cluster,minimum_percent_present,bin_space_tolerance)