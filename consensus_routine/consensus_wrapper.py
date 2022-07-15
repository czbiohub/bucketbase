import sqlalchemy
import sys
sys.path.insert(0, '../utils/')
from utils import execute_query

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
   


def guide_consensus_routine(database_address,spectrum_cutoff):
    bins_without_autocuration_validity=aquire_bin_ids_without_autocuration_status(database_address)
    bins_without_autocuration_non_zero_spectrum_count=aquire_bin_ids_and_spectrum_count_without_autocuration(database_address)
    bins_enough_spectra_for_auto_curate=[
        element[0] for element in bins_without_autocuration_non_zero_spectrum_count if (element[1]>spectrum_cutoff)
    ]
    bins_that_wont_be_used_for_autocuration=set(bins_without_autocuration_validity).difference(
        set(bins_enough_spectra_for_auto_curate)
    )
    bins_mzrt_only=set(bins_without_autocuration_validity).difference(
         #do the comprehension to extract the bin ids on the fly
         {element[0] for element in bins_without_autocuration_non_zero_spectrum_count}
    )
    set_autocurate_valid(bins_that_wont_be_used_for_autocuration,0)
    #we do not yet set the remaining bins to autocurate valid=1, because they need to undergo additional tests
    #like entropic check
    set_mzrt_only(bins_without_autocuration_non_zero_spectrum_count,0)
    set_mzrt_only(bins_mzrt_only,1)

    print(len(bins_without_autocuration_validity))
    print(len(bins_without_autocuration_non_zero_spectrum_count))
    print(len(bins_mzrt_only))
    print(len(bins_enough_spectra_for_auto_curate))

    set_autocurate_valid(bins_enough_spectra_for_autocurate,1)
    

if __name__=="__main__":
    minimum_count_for_auto_curation_possible=20
    
    database_address='../../data/database/bucketbase.db'
    guide_consensus_routine(database_address,minimum_count_for_auto_curation_possible)