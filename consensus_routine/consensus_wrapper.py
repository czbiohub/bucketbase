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

def guide_consensus_routine(database_address):
    pass
    bins_without_autocuration_validity=aquire_bin_ids_without_autocuration_status(database_address)
    print(bins_without_autocuration_validity)
    #assign_autocuration_validity()

if __name__=="__main__":
    database_address='../../data/database/bucketbase.db'
    guide_consensus_routine(database_address)