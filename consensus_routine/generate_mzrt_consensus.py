import pandas as pd
from utils import execute_query_connection_established

def select_rts_for_bin(database_connection,bin_id):
    query=f'''
    select retention_time 
    from annotations
    inner join
    runs
    on annotations.run_id=runs.run_id
    where (annotations.bin_id={bin_id}) and (annotations.retention_time is not null)
    '''
    return [element[0] for element in execute_query_connection_established(database_connection,query)]

def select_mzs_for_bin(database_connection,bin_id):
    query=f'''
    select precursor_mz 
    from annotations
    inner join
    runs
    on annotations.run_id=runs.run_id
    where (annotations.bin_id={bin_id}) and (annotations.precursor_mz is not null)
    '''
    return [element[0] for element in execute_query_connection_established(database_connection,query)]


def generate_mzrt_wrapper(
        database_connection,
        bin_list,

    ):
    '''
    for each cluster, make a consensus spectra in text format,
    then concatentate them and return that
    '''
    result_dict={
        'bin_id':[],
        'consensus_mz':[],
        'consensus_rt':[]
    }


    for z,temp_bin in enumerate(bin_list):
        print(f'bin number {temp_bin} iteration number {z}')
        
        temp_mz_list=select_mzs_for_bin(database_connection,temp_bin)
        temp_rt_list=select_rts_for_bin(database_connection,temp_bin)

        result_dict['bin_id'].append(temp_bin)

        result_dict['consensus_mz'].append(sum(temp_mz_list)/len(temp_mz_list))
        result_dict['consensus_rt'].append(sum(temp_rt_list)/len(temp_rt_list))
    return pd.DataFrame.from_dict(result_dict)

if __name__=="__main__":
    temp_rts=select_rts_for_bin(database_connection,temp_bin)
    print(len(temp_rts))
    temp_mzs=select_mzs_for_bin(database_connection,temp_bin)
    print(len(temp_mzs))