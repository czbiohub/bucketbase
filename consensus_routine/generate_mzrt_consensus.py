def select_rts_for_bin(database_address,bin_id):
    query=f'''
    select retention_time 
    from annotations
    inner join
    runs
    on annotations.run_id=runs.run_id
    where (annotations.bin_id={bin_id}) and (annotations.retention_time is not null) and (runs.run_type='Sample')
    '''
    return [element[0] for element in execute_query(database_address,query)]

def select_mzs_for_bin(database_address,bin_id):
    query=f'''
    select precursor_mz 
    from annotations
    inner join
    runs
    on annotations.run_id=runs.run_id
    where (annotations.bin_id={bin_id}) and (annotations.precursor_mz is not null) and (runs.run_type='Sample')
    '''
    return [element[0] for element in execute_query(database_address,query)]


if __name__=="__main__":
    temp_rts=select_rts_for_bin(database_address,temp_bin)
    print(len(temp_rts))
    temp_mzs=select_mzs_for_bin(database_address,temp_bin)
    print(len(temp_mzs))