import numpy as np
import sqlalchemy

def execute_query(database_address,query_string,returns_rows=True):
    '''
    execute some query/operation on our DB
    '''
    engine=sqlalchemy.create_engine(f"sqlite:///{database_address}")
    connection=engine.connect()

    temp_cursor=connection.execute(
        query_string
    )

    if returns_rows==True:
        temp_result=temp_cursor.fetchall()
    
        connection.close()
        engine.dispose()

        if len(temp_result)==0:
            return -1
        elif len(temp_result)>1:
            return temp_result

    elif temp_cursor==False:
        connection.close()
        engine.dispose()

def parse_text_spectra_return_pairs(spectra_text):
    '''
    '''
    output_list=list()
    for spectrum in spectra_text:
        mz_int_pair_list=spectrum.split(' ')
        mz_list=[float(temp_pair.split(':')[0]) for temp_pair in mz_int_pair_list]
        intensity_list=[float(temp_pair.split(':')[1]) for temp_pair in mz_int_pair_list]   
        #we make it parallel and then pair for time saving during development
        temp_spectrum_parallel=np.array([mz_list,intensity_list])
        temp_spectrum_pair=np.column_stack(temp_spectrum_parallel)
        output_list.append(temp_spectrum_pair)
    return output_list


def parse_text_spectra_return_parallel(spectra_text):
    '''

    '''
    pass