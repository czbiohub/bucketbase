#import sqlite3
#from sqlite3 import Error
import sqlalchemy
import sys
import os


# def create_database(database_address):
#     '''
#     Create the databse. will crash if there is no connection made.
#     But hey at least you get the error
#     '''
#     try:
#         connection=sqlite3.connect(database_address)
#         #print('here')
#     except Error as error:
#         print(error)
#     return connection

def create_connection(database_address):
    ## sqlite://<nohostname>/<path>
    engine=sqlalchemy.create_engine(f"sqlite:///{my_database_location}")
    connection=engine.connect()
    return connection

if __name__ =="__main__":

    
    #main_db_location=sys.argv[1]
    study_name=sys.argv[1]
    to_transient_for_pycutter_pipeline=sys.argv[2]
    #main_db_exists=sys.argv[4]


    if to_transient_for_pycutter_pipeline=='transient':
        os.system(f'mkdir -p ../data/{study_name}/database/')

    my_database_location=f"../data/{study_name}/database/{to_transient_for_pycutter_pipeline}_bucketbase.db"



    db_type='sqlite'
    
    connection=create_connection(my_database_location)

    connection.execute(
        '''
        CREATE TABLE bins(
            bin_id INTEGER PRIMARY KEY,
            english_name TEXT,
            inchikey TEXT,
            adduct TEXT,
            group_id INTEGER,
            is_istd INTEGER,
            is_known INTEGER,
            comment TEXT,
            polarity TEXT,
            valid_for_autocuration INTEGER,
            consensus_rt REAL,
            consensus_mz REAL,
            consensus_spectrum TEXT,
            mzrt_only INTEGER
        )
        '''
    )

    connection.execute(
        '''
        CREATE TABLE annotations(
            annotation_id INTEGER PRIMARY KEY,
            precursor_mz REAL,
            retention_time REAL,
            intensity_height REAL,
            spectrum TEXT,
            bin_id INTEGER,
            adduct TEXT,
            member_of_consensus INTEGER,
            run_id TEXT,
            FOREIGN KEY (bin_id)
                REFERENCES bins (bin_id),
            FOREIGN KEY (run_id)
                REFERENCES runs (run_id)
        )
        '''
    )

    #method_id and sample_id should be foreign keys of 
    #tables that extend the work that sarah is doing with bulk loader
    connection.execute(
        '''
        CREATE TABLE runs(
            run_id TEXT PRIMARY KEY,
            method_id TEXT,
            sample_des_id TEXT,
            run_type TEXT
        )
        '''
    )

    # connection.execute(
    #     '''
    #     INSERT INTO bins
    #     VALUES (
    #         0, 'ABCDEFG-UOSOHZXAS-N', '[M+H]+', null, 0,1, 'this is my first row!'
    #     )
    #     '''
    # )
