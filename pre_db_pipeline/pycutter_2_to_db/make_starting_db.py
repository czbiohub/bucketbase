#import sqlite3
#from sqlite3 import Error
import sqlalchemy


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

    db_type='sqlite'
    my_database_location="../../../data/database/bucketbase.db"
    #create_database(my_database_location)


    #if db_type=='sqlite':
    connection=create_connection(my_database_location)

    # #example code for reference
    # elif db_type=='postgres':
    #     my_server='localhost'
    #     my_database='binvestigate_first'
    #     my_dialect='postgresql'
    #     my_driver='psycopg2'
    #     my_username='rictuar'
    #     my_password='password'
    #     my_port='5432'

    #     my_connection_string=f'{my_dialect}+{my_driver}://{my_username}:{my_password}@{my_server}:{my_port}/{my_database}'
    #     engine=create_engine(my_connection_string)#,echo=True)
    #     connection=engine.connect()

    connection.execute(
        '''
        CREATE TABLE bins(
            bin_id INTEGER PRIMARY KEY,
            inchikey TEXT,
            adduct TEXT,
            group_id INTEGER,
            is_istd INTEGER,
            is_known INTEGER,
            comment TEXT,
            polarity TEXT
        )
        '''
    )

    connection.execute(
        '''
        CREATE TABLE annotations(
            annotation_id INTEGER PRIMARY KEY,
            precursor_mass REAL,
            retention_time REAL,
            intensity_height REAL,
            spectrum TEXT,
            bin_id INTEGER,
            member_of_consenus INTEGER,
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
            method_id,
            sample_description_id,
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
