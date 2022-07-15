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