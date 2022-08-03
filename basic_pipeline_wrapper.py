import os
#there is a question of how to do this if we are doing the first study
#there is a computationally efficient approach where wedo something special like copy the db
#maybe that is true, but for the moment, i think that we ill take an approach wehre
#we 

#these parameters are for the user to decide
#there are other parameters at other places in the pipeline (like uploading to transient or main db)
#that are not "chosen" by the user, but instead are "chosen" based on where we are in the pipeline
parameter_dict={
    
    #'main_db_exists':'no', #yes or no
    #'main_db_location':'../../data/database/', #irrelevant if main_db_exists==no. could be implicity, but why not centralize location
    #                                           #in case it changes?
    #opted against because different scripts have different relative locations
    'study_name':'BRYU005',

}

#parameters decided by location
#to_transient_for_pycutter_pipeline 'transient' and 'main'
#pycutter's step 1 vs step 2

#run pycutter step 1
os.system(
    'python3 ./pycutter/PyCutterProcessing.py BRYU005 one pos'
)
os.system(
    'python3 ./pycutter/PyCutterProcessing.py BRYU005 one neg'
)

#make the initial database
os.system(
    'python3 ./create_database/make_starting_db.py BRYU005 transient'
)

#fill the initial database
os.system(
    'python3 ./create_database/upload_to_db.py BRYU005 transient pos'
)
os.system(
    'python3 ./create_database/upload_to_db.py BRYU005 transient neg'
)

