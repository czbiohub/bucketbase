# import numpy as np
# import pandas as pd
# import sqlalchemy
import os
import numpy as np

from prepare_bin_table_upload import * #find_lowest_bin,create_bin_table_upload

def find_lowest_annotation(database_address):
    '''
    if the annotation number is not present, then we must find the lowest valid annotation numebr by querying the db
    if the db just got created, infer -1 as the lowest annotation number
    so that the numbering starts at 0
    '''
    engine=sqlalchemy.create_engine(f"sqlite:///{database_address}")
    connection=engine.connect()

    temp_cursor=connection.execute(
        '''
        select annotation_id from annotations
        order by annotation_id desc
        limit 1
        '''
    )

    temp_result=temp_cursor.fetchall()

    connection.close()
    engine.dispose()

    if len(temp_result)==0:
        return -1
    elif len(temp_result)==1:
        return temp_result[0][0]


def get_alignment_id_bin_id_map(alignment_panda,bin_panda):
    '''
    returns a panda with one column alignment_id and one column bin_id

    '''
    
    #alignment_id_bin_id_panda=alignment_panda['alignment_id'].copy().to_frame()
    #alignment_id_bin_id_panda['bin_id']='junk'
    alignment_id_bin_id_panda=alignment_panda.loc[
        alignment_panda.bin_id.astype(str).str.isdigit()==False,
        [
            'alignment_id',
        ]
    ].copy()
    alignment_id_bin_id_panda.reset_index(inplace=True,drop=True)

    alignment_id_bin_id_panda['bin_id']=bin_panda['bin_id']
    print(alignment_panda)
    print(alignment_id_bin_id_panda)


    alignment_id_bin_id_panda_2=alignment_panda.loc[
        alignment_panda.bin_id.astype(str).str.isdigit()==True,
        [
            'alignment_id','bin_id'
        ]
    ].copy()
    alignment_id_bin_id_panda_2.reset_index(inplace=True,drop=True)

    print(alignment_id_bin_id_panda_2)
    # alignment_id_bin_id_panda=alignment_id_bin_id_panda.loc[
    #     :,
    #     ['alignment_id','bin_id']
    # ]

    # alignment_id_bin_id_panda=alignment_id_bin_id_panda.merge(
    #     bin_panda,
    #     how='inner',
    #     left_on='alignment_id',
    #     right_on='alignment_id'
    # )

    # print(alignment_panda)
    # print(bin_panda)
    # print(bin_panda.columns)

    # alignment_id_bin_id_panda=alignment_id_bin_id_panda.merge(
    #     bin_panda,
    #     how='inner',
    #     left_on='alignment_id',
    #     right_on='alignment_id'
    # )
    
    # #bin_panda['bin_id'].copy()
    # print(alignment_id_bin_id_panda)

    # temp=bin_panda_already_binned=alignment_panda.loc[
    #     alignment_panda.bin_id.astype(str).str.isdigit()==True,
    #     [
    #         'bin_id','inchikey','adduct','english_name',
    #         'comment','polarity'
    #     ]
    # ]
    total_alignment_id_bin_id_panda=pd.concat(
        [alignment_id_bin_id_panda,alignment_id_bin_id_panda_2],
        axis='index',
        ignore_index=True
    )
    print('assessing alignment id bin id panda')
    hold=input('hold')
    return total_alignment_id_bin_id_panda

def clean_mapping_panda(mapping_panda,alignment_id_bin_id_panda):
    '''
    we simplify the mapping file by
    0) swapping column names
    1) removing all alignment ids that are not in the output of the final pycutter step
    2) we remove all columns except the alignment ids and samples
    '''
    
    
    column_swap_dict={
        'Alignment ID':'alignment_id',
    }
    mapping_panda.rename(
        mapper=column_swap_dict,
        inplace=True,
        axis='columns'
    )    
    
    mapping_panda=mapping_panda.loc[
        mapping_panda.alignment_id.isin(set(alignment_id_bin_id_panda.alignment_id.tolist())),
        :
    ]
    mapping_panda=mapping_panda.reset_index(drop=True)

    index_of_msms_spectrum_column=mapping_panda.columns.tolist().index('MS/MS spectrum')
    #0 is the location of the alignment_id property. we want to keep that and the samples, only
    mapping_panda=mapping_panda.iloc[
        :,
        [0]+list(range(index_of_msms_spectrum_column+1,len(mapping_panda.columns)))
    ]


    mapping_panda=mapping_panda.merge(
        right=alignment_id_bin_id_panda,
        on='alignment_id',
        how='inner'
    )

    #put bin_id at the front
    cols = list(mapping_panda.columns)
    cols = [cols[-1]] + cols[:-1]
    mapping_panda = mapping_panda[cols]


    return mapping_panda

def create_annotation_table_wrapper(individual_files_directory,mapping_panda,database_address):
    '''
    This function, for each individual file, creates an annotation upload panda.
    It then combines them to one gigantic annotation upload panda
    then assigns the annotation_id by finding the smallest ID from the sqlite DB
    '''
    file_list=os.listdir(individual_files_directory)

    annotation_panda_upload_list=list()
    for i,temp_file in enumerate(file_list):
        #if i>5 and i<255:
        #    continue
        print(i)
        annotation_panda_upload_list.append(
            create_annotation_table_one_individual_file(individual_files_directory,temp_file,mapping_panda)
        )

    annotation_upload_panda=pd.concat(
        annotation_panda_upload_list,
        axis='index',
        ignore_index=True,
    ) 
    lowest_annotation_id=1+find_lowest_annotation(database_address)

    annotation_upload_panda['annotation_id']=list(range(
        lowest_annotation_id,(lowest_annotation_id+len(annotation_upload_panda.index))
    ))

    return annotation_upload_panda

def normalize_spectrum(spectrum):
    '''
    divides each intensity by the  max of the intensities
    '''
    spectrum[1]=spectrum[1]/(spectrum[1].max())
    return spectrum

def parse_one_ms_dial_spectrum(spectrum_text):
    '''
    takes an ms/ms spectrum as a string and returns an array
    '''

    mz_int_pair_list=spectrum_text.split(' ')
    mz_list=[float(temp_pair.split(':')[0]) for temp_pair in mz_int_pair_list]
    intensity_list=[float(temp_pair.split(':')[1]) for temp_pair in mz_int_pair_list]
    
    return np.array([
        mz_list,
        intensity_list
    ])

def convert_np_spectrum_to_text(spectrum_array):
    '''
    takes an ms/ms spectrum in our numpy array format and converts it to string
    '''
    string_rep=[
        str(spectrum_array[0][i])+':'+str(spectrum_array[1][i]) for i in range(len(spectrum_array[0]))
    ]
    return ' '.join(string_rep)


def create_annotation_table_one_individual_file(individual_files_directory,temp_file,mapping_panda):    
    '''
    The order of events is to 
    1) take subset of the mapping panda. keep only bin_id, alignment_id, and the particular sample that we want (peak ids)
    2) read in the individual file and clean up the columns a little bit. drop junk and make names easier ot work with
    3) merge the mapping panda and the individual file panda. the "invididual file column" in the mapping panda
    provides the PeakID that we need to access in the individual file
    4) remove unnecessary stuff like the peakID and annotation ID
    5) add the run_id
    '''
    run_id=temp_file[:-4]

    temp_annotation_upload_panda=mapping_panda.loc[
        :,
        #we want to remove the .txt from the file address
        #the last column is the peak ids
        ['bin_id','alignment_id',run_id]
    ]

    individual_file_panda=pd.read_csv(individual_files_directory+temp_file,sep='\t')
    
    column_swap_dict={
        'Height':'intensity_height',
        'RT (min)':'retention_time',
        'MSMS spectrum':'spectrum',
        'Adduct':'adduct',
        'Precursor m/z':'precursor_mz',
        'PeakID':'peak_id'
    }
    individual_file_panda=individual_file_panda.loc[
        :,
        column_swap_dict.keys()
    ]

    individual_file_panda.rename(
        mapper=column_swap_dict,
        inplace=True,
        axis='columns'
    ) 


    temp_annotation_upload_panda=temp_annotation_upload_panda.merge(
        right=individual_file_panda,
        left_on=run_id,
        right_on='peak_id',
        #i think we chose left so that "peak not prsent" goes into the databse
        how='left'
    )

    temp_annotation_upload_panda['run_id']=run_id
    temp_annotation_upload_panda['member_of_consensus']=0
    temp_annotation_upload_panda.drop(
        ['alignment_id','peak_id',run_id],
        inplace=True,
        axis='columns'
    )

    #map, like apply, is not true vectorizaiton
    #opted for map to avoid putting nan case handlers in called functions
    #the idea is to normalize spectra then return them to text for upload
    temp_annotation_upload_panda.loc[
        temp_annotation_upload_panda['spectrum'].notna(),'spectrum'
    ]=temp_annotation_upload_panda.loc[
        temp_annotation_upload_panda['spectrum'].notna(),'spectrum'
    ].apply(parse_one_ms_dial_spectrum)

    temp_annotation_upload_panda.loc[
        temp_annotation_upload_panda['spectrum'].notna(),'spectrum'
    ]=temp_annotation_upload_panda.loc[
        temp_annotation_upload_panda['spectrum'].notna(),'spectrum'
    ].apply(normalize_spectrum)

    temp_annotation_upload_panda.loc[
        temp_annotation_upload_panda['spectrum'].notna(),'spectrum'
    ]=temp_annotation_upload_panda.loc[
        temp_annotation_upload_panda['spectrum'].notna(),'spectrum'
    ].apply(convert_np_spectrum_to_text)

    return temp_annotation_upload_panda


if __name__=="__main__":
    '''
    the "secret insight" is realizing that the alignment_panda_upload gives
    bin_id:alignment_id
    and that the mapping file gives
    alignment_id:peak_id in individual files
    in this way, we connect bins with rows in the individual files
    '''
    
    final_alignment_address='../../../data/BRYU005_pipeline_test/step_2_final_alignment/BRYU005_CombineSubmit_June2022_pos.txt'
    database_address='../../../data/database/bucketbase.db'
    alignment_panda=pd.read_csv(final_alignment_address,sep='\t',skiprows=3)
    bin_panda=create_bin_table_upload(alignment_panda,database_address)
    alignment_id_bin_id_panda=get_alignment_id_bin_id_map(alignment_panda,bin_panda)

    mapping_file_address='../../../data/three_studies/alignment_individual_mappings/BRYU005_pos_mapping.txt'
    mapping_panda=pd.read_csv(mapping_file_address,sep='\t',skiprows=4)
    mapping_panda=clean_mapping_panda(mapping_panda,alignment_id_bin_id_panda)
    

    individual_files_directory='../../../data/three_studies/individual_sample_data_subset/unzipped/BRYU005_Bacterial_Supernatant/pos/'
    annotation_panda_for_upload=create_annotation_table_wrapper(individual_files_directory,mapping_panda,database_address)
    print(annotation_panda_for_upload)
