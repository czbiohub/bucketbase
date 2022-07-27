import numpy as np
import pandas as pd
import sys
sys.path.insert(0, '../utils/')
from utils import execute_query_connection_established
from utils import parse_text_spectra_return_pairs
import sqlalchemy
import spectral_entropy

'''
the purpose of this is to provide some level of automatic curation to the post-pycutter-step-1 file
at this point,we have already taken the pycutter step 1 output and put it into a transient databse
we go through every "bin" in the transient database, and attempt curations by comparing with the main databsae
at current, we select the entire bins table from the transient databse to make this faster

i am deliberating on whether to select the entire main databse bin table. i will go against that strategy for now
but i expect that it would be significantly faster
i compromise by making a persistent connection (as i always should have) and then selecting a subset
'''

'''
essentially, we go through every single bin in the transient databsae. for each bin, we attempt to match
with bins from the permanent database. matches are made based on same adduct, polarity, mz,rt within tolerance, etc

for each candidate, we take the dot product/reverse dot product. we take the average of those similarities.

#we then store the "candidates" as a dataframe as an element in the transient-to-be-identified dataframe. we also
store the bin_id of the top match for convenience in writing to excel (note, this could be a tie. morstly
for development at this point)
'''


def select_all_bins(a_database_connection,ion_mode):
    '''
    get the entire bins table from the specified databse
    '''
    query=f'''PRAGMA table_info(bins)'''
    bins_column_reply=execute_query_connection_established(a_database_connection,query)
    bin_columns=[element[1] for element in bins_column_reply]
    
    query=f'''
    select * from bins
    where
    (polarity = "{ion_mode}")
    order by
    bin_id'''
    #note, recently removed AND (valid_for_autocuration = 1)
    bins_reply=execute_query_connection_established(a_database_connection,query)
    bins_df=pd.DataFrame.from_records(
        data=bins_reply,
        columns=bin_columns
    )

    #print(bins_df)
    return bins_df
    #return [element[0] for element in execute_query(database_address,query)]

def select_candidate_bins(
    a_database_connection,
    #i guess we really dont need ion mode if we specify adduct
    ion_mode,
    transient_adduct,
    transient_mz,
    transient_rt,
    mz_tolerance,
    rt_tolerance
):

    query=f'''PRAGMA table_info(bins)'''
    bins_column_reply=execute_query_connection_established(a_database_connection,query)
    bin_columns=[element[1] for element in bins_column_reply]

    mz_min=transient_mz-mz_tolerance
    mz_max=transient_mz+mz_tolerance
    rt_min=transient_rt-rt_tolerance
    rt_max=transient_rt+rt_tolerance   

    query=f'''
    select * from bins
    where
    (valid_for_autocuration=1) AND
    (polarity="{ion_mode}") AND
    (adduct="{transient_adduct}") AND
    ((consensus_mz>{mz_min}) AND (consensus_mz<{mz_max})) AND
    ((consensus_rt>{rt_min}) AND (consensus_rt<{rt_max}))
    '''
    bins_reply=execute_query_connection_established(a_database_connection,query)
    if bins_reply=='no results found':
       return 'no results found'
    bins_df=pd.DataFrame.from_records(
        data=bins_reply,
        columns=bin_columns
    )

    return bins_df


#can refactor as an apply if desired once we see what columns we weant
#expect a 4x speedup compared to iterrows
def autocurate_transient_bins(
    transient_bins_panda,
    database_connection,
    ion_mode,
    weight_of_dot_product,
    weight_of_reverse_dot_product,
    mz_tolerance,
    rt_tolerance,
    ms2_tolerance,
    output_panda
):
    '''
    '''

    
    transient_bins_panda['top_match_bin_id']='no_match_found'
    transient_bins_panda['top_match_english_name']='no_match_found'
    transient_bins_panda['matching_bins_text']='no_match_found'
    transient_bins_panda['matching_bins_df']='no_match_found'

    #select_database_bins=
    for row,series in transient_bins_panda.iterrows():
        english_name=series['english_name']
        transient_bin_id=series['bin_id']
        print(f'we are trying to curate {english_name} which has transient bin_id {transient_bin_id}') 
        if series['valid_for_autocuration']==0:
            print('it is not valid for autocuration, so we continue')
            transient_bins_panda.at[transient_bin_id,'top_match_bin_id']='not valid for autocuration'
            transient_bins_panda.at[transient_bin_id,'top_match_english_name']='not valid for autocuration'
            transient_bins_panda.at[transient_bin_id,'matching_bins_df']='not valid for autocuration'
            transient_bins_panda.at[transient_bin_id,'matching_bins_text']='not valid for autocuration'
            continue
        
        candidate_bins_df=select_candidate_bins(
            database_connection,
            ion_mode,
            series['adduct'],
            series['consensus_mz'],
            series['consensus_rt'],
            mz_tolerance,
            rt_tolerance
        )
        english_name=series['english_name']
        transient_bin_id=series['bin_id']
        print(f'we are trying to curate {english_name} which has transient bin_id {transient_bin_id}')

        # candidates_bins_df['dot_product']=cadidates_bin_df.apply(
        #     get_dot_product_of_largest_clusters,axis=1
        # )
        if isinstance(candidate_bins_df,pd.DataFrame):
            #new_columns_panda=candidate_bins_df.apply(
            # for i in candidate_bins_df['consensus_spectrum']:
            #     candidate_bins_df[['dot_product','reverse_dot_product']]=candidate_bins_df.apply(
            #     lambda x: get_dot_product_of_largest_clusters(i,series['consensus_spectrum'],ms2_tolerance),
            #     axis=1,
            #     result_type='expand'
            #     )

            candidate_bins_df[['dot_product','reverse_dot_product']]=candidate_bins_df.apply(
                lambda x: get_dot_product_of_largest_clusters(x['consensus_spectrum'],series['consensus_spectrum'],ms2_tolerance),
                axis=1,
                result_type='expand'
            )

            #candidate_bins_df=pd.concat([candidate_bins_df,new_columns_panda],axis='columns')

            candidate_bins_df['avg_similarity']=weight_of_dot_product*candidate_bins_df['dot_product']+weight_of_reverse_dot_product*candidate_bins_df['reverse_dot_product']
            candidate_bins_df.sort_values(by='avg_similarity',inplace=True,ascending=False)
            
            
            
            
            print('we found:')
            print(candidate_bins_df)
            print(candidate_bins_df.columns)
            #hold=input('hold')

            transient_bins_panda.at[transient_bin_id,'top_match_bin_id']=candidate_bins_df.at[0,'bin_id']
            transient_bins_panda.at[transient_bin_id,'top_match_english_name']=candidate_bins_df.at[0,'english_name']
            transient_bins_panda.at[transient_bin_id,'matching_bins_df']=candidate_bins_df
            temp_json=candidate_bins_df.to_json(orient='records')
            #Microsoft Excel has a character limit of 32,767 characters in each cell.
            if len(temp_json)>32000:
                temp_json=temp_json[:32000]
            transient_bins_panda.at[transient_bin_id,'matching_bins_text']=temp_json


        else:
            continue

    return transient_bins_panda

def guide_auto_curation_wrapper(
        database_connection,
        transient_database_connection,
        ion_mode,
        weight_of_dot_product,
        weight_of_reverse_dot_product,
        mz_tolerance,
        rt_tolerance,
        ms2_tolerance,
        output_panda
    ):
    '''
    we hae this lil mini function so that we can pseudo-vectrize (as apply) and/or multip
    rocess autocurate_transient_bins if we desire
    '''
    transient_bins_panda=select_all_bins(transient_connection,ion_mode)
    #database_bins_panda=select_all_bins(database_connection,ion_mode)
    print(transient_bins_panda)
    print(transient_bins_panda.columns)
    transient_database_auto_curated_as_df=autocurate_transient_bins(
        transient_bins_panda,
        database_connection,
        ion_mode,
        weight_of_dot_product,
        weight_of_reverse_dot_product,
        mz_tolerance,
        rt_tolerance,
        ms2_tolerance,
        output_panda
    )
    return transient_database_auto_curated_as_df

def get_dot_product_of_largest_clusters(transient_spectrum_text,database_spectrum_text,ms2_tolerance):
    '''
    '''

    #print('-'*50)
    ##print(transient_spectrum_text)
    #print(database_spectrum_text)

    #we send a list because thats what the util function expects, whoops.
    candidate_spectrum=parse_text_spectra_return_pairs(
        [transient_spectrum_text.split('@')[0]]
    )
    database_spectrum=parse_text_spectra_return_pairs(
        [database_spectrum_text.split('@')[0]]
    )


    #print(database_spectrum[0])
    #print(candidate_spectrum[0])
    #eprint(ms2_tolerance)

    #print(database_spectrum)
    dot_product=spectral_entropy.similarity(
                    #we select an element because of the above problem. whoops.
                    database_spectrum[0], 
                    candidate_spectrum[0], 
                    method='dot_product',
                    ms2_da=ms2_tolerance,
                    need_clean_spectra=False,
                )
    #print(dot_product)

    reverse_dot_product=spectral_entropy.similarity(
                    database_spectrum[0],
                    candidate_spectrum[0],
                    method='dot_product_reverse',
                    ms2_da=ms2_tolerance,
                    need_clean_spectra=False,
                )
    #print(reverse_dot_product)

    return dot_product,reverse_dot_product

def read_pycutter_step_1_input(pycutter_autocurated_input_address,inserted_columns):
    input_panda=pd.read_csv(pycutter_autocurated_input_address,sep='\t')
    # Define 5th row as column header
    input_panda.columns = input_panda.iloc[3]
    #temporarily lose the top rows. we will reattach them at the end with a concatenation
    input_panda.drop(input_panda.index[0:4], inplace=True)
    #print(input_panda)
    #print(input_panda.columns)
    inchikey_location=input_panda.columns.get_loc('Adduct type')
    for i in range(len(inserted_columns)-1,-1,-1):
        input_panda.insert(loc=inchikey_location+1,column=inserted_columns[i],value=np.nan)
        #input_panda.insert(loc=0,column=inserted_columns[i],value=np.nan)
        ##input_panda.at[3,inserted_columns[i]]=inserted_columns[i]
    #print(inchikey_location)
    #print(input_panda.columns)
    #print(input_panda[inserted_columns])
    #hold=input('hold')
    input_panda.reset_index(drop=True,inplace=True)
    #input_panda.insert()
    return input_panda

def insert_autocurations_into_pycutter_input(pycutter_input,transient_database_auto_curated_as_df,pycutter_step_1_output_address,inserted_columns):
    '''
    attach the results form the transient database autocuration
    then, attach the top 3 rows of the original pycutter step 1 output

    '''
    pycutter_input['bin_id']=transient_database_auto_curated_as_df['top_match_bin_id']
    pycutter_input['english_name']=transient_database_auto_curated_as_df['top_match_english_name']
    pycutter_input['curation_text']=transient_database_auto_curated_as_df['matching_bins_text']

    temp=pd.read_csv(pycutter_step_1_output_address,sep='\t',nrows=5,header=None)
    print(temp)
    print(temp.iloc[4])
    temp.columns = temp.iloc[4]
    temp.drop(list(range(5,temp.index.stop)),inplace=True,axis='index')
    inchikey_location=temp.columns.get_loc('Adduct type')
    for i in range(len(inserted_columns)-1,-1,-1):
        temp.insert(loc=inchikey_location+1,column=inserted_columns[i],value=np.nan)
    temp.iloc[4]=temp.columns
    print(temp)
    print(pycutter_input)
    hold=input('just before conczt ')
    final_result=pd.concat([temp,pycutter_input],axis='index',ignore_index=True)

    return final_result

if __name__=="__main__":

    ion_mode='pos'
    weight_of_dot_product=0.5
    weight_of_reverse_dot_product=0.5
    mz_tolerance=0.01
    rt_tolerance=1
    ms2_tolerance=0.015


    transient_database_address='../../data/database/transient_bucketbase.db'
    database_address='../../data/database/bucketbase.db'
    
    pycutter_autocurated_input_address='../../data/BRYU005_pipeline_test/step_1_post_pycutter/py_cutter_step_1_output.tsv'
    pycutter_autocurated_output_address='../../data/BRYU005_pipeline_test/step_1_b_post_auto_curation/pycutter_step_1_autocurated.tsv'
    inserted_columns=['comment','bin_id','english_name','curation_text']
    pycutter_input=read_pycutter_step_1_input(pycutter_autocurated_input_address,inserted_columns)
    #recreate the same bin IDs that we see in the transient databse.
    #https://github.com/czbiohub/bucketbase/issues/30
    #print(pycutter_input.index)
    #print(list(range(0,pycutter_input.index.stop)))
    pycutter_input['bin_id']=list(range(0,pycutter_input.index.stop))
    print(pycutter_input.columns)
    hold=input('just after read in')

    transient_engine=sqlalchemy.create_engine(f"sqlite:///{transient_database_address}")
    transient_connection=transient_engine.connect()

    database_engine=sqlalchemy.create_engine(f"sqlite:///{database_address}")
    database_connection=database_engine.connect()

    transient_database_auto_curated_as_df=guide_auto_curation_wrapper(
        database_connection,
        transient_connection,
        ion_mode,
        weight_of_dot_product,
        weight_of_reverse_dot_product,
        mz_tolerance,
        rt_tolerance,
        ms2_tolerance,
        pycutter_input
        )

    print(transient_database_auto_curated_as_df)
    print(pycutter_input.columns)
    hold=input('hold')

    #print(transient_database_auto_curated_as_df)
    
    transient_connection.close()
    transient_engine.dispose()

    database_connection.close()
    database_engine.dispose()

    #print(pycutter_input)
    final_result=insert_autocurations_into_pycutter_input(pycutter_input,transient_database_auto_curated_as_df,pycutter_autocurated_input_address,inserted_columns)
    final_result.to_csv(pycutter_autocurated_output_address,sep='\t',index=False,header=False)
    print(final_result)
    #we now need to merge the autocurated transient DB panda with the pycutter input panda