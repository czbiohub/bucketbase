import numpy as np
import sys
sys.path.insert(0, '../utils/')
from utils import execute_query
from utils import parse_text_spectra_return_pairs
import spectral_entropy
from scipy.spatial import distance
from scipy.cluster import hierarchy

def select_spectra_for_bin(database_address,bin_id):
    query=f'''
    select spectrum 
    from annotations
    inner join
    runs
    on annotations.run_id=runs.run_id
    where (annotations.bin_id={bin_id}) and (annotations.spectrum is not null) and (runs.run_type='Sample')
    '''
    return [element[0] for element in execute_query(database_address,query)]

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

def make_distance_matrix(spectra,similarity_metric,ms2_tolerance):
    '''
    makes a distance matrix. currently very slow because of explicit for loops.
    see this github issue https://github.com/czbiohub/bucketbase/issues/16
    so, its clear that the dot product is assumed to be symmetric here even tho in ms
    the components of the dot product are calculated directionally
    '''
    similarity_matrix=np.zeros(
        shape=(len(spectra),len(spectra))
    )
    
    if similarity_metric=='dot_product':
        for i in range(len(spectra)):
            for j in range(i,len(spectra)):
                similarity_matrix[i][j]=spectral_entropy.similarity(
                    spectra[i], spectra[j], 
                    method='dot_product',
                    ms2_da=ms2_tolerance,
                    need_clean_spectra=True,
                )    
    
    similarity_matrix=np.triu(similarity_matrix)
    similarity_matrix=similarity_matrix+similarity_matrix.T-np.diag(np.diag(similarity_matrix))
    
    ############DANGER##############
    #getting some error where 1 was being rep'd as 0.9999999
    #so i just manually set
    np.fill_diagonal(similarity_matrix,1)

    distance_matrix=1-similarity_matrix
    #the squareform is done because that is what scipys hierarchial routine expects    
    distance_matrix_flattened=distance.squareform(distance_matrix)
    return similarity_matrix,distance_matrix_flattened


def perform_hierarchical_clustering_routine(spectra,similarity_metric,ms2_tolerance,mutual_distance_for_cluster):
    '''
    receives a list of spectra as mz-int pairs and outputs cluster assignments
    '''
    #generate the similarity matrix
    similarity_matrix,distance_matrix_flattened=make_distance_matrix(spectra,similarity_metric,ms2_tolerance)

    #generate the dendrogram. it is symmetric and we arbitrarily chose row
    #we choose average because it feels close to the notion of how we will do
    #consensus spectrum (as the average of spectra) (even tho thats not quite what we do
    row_linkage=hierarchy.linkage(
        distance_matrix_flattened,method='average'
    )


    cluster_identities=hierarchy.fcluster(
        Z=row_linkage,
        t=mutual_distance_for_cluster,
        criterion='distance',
    )

    print(cluster_identities)
    hold=input('cluster identities')




def valid_for_autocuration_test_wrapper(
    database_address,
    bin_list,
    similarity_metric,
    ms2_tolerance,
    noise_level,
    mutual_distance_for_cluster,
    max_entropy_parameter,
    min_mz_range_parameter,
    largest_cluster_membership_parameter,
    consensus_spectrum_routine_parameters
):

    print(bin_list)
    for temp_bin in bin_list:
        
        #note that mz and rt are decided whether or not there is an associated spectrum for 
        #that individual annotation
        #however, it will not be equal to the number of runs for that study
        #these spectra, mzs, and rts do not depend on blanks or qs
        #there is probably a more fine-grained way to do this, but we pass on it for now.
        #note that the annotaitons from blanks and qcs still go into the bin, but the bin is only
        #characterized based on samples
        temp_spectra_text=select_spectra_for_bin(database_address,temp_bin)
        print(len(temp_spectra_text))
        temp_rts=select_rts_for_bin(database_address,temp_bin)
        print(len(temp_rts))
        temp_mzs=select_mzs_for_bin(database_address,temp_bin)
        print(len(temp_mzs))
        
        
        #there was a problem in that certain modules expect spectra as mz/rt pairs and some (the ones i wrote)
        #expect parallel lists. fixing this would be an easy way to reduce code complexity
        #but the main goal for the moment is a working version by next week
        temp_spectra_paired=parse_text_spectra_return_pairs(temp_spectra_text)
        print(temp_spectra_paired)
        cluster_assignments=perform_hierarchical_clustering_routine(temp_spectra_paired,similarity_metric,ms2_tolerance,mutual_distance_for_cluster)
        
        hold=input('hold')

