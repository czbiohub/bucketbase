import numpy as np
import pandas as pd
import sys
sys.path.insert(0, '../utils/')
from utils import execute_query
from utils import parse_text_spectra_return_pairs
import spectral_entropy
from scipy.spatial import distance
from scipy.cluster import hierarchy
from collections import Counter
from scipy.stats import entropy
from generate_consensus_spectra import *

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
                    need_clean_spectra=False,
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

    #print(cluster_identities)
    #hold=input('cluster identities')
    return cluster_identities

def get_cluster_membership_ordering(cluster_assignments):
    '''
    '''
    total_number_of_elements=len(cluster_assignments)
    
    cluster_counts=Counter(cluster_assignments)
    cluster_counts=[(key,cluster_counts[key]) for key in cluster_counts]
    cluster_counts_ordered_by_population=[
        x for x in sorted(cluster_counts, key=lambda pair: pair[1],reverse=True)
    ]
    
    clusters_ordered_by_population=[temp[0] for temp in cluster_counts_ordered_by_population]

    biggest_percent=cluster_counts_ordered_by_population[0][1]/total_number_of_elements

    return clusters_ordered_by_population,biggest_percent

def get_mz_range_of_bin_spectra(temp_spectra_paired):
    '''
    
    '''
    temp_spectra_parallel=[np.stack(temp_spectrum,axis=1) for temp_spectrum in temp_spectra_paired]
    all_mz_list=[spectrum[0] for spectrum in temp_spectra_parallel]
    all_mz=np.concatenate(all_mz_list,dtype=object)
    return all_mz.max()-all_mz.min()

def get_cleaned_spectra(temp_spectra_paired,noise_level,ms2_tolerance):
    '''
    '''
    cleaned_spectra=[spectral_entropy.tools.clean_spectrum(
                        temp_spectrum,
                        noise_removal=noise_level,
                        ms2_da=ms2_tolerance
                    ) for temp_spectrum in temp_spectra_paired]
    return cleaned_spectra

def get_bin_entropy(temp_spectra_paired,use_ceiling,ms2_tolerance):
    '''
    first portion is same strategy as generate_consensus_spectrum
    the second half is from find_average_mz_and_intensity
    
    the strategy, unlike making a consensus spectrum, is that we sum a bin and then normlize by count

    if use_ceiling is true then the strategy doesnt matter
    '''
    
    spectra=[np.stack(temp_spectrum,axis=1) for temp_spectrum in temp_spectra_paired]
    spectrum_count=len(spectra)
    
    all_mz_list=[spectrum[0] for spectrum in spectra]
    all_intensity_list=[spectrum[1] for spectrum in spectra]
    
    all_mz=np.concatenate(all_mz_list,dtype=object)
    all_intensity=np.concatenate(all_intensity_list,dtype=object)
    
    bins=np.arange(all_mz.min(),
                   all_mz.max()+ms2_tolerance,
                  (all_mz.max()-all_mz.min())/100
                  )
    
    bin_identities=np.digitize(
        all_mz,
        bins
    ) 

    output_intensity_list=list()
    for i in range(len(bins)):
        interesting_indexes=[
            j for j, element in enumerate(bin_identities) if (element==i)
        ]
        interesting_intensities=all_intensity[interesting_indexes]
        ####output_intensity_list.append(interesting_intensities.mean())
        output_intensity_list.append((interesting_intensities.sum())/spectrum_count)

    output_intensity_list=np.nan_to_num(output_intensity_list)
    if use_ceiling==False:
        return np.exp(entropy(output_intensity_list))
    elif use_ceiling==True:
        return np.exp(entropy(np.ceil(output_intensity_list)))

def valid_for_autocuration_test_wrapper(
    database_address,
    bin_list,
    similarity_metric,
    ms2_tolerance,
    noise_level,
    mutual_distance_for_cluster,
    use_ceiling,
    max_entropy_parameter,
    min_mz_range_parameter,
    largest_cluster_membership_parameter_percent,
    bin_spectra_count_minimum_parameter,
    minimum_percent_present,
    bin_space_tolerance
):

    result_dict={
        'bin_id':[],
        'valid_for_autocuration':[],
        'consensus_spectrum':[]
    }

    for z,temp_bin in enumerate(bin_list):
        #print(result_dict)
        #hold=input('hold')
        #if temp_bin !=42:
        #    continue

        print(f'bin number {temp_bin} iteration number {z}')
        #note that mz and rt are decided whether or not there is an associated spectrum for 
        #that individual annotation
        #however, it will not be equal to the number of runs for that study
        #these spectra, mzs, and rts do not depend on blanks or qs
        #there is probably a more fine-grained way to do this, but we pass on it for now.
        #note that the annotaitons from blanks and qcs still go into the bin, but the bin is only
        #characterized based on samples
        temp_spectra_text=select_spectra_for_bin(database_address,temp_bin)
        print(f'{len(temp_spectra_text)} spectra') #len(temp_spectra_text))

        #there was a problem in that certain modules expect spectra as mz/rt pairs and some (the ones i wrote)
        #expect parallel lists. fixing this would be an easy way to reduce code complexity
        #but the main goal for the moment is a working version by next week
        temp_spectra_paired=parse_text_spectra_return_pairs(temp_spectra_text)
        #we need to clean to make sure that subsequent stuff is homogenous
        temp_spectra_paired_cleaned=get_cleaned_spectra(temp_spectra_paired,noise_level,ms2_tolerance)
        
        #now, the general logic is to perform a series of tests
        #the tests are all kind of different. at each place, if a test is failed, the procedure is the same.
        #get the consensus spectrum for each cluster and set the auto_curate_valid to False
        #the tests (order matters, for example, it doesnt make sense to do entropy test if mz range is very small
        #1) membership percent of largest cluster
        #2) total members of cluster
        #3) spread of mz
        #4) entropy

        #1)
        #get the cluster assignments. if there is only one spectrum, then cheese it and assign the cluster membership manually
        if len(temp_spectra_paired_cleaned)==1:
            cluster_assignments=np.array([1])
            cluster_assignments_sorted_by_membership=np.array([1])
            biggest_cluster_percent=1.0
        else:
            cluster_assignments=perform_hierarchical_clustering_routine(temp_spectra_paired_cleaned,similarity_metric,ms2_tolerance,mutual_distance_for_cluster)
            #count membership and get cluster percent
            cluster_assignments_sorted_by_membership,biggest_cluster_percent=get_cluster_membership_ordering(cluster_assignments)       
        if biggest_cluster_percent<largest_cluster_membership_parameter_percent:
            consensus_spectra_text=generate_consensus_spectra_text_wrapper(
                temp_spectra_paired_cleaned,
                cluster_assignments,
                cluster_assignments_sorted_by_membership,
                ms2_tolerance,
                minimum_percent_present,
                bin_space_tolerance
            )
            result_dict['bin_id'].append(temp_bin)
            result_dict['valid_for_autocuration'].append(0)
            result_dict['consensus_spectrum'].append(consensus_spectra_text)
            continue

        #2)
        if len(cluster_assignments)<bin_spectra_count_minimum_parameter:
            consensus_spectra_text=generate_consensus_spectra_text_wrapper(
                temp_spectra_paired_cleaned,
                cluster_assignments,
                cluster_assignments_sorted_by_membership,
                ms2_tolerance,
                minimum_percent_present,
                bin_space_tolerance
            )
            result_dict['bin_id'].append(temp_bin)
            result_dict['valid_for_autocuration'].append(0)
            result_dict['consensus_spectrum'].append(consensus_spectra_text)       
            continue    

        #3
        mz_range=get_mz_range_of_bin_spectra(temp_spectra_paired_cleaned)
        if mz_range<min_mz_range_parameter:
            consensus_spectra_text=generate_consensus_spectra_text_wrapper(
                temp_spectra_paired_cleaned,
                cluster_assignments,
                cluster_assignments_sorted_by_membership,
                ms2_tolerance,
                minimum_percent_present,
                bin_space_tolerance
            )
            result_dict['bin_id'].append(temp_bin)
            #Note that a fail here still means we curate with it
            result_dict['valid_for_autocuration'].append(1)
            result_dict['consensus_spectrum'].append(consensus_spectra_text)       
            continue    

        #4
        bin_entropy=get_bin_entropy(temp_spectra_paired_cleaned,False,ms2_tolerance)
        if bin_entropy>max_entropy_parameter:
            consensus_spectra_text=generate_consensus_spectra_text_wrapper(
                temp_spectra_paired_cleaned,
                cluster_assignments,
                cluster_assignments_sorted_by_membership,
                ms2_tolerance,
                minimum_percent_present,
                bin_space_tolerance
            )
            result_dict['bin_id'].append(temp_bin)
            result_dict['valid_for_autocuration'].append(0)
            result_dict['consensus_spectrum'].append(consensus_spectra_text)       
            continue  
        elif bin_entropy<max_entropy_parameter:
            consensus_spectra_text=generate_consensus_spectra_text_wrapper(
                temp_spectra_paired_cleaned,
                cluster_assignments,
                cluster_assignments_sorted_by_membership,
                ms2_tolerance,
                minimum_percent_present,
                bin_space_tolerance
            )
            result_dict['bin_id'].append(temp_bin)
            result_dict['valid_for_autocuration'].append(1)
            result_dict['consensus_spectrum'].append(consensus_spectra_text)       
            continue  



    #after going through all bins, we convert to panda and return
    return pd.DataFrame.from_dict(result_dict)