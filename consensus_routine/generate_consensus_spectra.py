import numpy as np
import math
from collections import Counter

def generate_bin_groups(bins, bin_space_tolerance):
    '''
    the idea is that we want to group together small clusters of bins (in case there)
    is a spread of mz. to do this, we iterate over all of the bins. for each bin,
    if it is very lcose to an already encountered bin, we group them. else, they are
    put in a new, separate bin.
    
    #an m/z can be assigned to multiple m/z groups
    '''
    output_dict=dict()
    for temp_bin in bins:
        found_match=False
        for temp_key in output_dict.keys():
            for temp_value in output_dict[temp_key].copy():
                if (abs(temp_bin-temp_value)<=bin_space_tolerance):
                    output_dict[temp_key].add(temp_bin)
                    found_match=True
        if found_match==False:
            output_dict[temp_bin]={temp_bin}
    
    for temp_key in output_dict.keys():
        for temp_value in output_dict[temp_key].copy():
            output_dict[temp_key].add(temp_value+1)
            output_dict[temp_key].add(temp_value-1)
    
    return output_dict

def find_bin_identities_with_x_percent_present(bin_identities,percent_present,spectrum_count):
    '''
    we want to find bins where the ion is found "at least 80 percent of the time"
    or whatever
    so, we calculate a minimum count from percent_present and spectrum_count
    then find the identity of each bin posseing that count
    '''
    minimum_count=math.floor(spectrum_count*percent_present)

    my_counter=Counter(bin_identities)
    bins_meeting_count=[
        temp_bin for temp_bin in my_counter.keys() if my_counter[temp_bin]>minimum_count
    ]
    return bins_meeting_count

def find_average_mz_and_intensity(
    bin_identities,
    all_mz,
    all_intensities,
    meaningful_bin_groupings
):
    '''
    the idea is that we 
    '''
    output_mz_list=list()
    output_intensity_list=list()

    for temp_key in meaningful_bin_groupings.keys():
        #we get a set of bin identities that we want
        temp_bins_of_interest=meaningful_bin_groupings[temp_key]
        #we map those bin identities to the indexes of bin_identities
        interesting_indexes=[
            i for i,element in enumerate(bin_identities) if (element in temp_bins_of_interest)
        ]
        interesting_mz=all_mz[interesting_indexes]
        interesting_intensities=all_intensities[interesting_indexes]
        output_mz_list.append(interesting_mz.mean())
        output_intensity_list.append(interesting_intensities.mean())
        
    return output_mz_list,output_intensity_list

def convert_paired_spectrum_to_text(spectrum_array):
    '''
    takes an ms/ms spectrum in our numpy array format and converts it to string
    '''
    #print(spectrum_array)
    #hold=input('spectrum_array')
    string_rep=[
        str(spectrum_array[i][0])+':'+str(spectrum_array[i][1]) for i in range(len(spectrum_array))
    ]
    return ' '.join(string_rep)

def normalize_spectrum(spectrum):
    '''
    divides each intensity by the  max of the intensities
    '''
    spectrum[1]=spectrum[1]/(spectrum[1].max())
    return spectrum

def generate_consensus_spectrum(
    temp_spectra_paired_cleaned,
    ms2_tolerance,
    minimum_percent_present,
    bin_space_tolerance
):
    spectra=[np.stack(temp_spectrum,axis=1) for temp_spectrum in temp_spectra_paired_cleaned]
    
    spectrum_count=len(spectra)
    all_mz_list=[spectrum[0] for spectrum in spectra]
    all_intensity_list=[spectrum[1] for spectrum in spectra]
    all_mz=np.concatenate(all_mz_list,dtype=object)
    all_intensity=np.concatenate(all_intensity_list,dtype=object)
    
    bins=np.arange(all_mz.min(),all_mz.max()+ms2_tolerance,ms2_tolerance)
    
    bin_identities=np.digitize(
        all_mz,
        bins
    )
    
    meaningful_bin_identities=find_bin_identities_with_x_percent_present(bin_identities,minimum_percent_present,spectrum_count)
    meaningful_bin_groupings=generate_bin_groups(meaningful_bin_identities,bin_space_tolerance)

    average_mz_list,average_intensity_list=find_average_mz_and_intensity(
        bin_identities,
        all_mz,
        all_intensity,
        meaningful_bin_groupings
    )

    consensus_spectrum_rennormalized=normalize_spectrum(np.array([average_mz_list,average_intensity_list]))
    consensus_spectrum_paired=np.stack(
        consensus_spectrum_rennormalized,
        axis=1
    )

    consensus_spectrum_paired_as_text=convert_paired_spectrum_to_text(consensus_spectrum_paired)
    return consensus_spectrum_paired_as_text

def generate_consensus_spectra_text_wrapper(
        temp_spectra_paired_cleaned,
        cluster_assignments,
        cluster_assignments_sorted_by_membership,
        ms2_tolerance,
        minimum_percent_present,
        bin_space_tolerance
    ):
    '''
    for each cluster, make a consensus spectra in text format,
    then concatentate them and return that
    '''

    consensus_spectra_text_list=list()

    for temp_cluster_identity in cluster_assignments_sorted_by_membership:
        
        #subset spectra according to their cluster assignments
        #we need only those whose assignment is equal to the current "cluster number"
        temp_spectra_for_consensus=[
            temp_spectrum for i,temp_spectrum in enumerate(temp_spectra_paired_cleaned) if cluster_assignments[i]==temp_cluster_identity 
        ]
        
        temp_consensus_spectrum_text=generate_consensus_spectrum(
            temp_spectra_for_consensus,
            ms2_tolerance,
            minimum_percent_present,
            bin_space_tolerance
        )

        consensus_spectra_text_list.append(temp_consensus_spectrum_text)

    return '@'.join(consensus_spectra_text_list)
