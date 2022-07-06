import pandas as pd


def find_lowest_bin(self,bin):
    '''
    if the bin number is not present, then we must find the lowest valid bin numebr by querying the db
    '''



if __name__=="__main__":
    final_alignment_address='../../../data/BRYU005_pipeline_test/step_0_raw_from_ms_dial/BRYU005_pos_alignment_raw.txt'

    alignment_panda=pd.read_csv(final_alignment_address,sep='\t',skiprows=4)