import os, sys, subprocess, time
import pandas as pd
import numpy as np

class PyCutterProcessing:

    def __init__(self):

        self.deleted = 0
        self.log = []


    def process_alignment(self, raw_alignment_file, metadata_file):

        """
        Processes MS-DIAL alignment file for metabolomics data curation readiness
        """

        # Convert files into pandas DataFrames
        df = pd.read_table(raw_alignment_file, header=None, index_col=False)
        self.write_to_log('Initiated PyCutter processing for ' + os.path.basename(raw_alignment_file))

        # Define 5th row as column header
        df.columns = df.iloc[4]
        df = df.drop(df.index[4])

        # Check whether data is positive or negative mode
        if '[M+H]+' in df['Adduct type'].astype(str).values:
            _mode = 'Positive'
        elif '[M-H]-' in df['Adduct type'].astype(str).values:
            _mode = 'Negative'
        else:
            _mode = 'Error'

        # Add Algorithm column
        df['Algorithm'] = None

        # Add MSI column
        df['MSI'] = None

        # Add Delta RT column
        df['deltaRT'] = df['Reference RT'].astype(float) - df['Average Rt(min)'].astype(float)

        # Add MS2Score column
        df['Dot product'] = df['Dot product'].fillna(0)
        df['Reverse dot product'] = df['Reverse dot product'].fillna(0)
        df.loc[0:4, 'Dot product'] = np.nan
        df.loc[0:4, 'Reverse dot product'] = np.nan
        df['MS2Score'] = (df['Dot product'].astype(float) + df['Reverse dot product'].astype(float)) / 2

        # Fill NaN values in Fragment presence column
        df['Fragment presence %'] = df['Fragment presence %'].fillna(0)
        df.loc[0:4, 'Fragment presence %'] = np.nan

        # Add SAmax/BKavg column
        sample_cols = np.array(df.iloc[1] == 'Sample')  # Indices of columns of samples
        blank_cols = np.array(df.iloc[1] == 'Blank')  # Indices of columns of blanks

        sample_max = df.iloc[4:, sample_cols].astype(float).max(axis=1)  # Maximum value over sample columns
        blank_average = df.iloc[4:, blank_cols].astype(float).mean(axis=1)  # Average value over blank columns

        df['Sample Max'] = sample_max
        df['SAmax/BKavg'] = sample_max / blank_average
        df['SAmax/BKavg'] = np.where(df['SAmax/BKavg'] == np.inf, df['Sample Max'], df['SAmax/BKavg'])

        # Add Percent CV column for duplicates
        pool_cols = np.array(df.iloc[1] == 'QC')  # Indices of columns of samples
        pool_std = df.iloc[4:, pool_cols].astype(float).std(axis=1)  # Standard deviation over sample columns
        pool_average = df.iloc[4:, pool_cols].astype(float).mean(axis=1)  # Average value over sample columns
        df['Percent CV'] = (pool_std / pool_average) * 100

        # Add Sample Average column for duplicates
        sample_cols = np.array(df.iloc[1] == 'Sample')  # Indices of columns of samples
        sample_average = df.iloc[4:, sample_cols].astype(float).mean(axis=1)  # Average value over sample columns
        df['Sample Average'] = sample_average

        # Calculate weighted MS2 score for duplicate handling
        r = 0.5 * df['Reverse dot product'].astype(float)
        d = 0.3 * df['Dot product'].astype(float)
        f = 0.2 * df['Fragment presence %'].astype(float)
        df['Weighted MS2 Score'] = r + d + f

        self.write_to_log(
            'Added columns: MSI, deltaRT, MS2Score, Weighted MS2Score, SAmax/BKavg, Percent CV, Sample Average, Sample Max')

        # Rearrange columns
        columns_to_rearrange = ["MSI", "Alignment ID", "Average Rt(min)", "Average Mz", "Reference RT", "deltaRT",
                                "Metabolite name", "INCHIKEY", "Adduct type", "Algorithm", "Weighted MS2 Score",
                                "Sample Average", "MS2Score", "Dot product", "Reverse dot product", "SAmax/BKavg",
                                "Percent CV", "Spectrum reference file name", "Post curation result", "Fill %"]

        df = self.movecol(df, cols_to_move=columns_to_rearrange, ref_col='MS/MS assigned', place='Before')

        self.write_to_log('Rearranged columns successfully.')

        # Fill sample ID's from BulkLoader metadata into alignment
        if metadata_file != '':

            mdf = pd.read_csv(metadata_file)

            if _mode == 'Positive':

                negatives = mdf.loc[mdf['Filename'].str.contains(r'Neg', na=False)]
                mdf['Filename'].drop(negatives.index, inplace=True)

            elif _mode == 'Negative':

                positives = mdf.loc[mdf['Filename'].str.contains(r'Pos', na=False)]
                mdf['Filename'].drop(positives.index, inplace=True)

            mdf = mdf.transpose()
            sample_index = 40 + np.count_nonzero(sample_cols)
            df.iloc[3:4, 40:sample_index] = mdf[1:2]

        self.write_to_log('Filled sample ID''s from BulkLoader metadata file successfully')

        # Copy data to second sheet
        df2 = df.copy()
        df2 = df2.reset_index()
        df2.drop('index', 1, inplace=True)

        self.write_to_log('Created Reduced sheet.')

        # Identify annotated features, unknowns, and features without MS2
        df2.loc[(df2['Metabolite name'] != 'Unknown'), 'Feature Type'] = 'Metabolite'
        df2.loc[(df2['INCHIKEY'] == 'Internal Standard'), 'Feature Type'] = 'Internal Standard'
        df2.loc[(df2['Metabolite name'].str.contains(r'w/o MS2', na=False)) & (
                df2['INCHIKEY'] != 'Internal Standard'), 'Metabolite name'] = 'Unknown'
        df2.loc[(df2['Metabolite name'] == 'Unknown') & (
                df2['INCHIKEY'] != 'Internal Standard'), 'Feature Type'] = 'Unknown'
        df2.loc[(df2['Metabolite name'].isnull()), 'Feature Type'] = "Header"

        # Categorize all features by type
        categories = ['Header', 'Feature Type', 'Internal Standard', 'Metabolite', 'Unknown']
        feature_map = {categories[i]: i for i in range(len(categories))}
        df2['Feature Index'] = df2['Feature Type'].map(feature_map)

        self.write_to_log('Categorized internal standards, metabolites, and unknowns')

        # Remove annotated features with SAmax/BKavg < 3 (unless iSTD)
        about_to_drop = df2[(df2['Feature Index'] == 3) & (df2['SAmax/BKavg'] < 3)]
        df2 = df2.drop(about_to_drop.index)

        # LOG WRITING – Remove annotated features with SAmax/BKavg < 3 (unless iSTD)
        self.deleted = self.deleted + len(about_to_drop)
        about_to_drop['Alignment ID'] = about_to_drop['Alignment ID'].astype(str)
        dropped_list = about_to_drop[['Alignment ID', 'Metabolite name']].values.tolist()
        self.write_to_log('Deleting ' + str(len(dropped_list)) + ' annotated features with SAmax/BKavg < 3...')
        self.write_to_log('\n'.join(' '.join(metabo) for metabo in dropped_list))

        # LOG WRITING – Mark annotated features with MS2Score < 70 and no Reference RT as "Unknown"
        unknowned_list = df2[(df2['Feature Index'] == 3)
                             & (df2['Reference RT'].isnull())
                             & (df2['MS2Score'] < 70)]

        unknowned_list['Alignment ID'] = unknowned_list['Alignment ID'].astype(str)
        unknowned_list = unknowned_list[['Alignment ID', 'Metabolite name']].values.tolist()
        self.write_to_log('Marking ' + str(
            len(unknowned_list)) + ' annotated features with MS2Score < 70 and no Reference RT as "Unknown"...')
        self.write_to_log('\n'.join(' '.join(metabo) for metabo in unknowned_list))

        # Mark annotated features with MS2Score < 70 and no Reference RT as "Unknown"
        df2.loc[(df2['Feature Index'] == 3)
                & (df2['Reference RT'].isnull())
                & (df2['MS2Score'] < 70),
                ['Metabolite name', 'INCHIKEY', 'Feature Index']] = ['Unknown', 'null', 4]

        # LOG WRITING – Mark annotated features with MS2Score < 70 and no RT match as "Unknown"
        unknowned_list = df2[(df2['Feature Index'] == 3)
                             & (df2['deltaRT'].abs() > 0.4)
                             & (df2['MS2Score'] < 70)]

        unknowned_list['Alignment ID'] = unknowned_list['Alignment ID'].astype(str)
        unknowned_list = unknowned_list[['Alignment ID', 'Metabolite name']].values.tolist()
        self.write_to_log('Marking ' + str(
            len(unknowned_list)) + ' annotated features with MS2Score < 70 and no RT match as "Unknown"...')
        self.write_to_log('\n'.join(' '.join(metabo) for metabo in unknowned_list))

        # Mark annotated features with MS2Score < 70 and no RT match as "Unknown"
        df2.loc[(df2['Feature Index'] == 3)
                & (df2['deltaRT'].abs() > 0.4)
                & (df2['MS2Score'] < 70),
                ['Metabolite name', 'INCHIKEY', 'Feature Index']] = ['Unknown', 'null', 4]

        # LOG WRITING – Mark annotated features with MS2Score < 55 and RT match as "Unknown"
        unknowned_list = df2[(df2['Feature Index'] == 3)
                             & (df2['deltaRT'].abs() < 0.4)
                             & (df2['MS2Score'] < 55)]

        unknowned_list['Alignment ID'] = unknowned_list['Alignment ID'].astype(str)
        unknowned_list = unknowned_list[['Alignment ID', 'Metabolite name']].values.tolist()

        self.write_to_log(
            'Marking ' + str(len(unknowned_list)) + ' annotated features with MS2Score < 55 and RT match as "Unknown"...')
        self.write_to_log('\n'.join(' '.join(metabo) for metabo in unknowned_list))

        # Mark annotated features with MS2Score < 55 and RT match as "Unknown"
        df2.loc[(df2['Feature Index'] == 3)
                & (df2['deltaRT'].abs() < 0.4)
                & (df2['MS2Score'] < 55),
                ['Metabolite name', 'INCHIKEY', 'Feature Index']] = ['Unknown', 'null', 4]

        # Remove unknowns in Reduced sheet with SAmax/BKavg < 5
        about_to_drop = df2[(df2['Feature Index'] == 4) & (df2['SAmax/BKavg'] < 5)]
        df2 = df2.drop(about_to_drop.index)

        self.deleted = self.deleted + len(about_to_drop)
        about_to_drop['Alignment ID'] = about_to_drop['Alignment ID'].astype(str)
        dropped_list = about_to_drop[['Alignment ID', 'Metabolite name']].values.tolist()
        self.write_to_log('Deleting ' + str(len(about_to_drop)) + ' unknowns in Reduced sheet with SAmax/BKavg < 5...')
        self.write_to_log(', '.join(' '.join(metabo) for metabo in dropped_list))

        # Remove MS-DIAL summary stats
        df2 = df2.drop(columns=df2.columns[(df2 == 'Average').any()])
        df2 = df2.drop(columns=df2.columns[(df2 == 'Stdev').any()])

        self.write_to_log('Removed MS-DIAL summary stats in Reduced sheet')

        # Move all duplicate features to new DataFrame
        df3 = df2[0:4]
        df2_metabolites = df2[df2['Feature Index'] == 3]
        duplicates = df2_metabolites[df2_metabolites.duplicated(subset=['INCHIKEY'], keep=False)]

        duplicates_list = duplicates.copy()
        duplicates_list['Alignment ID'] = duplicates_list['Alignment ID'].astype(str)
        duplicates_list = duplicates_list[['Alignment ID', 'Metabolite name']].values.tolist()
        self.write_to_log('Capturing ' + str(len(duplicates)) + ' duplicates...')
        self.write_to_log('\n'.join(' '.join(metabo) for metabo in duplicates_list))

        # Sort duplicates
        duplicates.sort_values(by=['INCHIKEY', 'Metabolite name', 'Weighted MS2 Score', 'Sample Average'],
                               ascending=[True, True, False, False], inplace=True)

        # Grabs list of all unique InChIKeys
        entries = duplicates[~duplicates.duplicated(subset=['INCHIKEY'])]
        entries['INCHIKEY'] = entries['INCHIKEY'].fillna('null')
        entries = entries['INCHIKEY'].tolist()

        # For each unique metabolite,
        for entry in entries:

            loopdf = duplicates[duplicates['INCHIKEY'] == entry]

            if entry != 'null':

                # If the difference in Weighted MS2 Score < 5 between duplicate [M+H]+ adducts,
                adducts_only = loopdf[(loopdf['Adduct type'] == '[M+H]+') | (loopdf['Adduct type'] == '[M-H]-')]

                threshold_score = adducts_only['Weighted MS2 Score'].max() - 5

                entries_above_threshold = loopdf[(loopdf['Weighted MS2 Score'] >= threshold_score) &
                                                 ((loopdf['Adduct type'] == '[M+H]+') | (
                                                             loopdf['Adduct type'] == '[M-H]-'))]

                entries_below_threshold = loopdf[((loopdf['Weighted MS2 Score'] >= threshold_score) &
                                                  (loopdf['Adduct type'] != '[M+H]+') & (loopdf['Adduct type'] != '[M-H]-'))
                                                 | (loopdf['Weighted MS2 Score'] < threshold_score)]

                if not entries_above_threshold.empty:

                    # Choose the duplicate with the highest sample average
                    greatest_sample_average = entries_above_threshold['Sample Average'].max()
                    entries_above_threshold.loc[
                        (entries_above_threshold['Sample Average'] == greatest_sample_average), 'Algorithm'] = 'Correct'

                    # Mark duplicates for machine learning
                    entries_above_threshold.loc[
                        (entries_above_threshold['Sample Average'] != greatest_sample_average), 'Algorithm'] = 'Duplicate'
                    entries_below_threshold['Algorithm'] = 'Duplicate'

                    metabolite = entries_above_threshold.append(entries_below_threshold)

                    # Append set of duplicates to sheet
                    df3 = df3.append(metabolite, ignore_index=True)

                else:

                    loopdf['Algorithm'] = 'Duplicate'

                    # Append set of duplicates to sheet
                    df3 = df3.append(loopdf, ignore_index=True)

            else:

                df3 = df3.append(loopdf, ignore_index=True)

            # After duplicate has been rectified, drop it from DataFrame
            duplicates = duplicates.drop(duplicates[duplicates['INCHIKEY'] == entry].index)

        # If duplicate rectifying algorithm misses any duplicates for any reason,
        # they will still be re-appended to df3, the DataFrame containing rectified duplicates
        df3 = df3.append(duplicates, ignore_index=True)
        rectified_duplicates = df3.drop([0, 1, 2, 3])

        # Re-append duplicates to Reduced sheet
        reduced_header = df2[0:4]
        metabolites_no_duplicates = df2_metabolites[
            ~(df2_metabolites.duplicated(subset=['INCHIKEY'], keep=False))]
        internal_standards = df2[df2['Feature Index'] == 2]
        unknowns = df2[df2['Feature Index'] == 4]

        reduced_sheet = pd.concat([metabolites_no_duplicates, unknowns])
        reduced_sheet.sort_values(by=['Feature Index', 'SAmax/BKavg'], inplace=True)

        df2 = pd.concat([reduced_header, internal_standards, rectified_duplicates, reduced_sheet])

        self.write_to_log('Rectified ' + str(len(rectified_duplicates)) + ' duplicates')
        self.write_to_log('Sorted metabolites by SAmax/BKavg')

        # Report columns as integers
        for column in ["Sample Average", "SAmax/BKavg", "Percent CV", "Weighted MS2 Score", "MS2Score",
                       "Dot product", "Reverse dot product"]:
            df2[column] = df2[column].astype(float).round()

        # Reset index
        new_index = pd.Index(range(0, len(df2)))
        df2 = df2.set_index(new_index)

        # Header formatting
        df = self.insert_row(4, df, list(df.columns))
        df2 = self.insert_row(4, df2, list(df2.columns))

        expected = len(df) - self.deleted

        # Final save of formatted spreadsheet
        input_file_name = os.path.basename(raw_alignment_file)
        save_directory = raw_alignment_file.replace(input_file_name, '')

        # TODO – Name the alignment file using project ID
        filename = 'Project_ID.tsv'

        # writer = pd.ExcelWriter(save_directory + filename,
        #                         engine='xlsxwriter',
        #                         engine_kwargs={'options': {'strings_to_numbers': True}})

        # df.to_excel(writer, sheet_name='Raw', header=False, index=False)

        self.write_to_log('Raw sheet size: ' + str(len(df)))

        #df2.to_excel(writer, sheet_name='Reduced', header=False, index=False)
        df2.to_csv(save_directory+filename,sep='\t',header=False,index=False)

        self.write_to_log('Reduced sheet size: ' + str(len(df2)))
        self.write_to_log('Expected Reduced sheet size: ' + str(expected))

        #writer.save()
        self.write_to_log('Processing complete')

        # # Generate post-processing log
        # with open(save_directory + filename.replace('.xlsx', ' ') + 'Log.txt', 'w', encoding='utf-8') as logfile:

        #     for line in self.log:
        #         logfile.write(line)
        #         logfile.write('\n')
        #         logfile.write('\n')

        # self.open_folder(path=save_directory)
        return df2


    def write_to_log(self, text):

        """
        Logging function
        """

        current_time = time.strftime("%H:%M:%S")
        self.log.append(current_time + ' – ' + text)


    def open_folder(self, path):

        """
        Opens folder containing alignment/combined file after processing
        """

        if sys.platform == 'darwin':
            subprocess.check_call(['open', '--', path])
        elif sys.platform == 'linux2':
            subprocess.check_call(['gnome-open', '--', path])
        elif sys.platform == 'win32':
            subprocess.run(['explorer', os.path.realpath(path)])


    def movecol(self, df, cols_to_move, ref_col='', place='After'):

        """
        Rearranges columns in pandas DataFrame
        """

        cols = df.columns.tolist()
        if place == 'After':
            seg1 = cols[:list(cols).index(ref_col) + 1]
            seg2 = cols_to_move
        if place == 'Before':
            seg1 = cols[:list(cols).index(ref_col)]
            seg2 = cols_to_move + [ref_col]

        seg1 = [i for i in seg1 if i not in seg2]
        seg3 = [i for i in cols if i not in seg1 + seg2]

        return (df[seg1 + seg2 + seg3])


    def insert_row(self, row_number, df, row_value):

        """
        Inserts rows in pandas DataFrame
        """

        # Slice the upper half of the dataframe
        df1 = df[0:row_number]

        # Store the result of lower half of the dataframe
        df2 = df[row_number:]

        # Insert the row in the upper half dataframe
        df1.loc[row_number] = row_value

        # Concat the two dataframes
        df_result = pd.concat([df1, df2])

        # Reassign the index labels
        df_result.index = [*range(df_result.shape[0])]

        # Return the updated dataframe
        return df_result


# Use this to run a simple test of the script
if __name__ == "__main__":
    pycutter = PyCutterProcessing()
    pycutter.process_alignment(
        "../../../data/pipeline_test/step_0_raw_from_ms_dial/STQU002_pos_raw_alignment.txt", 
        ""
    )