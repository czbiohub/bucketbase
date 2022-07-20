import os, sys, subprocess, time
import pandas as pd
import numpy as np

class PyCutterProcessing:

    """
    Class for processing raw MS-DIAL alignment export for metabolomics data curation readiness
    """

    def __init__(self):

        self.deleted = 0
        self.log = []


    def process_alignment(self, raw_alignment_file, metadata_file=''):

        """
        Step 1: processes raw MS-DIAL alignment file for metabolomics data curation readiness
        """

        self.log = []
        self.write_to_log('Initiated PyCutter processing for ' + os.path.basename(raw_alignment_file))

        # Convert files into pandas DataFrames
        df = pd.read_table(raw_alignment_file, header=None, index_col=False)

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

        self.write_to_log('Added columns: MSI, deltaRT, MS2Score, Weighted MS2Score, SAmax/BKavg, Percent CV, Sample Average, Sample Max')

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
        self.write_to_log('Marking ' + str(len(unknowned_list)) + ' annotated features with MS2Score < 70 and no Reference RT as "Unknown"...')
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
        self.write_to_log('Marking ' + str(len(unknowned_list)) + ' annotated features with MS2Score < 70 and no RT match as "Unknown"...')
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

        self.write_to_log('Marking ' + str(len(unknowned_list)) + ' annotated features with MS2Score < 55 and RT match as "Unknown"...')
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

        # Final save of formatted spreadsheet
        filename = 'Project_ID.tsv'
        df2.to_csv(filename, sep='\t', header=False, index=False)

        # Final log updates
        expected = len(df) - self.deleted
        self.write_to_log('Raw sheet size: ' + str(len(df)))
        self.write_to_log('Reduced sheet size: ' + str(len(df2)))
        self.write_to_log('Expected Reduced sheet size: ' + str(expected))
        self.write_to_log('Processing complete')
        step_one_log = self.log.copy()

        return {
            "Raw": df,
            "Reduced": df2,
            "Log": step_one_log
        }


    def combine_annotations(self, pos_alignment, neg_alignment):

        """
        Step 2: combine annotations from positive/negative alignments into a single uniform dataset
        """

        self.log = []
        self.write_to_log("Initiating combine of annotations from curated positive/negative alignments.")

        # Convert files into pandas DataFrames
        df_pos = pd.read_excel(pos_alignment, sheet_name='Reduced', header=None, index_col=False)
        df_neg = pd.read_excel(neg_alignment, sheet_name='Reduced', header=None, index_col=False)

        # Define 5th row as column header
        df_pos.columns = df_pos.iloc[4]
        df_pos.drop(df_pos.index[4], inplace=True)

        df_neg.columns = df_neg.iloc[4]
        df_neg.drop(df_neg.index[4], inplace=True)

        # Remove columns generated by MS-DIAL
        columns_to_remove = ['Reference RT', 'deltaRT', 'Adduct type', 'Total score', 'MS2Score', 'RT similarity',
                             'Dot product', 'Reverse dot product', 'Spectrum reference file name',
                             'Post curation result', 'Fill %', 'MS/MS assigned', 'SAmax/BKavg', 'Percent CV', 'Reference m/z', 'Formula',
                             'Ontology', 'SMILES', 'Annotation tag (VS1.0)', 'RT matched', 'm/z matched', 'MS/MS matched',
                             'Comment', 'Manually modified for quantification', 'Manually modified for annotation',
                             'Algorithm', 'Isotope tracking parent ID', 'Isotope tracking weight number', 'Fragment presence %',
                             'S/N average', 'MS1 isotopic spectrum', 'MS/MS spectrum', 'Sample Average', 'Sample Max']

        for column in columns_to_remove:

            if column in df_pos.columns:
                df_pos.drop(column, inplace=True, axis=1)

            if column in df_neg.columns:
                df_neg.drop(column, inplace=True, axis=1)

        # Add polarity column
        df_pos['Polarity'] = 'pos'
        df_neg['Polarity'] = 'neg'
        df_pos.loc[0:4, 'Polarity'] = np.nan
        df_neg.loc[0:4, 'Polarity'] = np.nan

        # Add Pearson's correlation coefficient column
        df_pos['Pearson'] = 0
        df_neg['Pearson'] = 0
        df_pos.loc[0:4, 'Pearson'] = np.nan
        df_neg.loc[0:4, 'Pearson'] = np.nan

        # Rearrange Polarity
        rearranged_columns = ['Alignment ID', 'Average Rt(min)', 'Average Mz', 'Metabolite name', 'INCHIKEY',
                              'Polarity', 'Pearson', 'Feature Type', 'Feature Index']
        df_pos = self.movecol(df_pos, cols_to_move=rearranged_columns, ref_col='MSI', place='After')
        df_neg = self.movecol(df_neg, cols_to_move=rearranged_columns, ref_col='MSI', place='After')

        self.write_to_log("Removed unnecessary columns, added Polarity and Pearson columns")

        # Mark leftover unknowns
        df_pos['MSI'] = df_pos['MSI'].fillna('')
        df_pos[5:].loc[df_pos['MSI'] == '', ['Metabolite name', 'INCHIKEY',
                                             'Feature Type', 'Feature Index']] = ['Unknown', 'null', 'Unknown', 4]

        df_neg['MSI'] = df_neg['MSI'].fillna('')
        df_neg[5:].loc[df_neg['MSI'] == '', ['Metabolite name', 'INCHIKEY',
                                             'Feature Type', 'Feature Index']] = ['Unknown', 'null', 'Unknown', 4]

        self.write_to_log("Marked leftover unknowns as INCHIKEY = null and MSI = 4")

        # Save copies of pos/neg alignments
        df_pos_final = df_pos.copy()
        df_neg_final = df_neg.copy()

        # Move filenames from negative alignment into positive alignment
        sample_cols = np.array(df_neg.iloc[1] == 'Sample')
        sample_index = 10 + np.count_nonzero(sample_cols)
        df_pos.iloc[2:3, 10:sample_index] = df_neg.columns.values.tolist()[10:sample_index]

        blank_cols = np.array(df_neg.iloc[1] == 'Blank')
        blank_index = sample_index + np.count_nonzero(blank_cols)
        df_pos.iloc[2:3, sample_index:blank_index] = df_neg.columns.values.tolist()[sample_index:blank_index]

        qc_cols = np.array(df_neg.iloc[1] == 'QC')
        qc_index = blank_index + np.count_nonzero(qc_cols)
        df_pos.iloc[2:3, blank_index:qc_index] = df_neg.columns.values.tolist()[blank_index:qc_index]

        # Change columns in negative alignment for correct combining
        df_neg.columns = df_pos.columns

        # Create Combined Annotated sheet
        df_neg = df_neg[5:]
        df_combined = df_pos.append(df_neg, ignore_index=True)

        self.write_to_log("Size of Combined Annotated sheet: " + str(len(df_combined)))
        self.write_to_log("Expected size: " + str(len(df_pos) + len(df_neg)))

        # Initialize 'Unique Annotated' DataFrame
        df_unique = df_combined[0:4]

        # Remove unknowns
        df_combined.drop(df_combined[df_combined['Feature Index'] == 4].index, inplace=True)

        self.write_to_log("Removed unknowns from Combined Annotated sheet")

        # Get features (inchikeys) that appear in both positive and negative mode
        df_combined_features = df_combined[df_combined['Feature Index'] == 3]
        duplicate_inchikeys = df_combined_features[df_combined_features.duplicated(subset=['INCHIKEY'], keep=False)]

        # Get list of unique inchikeys
        inchikeys = duplicate_inchikeys[~duplicate_inchikeys.duplicated(subset=['INCHIKEY'])]
        inchikeys['INCHIKEY'] = inchikeys['INCHIKEY'].fillna('null')
        inchikeys = inchikeys['INCHIKEY'].tolist()

        self.write_to_log(
            "Capturing " + str(len(inchikeys)) + " features found in both positive and negative modes for merging")

        # For each inchikey,
        for index, inchikey in enumerate(inchikeys):

            # Put the duplicate pos/neg features into a single DataFrame
            loopdf = duplicate_inchikeys[duplicate_inchikeys['INCHIKEY'] == inchikey]

            if inchikey != 'null':

                # Remove the duplicate features from the Combined Annotated sheet for now
                df_combined.drop(df_combined[(df_combined['INCHIKEY'] == inchikey)
                                             & (df_combined['Feature Index'] == 3)].index, inplace=True)

                # Now separate pos/neg features into separate DataFrames
                pos_feature = loopdf.loc[loopdf['Polarity'] == 'pos']
                neg_feature = loopdf.loc[loopdf['Polarity'] == 'neg']

                if not neg_feature.empty and not pos_feature.empty:

                    pos_rt = pos_feature['Average Rt(min)'].astype(float).values
                    neg_rt = neg_feature['Average Rt(min)'].astype(float).values
                    rt_difference = abs(pos_rt - neg_rt)

                    # Calculate Pearson's correlation coefficient
                    sample_cols = np.array(df_combined.iloc[1] == 'Sample')
                    pos_feature_sample_intensity = pos_feature.iloc[:, sample_cols].astype(float)
                    neg_feature_sample_intensity = neg_feature.iloc[:, sample_cols].astype(float)

                    sample_cols = np.array(df_combined.iloc[1] == 'Sample')
                    pos_feature_sample_average = pos_feature.iloc[:, sample_cols].fillna(0).astype(
                        float).values.mean(axis=1)
                    neg_feature_sample_average = neg_feature.iloc[:, sample_cols].fillna(0).astype(
                        float).values.mean(axis=1)

                    x_diff = pos_feature_sample_intensity.sub(pos_feature_sample_average.squeeze())
                    y_diff = neg_feature_sample_intensity.sub(neg_feature_sample_average.squeeze())

                    x_diff_squared = x_diff ** 2
                    y_diff_squared = y_diff ** 2

                    numerator = x_diff.astype(float).values * y_diff.astype(float).values
                    numerator = numerator.sum()
                    denominator = (x_diff_squared.astype(float).values.sum() * y_diff_squared.astype(
                        float).values.sum()) ** 0.5

                    # Set Pearson correlation coefficient values
                    loopdf['Pearson'] = numerator / denominator
                    pos_feature['Pearson'] = numerator / denominator
                    neg_feature['Pearson'] = numerator / denominator

                    # Comma-separate values together
                    merged_feature = pos_feature.copy()
                    merged_feature['MSI'] = pos_feature['MSI'].astype(str).values + ', ' + neg_feature['MSI'].astype(str).values
                    merged_feature['Alignment ID'] = pos_feature['Alignment ID'].astype(str).values + ', ' + \
                                                     neg_feature['Alignment ID'].astype(str).values
                    merged_feature['Average Rt(min)'] = pos_feature['Average Rt(min)'].astype(str).values + ', ' + \
                                                        neg_feature['Average Rt(min)'].astype(str).values
                    merged_feature['Average Mz'] = pos_feature['Average Mz'].astype(str).values + ', ' + \
                                                   neg_feature['Average Mz'].astype(str).values
                    merged_feature['Polarity'] = 'both'

                    # Average sample intensities
                    sample_cols = np.array(df_combined.iloc[1] == 'Sample')
                    pos_feature_sample_intensity = pos_feature.iloc[:, sample_cols].astype(float).values
                    neg_feature_sample_intensity = neg_feature.iloc[:, sample_cols].astype(float).values
                    merged_feature.iloc[:, sample_cols] = (
                                                                      pos_feature_sample_intensity + neg_feature_sample_intensity) / 2
                    merged_feature.iloc[:, sample_cols] = merged_feature.iloc[:, sample_cols].round(0)

                    # Average blank intensities
                    blank_cols = np.array(df_combined.iloc[1] == 'Blank')
                    pos_feature_blank = pos_feature.iloc[:, blank_cols].astype(float).values
                    neg_feature_blank = neg_feature.iloc[:, blank_cols].astype(float).values
                    merged_feature.iloc[:, blank_cols] = (pos_feature_blank + neg_feature_blank) / 2
                    merged_feature.iloc[:, blank_cols] = merged_feature.iloc[:, blank_cols].round(0)

                    # Average pool intensities
                    qc_cols = np.array(df_combined.iloc[1] == 'QC')
                    pos_feature_qc = pos_feature.iloc[:, qc_cols].astype(float).values
                    neg_feature_qc = neg_feature.iloc[:, qc_cols].astype(float).values
                    merged_feature.iloc[:, qc_cols] = (pos_feature_qc + neg_feature_qc) / 2
                    merged_feature.iloc[:, qc_cols] = merged_feature.iloc[:, qc_cols].round(0)

                    # Add pos/neg features with Pearson correlation back to Combined Annotated sheet
                    df_combined = df_combined.append(loopdf, ignore_index=True)

                    # If RT match, add merged feature to Unique Annotated sheet... if not, add pos/neg features separately
                    if rt_difference <= 0.1:

                        df_unique = df_unique.append(merged_feature, ignore_index=True)

                    else:

                        loopdf['Pearson'] = 'No RT match'
                        df_unique = df_unique.append(loopdf, ignore_index=True)
                        self.write_to_log(
                            "No RT match for " + pos_feature['Metabolite name'].iloc[0] + ", leaving as-is")

                else:

                    # If there's some sort of mistake, add both features for human review
                    loopdf['Pearson'] = 'Error 1'
                    df_unique = df_unique.append(loopdf, ignore_index=True)


            else:

                # If there's some sort of mistake, add both features for human review
                loopdf['Pearson'] = 'Error 2'
                df_unique = df_unique.append(loopdf, ignore_index=True)

            # After duplicate inchikey has been rectified, drop it from DataFrame
            duplicate_inchikeys.drop(duplicate_inchikeys[
                                         duplicate_inchikeys['INCHIKEY'] == inchikey].index, inplace=True)

        # Append rest of duplicates DataFrame, if merging algorithm misses any duplicates for any reason
        df_unique = df_unique.append(duplicate_inchikeys, ignore_index=True)

        self.write_to_log("Did not merge " + str(len(df_unique) - len(
            inchikeys)) + " features found in both positive/negative mode due to significant RT difference")

        # Add remaining features to Unique Annotated
        internal_standards = df_combined[df_combined['Feature Index'] == 2]
        non_duplicate_inchikeys = df_combined_features.drop_duplicates(subset=['INCHIKEY'], keep=False)
        df_unique = df_unique.append([internal_standards, non_duplicate_inchikeys], ignore_index=True)

        self.write_to_log("Size of Unique Annotated sheet: " + str(len(df_unique)))

        # Sort by INCHIKEY
        df_unique.sort_values(by=['Feature Index', 'INCHIKEY'], inplace=True)
        df_unique.reset_index(drop=True, inplace=True)
        df_combined.sort_values(by=['Feature Index', 'INCHIKEY'], inplace=True)
        df_combined.reset_index(drop=True, inplace=True)

        # Clean up sheets
        df_unique.drop(columns=['Feature Index', 'Feature Type'], inplace=True)
        df_combined.drop(columns=['Feature Index', 'Feature Type'], inplace=True)
        df_pos_final.drop(columns=['Feature Index', 'Feature Type', 'Pearson'], inplace=True)
        df_neg_final.drop(columns=['Feature Index', 'Feature Type', 'Pearson'], inplace=True)

        # Header formatting
        df_pos_final = self.insert_row(4, df_pos_final, list(df_pos_final.columns))
        df_neg_final = self.insert_row(4, df_neg_final, list(df_neg_final.columns))
        df_combined = self.insert_row(4, df_combined, list(df_combined.columns))
        df_unique = self.insert_row(4, df_unique, list(df_unique.columns))

        self.write_to_log("Sorting by INCHIKEY and cleaning up sheets...")

        # Final save of formatted spreadsheet
        filename = 'Project_ID_Combined.tsv'

        self.write_to_log("Writing Unique Annotated sheet to Excel file")
        df_unique.to_csv(filename, sep='\t', header=False, index=False)
        self.write_to_log("Combining complete")

        step_two_log = self.log.copy()

        return {
            "Unique Annotated": df_unique,
            "Combined Annotated": df_combined,
            "Pos All Features": df_pos_final,
            "Neg All Features": df_neg_final,
            "Log": step_two_log
        }


    def write_to_log(self, text):

        """
        Logging function
        """

        current_time = time.strftime("%H:%M:%S")
        self.log.append(current_time + ' – ' + text)


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

    #step_one_files = pycutter.process_alignment(
    #    "../../../data/pipeline_test/step_0_raw_from_ms_dial/STQU002_pos_raw_alignment.txt")
    step_one_files=pycutter.process_alignment(
        '../../../data/BRYU005_pipeline_test/step_0_raw_from_ms_dial/BRYU005_pos_alignment_raw.txt',
        '../../../data/BRYU005_pipeline_test/step_0_raw_from_ms_dial/BRYU005_seq_MetaData.csv'
    )
    step_one_files['Reduced'].to_csv(
        f'../../../data/BRYU005_pipeline_test/step_1_post_pycutter/py_cutter_step_1_output.tsv',
        sep='\t',
        index=None,
        header=None
    )


    # Step 1 – Processing raw MS-DIAL alignment
    # step_one_files = pycutter.process_alignment("alignment.txt", "metadata.csv")
    # df_raw = step_one_files["Raw"]
    # df_reduced = step_one_files["Reduced"]
    # processing_log = step_one_files["Log"]

    # Step 2 – Combining curated positive and negative mode alignments
    # step_two_files = pycutter.combine_annotations("pos_alignment.xlsx", "neg_alignment.xlsx")
    # df_unique = step_two_files["Unique Annotated"]
    # df_combined = step_two_files["Combined Annotated"]
    # df_pos_features = step_two_files["Pos All Features"]
    # df_neg_features = step_two_files["Neg All Features"]
    # combining_log = step_two_files["Log"]
