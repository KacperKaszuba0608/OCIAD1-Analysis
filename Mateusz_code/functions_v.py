import pandas as pd
import numpy as np
import scipy.stats as st
import re


def get_sets(df, *args):
    '''
    Extract column names that contain specified words and the corresponding columns

    Args:
        df (DataFrame): The DataFrame
        *args (str): Words to identify column names
    
    Returns:
        samples_all (list of str): All column names containing the words
        *samples (list of str): Column names containing each of the words, ordered the same as the input words
        *columns (list of DataFrame): Columns corresponding to the names, ordered the same as the input words
    '''
    samples = []
    columns = []

    columns_all = list(df.columns.values)
    samples_all = [columns_all[i] for i in range(0,len(columns_all))   if [re.search('|'.join(args), j) for j in columns_all][i]]

    for req in args:
        samples.append([columns_all[i] for i in range(0,len(columns_all))   if [re.search(req, j) for j in columns_all][i]])
        columns.append(df[samples[-1]])

    return (samples_all,) + (*samples,) + (*columns,)


def fatty_acids(df):
    '''
    Get information about the fatty acids of the lipids based on their identification and add it as columns of the DataFrame

    Args:
        df (DataFrame): The DataFrame

    Returns:
        df['Chain Lengths'] (list of int): Length of each identified chain
        df['Length Even'] (bool): Information if the chain length is even (True) or odd (False)
        df['Double Bonds'] (list of int): Number of double bonds of each identified chain
        df['Double Bonds Total'] (int): Total number of double bonds in the chains of the lipid
        df['Saturation'] (int): Type of saturation of the lipid. 
            0- Saturated
            1- Monounsaturated
            2- Polyunsaturated

    '''
    chain_lengths = []
    double_bonds = []

    chains = [re.findall(r'\d\d:\d?\d', id) for id in df['Identification']]

    ether_info = [re.findall(r'[PO]-', id) for id in df['Identification']]
    ether = [False if not x else True   for x in ether_info]

    for chain in chains:
        chains_bonds = [x.split(':') for x in chain]
        chain_lengths.append([int(item[0]) for item in chains_bonds])
        double_bonds.append([int(item[1]) for item in chains_bonds])
    
    is_even = [all(i % 2 == 0 for i in j) for j in chain_lengths]
    double_bonds_sums = [sum(i) for i in double_bonds]    
    saturation = [0 if db == 0 else 1 if db == 1 else 2   for db in double_bonds_sums]
    
    df.insert(2, 'Ether', ether, True)
    df.insert(2, 'Saturation', saturation, True)
    df.insert(2, 'Double Bonds Total', double_bonds_sums, True)
    df.insert(2, 'Double Bonds', double_bonds, True)
    df.insert(2, 'Length Even', is_even, True)
    df.insert(2, 'Chain Lengths', chain_lengths, True)

    return df


def normalise(df):
    '''
    Normalise the columns of the DataFrame using sum normalisation

    Args:
        df (DataFrame): The DataFrame

    Returns:
        df_norm (DataFrame): The normalised DataFrame
    '''
    df = 2 ** df
    df_norm = pd.DataFrame(columns = df.columns)
    sums = df.sum()

    for i in range(0,len(df.columns)):
        df_norm.iloc[:,i] = df.iloc[:,i] / sums.iloc[i]
    
    df_norm = np.log2(df_norm.astype('float64'))

    return df_norm


def fc(df):
    '''
    Calculate the fold change between WT and KO

    Args:
        df (DataFrame): The DataFrame

    Returns:
        df['FC'] (float): The fold change between WT and KO
        df['FC up/down'] (int): The direction of change
            1- change up
            0- change down
    '''

    df = df.copy()
    df['FC'] = df.iloc[:,0:3].mean(axis = 1) - df.iloc[:,3:6].mean(axis = 1)
    df['FC up/down'] = [1 if k > 0 else 0   for k in df['FC']]

    return df


def ttest(df):
    '''
    Perform a t-test on the WT and KO samples

    Args:
        df (DataFrame): The DataFrame

    Returns:
        df['pval'] (float): p-value of the t-test
        df['Significant'] (int): Statistical significance of the difference between WT and KO
            1- Statistically significant
            0- Statistically insignificant
    '''
    df = 2 ** df
    pval = [st.ttest_ind(df.iloc[i, 0:3].tolist(), df.iloc[i, 3:6].tolist(), equal_var = False).pvalue for i in range(0,len(df.index))]
    significant = [1 if j < 0.05 else 0   for j in pval]
    df = np.log2(df)

    df['pval'] = pval
    df['Significant'] = significant

    return df


def filter_classes(df, n):
    '''
    Extract lipid classes from the DataFrame that have at least n data points

    Args:
        df (DataFrame): The Data Frame
        n (int): Threshold of the number of data points to filter the classes by

    Returns:
        classes_filt (list of str): Names of the lipid classes satisfying the condition
        df_filt (DataFrame): Filtered DataFrame containing only the lipid classes satisfying the condition
    '''
    classes_filt = df['Lipid Class'].value_counts()[lambda cl : cl > n].index.to_list()
    df_filt = df[df['Lipid Class'].isin(classes_filt)].reset_index(drop = True)

    return classes_filt, df_filt


def compare_mean_median(df1, df2):
    '''
    Calculate and compare means and medians of two Series

    Args:
        df1 (Series): First Series
        df2 (Series): Second Series
    
    Returns:
        mean1 (float): mean of the first Series
        mean2 (florat): mean of the second Series
        median1 (float): median of the first Series
        median2 (florat): median of the second Series
        mean_diff (float): difference of both means
        median_diff (float): difference of both medians
    '''
    df1 = 2 ** df1
    df2 = 2 ** df2

    mean1 = df1.mean()
    mean2 = df2.mean()
    median1 = df1.median()
    median2 = df2.median()
    mean_diff = mean1 - mean2
    median_diff = median1 - median2

    return mean1, mean2, median1, median2, mean_diff, median_diff


def cv(df):
    '''
    Calculate the coefficient of variation (CV) for each row of the DataFrame and the median of those CVs

    Args:
        df (DataFrame): The DataFrame
    
    Returns:
        df['CV'] (float): CV of the row
        cv_median (float): median of the calculated CVs
    '''
    df = df.copy()
    df = 2 ** df
    cv = (df.std(axis=1) / df.mean(axis=1)) * 100
    df = np.log2(df)
    df['CV'] = cv
    cv_median = df['CV'].median()

    return df, cv_median


def group_wt_ko(df, samples_ko, samples_wt):
    '''
    Group the DataFrame into two columns containg all of the intensities for WT and KO

    Args:
        df (DataFrame): The DataFrame
        samples_ko (list of str): Names of the columns containing KO samples
        samples_wt (list of str): Names of the columns containing WT samples
    
    Returns:
        df_grouped (DataFrame): Two-column DataFrame containing all data for WT and KO in each column
    '''
    df_ko = pd.concat([df[x] for x in samples_ko], ignore_index = True)
    df_wt = pd.concat([df[x] for x in samples_wt], ignore_index = True)
    df_grouped = pd.concat([df_ko, df_wt], axis = 1)
    df_grouped.columns = ['KO', 'WT']

    return df_grouped


def significance_score(df):
    '''
    Calculate the significance scores and their directionality for the lipid head groups in the DataFrame

    Significance score: 
    Percentage of the statistically significant fold changes

    Significance score directionality:
    Percentage of positive significant fold changes, 
    percentage of negative significant fold changes

    Args:
        df (DataFrame): The DataFrame
    
    Returns:
        df['Significance Score'] (float): Significance score of the head group of the lipid in the row
        lipid_scores (dict of str: list of float): Significance scores and their directionalities for each head group present in the DataFrame
    '''
    lipid_sig_scores = []
    lipid_scores = {}
    for lipid in df['Lipid Class'].unique():
        df_lipid = df[df['Lipid Class'] == lipid]
        df_lipid_significant = df_lipid[df_lipid['Significant'] == 1]

        all = len(df_lipid.index)
        significants = len(df_lipid_significant.index)
        ups = df_lipid_significant['FC up/down'].sum()

        sig_score = (significants / all) * 100
        lipid_sig_scores.extend([sig_score] * all)
        lipid_scores[lipid] = [sig_score,
                               (ups / all) * 100,
                               ((significants - ups) / all) * 100
                               ]
    
    df['Significance Score'] = lipid_sig_scores

    return lipid_scores


def even_odd_score(df):
    '''
    Calcualate the percentage of the even- and odd-chain lipids with positive and negative fold changes

    Args:
        df (DataFrame): The DataFrame
    
    Returns:
        even_scrore_up (float): Percentage of even-chain lipids with positive fold changes
        even_scrore_down (float): Percentage of even-chain lipids with negative fold changes
        odd_scrore_up (float): Percentage of odd-chain lipids with positive fold changes
        odd_scrore_down (float): Percentage of odd-chain lipids with negative fold changes
    '''
    df_even = df[df['Length Even'] == True]
    df_even_up = df_even[df_even['FC up/down'] == 1]
    df_even_down = df_even[df_even['FC up/down'] == 0]
    even_score_up = len(df_even_up.index) / len(df_even.index) * 100
    even_score_down = len(df_even_down.index) / len(df_even.index) * 100

    df_odd = df[df['Length Even'] == False]
    df_odd_up = df_odd[df_odd['FC up/down'] == 1]
    df_odd_down = df_odd[df_odd['FC up/down'] == 0]
    odd_score_up = len(df_odd_up.index) / len(df_odd.index) * 100
    odd_score_down = len(df_odd_down.index) / len(df_odd.index) * 100

    return even_score_up, even_score_down, odd_score_up, odd_score_down


def compare_even_odd(df, classes):
    '''
    Perform t-tests to compare the fold changes of the even-chain and odd-chain lipids in the DataFrame for each specified head group

    Args:
        df (DataFrame): The DataFrame
        classes (list of str): specified lipid head groups
    
    Returns:
        df_comp (DataFrame): DataFrame containing the p-values and the statistical significance of the differences between fold changes of even-chain and odd-chain lipids for the specified head groups
            1- Statistically significant
            0- Statistically insignificant
    '''
    df_comp = pd.DataFrame(index = classes)
    for cls in classes:
        df_class = df[df['Lipid Class'] == cls]
        df_even = df_class[df_class['Length Even'] == True]
        df_odd = df_class[df_class['Length Even'] == False]
        df_comp.loc[cls, 'pval'] = st.ttest_ind(df_even['FC'].to_list(), df_odd['FC'].to_list(), equal_var = False).pvalue
        df_comp.loc[cls, 'Significant'] = 1 if df_comp.loc[cls, 'pval'] < 0.05 else 0

    return df_comp


def get_even_odd(df, even, samp):
    '''
    Save the DataFrame containing the even- or odd-chain lipids as a .csv file

    Args:
        df (DataFrame): The DataFrame
        even (bool): Whether to save even- or odd-chain lipids
        samp (string): Name of the sample to set as the .csv file name
    '''
    df_even_odd = df[df['Length Even'] == even].reset_index(drop = True)
    df_even_odd = df_even_odd[df_even_odd['Significant'] == 1].reset_index(drop = True)
    df_even_odd.to_csv(f'even-chains-{samp}.csv', decimal = ',')


def common_classes(*args):
    '''
    Get a list of common lipid head-groups present in all DataFrames

    Args:
        *args (DataFrame): The DataFrames
    
    Returns:
        common_classes (list of str): Lipid head-groups present in all DataFrames
    '''
    common_classes = set(args[0]['Lipid Class'].value_counts().index)
    for df in args[1:]:
        common_classes = common_classes.intersection(set(df['Lipid Class'].value_counts().index))
    common_classes = list(common_classes)
    
    return common_classes


def median_fc(df, classes, samp):
    '''
    Calculate the median fold change for each specified lipid head-group in the DataFrame

    Args:
        df (DataFrame): The DataFrame
        classes (list of str): specified lipid head-groups
        samp (str): name of the sample, for naming purpouses

    Returns:
        df_m (DataFrame): DataFrame contiaing median fold changes for the specified lipid head-groups
    '''
    df_m = pd.DataFrame(columns = classes, index = [f'{samp}'])
    for cls in classes:
        df_class = df[df['Lipid Class'] == cls]
        df_m[cls] = df_class['FC'].median()

    return df_m
