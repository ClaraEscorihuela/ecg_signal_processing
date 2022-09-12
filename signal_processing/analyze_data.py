"""
Library for reading and analysing demographic information from phisionet dataset
"""
import os
import sys
curr_dir = os.getcwd()+'\signal_processing'
sys.path.append(curr_dir)

from IPython.display import display
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import difflib
import wfdb
import glob
import os

"""
@brief Upload the data and names from ecg signals in phisionet dataset
@param path_name: patient directory with the ecg signals
@param derivations: signal data
@param derivations_names: signal names
"""
def upload_data(path_name):
    record = wfdb.rdrecord(f'{path_name}')  #f'{path_name}\\s0015lre'
    derivations = record.p_signal
    derivations_names = record.sig_name
    return derivations, derivations_names

"""
@brief Select one signal by its derivation name
@param derivations: all signals
@param derivations_names: all names from the signals
@param ders: selected derivation
@return signal data from the selected derivation
"""
def select_derivation(derivations, derivations_names, der = 'ii'):
    ind = derivations_names.index(der)
    return derivations[:,ind]

"""
@brief Create dataframe with the demographic information from the physionet dataframe
@param path_name: directory with all patient folders
@param by_patient: boolean, true if the information from the last ecg of a patient 
                   with multiple ecgs is saved, false whether all the ecgs information
                   are saved
@param df_glob: dataframe with all demographic information
"""
def create_dataframe(path_name, by_patient = True):

    def listToDict(lst):
        op = {i.split(':')[0].replace(' ', ''): i.split(':')[1].replace(' ', '') for i in lst}
        return op

    # Create global dataframe
    df_glob = pd.DataFrame()

    # Add patient information per dataframe
    for pat in os.listdir(path_name):
        if 'patient' not in pat:
            continue

        files = glob.glob(f"{path_name}/{pat}/*.hea")
        df = pd.DataFrame()

        # Collect demographic info
        for file in files:
            record = wfdb.rdheader(file[:-4])
            dic = listToDict(record.comments)
            df_now = pd.DataFrame.from_records([dic])

            # Do not select patients without ECG date
            if df_now['ECGdate'].values == 'n/a' and by_patient:
                continue

            # Add subject ID to the dataframe
            df['SubjID'] = pat[-3:]
            df = pd.concat([df, df_now], axis=0)

        # Make sure there is information
        if len(df) == 0:
            continue

        # Select info from last ECG recording
        if by_patient:
            df['ECGdate'] = pd.to_datetime(df['ECGdate'], format="%d/%m/%Y")
            df = df.sort_values(by='ECGdate')

        df_glob = pd.concat([df_glob, df.iloc[-1:]])

        #print(f'patient {pat[-3:]} done')
    #print('yay, we have it!')

    return df_glob

"""
@brief Fix dataframe by cleaning data, including transforming strings to numeric data, 
                                                 solving nan 
                                                 equalizing categories
                                                 ...
@param df: dataframe to fix
@return df: df fix
"""
def fix_dataframe(df):

    def df1(str):
        return int(str.split('/')[0])
    def df2(str):
        return int(str.split('/')[1][:-4])
    def fuzzy_replace(x, names):
        aliases = difflib.get_close_matches(x, names, cutoff=0.9)
        closest = pd.Series(aliases).mode()
        closest = aliases[0] if closest.empty else closest[0]
        return closest

    # Fill nan values
    df = df.fillna('n/a')
    df['Smoker'] = df['Smoker'].replace('unknown','n/a')

    # Transform numerical columns
    df['age'] = df['age'].apply(pd.to_numeric, errors='coerce')

    # Make sure string columns present the same name
    # Add all the name columns you would like to check
    for col in ['Reasonforadmission']:
        df[col] = df[col].str.lower()
        # Apply fuzzy matching strategy to make a column of equal categories
        df[col].apply(lambda x: fuzzy_replace(x, df[col].unique()))

    # Create systole and diastole pressure
    df['systole'] = df['PeripheralbloodPressure(syst/diast)'].apply(lambda row: df1(row) if row != 'n/a' else float("nan"))
    df['diastole'] = df['PeripheralbloodPressure(syst/diast)'].apply(lambda row: df2(row) if row != 'n/a' else float("nan"))

    return df

"""
@brief Plot histogram with density funciton for contnious data
@param df: dataframe to plot
@return val: list with x,y column names 
"""
def histogram_plot(df,val):
    fig, ax = plt.subplots(figsize=(5,4))
    sns.histplot(df[val].dropna(),kde=True,line_kws={"linewidth":1.5}, ax = ax)
    sns.despine(bottom=False, left=False)
    plt.xlabel(val, fontsize=10)
    plt.xticks(rotation=45)
    plt.ylabel('Count', fontsize=10)
    plt.title(f'{val} frequency',loc='left', fontsize = 12)
    plt.tight_layout()
    plt.show()

"""
@brief Draw bar plot for categorical data
@param df: dataframe to plot
@return val: list with x,y column names 
@return title: title for the plot
"""
def bar_plot(df, values, title):

    fig, ax = plt.subplots(figsize=(5,4))

    if len(values) == 2:
        sns.countplot(x=values[0],
                           hue=values[1],
                           data=df,
                           palette= "Blues",
                           dodge = True, ax = ax)
        title = f'{title[0]} & {title[1]}'
        leg = plt.legend(loc='best', fontsize=9)
        leg.get_frame().set_linewidth(0.0)
    else:
        sns.countplot(x=values[0], data=df, palette= "Blues", ax= ax)
        title = f'{title[0]}'
    plt.xticks(rotation=45, ha='right')

    sns.despine(bottom=False, left=False)
    plt.tight_layout()
    plt.ylabel('Count', fontsize = 10)
    #plt.xlabel(title[], fontsize=10)
    plt.title(title, loc='left', fontsize=12)
    plt.tight_layout()
    plt.show()

"""
@brief Draw boxplot to relate continious and categorical data
@param df: dataframe to plot
@param val: list with x,y column names 
@param title: title for the plot
@param bins: for continious data, bins to classify it (usually age)
@param labels: for continious data, labels to locate the bins (usually age)
"""
def boxplot_plot(df,values, bins = None, labels = None):

    fig, ax = plt.subplots()
    if len(values) > 2:
        ax = sns.boxplot(x=df[values[0]], y=df[values[1]], hue = df[values[2]], palette="Blues")
        #plt.title(f'{title[1]} vs ({title[0]} & {title[2]})', loc='left', fontsize=12)
        leg = plt.legend()
        leg.get_frame().set_linewidth(0.0)

    else:
        ax = sns.boxplot(x=df[values[0]], y=df[values[1]],  palette="Blues")
        plt.title(f'{values[1]} vs {values[0]}', loc='left', fontsize=12)

    #ax = sns.stripplot(x=df[values[0]], y=df[values[1]], color="orange", jitter=0.2,
    #                   size=3, alpha=0.5)
    sns.despine(bottom=False, left=False)

    rot = 90
    plt.xlabel(values[0], fontsize = 10)
    plt.xticks(rotation=rot)
    plt.ylabel(values[1], fontsize=10)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":

    values = ['age','systole','diastole']
    values = [['age','Reasonforadmission']]
    title = 'Categorical data'

    wk_dir = os.getcwd()
    path_name = f'{wk_dir}\\physionet.org\\files\\ptbdb\\1.0.0'

    df = create_dataframe(path_name)
    df = fix_dataframe(df)
    for val in values:
        boxplot_plot(df, val)
    print('stop')
