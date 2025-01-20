from aflib import decipher
import sqlite3
import os
import pickle
import pandas as pd
import json
import itertools
from tqdm.notebook import tqdm
from multiprocessing import Pool


def identity_diff(base_path, previous_data, current_data):
    changed_items = {'added': [], 'removed': [], 'modified_files': {}}

    # Compare current data with previous data
    current_paths = set(current_data.keys())
    previous_paths = set(previous_data.keys())

    # Identify added and removed paths
    added_paths = current_paths - previous_paths
    removed_paths = previous_paths - current_paths

    for path_local in added_paths:
        changed_items['added'].append(path_local)
    for path_local in removed_paths:
        changed_items['removed'].append(path_local)

    # Check for modifications in existing paths
    for path_local in current_paths.intersection(previous_paths):
        current_info = current_data[path_local]
        previous_info = previous_data[path_local]
        if current_info != previous_info:
            changed_items[path_local] = current_info
            modified_files = []
            for file in current_info['files']:
                if file not in previous_info['files'] or \
                        os.path.getmtime(os.path.join(base_path, path_local, file)) != previous_info['last_modified']:
                    modified_files.append(file)
            if modified_files:
                changed_items['modified_files'][path_local] = modified_files

    return changed_items


def analyze_db(dir_path):
    # Check if previous data exists
    previous_data = {}
    if os.path.exists(os.path.join(dir_path, 'dir_data.pickle')):
        with open(os.path.join(dir_path, 'dir_data.pickle'), 'rb') as f:
            previous_data = pickle.load(f)

    # Get current directory structure
    current_data = {}
    for root, dirs, files in os.walk(dir_path):
        relative_path = os.path.relpath(root, dir_path)
        current_data[relative_path] = {
            'dirs': dirs,
            'files': files,
            'last_modified': os.path.getmtime(root)
        }

    # Identify differences
    changes = identity_diff(dir_path, previous_data, current_data)

    # Clean MSA folder from changes for the remvoed folder/
    changes['removed'] = [x for x in changes['removed'] if 'msas' not in x]

    # Clean MSA folder from changes for the added folder
    changes['added'] = [x for x in changes['added'] if 'msas' not in x]
    # Save current data for future comparison
    with open(os.path.join(dir_path, 'dir_data.pickle'), 'wb') as f:
        pickle.dump(current_data, f)
    return changes


def get_added_df(changes, library_path, select_intelligent=True):
    # Generate a list f paths to reload
    reloadable_paths_raw = changes['added'] + list(changes['modified_files'].keys())
    reloadable_paths = [os.path.join(library_path, path) for path in reloadable_paths_raw if '-' in path]

    if not reloadable_paths:
        return pd.DataFrame()
    # Load the library
    res = decipher.load_library_aslist(reloadable_paths)
    unique_df = decipher.uniquify(pd.DataFrame(res))
    # Remove rows with missing values
    unique_df.dropna(subset=['Mean pDockQ', 'pTM', 'ipTM'], inplace=True, ignore_index=True)
    # Remove the path prefix
    unique_df['Relative Directory'] = unique_df['Directory'].apply(lambda x: x.replace(library_path + os.path.sep, ''))

    return unique_df


def update_db(library_path):
    """
    Update the SQLite database by scanning the library for new folders and folders that have been modified.
    :param library_path: path to library
    :return:
    """
    all_potential_folders, reload_folders = scan_library(library_path)
    load_library(library_path, reload_folders)


# # This is the function can be activated if a pickle database is used (didn't work that well for me)
# def update_db(library_path, pickle_db=True):
#     # Identify the changes in the database
#     changes = analyze_db(library_path)
#     df_toadd = get_added_df(changes, library_path)
#
#     # Load current database
#     if not os.path.exists(os.path.join(library_path, 'PPIDB_full.db')):
#         print('Database not found – creating a new one')
#         cnx = sqlite3.connect(os.path.join(library_path, 'PPIDB_full.db'))
#
#     try:
#         if pickle_db:
#             df_in = pd.read_pickle(os.path.join(library_path, 'PPIDB_full.pkl'))
#         else:
#             df_in = pd.read_sql('AF2 full-lib', f"sqlite:///{os.path.join(library_path, 'PPIDB_full.db')}")
#
#         # Remove the entries that were removed from the database
#         for rem_path in changes['removed']:
#             components = os.path.basename(rem_path).split('-')
#             # Show row where 'Job name' contains both components
#             index = df_in[
#                 df_in['Job name'].str.contains(components[0]) & df_in['Job name'].str.contains(components[1])].index
#             df_in.drop(index, inplace=True)
#     except Exception as _:
#         print('Database not found – considering it do be empty')
#         df_in = pd.DataFrame()
#
#     # Add the new entries
#     df = pd.concat([df_in, df_toadd], ignore_index=True)
#     # print(df['Job name'].head())
#
#     # Save and replace current database with updated one
#     if pickle_db:
#         df.to_pickle(os.path.join(library_path, 'PPIDB_full.pkl'))
#     else:
#         cnx = sqlite3.connect(os.path.join(library_path, 'PPIDB_full.db'))
#         df.to_sql(name='AF2 full-lib', con=cnx, if_exists='replace', index=False)


def get_prediction_features(path):
    with open(path, 'rb') as f:
        return json.load(f)


def get_data(library_path, verbose=False):
    """
    Get the data from the library by iterating over all the json files in the library
    :param library_path: path to the library
    :param verbose: Should the code be chatty
    :return: dataframe containing unique combinations of proteins from the database
    """

    # *** The following section is only meant to be used if the database is not working for some reason ***
    all_json_files = []
    for root, dirs, files in os.walk(library_path):
        for file in files:
            if file.lower() == 'docking_results.json':
                all_json_files.append(os.path.join(root, file))

    all_data = list(map(get_prediction_features, all_json_files))
    # Get unique protiens from the database
    unique_genes = set(x['Chain A'] for x in all_data) | set(x['Chain B'] for x in all_data)
    # Get all possible combinations of the unique proteins
    combinations = list(itertools.combinations_with_replacement(unique_genes, 2))
    # Make DataFrame of docking data
    df = pd.DataFrame(all_data)

    # Create a list where the unique combination of the proteins are stored
    list_of_unique_combs = []

    # For each combination, check if the combination is in the DataFrame
    for comb_a, comb_b in combinations:
        temp_dict = {
            'Chain A': comb_a,
            'Chain B': comb_b,
            'pTM': None,
            'ipTM': None,
            'pLDDT': None,
            'pDockQ': None
        }
        if comb_a.lower() == 'nanos3' and comb_b.lower() == 'dnd1' or comb_a.lower() == 'dnd1' and comb_b.lower() == 'nanos3':
            temp_dict['Condition'] = 'Both'
            temp_dict['Target'] = 'None'
        elif comb_a.lower() == 'dnd1' or comb_b.lower() == 'dnd1':
            temp_dict['Condition'] = 'dnd1'
            if comb_a.lower() == 'dnd1':
                temp_dict['Target'] = comb_b
            else:
                temp_dict['Target'] = comb_a
        elif comb_a.lower() == 'nanos3' or comb_b.lower() == 'nanos3':
            temp_dict['Condition'] = 'nanos3'
            if comb_a.lower() == 'nanos3':
                temp_dict['Target'] = comb_b
            else:
                temp_dict['Target'] = comb_a
        subframe = df.query(
            '(`Chain A` == @comb_a) & (`Chain B` == @comb_b) | (`Chain A` == @comb_b) & (`Chain B` == @comb_a)')
        if len(subframe) > 1:
            max_iptm = subframe['Max ipTM'].max()
            temp_dict['ipTM'] = max_iptm
            temp_dict['pTM'] = subframe.query('`Max ipTM` == @max_iptm')['Max pTM'].values[0]
            temp_dict['pLDDT'] = subframe.query('`Max ipTM` == @max_iptm')['Average pLDDT'].values[0]
            temp_dict['pDockQ'] = subframe.query('`Max ipTM` == @max_iptm')['Max pDockQ'].values[0]
            list_of_unique_combs.append(temp_dict)
            if verbose:
                print('Found multiple entries for', comb_a, comb_b)
        elif len(subframe) == 1:
            temp_dict['ipTM'] = subframe['Max ipTM'].values[0]
            temp_dict['pTM'] = subframe['Max pTM'].values[0]
            temp_dict['pLDDT'] = subframe['Average pLDDT'].values[0]
            temp_dict['pDockQ'] = subframe['Max pDockQ'].values[0]
            list_of_unique_combs.append(temp_dict)

    # Create a DataFrame from the list
    return pd.DataFrame(list_of_unique_combs)


def scan_library(lib_path, reload_all=False):
    """
    Scans the library for new folders and folders that have been modified.
    :param lib_path:
    :param reload_all:
    :return: List [All potential folders, Folders to reload]
    """
    updated_dir_data = {}

    try:
        cnx = sqlite3.connect(os.path.join(lib_path, 'af2db.db'))
        cur = cnx.cursor()
        cur.execute('SELECT Directory FROM AF2DB_DETAILED')
        outdated_dir_data = {x[0]: os.path.getmtime(x[0]) for x in cur.fetchall()}
        # with open(os.path.join(lib_path, 'dir_data.json'), 'r') as f:
        #     outdated_dir_data = json.load(f)
    except FileNotFoundError:
        outdated_dir_data = {}

    if reload_all:
        outdated_dir_data = {}

    all_potential_folders = []
    reload_folders = []
    for root, dirs, files in os.walk(lib_path):
        for pred_dir in dirs:
            full_pred_dir = os.path.join(root, pred_dir)
            if pred_dir.find('-') != -1 and not pred_dir.startswith('.'):
                updated_dir_data[full_pred_dir] = os.path.getmtime(full_pred_dir)
                all_potential_folders.append(full_pred_dir)
                if outdated_dir_data.get(full_pred_dir) is None:
                    reload_folders.append(full_pred_dir)
                elif outdated_dir_data.get(full_pred_dir) != os.path.getmtime(full_pred_dir):
                    reload_folders.append(full_pred_dir)
                    print('Folder modified:', pred_dir)

    with open(os.path.join(lib_path, 'dir_data.json'), 'w') as f:
        json.dump(updated_dir_data, f)

    return all_potential_folders, reload_folders


def process_combination(folder):
    # print(f'Processing {folder}')
    calc_new = False if os.path.basename(folder).find('%') != -1 else True

    if calc_new:
        comb_pred = decipher.CombinationPrediction(folder)
        return comb_pred.main()
    else:
        with open(os.path.join(folder, 'docking_results.json'), 'r') as f:
            return json.load(f)


def load_library(lib_path, reload_folders):
    """
    Loads the library into the database.
    :param lib_path:
    :param reload_folders:
    :return:
    """
    db_path = os.path.join(lib_path, 'af2db.db')

    if not os.path.exists(db_path):

        # Initialize the database
        cnx = sqlite3.connect(db_path)

        # Create a cursor
        cur = cnx.cursor()

        # Create the detailed table
        cur.execute("""CREATE TABLE AF2DB_DETAILED (
                    ID INTEGER UNIQUE ,
                    JobName TEXT,
                    ChainA TEXT,
                    ChainB TEXT,
                    GeneA TEXT,
                    GeneB TEXT,
                    Directory TEXT,
                    MAXpDockQ REAL,
                    AVGpDockQ REAL,
                    WEIGHTEDpDockQ REAL,
                    MAXipTM REAL,
                    MAXpTM REAL,
                    ChainAseq TEXT,
                    ChainBseq TEXT,
                    ChainAseqLen INTEGER,
                    ChainBseqLen INTEGER,
                    AVGpLDDT REAL)""")

        cur.execute("""CREATE TABLE AF2DB_SUMMARY (
                    ChainA TEXT,
                    ChainB TEXT,
                    ChainAB TEXT UNIQUE,
                    GeneA TEXT,
                    GeneB TEXT,
                    pDockQ REAL,
                    ipTM REAL,
                    pTM REAL,
                    pLDDT REAL,
                    ChainAseq TEXT,
                    ChainBseq TEXT,
                    ChainAseqLen INTEGER,
                    ChainBseqLen INTEGER)""")
        cnx.commit()
    else:
        cnx = sqlite3.connect(db_path)
        cur = cnx.cursor()

    # Multiprocessing of the to be reloaded folders
    with Pool(processes=6) as p:
        processed_combinations = list(
            tqdm(p.imap_unordered(process_combination, reload_folders), total=len(reload_folders)))

    # Iterate though the processed combinations and load them into the database
    for docking_results in processed_combinations:
        # print(docking_results['Job name'])
        # Translate the docking results
        docking_translated = {
            'ID': docking_results['ID'],
            'Jobname': docking_results['Job name'],
            'ChainA': docking_results['Chain A'],
            'ChainB': docking_results['Chain B'],
            'GeneA': docking_results['Gene A'],
            'GeneB': docking_results['Gene B'],
            'Directory': docking_results['Directory'],
            'MAXpDockQ': docking_results['Max pDockQ'],
            'AVGpDock': docking_results['Avg pDockQ'],
            'WEIGHTEDpDock': docking_results['Weighted pDockQ'],
            'MAXipTM': docking_results['Max ipTM'],
            'MAXpTM': docking_results['Max pTM'],
            'ChainAseq': docking_results['Chain A sequence'],
            'ChainBseq': docking_results['Chain B sequence'],
            'ChainAseqLen': docking_results['Chain A length'],
            'ChainBseqLen': docking_results['Chain B length'],
            'AVGpLDDT': docking_results['Average pLDDT'],
            # 'DistanceMap': docking['Distance map'],
            # 'pLDDTmap': docking['pLDDT map']
        }

        try:
            with cnx:
                cur.execute(
                    "INSERT INTO AF2DB_DETAILED VALUES (:ID,:Jobname,:ChainA,:ChainB,:GeneA,:GeneB,:Directory,:MAXpDockQ,:AVGpDock,:WEIGHTEDpDock,:MAXipTM,:MAXpTM,:ChainAseq,:ChainBseq,:ChainAseqLen,:ChainBseqLen,:AVGpLDDT)",
                    docking_translated)
        except sqlite3.IntegrityError:
            print(
                f'Integrity error: {docking_translated["Jobname"]} with ID <{docking_translated["ID"]}> was already in the database.')
            continue


def summarize_detailed_db(lib_path):
    """
    Summarizes the detailed database into a summary database.
    :param lib_path:
    :return:
    """
    db_path = os.path.join(lib_path, 'af2db.db')

    if not os.path.exists(db_path):
        raise FileExistsError('Database not found!')

    cnx = sqlite3.connect(db_path)

    with cnx:
        # Connect to the database
        cur = cnx.cursor()

        # Fetch all the distinct chains
        cur.execute("SELECT DISTINCT ChainA, ChainB FROM AF2DB_DETAILED")
        data = cur.fetchall()
        unique_chains = set()
        for entry in data:
            for chain in entry:
                unique_chains.add(chain)

    # Get all possible combinations that can be made with the present chains
    combinations = list(itertools.combinations_with_replacement(unique_chains, 2))

    # Empty list to store all the   results
    results = []

    # For each combination, check if it is present in the database
    for comb in combinations:
        with cnx:
            df_comb = pd.read_sql_query(
                "SELECT * FROM AF2DB_DETAILED WHERE (ChainA = ? AND ChainB = ?) OR (ChainA = ? AND ChainB = ?)", cnx,
                params=(comb[0], comb[1], comb[1], comb[0]))
        if not len(df_comb) == 0:
            data = {
                "ChainA": comb[0],
                "ChainB": comb[1],
                "ChainAB": f'{comb[0]}-{comb[1]}',
                "GeneA": df_comb['GeneA'].iloc[0],
                "GeneB": df_comb['GeneB'].iloc[0],
                "pDockQ": float(df_comb['MAXpDockQ'].max()),
                "ipTM": float(df_comb['MAXipTM'].max()),
                "pTM": float(df_comb['MAXpTM'].max()),
                "pLDDT": float(df_comb['AVGpLDDT'].max()),
                "ChainAseq": df_comb['ChainAseq'].iloc[0],
                "ChainBseq": df_comb['ChainBseq'].iloc[0],
                "ChainAseqLen": int(df_comb['ChainAseqLen'].iloc[0]),
                "ChainBseqLen": int(df_comb['ChainBseqLen'].iloc[0]),

            }

            results.append(data)

    for result in results:
        with cnx:
            cur.execute(
                "INSERT INTO AF2DB_SUMMARY VALUES (:ChainA,:ChainB, :ChainAB, :GeneA,:GeneB,:pDockQ,:ipTM,:pTM,:pLDDT,:ChainAseq,:ChainBseq,:ChainAseqLen,:ChainBseqLen)",
                result)


def get_table_asdf(library_path, table_name):
    if not library_path.endswith('.db'):
        db_path = os.path.join(library_path, 'af2db.db')
    else:
        db_path = library_path
    if not os.path.exists(db_path):
        raise FileExistsError('Database not found!')

    # Connect to the database
    con = sqlite3.connect(db_path)

    # Check if the table exists
    try:
        con.execute(f"SELECT * FROM {table_name}")
    except sqlite3.OperationalError:
        print(f"Table {table_name} does not exist. Please choose a valid table from the options below:")
        sql_query = """SELECT name FROM sqlite_master WHERE type='table';"""
        for table in con.execute(sql_query).fetchall():
            print('-', table[0])

    return pd.read_sql_query(f'SELECT * FROM {table_name}', con)
