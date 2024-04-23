from aflib.decipher import load_library_aslist, uniquify
import sqlite3
import os
import pickle
import pandas as pd


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


def get_added_df(changes, library_path):
    # Generate a list f paths to reload
    reloadable_paths_raw = changes['added'] + list(changes['modified_files'].keys())
    reloadable_paths = [os.path.join(library_path, path) for path in reloadable_paths_raw if '-' in path]

    if not reloadable_paths:
        return pd.DataFrame()
    # Load the library
    res = load_library_aslist(reloadable_paths)
    unique_df = uniquify(pd.DataFrame(res))
    # Remove rows with missing values
    unique_df.dropna(subset=['Mean pDockQ', 'pTM', 'ipTM'], inplace=True, ignore_index=True)
    # Remove the path prefix
    unique_df['Relative Directory'] = unique_df['Directory'].apply(lambda x: x.replace(library_path + os.path.sep, ''))

    return unique_df


def update_db(library_path, pickle_db=True):
    # Identify the changes in the database
    changes = analyze_db(library_path)
    df_toadd = get_added_df(changes, library_path)


    # Load current database
    if not os.path.exists(os.path.join(library_path, 'PPIDB_full.db')):
        print('Database not found – creating a new one')
        cnx = sqlite3.connect(os.path.join(library_path, 'PPIDB_full.db'))

    try:
        if pickle_db:
            df_in = pd.read_pickle(os.path.join(library_path, 'PPIDB_full.pkl'))
        else:
            df_in = pd.read_sql('AF2 full-lib', f"sqlite:///{os.path.join(library_path, 'PPIDB_full.db')}")

        # Remove the entries that were removed from the database
        for rem_path in changes['removed']:
            components = os.path.basename(rem_path).split('-')
            # Show row where 'Job name' contains both components
            index = df_in[
                df_in['Job name'].str.contains(components[0]) & df_in['Job name'].str.contains(components[1])].index
            df_in.drop(index, inplace=True)
    except Exception as _:
        print('Database not found – considering it do be empty')
        df_in = pd.DataFrame()

    # Add the new entries
    df = pd.concat([df_in, df_toadd], ignore_index=True)
    # print(df['Job name'].head())

    # Save and replace current database with updated one
    if pickle_db:
        df.to_pickle(os.path.join(library_path, 'PPIDB_full.pkl'))
    else:
        cnx = sqlite3.connect(os.path.join(library_path, 'PPIDB_full.db'))
        df.to_sql(name='AF2 full-lib', con=cnx, if_exists='replace', index=False)
