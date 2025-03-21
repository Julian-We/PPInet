import numpy as np
from Bio import SeqIO
from aflib import basebuilder
from Bio.Seq import Seq
import re
import os
from tqdm.notebook import tqdm
# import requests # Keep requests in case one has to use a ZFIN search in the "search engine"
import requests
import json
import pickle
from itertools import combinations

# --- Variables to be set
generate_proteome_freshly = False
generate_canonical_lib = False

# ---

dir_path = os.path.dirname(os.path.realpath(__file__))

record = list(SeqIO.parse(os.path.join(dir_path, "fastas/Danio_rerio.GRCz11.pep.all.fa"), "fasta"))


def get_protein_name(protein):
    match = re.search('gene_symbol:(.*?)(\.|\s|$)', protein.description)
    if match:
        return {match.group(1).lower(): protein}
    else:
        # print(protein.description)
        return {protein.id: protein}


def all_equal(iterator):
    iterator = iter(iterator)
    try:
        first = next(iterator)
    except StopIteration:
        return True
    return all(first == x for x in iterator)


# os.path.join(dir_path, 'fastas/ZFIN_1.0.1.4_basicGeneInformation.json')
with open(os.path.join(dir_path, 'fastas/ZFIN_1.0.1.4_basicGeneInformation.json')) as jsn_fle:
    zfin_db = json.load(jsn_fle)
    # print(type(zfin_db))


def get_uniprot_seq(protein_id):
    base_url = "https://www.uniprot.org/uniprot/"
    query_url = f"{base_url}{protein_id}.fasta"

    res = requests.get(query_url)
    if res.status_code == 200:
        return ''.join(res.text.splitlines()[1:])
    else:
        print('Protein sequence not found!')


def get_unique_items(input_list):
    unique_items = []
    seen_items = set()

    for item in input_list:
        item_seq = item.seq
        if item_seq not in seen_items:
            unique_items.append(item)
            seen_items.add(item_seq)

    return unique_items


def get_update_sequence(q_gene, s_gene):
    """

    :param q_gene:  query gene
    :param s_gene: searched gene
    :return:
    """
    multi = []
    if list(q_gene.keys())[0] == s_gene:
        multi.append(q_gene.get(list(q_gene.keys())[0]).seq)
    for var_num, protein_seq in enumerate(get_unique_items(multi)):
        if var_num == 0 and protein_seq is not None and len(protein_seq) > 1:
            return {s_gene: SeqIO.SeqRecord(Seq(protein_seq), id=s_gene, description='')}
        elif var_num != 0 and protein_seq is not None and len(protein_seq) > 1:
            return {f'{s_gene}_{var_num}': SeqIO.SeqRecord(Seq(protein_seq), id=s_gene, description='')}


def dict_get(di, key):
    if di.get(key) is not None:
        return di.get(key)  # .seq


def phoenix(record_dict, gene_list):
    gnnme = list(record_dict.keys())[0]

    all_multiplcates = list(map(dict_get, gene_list, len(gene_list) * [gnnme]))
    all_multiplcates = [i for i in all_multiplcates if i is not None]

    unique_multis = get_unique_items(all_multiplcates)

    unique_multis = [q for q in unique_multis if q is not None]

    out_dict = {}
    for var_num, mult in enumerate(unique_multis):
        if var_num == 0 and mult is not None and len(mult) > 1:
            # out_dict.update({gnnme: SeqIO.SeqRecord(Seq(mult), id=gnnme, description='')})
            out_dict.update({gnnme: mult})
        elif var_num != 0 and mult is not None and len(mult) > 1:
            # out_dict.update({f'{gnnme}_{var_num}': SeqIO.SeqRecord(Seq(mult), id=gnnme, description='')})
            out_dict.update({f'{gnnme}_{var_num}': mult})
    return out_dict


DICT_GENES = {}  # This stays no matter what
list_genes = list(map(get_protein_name, record))

# # Version 1 of how to generate the the DICT_GENES via data from a proteome fasta
# for d in tqdm(list_genes):
#     registered = False
#     turn = 1
#     while not registered:
#         local_gene_name = list((d.keys()))[0]
#         if DICT_GENES.get(local_gene_name) is None:
#             DICT_GENES.update(d)
#             registered = True
#         else:
#             if DICT_GENES.get(f'{local_gene_name}_{turn}') is None:
#                 DICT_GENES.update({f'{local_gene_name}_{turn}': d.get(local_gene_name)})
#                 registered = True
#             else:
#                 turn +=1

# # Version 2 of how to create the DICT_GENES via ZFIN DB and UniProt Pulls
# for zfin_entry in tqdm(zfin_db['data']):
#     gene_name = zfin_entry['symbol']
#     # print(zfin_entry['basicGeneticEntity']['crossReferences'])
#     multiple_entries = []
#     for uniprot_entry in zfin_entry['basicGeneticEntity']['crossReferences']:
#         if 'uniprot' in uniprot_entry['id'].lower():
#             matcha = re.match( '(UniProt.{2}):(.{6})', uniprot_entry['id'])
#             sequence = get_uniprot_seq(matcha.group(2))
#             multiple_entries.append(sequence)
#     for var_num, protein_seq in enumerate(get_unique_items(multiple_entries)):
#         if var_num == 0 and protein_seq is not None and len(protein_seq) > 1:
#             DICT_GENES.update({gene_name: SeqIO.SeqRecord(Seq(protein_seq), id=gene_name, description='')})
#         elif var_num != 0 and protein_seq is not None and len(protein_seq) > 1:
#             DICT_GENES.update({f'{gene_name}_{var_num}': SeqIO.SeqRecord(Seq(protein_seq), id=gene_name, description='')})


# Version 1.5 via proteome fasta- filtering extra sequences
if generate_proteome_freshly:
    len_list_genes = len(list_genes)
    pre_dicts = list(map(phoenix, tqdm(list_genes), len_list_genes * [list_genes]))

    for pre_dict in pre_dicts:
        DICT_GENES.update(pre_dict)

    # Dump into json

    # Writing to sample.json
    with open(os.path.join(dir_path, "data", "DICT_GENES.pkl"), "wb") as out_jsn_fle:
        pickle.dump(DICT_GENES, out_jsn_fle)

with open(os.path.join(dir_path, "data", "DICT_GENES.pkl"), 'rb') as in_jsn_fle:
    # Reading from json file
    DICT_GENES = pickle.load(in_jsn_fle)
    if type(DICT_GENES) == type(dict):
        raise TypeError(f'DICT_GENES does not have the right type. Instead it has type {type(DICT_GENES)}')


# This section generates a list of canonical transcripts

def check_canonical_label(transcript_ids_set):
    transcript_ids = list(transcript_ids_set)
    # Ensembl REST API endpoint for fetching transcript information
    endpoint = f"https://rest.ensembl.org/lookup/id/"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    r = requests.post(endpoint, headers=headers, data='{"ids" : ' + '[' + ', '.join(
        f'"{elem}"' for elem in transcript_ids) + ']' + '}')  # str({ "ids" : transcript_ids })
    # Make a GET request to the endpoint
    # response = requests.get(endpoint, headers={ "Content-Type" : "application/json"})
    # Check if the request was successful (status code 200)
    if r.status_code == 200:
        out_dict = {}
        print(r.status_code)
        transcript_info = r.json()
        for key, value in transcript_info.items():
            out_dict.update({key: value.get('canonical_transcript')})
        # print(transcript_info)#.get('canonical_transcript'))
        return out_dict
    else:
        print(r.status_code)
        return None


if generate_canonical_lib:
    ensemble_names = []
    for entry in list(DICT_GENES.values()):
        desc = entry.description
        desc_cut_1 = desc[desc.find('gene:') + 5:]
        desc_cut_2 = desc_cut_1[:desc_cut_1.find('.')]
        ensemble_names.append(desc_cut_2)

    # Create subsets
    list_of_ensmbl_subsets = []
    counter = 0
    subset = set()
    check_sum = 0
    neg_sum = 0
    for gene_name in ensemble_names:
        if counter < 50:
            subset.add(gene_name)
            counter += 1
        elif len(gene_name) != 18:
            neg_sum += 1
        else:
            list_of_ensmbl_subsets.append(subset)
            check_sum += len(subset)
            # print(len(subset))
            subset = {gene_name}
            counter = 1
    if len(subset) != 0:
        list_of_ensmbl_subsets.append(subset)

    list_canonical_dicts = list(map(check_canonical_label, tqdm(list_of_ensmbl_subsets)))

    CANONICAL_DICT = {}
    for canonical_gene in list_canonical_dicts:
        CANONICAL_DICT.update(canonical_gene)

    with open(os.path.join(dir_path, "data", "canonical_transcripts.pkl"), "wb") as out_jsn_fle_canon:
        pickle.dump(CANONICAL_DICT, out_jsn_fle_canon)

with open(os.path.join(dir_path, "data", "canonical_transcripts.pkl"), "rb") as in_jsn_fle_canon:
    CANONICAL_DICT = pickle.load(in_jsn_fle_canon)


def search_engine(search_query):
    foundit = False
    search_query = search_query.lower()
    while not foundit:
        if DICT_GENES.get(search_query) is not None:
            # print(f'found {search_query}')
            regex_search = f'{search_query}' + '_\d{1,3}'
            hit_list = [search_query]
            for some_gene_name in DICT_GENES.keys():
                if re.match(regex_search, some_gene_name) is not None:
                    hit_list.append(some_gene_name)
            foundit = True
            return hit_list

        response_selection = []
        for gene_name, gene_seqio in DICT_GENES.items():
            if search_query in gene_seqio.description.lower():
                if re.search('(.*)(_\d*)', gene_name) is None:
                    response_selection.append(gene_name)
            elif search_query.replace(" ", "") in gene_seqio.description.lower().replace(" ", ""):
                response_selection.append(gene_name)
            elif search_query in gene_name.lower():
                response_selection.append(gene_name)
        variants_dict = {}
        # print(response_selection)
        for variant in response_selection:
            does_match = re.match(r'(.*?)(_\d*)', variant)
            if does_match:
                prefix_name = does_match.group(1)
                variants_dict.setdefault(prefix_name, []).append(variant)
            else:
                variants_dict.update({variant: [variant]})
        if list(variants_dict.keys()):
            precise_name = input(
                f'Search query: {search_query}\n If your gene is in this list type its exact name:{list(variants_dict.keys())} \n If not, press "ok".')
        else:
            print(f'No result for {search_query} was found, thus its not considered')
            return [None]

        if precise_name == '':
            search_query = input('Please enter an alternative search query or give up by pressing "Ok"')
            if search_query == '':
                foundit = True
        else:
            if DICT_GENES.get(precise_name) is not None:
                return variants_dict.get(precise_name)
            else:
                search_query = input(f'The name was not found in the databank. Perhaps a typo in {precise_name}?')


def searchengine_master(candidate_list, canonical):
    clean_candidate_list = []
    for candidate in candidate_list:
        se_result = search_engine(candidate)
        for ses_result in se_result:
            if ses_result is not None:
                se_desc = DICT_GENES.get(ses_result).description
                # print(se_desc)
                if se_desc.startswith('**override**'):
                    clean_candidate_list.append(ses_result)
                    print('#################', ses_result)
                desc_cut_1 = se_desc[se_desc.find('gene:') + 5:]
                gene_name = desc_cut_1[:desc_cut_1.find('.')]
                se_tr_cut1 = se_desc[se_desc.find('transcript:') + 11:]
                transcript_name = se_tr_cut1[:se_tr_cut1.find(' ')]
                if ses_result is not None and canonical and CANONICAL_DICT.get(gene_name) == transcript_name:
                    clean_candidate_list.append(ses_result)
                elif ses_result is not None and not canonical:
                    clean_candidate_list.append(ses_result)
                else:
                    print(f'No Gene named <{ses_result}> was found in the canonical list. Please add yourself via palmatyalation.generate_pair()')

    return clean_candidate_list


def af_computetime(len_seq):
    # print(len_seq)
    range_to_m = {
        (0, 100): 4.9,
        (100, 200): 7.7,
        (200, 300): 13,
        (300, 400): 18,
        (400, 500): 29,
        (500, 600): 39,
        (600, 700): 53,
        (700, 800): 60,
        (800, 900): 91,
        (900, 1000): 96,
        (1000, 1100): 140,
        (1100, 1500): 280,
        (1500, 2000): 450,
        (2000, 2500): 969,
        (2500, 3000): 1240,
        (3000, 3500): 2465,
        (3500, 4000): 5660,
        (4000, 4500): 12475,
        (4500, 5000): 18824,
        (5000, 12000): 50000
    }

    for (low, high), value in range_to_m.items():
        if low < len_seq < high:
            ### Outdated version form AF2 website
            # len_prediction = range_to_m.get((low, high)) * 5 * 3
            # # Taking 0.2min per residue for the MSA as an estimation

            len_prediction = 17.2 * len_seq

            msa_prediction = 4.2 * len_seq

            # pred_time_min = (-300*np.exp(-0.006 * len_seq)) + (0.24 * len_seq) + 300
            # pred_time_min = (-300*np.exp(-0.006 * len_seq)) + (0.4 * len_seq) + 320
            pred_time_min = 0.01 * len_seq ** 1.625 + 200
            pred_time_min = 10000 if pred_time_min > 10000 else pred_time_min
            pred_time_sec = pred_time_min * 60

            # msa_prediction = 20 * len_seq + (2*10**-5 * len_seq**2)
            # predicted time +10%  tolerance
            # return (len_prediction + msa_prediction) * 1.1
            return pred_time_sec
    return 0


def target_system(target_gene, missile_gene, separation='entry'):
    """
    Function to generate a fasta file for a target and missile gene
    :param target_gene: One part of a pairwise prediction
    :param missile_gene:
    :param separation:
    :return:
    """
    try:
        missile_sequence = str(DICT_GENES.get(missile_gene).seq)
        target_sequence = str(DICT_GENES.get(target_gene).seq)
        if separation == ':':
            return {''.join([missile_gene, '-', target_gene]): ''.join([missile_sequence, ':', target_sequence])}
        else:
            missile_object = SeqIO.SeqRecord(Seq(missile_sequence), id=missile_gene, description='')
            target_object = SeqIO.SeqRecord(Seq(target_sequence), id=target_gene, description='')
            return {''.join([missile_gene, '-', target_gene]): [missile_object, target_object]}
    except AttributeError as e:
        print(f'Problem with mseq: {missile_gene}; tseq: {target_gene}')


def seconds_to_hms(secs):
    hours, remainder = divmod(secs, 3600)
    minutes, seconds = divmod(remainder, 60)
    return "{:02d}:{:02d}:{:02d}".format(int(hours), int(minutes), int(seconds))


def af_statistics(seq_dict, buffer_temp=0.075, **kwargs):
    """
    :param seq_dict:
    :param buffer_temp:
    :return:
    """
    str_total_len = f'The are a total of {len(seq_dict)} fasta files being folded'
    # print(str_total_len)
    total_legnth = 0
    total_runtime = 0
    for sequence in seq_dict.values():
        if type(sequence) != type(list()):
            tmp_len = len(sequence) - 1
        else:
            tmp_len = 0
            for sequ in sequence:
                tmp_len += len(sequ)
        total_legnth += tmp_len
        total_runtime += af_computetime(tmp_len)
        # print(seconds_to_hms(round(af_computetime(tmp_len), -2)))
    total_runtime += total_runtime * buffer_temp  # Plus 12% buffer in time since a few jobs timed out
    str_avg_len = f'With an average lenth of {total_legnth / len(seq_dict)}'
    str_time = f'And a predicted total runtime of {seconds_to_hms(total_runtime)}'
    l_master_str = [str_total_len, '\n', str_avg_len, '\n', str_time]
    # return ''.join(l_master_str)
    return [seconds_to_hms(total_runtime), total_runtime, total_legnth]


# THIS IS THE OG create_bash_script FUNCTION
def create_bash_script(experiment_name, comb_name, run_time, path, partition="gpu2080", user="jwegner1@uni-muenster.de",
                       nodes=1, cores=10, gres=1, memory=60, precomputed_msas=False, database='full_dbs'):
    """

    :param database:
    :param experiment_name:
    :param comb_name:
    :param run_time:
    :param path:
    :param partition:
    :param user:
    :param nodes:
    :param cores:
    :param gres:
    :param memory:
    :param precomputed_msas:
    :return:
    """

    job_name = f"job_{comb_name}"
    pre_msas = str(precomputed_msas).lower()

    # bash_script_lines = ["#!/bin/bash \n",
    #                      f"#SBATCH --partition={partition} \n",
    #                      f"#SBATCH --nodes={nodes} \n",
    #                      f"#SBATCH --gres=gpu:{gres} \n",
    #                      f"#SBATCH --cpus-per-task={cores} \n",
    #                      f"#SBATCH --mem={memory}G \n",
    #                      f"#SBATCH --time={run_time} \n",
    #                      f"#SBATCH --job-name={job_name} \n",
    #                      "#SBATCH --account=uni \n",
    #                      "#SBATCH --mail-type=ALL \n",
    #                      f"#SBATCH --mail-user={user} \n",
    #                      " \n",
    #                      "module load palma/2021a \n",
    #                      "module load foss/2021a \n",
    #                      "module load AlphaFold/2.1.2 \n",
    #                      "wait \n",
    #                      "export ALPHAFOLD_DATA_DIR=/Applic.HPC/data/alphafold \n",
    #                      "\n",
    #                      "alphafold \\n",
    #                      f"    --fasta_paths=/scratch/tmp/jwegner1/{experiment_name}/fasta/{comb_name}.fasta \\n",
    #                      "    --model_preset=multimer \\n",
    #                      f"    --output_dir=/scratch/tmp/jwegner1/{experiment_name}/xprt \\n",
    #                      f"    --use_precomputed_msas={pre_msas} \\n",
    #                      "    --max_template_date=2021-11-25 \\n",
    #                      "    --is_prokaryote_list=false \\n",
    #                      "    --db_preset=reduced_dbs \\n",
    #                      "    --data_dir=/Applic.HPC/data/alphafold\n"
    #                      ]




    bash_script_text = f"""#!/bin/bash
#SBATCH --partition={partition}
#SBATCH --nodes={nodes}
#SBATCH --gres=gpu:{gres}
#SBATCH --cpus-per-task={cores}
#SBATCH --mem={memory}G
#SBATCH --time={run_time}
#SBATCH --job-name={job_name}
#SBATCH --account=uni
#SBATCH --mail-type=ALL
#SBATCH --mail-user={user}

module load palma/2022a
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load AlphaFold/2.3.1-CUDA-11.7.0
wait
export ALPHAFOLD_DATA_DIR=/Applic.HPC/data/alphafold-2.3.0

alphafold \\
    --fasta_paths=/scratch/tmp/jwegner1/{experiment_name}/fasta/{comb_name}.fasta \\
    --model_preset=multimer \\
    --output_dir=/scratch/tmp/jwegner1/{experiment_name}/xprt \\
    --max_template_date=2021-11-25 \\
    --use_precomputed_msas={pre_msas} \\
    --num_multimer_predictions_per_model=3 \\
    --db_preset={database} \\
    --data_dir=/Applic.HPC/data/alphafold-2.3.0
"""

    with open(os.path.join(path, f'{job_name}.sh'), 'w+') as bsh:
        # bsh.writelines(bash_script_lines)
        bsh.write(bash_script_text)
    return os.path.join(f'/scratch/tmp/jwegner1/{experiment_name}/bash_scripts/', f'{job_name}.sh')



def export_fastas(pre_fasta_dict, output_directory, time_buffer_factor=1.5, t_min=6000, **kwargs):
    try:
        fasta_path = os.path.join(output_directory, 'fasta')
        os.makedirs(fasta_path)
    except FileExistsError:
        pass
    for fa_name, fa_sequence in pre_fasta_dict.items():
        fa_path = os.path.join(fasta_path, f'{fa_name}.fasta')
        if type(fa_sequence) == type(str()):
            fa_record = SeqIO.SeqRecord(Seq(fa_sequence), id=fa_name, description='')
            SeqIO.write(fa_record, fa_path, 'fasta-2line')
        elif type(fa_sequence) == type(list()):
            # print(fa_sequence)
            SeqIO.write(fa_sequence, fa_path, 'fasta-2line')
    # msg = af_statistics(pre_fasta_dict)
    try:
        bash_path = os.path.join(output_directory, 'bash_scripts')
        os.makedirs(bash_path)
        os.makedirs(os.path.join(output_directory, 'xprt'))
    except FileExistsError:
        pass
    master_bash_lines = ['#!/bin/bash \n']
    # Quick and dirty time to assesment
    tot_runtime = 0
    for k, v in pre_fasta_dict.items():
        msg = af_statistics({k: v}, **kwargs)
        if msg[1] < t_min:
            print(f'----- {k} underbet t_min with : \t {msg[0]}')
            msg = [seconds_to_hms(t_min), t_min]
            bash_mem = 60
        print(f'{k} needs: \t {msg[0]}')
        try:
            if msg[2] > 2500:
                bash_mem = 140
            elif 2500 > msg[2] > 2250:
                bash_mem = 120
            elif 2250 > msg[2] > 1750:
                bash_mem = 100
            elif 1750 > msg[2]:
                bash_mem = 80
            else:
                bash_mem = 120
            # bash_mem = 120 if msg[2] > 1600 else 60
        except IndexError:
            pass
        sh_pth = create_bash_script(os.path.basename(output_directory), k, msg[0], bash_path,
                                    memory=bash_mem)
        master_bash_lines.append('sbatch ' + sh_pth + ' \n')
        # print(msg)
        tot_runtime += msg[1]

    with open(os.path.join(output_directory, f'dealer.sh'), 'w+') as bsh:
        bsh.writelines(master_bash_lines)

    for x in range(3):
        print('')
    print(f'The job will take approximately {seconds_to_hms(tot_runtime)}')
    print(f'Folding {len(pre_fasta_dict)} pairs')

    # with open(os.path.join(output_directory,'computing_prediction.txt'), "w") as cp_file:
    #     print(msg)
    #     cp_file.write(msg)
    #     cp_file.close()


def some_vs_some(output_directory: str, genes_group1: list, genes_group2: list, only_canonical: bool,
                 both_ways=False, db_path=None,
                 **kwargs):
    """
    Generate fasta files to feed into alphafold 2.
    :param both_ways:
    :param only_canonical:
    :param output_directory:
    :param genes_group1:
    :param genes_group2:
    :param seq_mode:
    :param db_path: Path to the library
    :return:
    """
    clean_gene_group1 = searchengine_master(genes_group1, only_canonical)
    clean_gene_group2 = searchengine_master(genes_group2, only_canonical)
    print(clean_gene_group1, clean_gene_group2)

    if db_path is not None:
        if db_path.endswith('.db'):
            df = basebuilder.get_table_asdf(db_path, 'AF2DB_DETAILED')
            root = os.path.dirname(db_path) + os.path.sep
        else:
            df = basebuilder.get_table_asdf(os.path.join(db_path, 'af2db.db'), 'AF2DB_DETAILED')
            root = db_path + os.path.sep
        # with open(os.path.join(db_path, 'PPIDB_full.pkl'), 'rb') as db_fle:
        #     df = pickle.load(db_fle)
    else:
        raise ValueError('No database path provided')
    df['Relative Directory'] = df['Directory'].apply(lambda x: x.replace(root, ''))
    # For each gene in this list generate all interacting
    pre_fasta_dict = {}
    for gene1 in clean_gene_group1:
        # print('tick')
        if db_path is not None:
            exclusion_list = []
            screen_gene = gene1
            for candidate in genes_group2:
                for combination in df['Relative Directory']:
                    combination = combination.split('/')[-1]
                    if candidate in combination and screen_gene in combination and candidate != screen_gene:
                        print(f'Excluding combination: {candidate}-{screen_gene}')
                        exclusion_list.append(candidate)

            superclean_gene_group2 = [gene for gene in clean_gene_group2 if gene not in exclusion_list]
        else:
            superclean_gene_group2 = clean_gene_group2

        combination_list = list(map(target_system, superclean_gene_group2, len(superclean_gene_group2) * [gene1]))
        for combi_dict in combination_list:
            pre_fasta_dict.update(combi_dict)
    if both_ways:
        for gene2r in clean_gene_group2:
            # print('tick')
            combination_list = list(map(target_system, clean_gene_group1, len(clean_gene_group1) * [gene2r]))
            for combi_dict in combination_list:
                pre_fasta_dict.update(combi_dict)

    export_fastas(pre_fasta_dict, output_directory, **kwargs)


def all_vs_all(output_directory: str, gene_group: list, only_canonical: bool, both_ways=False):
    """
    Searches for transcripts of listed genes within the ensembl zebrafish DB (optionally discards non-canonical).
    Creates in the output directory: a "dealer" shell script,
    that can be executed on PALMA II (calls shell scripts from subdir "bash_scripts").
    Creates Bash scripts that start a job on PALMA II executing alphafold.
    Creates FASTA files with each combination
    :param both_ways:
    :param output_directory:
    :param gene_group:
    :param only_canonical:
    :return:
    """
    clean_gene_group = searchengine_master(gene_group, only_canonical)
    # print(clean_gene_group)

    # For each gene in this list generate all interacting
    pre_fasta_dict = {}

    for gene1, gene2 in combinations(clean_gene_group, 2):
        pre_fasta_dict.update(target_system(gene1, gene2))

    # for gene1 in clean_gene_group:
    #     # print('tick')
    #     combination_list = list(map(target_system, clean_gene_group, len(clean_gene_group) * [gene1]))
    #     for combi_dict in combination_list:
    #         pre_fasta_dict.update(combi_dict)

    export_fastas(pre_fasta_dict, output_directory)


def generate_pair(output_diectory, gene_name1, gene_name2, sequence1=None, sequence2=None, only_canonical=True):
    """
    Generates a fasta file for a pair of genes
    :param output_diectory: output directory
    :param gene_name1: gene name of sequnce 1
    :param gene_name2: gene name of sequnce 2
    :param sequence1: sequence of gene 1
    :param sequence2: sequence of gene 2
    :return:
    """
    if sequence1 is None:
        search_name1 = searchengine_master([gene_name1], only_canonical)
        sequence1 = str(dict_get(DICT_GENES, search_name1[0]).seq)

    if sequence2 is None:
        search_name2 = searchengine_master([gene_name2], only_canonical)
        sequence2 = str(dict_get(DICT_GENES, search_name2[0]).seq)

    pre_fasta_dict = {}
    if sequence1 is not None and sequence2 is not None:
        missile_object = SeqIO.SeqRecord(Seq(sequence1), id=gene_name1, description='')
        target_object = SeqIO.SeqRecord(Seq(sequence2), id=gene_name2, description='')
        pre_fasta_dict.update({''.join([gene_name1, '-', gene_name2]): [missile_object, target_object]})
    else:
        pre_fasta_dict.update(target_system(gene_name1, gene_name2))

    export_fastas(pre_fasta_dict, output_diectory)


# generate multiple pairs from a list of genes and one gene with one sequence
def generate_multiple_pairs(output_directory, gene_list, gene_name, sequence, only_canonical=True):
    """
    Generates a fasta file for a pair of genes
    :param output_directory: output directory
    :param gene_list: list of gene names
    :param gene_name: gene name of sequnce 2
    :param sequence: sequence of gene 2
    :return:
    """
    if sequence is None:
        search_name = searchengine_master([gene_name], only_canonical)
        sequence = str(dict_get(DICT_GENES, search_name[0]).seq)

    pre_fasta_dict = {}
    if sequence is not None:
        target_object = SeqIO.SeqRecord(Seq(sequence), id=gene_name, description='')
        clean_list = searchengine_master(gene_list, only_canonical)
        for gene in clean_list:
            missile_object = SeqIO.SeqRecord(Seq(str(dict_get(DICT_GENES, gene).seq)), id=gene, description='')
            pre_fasta_dict.update({''.join([gene, '-', gene_name]): [missile_object, target_object]})

    export_fastas(pre_fasta_dict, output_directory)
