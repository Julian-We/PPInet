import pandas as pd
import os
import json
from aflib.docking import compute_docking_quality
import aflib.data as data
from aflib.basebuilder import update_db
import plotly.express as px
import seaborn as sns
from itertools import combinations
from tqdm import tqdm
import networkx as nx
import matplotlib.pyplot as plt
from pycirclize import Circos
from pycirclize.parser import Matrix
import pickle



dir_path = os.path.dirname(os.path.realpath(__file__))

with open(os.path.join(dir_path, "data", "DICT_GENES.pkl"), 'rb') as in_jsn_fle:
    # Reading from json file
    DICT_GENES = pickle.load(in_jsn_fle)
    if type(DICT_GENES) == type(dict):
        raise TypeError(f'DICT_GENES does not have the right type. Instead it has type {type(DICT_GENES)}')


def get_database(db_path):
    update_db(db_path, pickle_db=True)
    return pd.read_pickle(db_path)


class CombinationPrediction:
    def __init__(self, path):
        self.chainA_name, self.chainB_name = os.path.basename(path).split('-') if len(
            os.path.basename(path).split('-')) == 2 else [None, None]
        self.path = path
        with open(os.path.join(self.path, 'msas', 'chain_id_map.json')) as seq_map:
            seq_map_json = json.load(seq_map)
            self.seqA = seq_map_json["A"]["sequence"]
            self.seqB = seq_map_json["B"]["sequence"]

    def process_docking(self):
        with open(os.path.join(self.path, 'ranking_debug.json')) as json_file:
            ranking_json = json.load(json_file)

        pdockq_scores = []
        model_ranks = []
        pdock_facs = []
        for num, rank in enumerate(ranking_json['order']):
            scores = compute_docking_quality(os.path.join(self.path, f'relaxed_{rank}.pdb'),
                                             os.path.join(self.path, f'result_{rank}.pkl'))
            pdockq_score = scores['pdockq']
            pdockq_scores.append(pdockq_score)
            model_rank = ranking_json['iptm+ptm'][rank]
            model_ranks.append(model_rank)
            pdock_fac = pdockq_score * model_rank
            pdock_facs.append(pdock_fac)

            if num == 0:
                highest_rank_pdockq = pdockq_score
                highest_iptm_score = scores['iptm']
                highest_ptm_score = scores['ptm']
                avg_plddt = scores['avg plddt']

        average_pdockq = sum(pdockq_scores) / len(pdockq_scores)
        ranked_average = sum(pdock_facs) / sum(model_ranks)

        dict_out = {
            "Job name": os.path.basename(self.path),
            "Directory": self.path,
            "Chain A": self.chainA_name,
            "Chain B": self.chainB_name,
            "Max pDockQ": float(highest_rank_pdockq),
            "Avg pDockQ": float(average_pdockq),
            "weighted pDockQ": float(ranked_average),
            "Max ipTM": float(highest_iptm_score),
            "Max pTM": float(highest_ptm_score),
            "Chain A sequence": self.seqA,
            "Chain A length": len(self.seqA),
            "Chain B sequence": self.seqB,
            "Chain B length": len(self.seqB),
            "Average pLDDT": avg_plddt,
        }

        with open(os.path.join(self.path, 'docking_results.json'), 'w') as json_file:
            json.dump(dict_out, json_file, indent=4)
        return {
            "Job name": os.path.basename(self.path),
            "Directory": self.path,
            "Chain A": self.chainA_name,
            "Chain B": self.chainB_name,
            "Max pDockQ": highest_rank_pdockq,
            "Avg pDockQ": average_pdockq,
            "weighted pDockQ": ranked_average,
            "Max ipTM": highest_iptm_score,
            "Max pTM": highest_ptm_score,
            "Chain A sequence": self.seqA,
            "Chain A length": len(self.seqA),
            "Chain B sequence": self.seqB,
            "Chain B length": len(self.seqB),
            "Average pLDDT": avg_plddt
        }


def helper_generate(path):
    try:
        exp = CombinationPrediction(path)
        return exp.process_docking()
    except FileNotFoundError as fileerr:
        print(f'The combination of the following sequences {os.path.basename(path)} have not been able to be processed'
              f' to the the following error')
        print(fileerr)
    except pickle.UnpicklingError as pickleerr:
        print(f'The combination of the following sequences {os.path.basename(path)} have not been able to be processed'
              f' to the the following error')
        print(pickleerr)


def load_library_aslist(lib_path):
    if isinstance(lib_path, str):
        list_exp_path = [os.path.join(lib_path, exp_pth) for exp_pth in os.listdir(lib_path) if
                         os.path.isdir(os.path.join(lib_path, exp_pth)) and not exp_pth.startswith('.')]
    elif isinstance(lib_path, list):
        list_exp_path = lib_path
    else:
        raise ValueError('lib_path must be a string or a list of strings')

    none_list_results = list(tqdm(map(helper_generate, list_exp_path), total=len(list_exp_path)))

    # Clean none_list_results from None values where helper function failed
    list_results = [result for result in none_list_results if result is not None]

    # # Legacy: As for loop in the (slower than mapping)
    # list_results = []
    # for exp_pth in tqdm(list_exp_path):
    #     exp = CombinationPrediction(exp_pth)
    #     list_results.append(exp.process_docking())

    return list_results


def uniquify(df):
    candidate_set = set()
    for chain in ['Chain A', 'Chain B']:
        candidate_set.update(df[chain].tolist())

    unique_combinations = list(combinations(list(candidate_set), 2))
    dict_list = []
    for chaina, chainb in unique_combinations:
        out_dict = {}
        # df_comb = df[df['Job name'].str.contains(chaina) & df['Job name'].str.contains(chainb)]
        df_comb = df.copy()
        df_comb['split_string'] = df_comb['Job name'].str.split('-')
        df_comb = df_comb[df_comb['split_string'].apply(lambda x: chaina in x and chainb in x)]
        parta_split = chaina.split('_')[0]
        partb_split = chainb.split('_')[0]
        if not df_comb.empty and not parta_split == partb_split:
            # print(chaina, chainb)
            out_dict["Chain A"] = chaina
            out_dict["Chain B"] = chainb
            out_dict["Gene A"] = chaina.split('_')[0]
            out_dict["Gene B"] = chainb.split('_')[0]
            out_dict["Job name"] = f'{chaina}-{chainb}'
            out_dict["Directory"] = df_comb["Directory"].values[0]
            out_dict["Mean pDockQ"] = round(df_comb["Max pDockQ"].mean(), 3)
            out_dict["pTM"] = round(df_comb["Max pTM"].mean(), 3)
            out_dict["ipTM"] = round(df_comb["Max ipTM"].mean(), 3)
            if len(df_comb['Max pDockQ'].values) == 2:
                fwd, rev = df_comb['Max pDockQ'].values
                # print(f'Difference: {round(abs(fwd - rev), 3)}')
                out_dict["delta pDockQ"] = round(abs(fwd - rev), 3)
                out_dict["pDockQ values"] = [fwd, rev]
                out_dict["Max pTM values"] = df_comb["Max pTM"].values
                out_dict["Max ipTM values"] = df_comb["Max ipTM"].values
                out_dict["delta pTM"] = abs(df_comb["Max pTM"].values[0] - df_comb["Max pTM"].values[1])
                out_dict["delta ipTM"] = abs(df_comb["Max ipTM"].values[0] - df_comb["Max ipTM"].values[1])
            else:
                out_dict["delta pDockQ"] = None
                out_dict["pDockQ values"] = None
            # print(f'Average pDockQ: {round(df_comb["Highest rank pDockQ score"].mean(), 3)}')
            # print('\n')
            if chaina == df_comb["Chain A"].values[0]:
                out_dict["Chain A sequence"] = df_comb["Chain A sequence"].values[0]
                out_dict["Chain B sequence"] = df_comb["Chain B sequence"].values[0]
                out_dict["Chain A length"] = int(df_comb["Chain A length"].values[0])
                out_dict["Chain B length"] = int(df_comb["Chain B length"].values[0])
            else:
                out_dict["Chain A sequence"] = df_comb["Chain B sequence"].values[0]
                out_dict["Chain B sequence"] = df_comb["Chain A sequence"].values[0]
                out_dict["Chain A length"] = int(df_comb["Chain B length"].values[0])
                out_dict["Chain B length"] = int(df_comb["Chain A length"].values[0])
        dict_list.append(out_dict)
    return pd.DataFrame(dict_list)


def search_gene_groups(row, static_gene):
    query_gene = None
    if row['Chain A'].lower() != static_gene.lower():
        query_gene = row['Chain A']
    elif row['Chain B'].lower() != static_gene.lower():
        query_gene = row['Chain B']
    else:
        raise ValueError(f'No Gene found in row, that does not match {static_gene}')

    for group, genes_list in data.gene_groups.items():
        if query_gene.lower().split('_')[0] in [x.lower() for x in genes_list]:
            # print(query_gene, 'is in', group)
            return group
    print(query_gene, 'is not in any group')
    return 'None'


pyfile_path = os.path.dirname(os.path.realpath(__file__))

with open(os.path.join(pyfile_path, "data", "gene_groups.json"), "r") as read_file:
    gene_groups = json.load(read_file)


class Interactions:
    def __init__(self, *args, mode=None, **kwargs):
        """
        Class to analyse the output of AF2 on PALMA
        :param mode: Select from a sereies of modes to analyse the output. If none assumptions are made

        'screen': Input one gene as a string. Searches the entire library for interactions with that gene. Additionally,
        if a list is given, the function will search for all interactions between the genes in the list.
        'network': Input one list of genes. Searches the entire library for interactions between the genes in the list.
        'reload': Input a single path. Reloads the library from the path.
        :param path: path to directory containing the output-folder from AF2 on PALMA
        """
        self.mode = mode
        self.path = None
        self.unique_df = None
        self.library_list = None
        if mode is None:
            if len(args) == 1 and os.path.isdir(args[0]):
                self.path = args[0]
                self.mode = 'reload'
            elif len(args) == 1 and not os.path.isdir(args[0]):
                mode = 'screen'
            elif len(args) == 1 and isinstance(args[0], list):
                mode = 'network'
            elif len(args) == 2 and isinstance(args[0], str) and isinstance(args[1], list):
                mode = 'screen'
            else:
                raise ValueError('No valid input given. Please provide a mode of operation: screen, network or reload')

        if mode == 'reload':
            self.both_ways()
        elif mode == 'screen':
            self.screen(*args, **kwargs)
        elif mode == 'network':
            self.network(*args, **kwargs)

    def both_ways(self):
        """
        This function will generate interaction libraries for both the forward and reverse chains.
        :param path: path to the directory containing all export folders from a
        :return: dataframe containing all unique combinations of chains and their respective pDockQ scores and (i)pTM scores
        """
        res = load_library_aslist(self.path)
        unique_df = uniquify(pd.DataFrame(res))

        self.library_list = res
        unique_df.dropna(subset=['Mean pDockQ', 'pTM', 'ipTM'], inplace=True, ignore_index=True)
        self.unique_df = unique_df

    def network(self, gene_list, db_path=None):
        if db_path is None:
            raise ValueError('No database path given. Use db_path=... to specify the path to the database')
        df = get_database(db_path)
        self.unique_df = df[(df['Gene A'].isin(gene_list) & df['Gene B'].isin(gene_list))]

        # Check if all possible combinations are present in unique_df
        for x, y in combinations(gene_list, 2):
            present = False
            if ((self.unique_df['Gene A'] == x) & (self.unique_df['Gene B'] == y)).any():
                continue
            elif ((self.unique_df['Gene A'] == y) & (self.unique_df['Gene B'] == x)).any():
                continue
            else:
                print(f'No interaction found between {y} and {x}')

    def screen(self, *args, db_path=None, **kwargs):
        if db_path is None:
            raise ValueError('No database path given. Use db_path=... to specify the path to the database')
        df = get_database(db_path)
        if len(args) == 1:
            gene = args[0]
            if isinstance(gene, str):
                self.unique_df = df[df['Gene A'].str.contains(gene) | df['Gene B'].str.contains(gene)]
            else:
                raise ValueError('Gene must be a string')
        elif len(args) == 2:
            gene = args[0]
            candidate_genes = args[1]
            if isinstance(gene, str) and isinstance(candidate_genes, list):
                self.unique_df = df[(df['Gene A'].str.contains(gene) & df['Gene B'].isin(candidate_genes)) |
                                    (df['Gene B'].str.contains(gene) & df['Gene A'].isin(candidate_genes))]
            else:
                raise ValueError('Gene must be a string and candidate_genes must be a list')
        else:
            raise ValueError('No valid input given. Please provide a gene as a string '
                             'or a gene and a list of candidate genes')

    def pTM_plot(self, axis, experiment_name="Experiment", **kwargs):
        """
        Plots
        :param axis: Matplotlib axis, the sns.heatmap is plotted on
        :param experiment_name: name of the experiment which is then written in the title
        :return:
        """
        sns.scatterplot(data=self.unique_df, x='pTM', y='ipTM', ax=axis, **kwargs)
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1)
        axis.set_title(f"pTM/ipTM scores for {experiment_name}")

    def interactive_pTM_plot(self, html_output_path, color_by_group=False,
                             static_gene=None,
                             color_map=data.gene_group_colors, **kwargs):
        """
        Plots
        :param color_map: dictionary containing the color map for the gene groups
        :param html_output_path: Output directory for the html file
        :param color_by_group: If True, the scatter plot will be colored by gene group
        :param static_gene: Gene, that is in every combination. Needs to be set if color_by_group is True
        :return:
        """

        if color_by_group and static_gene is None:
            raise ValueError('If color_by_group is True, static_gene must be a string')

        if color_by_group:
            self.unique_df['gene group'] = self.unique_df.apply(search_gene_groups, axis=1, static_gene=static_gene)

        fig = px.scatter(self.unique_df,
                         x="pTM",
                         y="ipTM",
                         hover_name='Job name',
                         range_x=[0, 1],
                         range_y=[0, 1],
                         template="simple_white",
                         color='gene group' if color_by_group else None,
                         color_discrete_map=color_map,
                         opacity=0.75,
                         **kwargs)

        fig.update_layout(autosize=False, width=800, height=800)
        fig.write_html(html_output_path)
        fig.show()

    def matrix_plot(self, axis, parameter):
        try:
            df_pivoted = pd.DataFrame(self.library_list).pivot(index="Chain A",
                                                               columns="Chain B",
                                                               values=parameter).astype('float')
        except KeyError as keyerr:
            print(f'KeyError: {keyerr}')
            df_pivoted = self.unique_df.pivot_table(index="Chain A", columns="Chain B", values=parameter).astype('float')
        sns.heatmap(df_pivoted, annot=True, ax=axis)
        axis.set_xlabel('Chain B – second entry')
        axis.set_ylabel('Chain A – first entry')

    def plot_networkx(self, axis, cmap_pTM='gray_r'):
        g = nx.from_pandas_edgelist(self.unique_df, source='Chain A', target='Chain B',
                                    edge_attr=['pTM', 'ipTM'])
        wlabels = [i['ipTM'] * 10 for i in dict(g.edges).values()]
        wcolor = [i['pTM'] for i in dict(g.edges).values()]
        nx.draw_circular(g,
                         ax=axis,
                         with_labels=True,
                         width=wlabels,
                         node_size=3500,
                         node_color='#e3aaaa',
                         node_shape='o',
                         edge_color=wcolor,
                         edge_cmap=plt.colormaps.get_cmap(cmap_pTM),
                         font_family='times')

    def plot_relationship_circle(self, output_path):
        matrix = Matrix.parse_fromto_table(self.unique_df[['Chain A', 'Chain B', 'ipTM']])

        circos = Circos.initialize_from_matrix(
            matrix,
            space=3,
            cmap="viridis",
            ticks_interval=5,
            label_kws=dict(size=12, r=110),
            link_kws=dict(direction=0, ec="black", lw=0.4),
        )

        circos.savefig(os.path.join(output_path, "circular_relationships.png"))
        circos.savefig(os.path.join(output_path, "circular_relationships.pdf"))
