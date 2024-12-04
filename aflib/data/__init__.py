import json
import os
import pickle


with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'gene_groups.json'), "r") as read_file:
    # print('Reading gene_groups.json...')
    gene_groups = json.load(read_file)

with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'gene_group_colors.json'), "r") as read_color_file:
    gene_group_colors = json.load(read_color_file)

with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'DICT_GENES.pkl'), "rb") as read_genes_file:
    DICT_GENES = pickle.load(read_genes_file)
