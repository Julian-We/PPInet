import json
import os


with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'gene_groups.json'), "r") as read_file:
    # print('Reading gene_groups.json...')
    gene_groups = json.load(read_file)

with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'gene_group_colors.json'), "r") as read_color_file:
    gene_group_colors = json.load(read_color_file)
