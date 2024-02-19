import networkx as nx
import pandas as pd
from pyvis.network import Network

df = pd.read_csv('/Users/izbuser/Downloads/data_from_SpeedPPI/julian_2/output/all_ppis_unfiltered.csv')
df[['interactor 1', 'interactor 2']] = df['ID'].str.extract(r'(\w+)-\d+-([\w]+)')
df["weight"] = df['pdockq'] * 10
g = nx.from_pandas_edgelist(df, source='interactor 1', target='interactor 2', edge_attr="weight")
net = Network(notebook=False)
net.from_nx(g)
net.show("Interactionnetwork_SpeedPPI.html")
