# Kontsioti, Maskell, Anderson & Pirmohamed, Identifying drug-drug interactions in spontaneous reports utilizing signal detection and biological plausibility aspects (2024)
# This script builds the biological attribute network by assembling data from DrugBank, Open Targets, and STRING DB.
# The network is then used to calculate the shortest path measure for: (i) the CRESCENDDI controls; and (ii) the drug pairs screened for DDI signals related to QT interval prolongation.

import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import psycopg2
from sqlalchemy import create_engine
import networkx as nx
from networkx_query import search_nodes, search_edges
import itertools
from itertools import permutations, product

tqdm.pandas()

# A. STRING data (PPI) Pre-processing 

# STRING data
# Read protein links file
protein_links_df = pd.read_csv('data/raw/string_db/9606.protein.links.v11.5.txt', sep=" ", header=[0])
# Read protein aliases file
protein_aliases_df = pd.read_csv('data/raw/string_db/9606.protein.aliases.v11.5.txt', sep='\t', header=[0])

ensp_to_ensg_dict = protein_aliases_df[protein_aliases_df['source'] == 'Ensembl_gene'].set_index('#string_protein_id').to_dict()['alias']

protein_links_df['protein1_alias'] = protein_links_df['protein1'].map(ensp_to_ensg_dict)
protein_links_df['protein2_alias'] = protein_links_df['protein2'].map(ensp_to_ensg_dict)
# Keep only PPIs of higher confidence (score >= 700)
conf_protein_links_df = protein_links_df[protein_links_df.combined_score >= 700]

# Get all nodes from PPI network
ppi_nodes = pd.unique(conf_protein_links_df[['protein1_alias', 'protein2_alias']].values.ravel('K'))

##############################################################################
# B. Read all types of links to be added to the network

# Read Open Targets data (target-AE)
target_AE_links_df = pd.read_csv('data/processed/open_targets/target_event_pairs.csv')
# Get AE nodes
ae_nodes = list(set(target_AE_links_df.PT))

# Read Open Targets data (drug-target)

drug_target_links_OT_df = pd.read_csv('data/processed/drug_target_pairs.csv')
# Get drug nodes
drug_OT_nodes = list(set(drug_target_links_OT_df.RxNorm_ingredient))

 # Read DrugBank data (drug-target)
drug_target_links_DB_df = pd.read_csv('data/processed/drugbank_drug_targets.csv')

# HGNC IDs to Ensembl Gene IDs mapping dict
hgnc_to_ensg_dict = protein_aliases_df[protein_aliases_df['source'] == 'Ensembl_HGNC_HGNC_ID']
hgnc_to_ensg_dict['Ensembl_Gene'] = hgnc_to_ensg_dict['#string_protein_id'].map(ensp_to_ensg_dict)
hgnc_to_ensg_dict = hgnc_to_ensg_dict[['alias', 'Ensembl_Gene']]

# Map target HGNC IDs to Ensembl Gene IDs in DrugBank data 
drug_target_links_DB_df_mapped = pd.merge(drug_target_links_DB_df, hgnc_to_ensg_dict,
                                          how = 'left', left_on = 'HGNC ID',
                                          right_on = 'alias')
# Drop unmapped rows (no Ensembl Gene ID available) - MIGHT NEED TO FIX
drug_target_links_DB_df_mapped.dropna(subset = ['Ensembl_Gene'], inplace = True)
# Drop unnecessary columns
drug_target_links_DB_df_mapped.drop('alias', axis = 1, inplace = True)
# Get drug nodes
drug_DB_nodes = list(set(drug_target_links_DB_df_mapped.RxNorm_ingredient))

###############################################################################
# C. Construct the network

G = nx.Graph()
G.add_nodes_from(ppi_nodes, att = 'protein')
G.add_edges_from(list(conf_protein_links_df[['protein1_alias', \
                                             'protein2_alias']].itertuples(index=False)), 
                 link = 'protein-protein')

G.add_nodes_from(drug_OT_nodes, att = 'drug')
G.add_edges_from(list(drug_target_links_OT_df.itertuples(index=False)),
                 link = 'drug-protein-OT')

G.add_nodes_from(ae_nodes, att = 'ae')
G.add_edges_from(list(target_AE_links_df[['id', 'PT']].itertuples(index=False)),
                 link = 'protein-ae')

# Filter by attribute/node type - Only PPIs
selected_nodes = [n for n, attrdict in G.nodes.items() if attrdict.get('att') == 'protein']
H = nx.subgraph(G, selected_nodes)

##############################################################################
# D. Calculate the drugs_to_event_path_node_union metric

# Read .csv file - SDA scores for CRESCENDDI Controls
crescenddi_ctls = pd.read_csv('output/crescenddi_sda_scores.csv', index_col = 0)
# Remove parentheses (first and last character)
crescenddi_ctls['DDE_TUPLE_COPY'] = crescenddi_ctls['DDE_TUPLE'].str[1:-1]
# Split ids to separate columns
crescenddi_ctls[['drug_1_id', 'drug_2_id', 'ae_id']] = crescenddi_ctls['DDE_TUPLE_COPY'].str.split(',', expand=True).astype('int')
# Remove DDE_TUPLE_COPY column
crescenddi_ctls.drop('DDE_TUPLE_COPY', axis = 1, inplace = True)
##############################################################################
# Read .csv file - SDA scores for screening signals of QT prolongation
qt_prolongation_scores = pd.read_csv('output/qt_prolongation_sda_scores.csv',
                                    index_col = 0)
# Remove parentheses (first and last character)
qt_prolongation_scores['DDE_TUPLE_COPY'] = qt_prolongation_scores['dde_tuple'].str[1:-1]
# Split ids to separate columns
qt_prolongation_scores[['drug_1_id', 'drug_2_id', 'ae_id']] = qt_prolongation_scores['DDE_TUPLE_COPY'].str.split(',', expand=True).astype('int')
# Remove DDE_TUPLE_COPY column
qt_prolongation_scores.drop('DDE_TUPLE_COPY', axis = 1, inplace = True)

# Function to calculate the number of nodes in the union of the shortest
# paths of each drug to the AE via connected targets
def drugs_to_event_path_node_union(d1_id, d2_id, ae_id, full_network, ppi_network):
    # Get all proteins that are directly connected to the drug_1
    drug_1_connectors = [u for u,v,e in full_network.edges(data=True)
                         if e['link'] == 'drug-protein-OT'
                         and v == d1_id]
    # Get all proteins that are directly connected to the drug_2
    drug_2_connectors = [u for u,v,e in full_network.edges(data=True)
                         if e['link'] == 'drug-protein-OT'
                         and v == d2_id]
    # Get all proteins that are directly connected to the AE (seed proteins)
    ae_connectors = [u for u,v,e in full_network.edges(data=True) 
                     if e['link'] == 'protein-ae' and v == ae_id]
    # An empty list to store the number of protein nodes that correspond to 
    # each subnetwork that connects the drugs to the AE
    union_nodes_no = []
    # If there are direct links from the drugs and AE to at least a protein in the network
    if len(ae_connectors) > 0 and len(drug_1_connectors) > 0 and len(drug_2_connectors) > 0:
        # Get all possible pairs of starting protein nodes (coming from drugs)
        c = product(drug_1_connectors, drug_2_connectors)
        # For each possible starting protein nodes (coming from AEs)
        for ae_conn in ae_connectors:
            # For each pair of starting protein nodes (coming from drugs)
            for pair in c:
                try:
                    # Get shortest path from drug1-protein node to ae-protein node
                    sp_d1_ae = nx.all_shortest_paths(ppi_network, source=pair[0],
                                                target=ae_conn)
                    # Get shortest path from drug2-protein node to ae-protein node
                    sp_d2_ae = nx.all_shortest_paths(ppi_network, source=pair[1],
                                                target=ae_conn)
                    # Create an empty list to store the combinations
                    unique_combinations = []
                    # Get all permutations of sp_d1_ae with length of sp_d2_ae
                    sp_prod = product(sp_d1_ae, sp_d2_ae)
                    for sp_pair in sp_prod:
                        # Get a list of the union of nodes contained in 
                        # the individual shortest paths
                        union_nodes_sp = tuple(set(sp_pair[0] + sp_pair[1]))
                        unique_combinations.append(union_nodes_sp)
                    union_nodes_no.append(min(map(len, unique_combinations)))
                except (nx.NodeNotFound, nx.NetworkXNoPath) as e:
                    pass
    if len(union_nodes_no) > 0:
        value = min(union_nodes_no)
    else:
        value = list()
        
    return value

## i. For CRESCENDDI controls
crescenddi_ctls['drug_to_ae_shortest_paths_node_no'] = \
    crescenddi_ctls.progress_apply(
        lambda x: drugs_to_event_path_node_union(x['drug_1_id'], 
                                                 x['drug_2_id'], 
                                                 x['ae_id'], G, H), axis = 1)
# Write to a .csv file
crescenddi_ctls.to_csv(output/crescenddi_banet_scores.csv')

# ii. For QT prolongation drug pairs
qt_prolongation_scores['drug_to_ae_shortest_paths_node_no'] = \
qt_prolongation_scores.progress_apply(
    lambda x: drugs_to_event_path_node_union(x['drug_1_id'], 
                                             x['drug_2_id'], 
                                             x['ae_id'], G, H), axis = 1)

# Connect to PostgreSQL
# Replace with your PostgreSQL credentials
engine = create_engine(
        'postgresql://username:password@host:port/database')

conn = psycopg2.connect("dbname=database user=username")
cur = conn.cursor()

# Map OHDSI concept ids to names
# Get RxNorm/RxNorm Extension mappings
cur.execute("set search_path = cdmv5; \
            select c.concept_id, c.concept_name \
            from concept c \
            where c.vocabulary_id in ('RxNorm', 'RxNorm Extension');") 
row3 = cur.fetchall() 
drug_mapping_df = pd.DataFrame(row3)
drug_mapping_df[0] = drug_mapping_df[0].astype(np.int32)
drug_mapping_df[1] = drug_mapping_df[1].str.upper()
drug_mapping_df.columns = ['OHDSI_concept_id', 'concept_name']
conn.close()

drug_mapping_dict = dict(zip(drug_mapping_df.OHDSI_concept_id, drug_mapping_df.concept_name))

# QT prolongation
# Map drug concept ids
qt_prolongation_scores['d1_name'] = qt_prolongation_scores['drug_1_id'].map(drug_mapping_dict)
qt_prolongation_scores['d2_name'] = qt_prolongation_scores['drug_2_id'].map(drug_mapping_dict)
# Write to a .csv file
qt_prolongation_scores.to_csv('qt_prolongation_banet_scrores.csv',
                                      index = False)
