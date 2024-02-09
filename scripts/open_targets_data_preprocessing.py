# Kontsioti, Maskell, Anderson & Pirmohamed, Identifying drug-drug interactions in spontaneous reports utilizing signal detection and biological plausibility aspects (2024)
# This script performs preprocessing of the data downloaded from Open Targets

from pyspark import SparkConf
from pyspark.sql import SparkSession
import pyspark.sql.functions as F
from itertools import combinations
import pandas as pd
import numpy as np
from sqlalchemy import create_engine
import psycopg2
import networkx as nx
from matplotlib import pylab
import matplotlib.pyplot as plt
from nxviz import CircosPlot
from collections import Counter
import requests

# A. Get drug-target associations from the Drug dataset

# Replace with the path to the molecule directory
drugs_data_path = "path/to/molecule/folder"

# establish spark connection
spark = (
    SparkSession.builder
    .master('local[*]')
    .getOrCreate()
)

# read Parquet files for drug data
drugs_data = spark.read.parquet(drugs_data_path)

# print drugs dataset schema
drugs_data.printSchema()

# convert to Pandas dataframe
drugs_df = drugs_data.toPandas()
# End spark session
spark.stop()
# Remove drugs with no associated targets
drugs_df_2 = drugs_df[drugs_df['linkedTargets'].notna()]
# Extract lists from pySpark rows
drugs_df_2.loc[:,'linkedTargets_rows'] = drugs_df_2.loc[:,'linkedTargets'].apply(lambda x: x.asDict()['rows'])
# Extract DrugBank IDs (if available)
drugs_df_2.loc[:,'DrugbankID'] = drugs_df_2.loc[:,'crossReferences'].apply(lambda x: x.get('drugbank', None) if x is not None else None)
# Keep only selected columns
drugs_df_clean = drugs_df_2[['name', 'synonyms', 'linkedTargets_rows']]
# Explode lists to individual rows
drug_target_df = drugs_df_clean.explode('linkedTargets_rows')
# Group by target and count associated drugs
targets_drug_counts = drug_target_df.groupby('linkedTargets_rows')['name'].nunique()

# B. Perform exact string matching to map drug names to RxNorm concepts

# Connect to PostgreSQL
# Replace with your PostgreSQL credentials
engine = create_engine(
        'postgresql://username:password@host:port/database')

conn = psycopg2.connect("dbname=database user=username")
cur = conn.cursor()

# Exact string matching of drug name to RxNorm
cur.execute("set search_path = cdmv5; \
            select c.concept_id, c.concept_name, c.concept_class_id, c.vocabulary_id \
            from concept c \
            where c.vocabulary_id in ('RxNorm', 'RxNorm Extension') \
            and c.invalid_reason is NULL;") 
row_1 = cur.fetchall() 
rxnorm_mapping_df = pd.DataFrame(row_1).astype(str)
rxnorm_mapping_df.columns = ['concept_id', 'concept_name', \
                             'concept_class_id', 'vocabulary_id']
rxnorm_mapping_df['concept_name'] = rxnorm_mapping_df['concept_name'].str.upper()

drugs_df_clean_with_rxnorm = drugs_df_clean.merge(rxnorm_mapping_df,
                                                       how = 'left',
                                                       left_on = 'name',
                                                       right_on = 'concept_name')

# C. Normalise all RxNorm concepts to RxNorm Ingredients 

# Use RxNorm relationships to get the ingredients
cur.execute("set search_path = cdmv5; \
            select cr.concept_id_1, c1.concept_name, c1.concept_class_id, \
            cr.concept_id_2, c2.concept_name, c2.concept_class_id, \
            cr.relationship_id \
            from concept_relationship cr\
            inner join concept c1 \
            on c1.concept_id = cr.concept_id_1 \
            inner join concept c2 \
            on c2.concept_id = cr.concept_id_2 \
            where cr.relationship_id in ('Tradename of', 'RxNorm has ing', \
                                         'Form of') \
            and c1.vocabulary_id in ('RxNorm', 'RxNorm Extension') \
            and c2.vocabulary_id in ('RxNorm', 'RxNorm Extension') \
            and c2.concept_class_id = 'Ingredient' \
            and c1.invalid_reason is NULL and c2.invalid_reason is NULL;") 
row_2 = cur.fetchall() 
rxnorm_relationship_df = pd.DataFrame(row_2).astype(str)
rxnorm_relationship_df.columns = ['concept_id_1', 'concept_name_1', \
                                  'concept_class_1', 'concept_id_2', \
                                  'concept_name_2', 'concept_class_2', \
                                  'relationship']
# Convert drug names to uppercase
rxnorm_relationship_df['concept_name_1'] = rxnorm_relationship_df['concept_name_1'].str.upper()
rxnorm_relationship_df['concept_name_2'] = rxnorm_relationship_df['concept_name_2'].str.upper()

rxnorm_rel_dict = dict(zip(rxnorm_relationship_df.concept_id_1, 
                           rxnorm_relationship_df.concept_id_2))

drugs_df_clean_with_rxnorm['RxNorm_ingredient'] = drugs_df_clean_with_rxnorm['concept_id'].map(rxnorm_rel_dict)

# For RxNorm ingredients, copy over the value from the concept_id column
cond = (drugs_df_clean_with_rxnorm['RxNorm_ingredient'].isnull()) & \
             (drugs_df_clean_with_rxnorm['concept_class_id'] == 'Ingredient')
drugs_df_clean_with_rxnorm['RxNorm_ingredient'] = \
    np.where(cond, drugs_df_clean_with_rxnorm['concept_id'],
             drugs_df_clean_with_rxnorm['RxNorm_ingredient'])

# D. Export a .csv file with the curated drug-target associations

# Construct the final df
drug_target_final_df = drugs_df_clean_with_rxnorm[['name','linkedTargets_rows',
                                                   'RxNorm_ingredient']] 
drug_target_final_df = drug_target_final_df.explode('linkedTargets_rows')
drug_target_final_df.dropna(inplace=True)
# Generate drug - target pairs
drug_target_pairs = drug_target_final_df[['RxNorm_ingredient',
                                          'linkedTargets_rows']].drop_duplicates().reset_index(drop = True)
drug_target_pairs.to_csv('data/open_targets/drug_target_pairs.csv',
                         index = False)

# E. Get target-AE associations from the Target dataset 

# Replace with the path to the targets directory
targets_data_path = "path/to/targets/folder"

# establish spark connection
spark = (
    SparkSession.builder
    .master('local[*]')
    .getOrCreate()
)

# read Parquet files for targets data
targets_data = spark.read.parquet(targets_data_path)

# print targets dataset schema
targets_data.printSchema()

# convert to Pandas dataframe
targets_df = targets_data.toPandas()
# End spark session
spark.stop()

# Extract lists from pySpark rows
targets_df = targets_df.explode('proteinIds')
targets_df.reset_index(inplace=True)

# Remove drugs with no safety liability data
target_safety_df = targets_df[targets_df['safetyLiabilities'].notna()]
# Extract lists from pySpark rows
target_safety_df = target_safety_df.explode('safetyLiabilities')
target_safety_df.reset_index(inplace=True)

# Extract events from pySpark rows and store them in a dataframe
events_list = []
event_ids_list = []
for i in range(0,len(target_safety_df)):
    event = target_safety_df.loc[i,'safetyLiabilities'].event
    eventId = target_safety_df.loc[i,'safetyLiabilities'].eventId
    events_list.append(event)
    event_ids_list.append(eventId)
events_df = pd.DataFrame(list(zip(events_list, event_ids_list)))

# Extract event names from PySpark rows
target_safety_df['safetyLiabilities_event'] = target_safety_df['safetyLiabilities'].apply(lambda x: x.asDict()['event'])
# Extract event IDs from PySpark rows
target_safety_df['safetyLiabilities_eventId'] = target_safety_df['safetyLiabilities'].apply(lambda x: x.asDict()['eventId'])
# Extract effects from PySpark rows
target_safety_df['safetyLiabilities_effects'] = target_safety_df['safetyLiabilities'].apply(lambda x: x.asDict()['effects'])
# Extract related biosystems from PySpark rows
target_safety_df['safetyLiabilities_biosystems'] = target_safety_df['safetyLiabilities'].apply(lambda x: x.asDict()['biosample'])

# Extract protein IDs from PySpark rows
target_safety_df['proteinIds_id'] = target_safety_df['proteinIds'].apply(lambda x: x.asDict()['id'] if x is not None else None)
# Extract protein ID source from PySpark rows
target_safety_df['proteinIds_source'] = target_safety_df['proteinIds'].apply(lambda x: x.asDict()['source'] if x is not None else None)

# Keep only selected columns
target_safety_df_clean = target_safety_df[['id', 'approvedSymbol', \
                                        'safetyLiabilities_event', \
                                        'safetyLiabilities_eventId', \
                                        'safetyLiabilities_effects', \
                                        'safetyLiabilities_biosystems', \
                                        'proteinIds_id', 'proteinIds_source']]
# Keep only mappings to Uniprot Swissprot
target_safety_df_clean = target_safety_df_clean[target_safety_df_clean.proteinIds_source == 'uniprot_swissprot']
# Reset index
target_safety_df_clean.reset_index(drop=True, inplace=True)

# F. Use the OLS API to get MedDRA synonyms of the HPO/EFO/Orphanet terms 

# Use OLS API to get mappings of event ID to MedDRA
base_url = 'https://www.ebi.ac.uk/ols/api/terms?id='

target_safety_df_clean['MedDRA_synonyms'] = ''
for index, row in target_safety_df_clean.iterrows():
    print(index)
    if row['safetyLiabilities_eventId'] is not None:
        code = row['safetyLiabilities_eventId'].replace('_', ':')
        first_response = requests.get(base_url+code)
        response_list = first_response.json()['_embedded']['terms']
        synonyms = []
        for item in response_list:
            if 'annotation' in item:
                synonym_term_list = item['annotation'].get('database_cross_reference', list())
                synonyms.append(synonym_term_list)
        # Flatten the list and convert it to a set
        flat_synonym_list = set([i for sublist in synonyms for i in sublist])  
        # Keep only MedDRA terms
        meddra_terms = filter(lambda x: x.startswith('MedDRA:'), flat_synonym_list)
        target_safety_df_clean.at[index, 'MedDRA_synonyms'] = list(meddra_terms)

events_to_meddra = target_safety_df_clean[['safetyLiabilities_event', \
                                           'safetyLiabilities_eventId', \
                                            'MedDRA_synonyms']].explode('MedDRA_synonyms')
# Remove duplicated rows
events_to_meddra.drop_duplicates(inplace=True)
# Find the unique number of MedDRA terms mapped
meddra_counts = events_to_meddra['MedDRA_synonyms'].value_counts()

# G. Use of external mapping file to get additional MedDRA-HPO mappings

# Use MedDRA to HPO mappings file from the HPO GitHub repository
meddra_to_hpo_mapping = pd.read_csv('data/mappings/hpo_to_meddra_Mapping.tsv',
                                    sep='\t')
# Keep only rows that are marked as 'checked' (Comment 1 column)
meddra_to_hpo_mapping = meddra_to_hpo_mapping[(meddra_to_hpo_mapping['Comment_1'] == 'checked')]
# Keep only rows that are marked as 'OK' or (blank) (Comment 2 column)
meddra_to_hpo_mapping = meddra_to_hpo_mapping[(meddra_to_hpo_mapping['Comment_2'] == 'OK')|(meddra_to_hpo_mapping['Comment_2'].isna())]
# Replace substring to match HPO code string
meddra_to_hpo_mapping['HPO_code'] = meddra_to_hpo_mapping['HPO_code'].apply(lambda x: x.replace(':','_')) 
# Second mapping
events_to_meddra = pd.merge(events_to_meddra, meddra_to_hpo_mapping[['MedDRA_name', 'HPO_code', 'Matching_type', 'Comment_1', 'Comment_2', 'Comment_3']], 
         how = 'left', left_on = 'safetyLiabilities_eventId', right_on = 'HPO_code')

events_to_meddra['safetyLiabilities_eventId'].fillna("None", inplace=True)
# Get mapping dataframes by vocabulary
# EFO terms
efo_to_meddra = events_to_meddra[events_to_meddra['safetyLiabilities_eventId'].str.startswith('EFO')]
# HPO terms
hpo_to_meddra = events_to_meddra[events_to_meddra['safetyLiabilities_eventId'].str.startswith('HP')]
# Orphanet terms
orphanet_to_meddra = events_to_meddra[events_to_meddra['safetyLiabilities_eventId'].str.startswith('Orphanet')]
# No mapping available
no_mapping = events_to_meddra[events_to_meddra['safetyLiabilities_eventId'] == 'None']

##############################################################################
# H. Count values of targets by AE and AEs by target

# Group by event and count associated targets
events_target_counts = target_safety_df_clean.groupby('safetyLiabilities_event')['approvedSymbol'].nunique()
# Group by target and count associated events
targets_event_counts = target_safety_df_clean.groupby('approvedSymbol')['safetyLiabilities_event'].nunique()

##############################################################################
# I. Map MedDRA terms to OHDSI concepts (PT level)

# Get MedDRA mappings (concept_name, concept_code and concept_id)
# Connect to PostgreSQL
# Replace with your PostgreSQL credentials
engine = create_engine(
        'postgresql://username:password@host:port/database')

conn = psycopg2.connect("dbname=database user=username")
cur = conn.cursor()

cur.execute("set search_path = cdmv5; \
            select c.concept_id, c.concept_name, \
            c.concept_code, c.concept_class_id \
            from concept c \
            where c.vocabulary_id = 'MedDRA';") 
row2 = cur.fetchall() 
event_mapping_df = pd.DataFrame(row2).astype(str)
event_mapping_df[1] = event_mapping_df[1].str.upper()
event_mapping_df.columns = ['OHDSI_concept_id', 'concept_name', \
                            'concept_code', 'MedDRA_type']
conn.close()

# Map MedDRA names (from MedDRA_name column) to MedDRA OHDSI concept ids
events_to_meddra['mapped_concept_id'] = events_to_meddra['MedDRA_name'].str.upper().replace(event_mapping_df.set_index('concept_name')['OHDSI_concept_id'])

# Remove common substring from MedDRA_synonyms column and map MedDRA ID
# to its respective OHDSI concept id
events_to_meddra.loc[events_to_meddra.mapped_concept_id.isna(), 'mapped_concept_id'] = events_to_meddra['MedDRA_synonyms'].str[7:].replace(event_mapping_df.set_index('concept_code')['OHDSI_concept_id'])
# Add a column with the concept_class
events_to_meddra['concept_class'] = events_to_meddra['mapped_concept_id'].replace(event_mapping_df.set_index('OHDSI_concept_id')['MedDRA_type'])

# Generate MedDRA mapping tables between different levels of hierarchy to PTs
# 1. MedDRA LLT to PT mapping
conn = psycopg2.connect("dbname=database user=username")
cur = conn.cursor()

cur.execute("set search_path = cdmv5; \
            select ca.ancestor_concept_id, ca.descendant_concept_id \
            from concept_ancestor ca \
            inner join concept c   \
            on ca.ancestor_concept_id = c.concept_id \
            where ca.min_levels_of_separation = 1 and \
            ca.max_levels_of_separation = 1 \
            and c.concept_class_id = 'PT' \
            and c.vocabulary_id = 'MedDRA';") 
row3 = cur.fetchall() 
conn.close()
llt_to_pt_df = pd.DataFrame(row3).astype(str)
llt_to_pt_df.set_index(1, inplace=True)
llt_to_pt_df.columns = ['PT']
llt_to_pt_df.index.name = 'LLT'

# 2. MedDRA HLT/HLGT to PT mapping
conn = psycopg2.connect("dbname=database user=username")
cur = conn.cursor()

cur.execute("set search_path = cdmv5; \
            select ca.ancestor_concept_id, ca.descendant_concept_id, c2.concept_name, c1.concept_name \
            from concept_ancestor ca \
            inner join concept c1   \
            on ca.descendant_concept_id = c1.concept_id \
            inner join concept c2    \
            on ca.ancestor_concept_id = c2.concept_id  \
            where ca.min_levels_of_separation in (1, 2) and \
            ca.max_levels_of_separation in (1, 2) \
            and c1.concept_class_id = 'PT' \
            and c1.vocabulary_id = 'MedDRA';") 
row4 = cur.fetchall() 
conn.close()
hlt_hlgt_to_pt_df = pd.DataFrame(row4).astype(str)
hlt_hlgt_to_pt_df.columns = ['HLT/HLGT_code', 'PT_code', 'HLT/HLGT_name', 'PT_name']


# Use dict to map LLTs to PTs in the events_to_meddra dataframe
events_to_meddra = events_to_meddra.merge(
    llt_to_pt_df, left_on = 'mapped_concept_id', right_on = llt_to_pt_df.index,
    how = 'left')
# Copy over the PTs to the new column
events_to_meddra.loc[events_to_meddra.concept_class == 'PT', 'PT'] = events_to_meddra['mapped_concept_id']

##############################################################################
# Use dict to map HLTs/HLGTs to PTs in the events_to_meddra dataframe
events_to_meddra = events_to_meddra.merge(
    hlt_hlgt_to_pt_df[['PT_code', 'HLT/HLGT_code']], left_on = 'mapped_concept_id', 
    right_on = 'HLT/HLGT_code', how = 'left')
# Copy over the values from PT_code column to the PT column
events_to_meddra['PT'].fillna(events_to_meddra['PT_code'], inplace = True)
# Remove the previously joined columns (PT_code and HLT/HLGT_code)
events_to_meddra.drop(['PT_code', 'HLT/HLGT_code'], axis=1, inplace=True)

# Create a final mapping df with event names mapped to PTs
event_to_pt_mapping = events_to_meddra[['safetyLiabilities_event', 'PT']]
event_to_pt_mapping['PT'].fillna('None', inplace = True)

##############################################################################
# J. Generate the final df of target safety information with mapped AEs to 
# OHDSI MedDRA concepts (if available) 

# Copy selected columns to a new dataframe
target_safety_final = target_safety_df_clean[['id', 'approvedSymbol', 
                                              'safetyLiabilities_event', 
                                              'safetyLiabilities_eventId', 
                                              'safetyLiabilities_effects',
                                              'proteinIds_id']]
# Map AEs using the event_to_pt_mapping df
target_safety_final = target_safety_final.merge(event_to_pt_mapping, 
                                                on = 'safetyLiabilities_event',
                                                how = 'left')
# Copy selected columns to a new dataframe
target_event_pairs = target_safety_final[['id', 'proteinIds_id', 'PT']]
# Remove rows that were not mapped to a PT
target_event_pairs = target_event_pairs[target_event_pairs.PT != 'None']
# Remove duplicated rows
target_event_pairs.drop_duplicates(inplace=True)
# Export df to a .csv file
target_event_pairs.to_csv('data/open_targets/target_event_pairs.csv', index = False)

# K. Merge drug-target and target-event data to generate drug-target-event triplets
# Merge dataframes using target IDs as a common column
d_t_e_triplets = pd.merge(drug_target_pairs, target_event_pairs, 
                          left_on = 'linkedTargets_rows',
                          right_on = 'id')
# Remove duplicated column with the target ID
d_t_e_triplets.drop('id', axis=1, inplace=True)
# Remove duplicated rows
d_t_e_triplets.drop_duplicates(inplace=True)
# Export df to a .csv file
d_t_e_triplets.to_csv('data/open_targets/drug_target_event_triplets.csv', index = False)
