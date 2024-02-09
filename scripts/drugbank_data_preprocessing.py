# Kontsioti, Maskell, Anderson & Pirmohamed, Identifying drug-drug interactions in spontaneous reports utilizing signal detection and biological plausibility aspects (2024)
# This script performs preprocessing of the data downloaded from DrugBank
  
## Adapted from: https://caseolap.github.io/docs/drug/drugbank/#parsing-drugs-data
import re
import itertools
import json
import sys
import os
import time
import traceback
from lxml import etree
import xml.etree.ElementTree as ET
import json as json
import pandas as pd
import numpy as np
from sqlalchemy import create_engine
import psycopg2

# Function to extract text from a tag
def get_text(element, tag):
    e = element.find(tag)
    if e is not None:
        return e.text
    else:
        return '' 
    
# Load xml file
database = open("data/raw/drugbank/full database.xml", 'r', encoding="utf8")
tree = ET.parse(database)
root = tree.getroot()

k = 0
f = open("data/processed/alldrugs.txt", 'w', encoding="utf8")
Data = []
name = None
for drug in root:
    k = k + 1
    name = drug.find("{http://www.drugbank.ca}name")
    if name is not None:
        d_name = name.text  
        line  = name.text

    #-----------------Targets ------------------------------------------------
    targets_element = drug.find("{http://www.drugbank.ca}targets")
    T1 = []
    if targets_element is not None:
        targets = targets_element.findall("{http://www.drugbank.ca}target")
        for target in targets:
            polypeptides = target.findall("{http://www.drugbank.ca}polypeptide")
            if polypeptides is None:
                continue
            for polypeptide in polypeptides:
                t = {}
                organism = polypeptide.find("{http://www.drugbank.ca}organism")
                uniprotid = polypeptide.attrib["id"]
                source = polypeptide.attrib["source"]
                t.update({"uniprot_id": uniprotid, \
                          "source":source, "organism":organism.text})
                T1.append(t)  
    #-----------------Enzymes ------------------------------------------------
    enzymes = drug.find("{http://www.drugbank.ca}enzymes")
    E = []
    for enzyme in enzymes:
        actions = enzyme.find("{http://www.drugbank.ca}actions")
        for action in iter(actions):
            e = {}
            for item,n in zip(enzyme,["id","name"]):
                e.update({n:item.text})
                e.update({"action":action.text})
            E.append(e)  
    #-----------------Transporters -------------------------------------------
    transporters = drug.find("{http://www.drugbank.ca}transporters")
    T2 = []
    for transporter in transporters:
        actions = transporter.find("{http://www.drugbank.ca}actions")
        for action in iter(actions):
            t = {}
            for item,n in zip(transporter,["id","name"]):
                t.update({n:item.text})
                t.update({"action":action.text})
            T2.append(t) 

    Data.append({"name":d_name,\
                #"description":d_description,\
                "targets":T1,\
                "transporters":T2,\
                "enzymes":E
                }
                )
    f.write(line)
    f.write("\n")

# Write drugs data to a json file
with open("data/processed/Drugs.json", "w") as f:
    json.dump(Data,f)

# Convert to a dataframe
df = pd.DataFrame(Data)
# Rename columns
df.columns = ['drug_name', 'targets', 'transporters', 'enzymes']

# Remove rows that contain only empty lists (i.e. no information found)
df_not_nan = df[(df['targets'].str.len() != 0) | (df['transporters'].str.len() != 0) |
                (df['enzymes'].str.len() != 0)]
# Convert drug name to uppercase
df_not_nan['drug_name'] = df_not_nan['drug_name'].str.upper()

## Connect to PostgreSQL database 
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

drugbank_with_rxnorm = df_not_nan.merge(rxnorm_mapping_df, how = 'left',
                                        left_on = 'drug_name', 
                                        right_on = 'concept_name')

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

drugbank_with_rxnorm['RxNorm_ingredient'] = drugbank_with_rxnorm['concept_id'].map(rxnorm_rel_dict)

# For RxNorm ingredients, copy over the value from the concept_id column
cond = (drugbank_with_rxnorm['RxNorm_ingredient'].isnull()) & \
             (drugbank_with_rxnorm['concept_class_id'] == 'Ingredient')
drugbank_with_rxnorm['RxNorm_ingredient'] = \
    np.where(cond, drugbank_with_rxnorm['concept_id'],
             drugbank_with_rxnorm['RxNorm_ingredient'])
    
# =============================================================================
# Enzyme dataframe
df1 = pd.concat([pd.DataFrame(i) for i in drugbank_with_rxnorm['enzymes']], 
                keys = drugbank_with_rxnorm.index).reset_index(level=1,drop=True)
df1.columns = ['enzyme_id','enzyme_action','enzyme_name']
enzyme_df = drugbank_with_rxnorm.drop(['enzymes','transporters','targets'], 
                                      axis=1).join(df1).reset_index(drop=True)
# Remove nan rows
enzyme_df.dropna(inplace=True)

# Read mapping data from DrugBank (csv downloaded directly from: 
# https://go.drugbank.com/releases/latest#protein-identifiers)
enzyme_mapping_table = pd.read_csv('data/raw/drugbank/drugbank_enzyme_mappings.csv')
enzyme_mapping_table.drop(['GenBank Protein ID', 'GenBank Gene ID', 'PDB ID', \
                           'GeneCard ID', 'GenAtlas ID', 'Drug IDs'], axis = 1,
                          inplace = True)
# Drop duplicate rows
enzyme_mapping_table.drop_duplicates(inplace = True)
# Add UniProt IDs and HGNC IDs
enzyme_df = enzyme_df.merge(enzyme_mapping_table, how = 'left',
                           left_on = 'enzyme_name', right_on = 'Name').reset_index(drop=True)
# Remove nan (unmapped) rows
enzyme_df.dropna(inplace=True)

# Transporter dataframe
df2 = pd.concat([pd.DataFrame(i) for i in drugbank_with_rxnorm['transporters']], 
                keys = drugbank_with_rxnorm.index).reset_index(level=1,drop=True)
df2.columns = ['transporter_id','transporter_action','transporter_name']
transporter_df = drugbank_with_rxnorm.drop(['enzymes', 'transporters', 'targets'], 
                                           axis=1).join(df2).reset_index(drop=True)
# Remove nan rows
transporter_df.dropna(inplace=True)

# Read mapping data from DrugBank (csv downloaded directly from: 
# https://go.drugbank.com/releases/latest#protein-identifiers)
transporter_mapping_table = pd.read_csv('data/raw/drugbank/drugbank_transporter_mappings.csv')
transporter_mapping_table.drop(['GenBank Protein ID', 'GenBank Gene ID', 'PDB ID', \
                           'GeneCard ID', 'GenAtlas ID', 'Drug IDs'], \
                               axis = 1, inplace = True)
# Drop duplicate rows
transporter_mapping_table.drop_duplicates(inplace = True)
# Add UniProt IDs and HGNC IDs
transporter_df = transporter_df.merge(transporter_mapping_table, how = 'left',
                           left_on = 'transporter_name', right_on = 'Name').reset_index(drop=True)
# Remove nan (unmapped) rows
transporter_df.dropna(inplace=True)

# Target dataframe
df3 = pd.concat([pd.DataFrame(i) for i in drugbank_with_rxnorm['targets']], 
                keys = drugbank_with_rxnorm.index).reset_index(level=1,drop=True)
df3.columns = ['target_id', 'source', 'organism']
# Restrict to only human targets
df3 = df3.loc[df3.organism == 'Humans']
# Join human target information to drug-target data
target_df = drugbank_with_rxnorm.drop(['enzymes','transporters','targets'], 
                                      axis=1).join(df3).reset_index(drop=True)
# Remove nan rows
target_df.dropna(inplace=True)

# Read mapping data from DrugBank (csv downloaded directly from: 
# https://go.drugbank.com/releases/latest#protein-identifiers)
target_mapping_table = pd.read_csv('data/raw/drugbank/drugbank_target_mappings.csv')
target_mapping_table = target_mapping_table[['ID', 'Name', 'UniProt ID', 'HGNC ID']]
                              
# Drop duplicate rows
target_mapping_table.drop_duplicates(inplace = True)
# Add UniProt IDs and HGNC IDs
target_df = target_df.merge(target_mapping_table, how = 'left',
                           left_on = 'target_id', right_on = 'UniProt ID').reset_index(drop=True)
# Remove nan (unmapped) rows
target_df.dropna(inplace=True)

##############################################################################
# Export dfs to csv files
enzyme_df.to_csv('data/processed/drugbank/drugbank_drug_enzymes.csv', index=False)
transporter_df.to_csv('data/processed/drugbank/drugbank_drug_transporters.csv', index=False)
target_df.to_csv('data/processed/drugbank/drugbank_drug_targets.csv', index=False)
