#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 11:21:37 2026

@author: galen2
"""

import pandas as pd
genes = pd.read_csv("curated_gene_sets_populated.csv")["gene_symbol"].unique()
expr  = pd.read_csv("DESeq2_full_combined_ssRNAseq.csv", index_col=0)
expr[expr.index.isin(genes)].to_csv("expr_pathway_genes.csv")
