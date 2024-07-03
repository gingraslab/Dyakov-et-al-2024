import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# # Read the data
# domain_df = pd.read_csv('pfam_dataset_gene_names.tsv', delimiter='\t')
# domain_df = domain_df[domain_df['type'] == 'Domain']

# cluster_df = pd.read_csv('7013_cleaned_v2_NMF_19_top10-07-005_clusters.csv')


# # Convert domain data to a dictionary mapping gene names to domains
# gene_to_domain = domain_df.groupby('gene_name')['hmm name'].apply(set).to_dict()

# # Prepare the background set (all genes with domain information)
# background_genes = set(gene_to_domain.keys())

# # Prepare clusters
# clusters = {col: set(cluster_df[col].dropna()) for col in cluster_df.columns}

# # Initialize results dataframe
# results = []

# # Perform analysis for each cluster
# for rank, genes in clusters.items():
#     # Filter genes that have domain information
#     cluster_genes_with_domain = genes.intersection(background_genes)
#     query_size = len(cluster_genes_with_domain)

#     # Count domains in this cluster
#     domain_counts = {}
#     for gene in cluster_genes_with_domain:
#         for domain in gene_to_domain[gene]:
#             if domain in domain_counts:
#                 domain_counts[domain]['overlap_size'] += 1
#             else:
#                 domain_counts[domain] = {'term_size': 0, 'overlap_size': 1}
    
#     # Count domains in the background
#     for gene, domains in gene_to_domain.items():
#         for domain in domains:
#             if domain in domain_counts:
#                 domain_counts[domain]['term_size'] += 1
#             else:
#                 domain_counts[domain] = {'term_size': 1, 'overlap_size': 0}
    
#     # Calculate enrichment and p-values
#     for domain, counts in domain_counts.items():
#         # Only consider domains that appear in the cluster
#         if counts['overlap_size'] > 0:
#             term_size = counts['term_size']
#             overlap_size = counts['overlap_size']
#             oddsratio, p_value = fisher_exact([
#                 [overlap_size, term_size - overlap_size],
#                 [query_size - overlap_size, len(background_genes) - (query_size + (term_size - overlap_size))]
#             ])
#             # Store the results
#             results.append({
#                 'NMF rank': rank,
#                 'protein domain': domain,
#                 'term_size': term_size,
#                 'query size': query_size,
#                 'overlap size': overlap_size,
#                 'fold enrichment': (overlap_size / term_size) / (query_size / len(background_genes)),
#                 'p-value': p_value,
#                 'intersection': ','.join([gene for gene in cluster_genes_with_domain if domain in gene_to_domain[gene]])
#             })

# # Convert results to DataFrame
# results_df = pd.DataFrame(results)

# # Function to apply BH correction within each group
# def apply_bh_correction(group):
#     # Apply Benjamini-Hochberg correction to the p-values within the group
#     corrected_pvals = multipletests(group['p-value'], method='fdr_bh')[1]
#     group['adj. p-value'] = corrected_pvals
#     return group

# # Group by 'NMF rank' and apply the correction
# results_df = results_df.groupby('NMF rank').apply(apply_bh_correction)



# # print(results_df[['NMF rank', 'protein domain', 'term_size', 'query size', 'overlap size', 'fold enrichment', 'p-value', 'adj. p-value', 'intersection']])

# results_df.to_csv('domain_enrichment_results_7013_2.csv', index=False)

# # Renaming columns according to the new specification
# results_df.rename(columns={
#     'protein domain': 'term_name',
#     'term_size': 'term_size',
#     'query size': 'query_size',
#     'overlap size': 'intersection_size',
#     'fold enrichment': 'fold enrichment',
#     'p-value': 'p-value',
#     'intersection': 'intersection',
#     'adj. p-value': 'adjusted_p_value'
# }, inplace=True)

# # #Using ExcelWriter to write each cluster's results to a separate sheet
# # with pd.ExcelWriter('domain_results_by_rank_2.xlsx') as writer:
# #     for rank, group in results_df.groupby('NMF rank'):
# #         # Write each group to a separate sheet named by its rank
# #         group.to_excel(writer, sheet_name=f'Rank {rank}', index=False)


# Read the data
domain_df = pd.read_csv('pfam_dataset_gene_names.tsv', delimiter='\t')
domain_df = domain_df[domain_df['type'] == 'Domain']

cluster_df = pd.read_csv('7013_cleaned_v2_NMF_19_top10-07-005_clusters.csv')

# Convert domain data to a dictionary mapping gene names to domains
gene_to_domain = domain_df.groupby('gene_name')['hmm name'].apply(set).to_dict()

# Prepare the background set (all genes with domain information)
background_genes = set(gene_to_domain.keys())

# Prepare clusters
clusters = {col: set(cluster_df[col].dropna()) for col in cluster_df.columns}

# Initialize results dataframe
results = []

# Perform analysis for each cluster
for rank, genes in clusters.items():
    # Filter genes that have domain information
    cluster_genes_with_domain = genes.intersection(background_genes)
    query_size = len(cluster_genes_with_domain)

    # Count domains in this cluster
    domain_counts = {}
    for gene in cluster_genes_with_domain:
        for domain in gene_to_domain[gene]:
            if domain in domain_counts:
                domain_counts[domain]['overlap_size'] += 1
            else:
                domain_counts[domain] = {'term_size': 0, 'overlap_size': 1}
    
    # Count domains in the background
    for gene, domains in gene_to_domain.items():
        for domain in domains:
            if domain in domain_counts:
                domain_counts[domain]['term_size'] += 1
            else:
                domain_counts[domain] = {'term_size': 1, 'overlap_size': 0}
    
    # Calculate enrichment and p-values
    for domain, counts in domain_counts.items():
        # Only consider domains that appear in the cluster
        if counts['overlap_size'] > 0:
            term_size = counts['term_size']
            overlap_size = counts['overlap_size']
            background_with_domain = term_size - overlap_size
            background_without_domain = len(background_genes) - term_size
            
            # Check if query_size or background_with_domain is zero
            if query_size == 0 or background_with_domain == 0:
                fold_enrichment = float('inf')  # Use infinity or another large number to indicate undefined enrichment
            else:
                fold_enrichment = (overlap_size / query_size) / (background_with_domain / len(background_genes))

            # Setting up the contingency table
            table = [
                [overlap_size, background_with_domain],  # Domain present in cluster, Domain present in background
                [query_size - overlap_size, background_without_domain]  # Domain absent in cluster, Domain absent in background
            ]

            # Perform Fisher's exact test
            oddsratio, p_value = fisher_exact(table)
            
            # Store the results
            results.append({
                'NMF rank': rank,
                'protein domain': domain,
                'term_size': term_size,
                'query_size': query_size,
                'overlap size': overlap_size,
                'fold enrichment': fold_enrichment,
                'p-value': p_value,
                'intersection': ','.join([gene for gene in cluster_genes_with_domain if domain in gene_to_domain[gene]])
            })
# Convert results to DataFrame
results_df = pd.DataFrame(results)

# Function to apply BH correction within each group
def apply_bh_correction(group):
    # Apply Benjamini-Hochberg correction to the p-values within the group
    corrected_pvals = multipletests(group['p-value'], method='fdr_bh')[1]
    group['adj. p-value'] = corrected_pvals
    return group

# Group by 'NMF rank' and apply the correction
results_df = results_df.groupby('NMF rank').apply(apply_bh_correction)

# Renaming columns according to the new specification
results_df.rename(columns={
    'protein domain': 'term_name',
    'term_size': 'term_size',
    'query_size': 'query_size',
    'overlap size': 'intersection_size',
    'fold enrichment': 'fold enrichment',
    'p-value': 'p-value',
    'intersection': 'intersection',
    'adj. p-value': 'adjusted_p_value'
}, inplace=True)

# Using ExcelWriter to write each cluster's results to a separate sheet
# with pd.ExcelWriter('domain_results_by_rank_3.xlsx') as writer:
#     for rank, group in results_df.groupby('NMF rank'):
#         # Write each group to a separate sheet named by its rank
#         group.to_excel(writer, sheet_name=f'Rank {rank}', index=False)

# # Save the dataframe to CSV
results_df.to_csv('domain_enrichment_results_7013_3.csv', index=False)