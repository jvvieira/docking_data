import pandas as pd


df_data = pd.read_csv('data/filtered_data.csv')
df_residues = pd.read_csv('data/ligand_residues.csv')

# Calculate statistics

# get the most frequent proteins with the mean value
most_frequent_proteins = df_data['protein'].value_counts().head(10)
mean_values = df_data.groupby('protein')['value'].mean().loc[most_frequent_proteins.index]

most_frequent_proteins_df = pd.DataFrame({
    'protein': most_frequent_proteins.index,
    'count': most_frequent_proteins.values,
    'mean_value': mean_values.values
})

print("Most Frequent Proteins with Mean Values:")
print(most_frequent_proteins_df)
print("\n")


# Get the most frequent ligands
most_frequent_residue = df_residues['residue'].value_counts().head(10)
print("Most Frequent Residues:")
print(most_frequent_residue.head(10))
print("\n")

# Where to find the top 1 residue
top_residue = most_frequent_residue.index[0]
top_residue_data = df_residues[df_residues['residue'] == top_residue]
print(f"Data for the most frequent residue ({top_residue}):")
print(top_residue_data.loc[:, ['protein', 'ligand']].drop_duplicates()) 