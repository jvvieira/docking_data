import os
from spicy import stats
import warnings
import pandas as pd
import shutil
    

percentile = 0.05
warnings.filterwarnings("ignore")

folders = ['DHFR','RAR','RXR','Ciclofilna']

files_found = []
for root, dirs, files in os.walk('./data/'):
    for name in files:
        if name.endswith('.sc'):
            files_found.append(os.path.join(root, name))
            # print(os.path.join(root, name))

rows_isc = []
rows_delta_interface = []

file_number = 1
for file in files_found:
    # print(f'Processing file: {file_number} / {len(files_found)}  {file}')
    file_number += 1
    with open(file,'r') as f:
        lines = f.readlines()
        for line in lines:
            if(line.find('total_score') != -1 or line.find('SEQUENCE') != -1):
                continue            
            
            if(lines[1].find('I_sc') != -1):
                rows_isc.append(line.split())
            
            elif(lines[1].find('interface_delta') != -1):
                rows_delta_interface.append(line.split())


## Load values from isc
df_isc = pd.DataFrame(rows_isc, columns=['SCORE','total_score','I_bsa','I_hb','I_pack','I_sc','I_unsat','dslf_fa13','fa_atr','fa_dun','fa_elec','fa_intra_rep','fa_intra_sol_xover4','fa_rep','fa_sol','hbond_bb_sc','hbond_lr_bb','hbond_sc','hbond_sr_bb','lk_ball_wtd','omega','p_aa_pp','pep_sc','pep_sc_noref','pro_close','rama_prepro','ref','reweighted_sc','rmsALL','rmsALL_CAPRI_if','rmsALL_allIF','rmsALL_if','rmsBB','rmsBB_CAPRI_if','rmsBB_allIF','rmsBB_if','rmsBB_intrf_lowres','rmsBB_lowres','rmsCA','rmsCA_if','rmsCA_intrf_lowres','rmsCA_lowres','rmsSC_CAPRI_if','rmsSC_allIF','score_lowres_opt','score_lowres_start','startRMSall','startRMSallif','startRMSbb','startRMSca','yhh_planarity','description'])
print(f'Total number of cases (ISC ): {len(df_isc)}')

## Load values from delta_interface
df_delta_interface = pd.DataFrame(rows_delta_interface, columns=['SCORE','total_score','angle_constraint','atom_pair_constraint','chainbreak','coordinate_constraint','dihedral_constraint','dslf_ca_dih','dslf_cs_ang','dslf_ss_dih','dslf_ss_dst','fa_atr','fa_dun','fa_elec','fa_intra_rep','fa_pair','fa_rep','fa_sol','frac_atoms_within_0.5','frac_atoms_within_1.0','frac_atoms_within_2.0','hbond_bb_sc','hbond_lr_bb','hbond_sc','hbond_sr_bb','if_angle_constraint','if_atom_pair_constraint','if_buried_sasa','if_buried_unsat_hbonds','if_chainbreak','if_coordinate_constraint','if_dihedral_constraint','if_dslf_ca_dih','if_dslf_cs_ang','if_dslf_ss_dih','if_dslf_ss_dst','if_fa_atr','if_fa_dun','if_fa_elec','if_fa_intra_rep','if_fa_pair','if_fa_rep','if_fa_sol','if_hbond_bb_sc','if_hbond_lr_bb','if_hbond_sc','if_hbond_sr_bb','if_omega','if_p_aa_pp','if_pro_close','if_rama','if_ref','interface_delta','ligand_auto_rms_no_super','ligand_auto_rms_with_super','ligand_centroid_travel','ligand_is_touching','ligand_num_chi','ligand_radius_of_gyration','omega','p_aa_pp','pro_close','rama','ref','description'])
print(f'Total number of cases (Delta Interface): {len(df_delta_interface)}')

df_isc['value'] = df_isc['I_sc'].astype(float)
df_delta_interface['value'] = df_delta_interface['interface_delta'].astype(float)

## Merge the dataframes
df_merged = pd.concat([df_isc.loc[:, ['value','description']], df_delta_interface.loc[:, ['value','description']]], ignore_index=True)
print(f'Total number of cases (ISC + Delta Interface): {len(df_merged)}')


## Filter non-binding cases (positive values)
df_merged = df_merged[df_merged['value'] < 0]
print(f'Number of binding cases: {len(df_merged)}')

## Normality test
df_merged['value'] = df_merged['value'].astype(float)
stat, p = stats.shapiro(df_merged['value'])
print(f'Shapiro-Wilk test statistic: {stat}, p-value: {p}')


percentile_value = df_merged['value'].quantile(percentile)
print(f'{100 - int(percentile * 100)}% percentile value: {percentile_value}')

## Filter only lower then the mean
df_merged = df_merged[df_merged['value'] < percentile_value]
print(f'Number of cases with value lower than the percentile: {len(df_merged)}')

df_merged['protein'] = df_merged['description'].apply(lambda x: x.split("_")[0])
df_merged['ligand'] = df_merged['description'].apply(lambda x: x.split("_")[1])
df_merged['model'] = df_merged['description'].apply(lambda x: x.split("_")[2])
df_merged['pose'] = df_merged['description'].apply(lambda x: x.split("_")[4])

## Best Comples 1
best_complex1 = df_merged.groupby('protein').count()
print(f'Best Complex 1: {best_complex1.sort_values("value", ascending=False).index[0]} with {best_complex1.sort_values("value", ascending=False).iloc[0]["value"]} cases')

## Best Complex 2
best_complex2 = df_merged[df_merged['protein'] == best_complex1.sort_values("value", ascending=False).index[0]].groupby('ligand').count()

print(f'Best Complex 2 with {best_complex1.sort_values("value", ascending=False).index[0]}: {best_complex2.sort_values("value", ascending=False).index[0]} with {best_complex2.sort_values("value", ascending=False).iloc[0]["value"]} cases')


## Clear any older results
if os.path.exists('./filtered_data/'):
    shutil.rmtree('./filtered_data/')

## Create a new folder for filtered data
os.makedirs('./filtered_data/', exist_ok=True)

## Copy the pdb to a new folder
for item in df_merged.itertuples():
    protein = item.protein
    ligand = item.ligand    
    model = item.model
    pose = item.pose
    description = item.description
    
    source_path = f'./data/{protein}/{protein}_{ligand}_{model}/{protein}_{ligand}_{model}_prep_{pose}.pdb'
    #verify if the source file exists
    if not os.path.exists(source_path):
        print(f'File not found: {source_path}')
        continue

    destination_path = f'./filtered_data/{description}.pdb'
    
    shutil.copy(source_path, destination_path)

print(df_merged.sort_values('value', ascending=True).head(10))