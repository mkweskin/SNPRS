#!/usr/bin/env python3

import os
import sys
import re
import pandas as pd
import argparse
import numpy as np
import uuid

def rebalance_first_level(level_df):
    
    level_df = level_df.copy()
    
    og_cols = [c for c in level_df.columns if c != "Adapted_Request"]
    og_cols = [c for c in og_cols if c != "Bases_Remaining"]
    backup_cols = [c for c in og_cols if c != "Allocated_Bases"]
    merge_cols = [c for c in backup_cols if c != "Requested_Bases"] 
        
    if (level_df["Base_Count"] >= level_df['Adapted_Request']).all():
        level_df['Allocated_Bases'] = level_df['Adapted_Request']
        return level_df[og_cols]
    
    # Check for second level
    second_level_df = level_df[level_df['Group_1'] != "NO_GROUP"].copy()
    has_second_level = not second_level_df.empty
    
    # Sample evenly if no second level
    if not has_second_level:
        backup_df = level_df[backup_cols].copy()
        level_df = (
            level_df.drop('Requested_Bases', axis=1)
                    .rename(columns={'Adapted_Request': 'Requested_Bases'})
        )
        
        even_level = (get_even_sampling(level_df)
                    .drop(columns=["Requested_Bases"],axis=1)
                    .merge(backup_df,on=merge_cols))
        
        return even_level[og_cols]
    
    # If there are second levels, is there an outgroup?
    outgroup_df = level_df[level_df['Group_1'] == "NO_GROUP"].copy()
    has_outgroup = not outgroup_df.empty
    
    # Are there enough reads to satisfy the second levels?
    second_level_requested = second_level_df['Requested_Bases'].sum()
    second_level_available = second_level_df['Base_Count'].sum()
    second_level_deficit = second_level_requested - second_level_available
    second_level_satisfied = second_level_available >= second_level_requested        
    
    if not has_outgroup:
        return_df = pd.DataFrame()
    else:
        backup_df = outgroup_df[backup_cols].copy()
        if second_level_satisfied:
            
            temp_outgroup_df = (
                outgroup_df.drop('Requested_Bases', axis=1)
                        .rename(columns={'Adapted_Request': 'Requested_Bases'})
            )
            temp_outgroup_df = (get_even_sampling(temp_outgroup_df)
                        .drop(columns=["Requested_Bases"])
                        .merge(backup_df,on=merge_cols,how="left")
            )
            
        else:
            outgroup_requested = outgroup_df['Adapted_Request'].sum()
            outgroup_count = outgroup_df.shape[0]
            
            new_outgroup_request = outgroup_requested + second_level_deficit
            temp_outgroup_df = outgroup_df.drop(['Adapted_Request','Requested_Bases'], axis=1)
            temp_outgroup_df['Requested_Bases'] = int(new_outgroup_request/outgroup_count)
            temp_outgroup_df = (get_even_sampling(temp_outgroup_df)
                        .drop(columns=["Requested_Bases"])
                        .merge(backup_df,on=merge_cols,how="left")
            )
            
        return_df = temp_outgroup_df[og_cols]
    
    # Handle second levels
    second_level_groups = second_level_df['Group_1'].unique()
    second_level_count = len(second_level_groups)
    bases_per_second_level = int(second_level_requested/second_level_count)
    second_level_deficit = 0
    
    for second_level in second_level_groups:
    
        second_level_df2 = second_level_df[second_level_df['Group_1'] == second_level].copy()
        second_level_bases = second_level_df2['Base_Count'].sum()
        
        if second_level_bases <= bases_per_second_level:
            second_level_deficit+= bases_per_second_level - second_level_bases
            second_level_df2['Allocated_Bases'] = second_level_df2['Base_Count']
            return_df = pd.concat([return_df,second_level_df2[og_cols]])
            second_level_df = second_level_df[second_level_df['Group_1'] != second_level].copy()
        
    while second_level_deficit > (0.01*bases_per_second_level) and second_level_df.shape[0]>0:
        
        remaining_second_levels = second_level_df['Group_1'].unique()
        remaining_count = len(remaining_second_levels)
        og_bases_per_second_level = bases_per_second_level
        bases_per_second_level = int(bases_per_second_level + int(second_level_deficit/remaining_count))
        
        for second_level in remaining_second_levels:
        
            second_level_df2 = second_level_df[second_level_df['Group_1'] == second_level].copy()
            second_level_bases = second_level_df2['Base_Count'].sum()
            
            if second_level_bases <= bases_per_second_level:
                second_level_deficit-= second_level_bases - og_bases_per_second_level
                second_level_df2['Allocated_Bases'] = second_level_df2['Base_Count']
                return_df = pd.concat([return_df,second_level_df2[og_cols]])
                second_level_df = second_level_df[second_level_df['Group_1'] != second_level].copy()
        
    if not second_level_df.empty:
        remaining_second_levels = second_level_df['Group_1'].unique()
        remaining_count = len(remaining_second_levels)
        
        for remaining in remaining_second_levels:
            second_level_df2 = second_level_df[second_level_df['Group_1']==remaining].copy()
            backup_df = second_level_df2[backup_cols].copy()

            second_level_requested = second_level_df2['Adapted_Request'].sum()
            second_level_count = second_level_df2.shape[0]
            new_request = second_level_requested + second_level_deficit
            second_level_df2 = second_level_df2.drop(['Adapted_Request','Requested_Bases'], axis=1)
            second_level_df2['Requested_Bases'] = int(new_request/second_level_count)
            second_level_df2 = (get_even_sampling(second_level_df2)
                        .drop(columns=["Requested_Bases"])
                        .merge(backup_df,on=merge_cols,how="left")
            )
            return_df = pd.concat([return_df,second_level_df2[og_cols]])
    
    return return_df
          
def get_even_sampling(request_df):
    
    total_available = request_df["Base_Count"].sum()
    total_requested = request_df["Requested_Bases"].sum()
    total_count = request_df.shape[0]
    requested_per_count = int(total_requested/total_count)
    
    # If there's not enough, sample everything
    if total_available <= total_requested:
        request_df["Allocated_Bases"] = request_df["Base_Count"]
        return request_df
    
    # Grab exhausted as return        
    return_df = request_df[request_df['Bases_Remaining'] == 0].copy()        
    high_datasets = request_df[request_df['Bases_Remaining'] > 0].copy()        
    total_shortage = return_df['Requested_Bases'].sum() - return_df['Base_Count'].sum()
            
    while total_shortage > (0.01*requested_per_count) and high_datasets.shape[0] > 0:
                
        high_count = high_datasets.shape[0]
        even_shortage = int(total_shortage/high_count)
        
        # If all the samples can upsample, do it and return
        if (high_datasets['Bases_Remaining'] >= even_shortage).all():
            
            high_datasets['Allocated_Bases'] = high_datasets['Allocated_Bases'] + even_shortage
            high_datasets['Bases_Remaining'] = high_datasets['Bases_Remaining'] - even_shortage
            
            return_df = pd.concat([return_df,high_datasets])
            return return_df

        # If some samples are exhausted while upsampling, keep going
        exhausted = high_datasets[high_datasets['Bases_Remaining'] <= even_shortage].copy()
        donated = (exhausted['Base_Count'] - exhausted['Allocated_Bases']).sum()
        exhausted['Allocated_Bases'] = exhausted['Base_Count']            
        total_shortage = total_shortage - donated

        # Append exhausted and return high datasets to loop again
        return_df = pd.concat([return_df,exhausted])            
        high_datasets = high_datasets[high_datasets['Bases_Remaining'] > even_shortage].copy()      
    
    return pd.concat([return_df,high_datasets])
    
def balance_bases(request_df):
    request_df = request_df.copy()
    request_df["Allocated_Bases"] = request_df["Requested_Bases"].clip(upper=request_df["Base_Count"])
    request_df["Bases_Remaining"] = request_df["Base_Count"] - request_df["Allocated_Bases"]

    # Sample everything
    if (request_df["Bases_Remaining"] == 0).all():
        return request_df.drop('Bases_Remaining',axis=1, errors="ignore")
    
    # Everything is good
    if (request_df["Bases_Remaining"] > 0).all():
        return request_df.drop('Bases_Remaining',axis=1, errors="ignore")
    
    # Handle shortages
    
    # Even sampling
    if (request_df['Final_Group'] == "Even_Sampling").all():
        even_df = get_even_sampling(request_df).drop('Bases_Remaining',axis=1, errors="ignore")
        return(even_df)
    else:

        # First check outgroup
        outgroup_deficit = 0       
        global_outgroup_df = request_df[request_df['Group_0'] == "NO_GROUP"].copy()
        if not global_outgroup_df.empty:
            if (global_outgroup_df['Bases_Remaining'] > 0).all():
                return_df = global_outgroup_df.copy()
            else:
                outgroup_requested = global_outgroup_df['Requested_Bases'].sum()
                outgroup_available = global_outgroup_df['Base_Count'].sum()
                outgroup_deficit = outgroup_requested - outgroup_available
                if outgroup_available <= outgroup_requested:
                    global_outgroup_df["Allocated_Bases"] = global_outgroup_df["Base_Count"]
                    return_df = global_outgroup_df.drop("Bases_Remaining",axis=1).copy()
                else:
                    return_df = get_even_sampling(global_outgroup_df).drop("Bases_Remaining",axis=1)
        else:
            return_df = pd.DataFrame()
        
        first_level_df = request_df[request_df['Group_0'] != "NO_GROUP"].copy()
        first_level_groups = first_level_df['Group_0'].unique() 
        first_level_count = len(first_level_groups)
        first_level_requested = first_level_df['Requested_Bases'].sum() + outgroup_deficit
        first_level_available = first_level_df['Base_Count'].sum() 
        
        # If there's not enough, add it all
        if first_level_available <= first_level_requested:
            first_level_df["Allocated_Bases"] = first_level_df["Base_Count"]
            return pd.concat([return_df,first_level_df])
        
        # If there is enough, add it proportionally as base counts allow
        bases_per_first_level = int(first_level_requested/first_level_count)
        first_level_deficit = 0
        
        for first_level in first_level_groups:
            
            level_df = first_level_df[first_level_df['Group_0'] == first_level].copy()
            first_level_bases = level_df['Base_Count'].sum()
            
            if first_level_bases <= bases_per_first_level:
                first_level_deficit+= bases_per_first_level - first_level_bases
                level_df['Allocated_Bases'] = level_df['Base_Count']
                return_df = pd.concat([return_df,level_df])
                first_level_df = first_level_df[first_level_df['Group_0'] != first_level].copy()
        
        while first_level_deficit > (0.01*bases_per_first_level) and first_level_df.shape[0]>0:
            
            remaining_first_levels = first_level_df['Group_0'].unique()
            remaining_count = len(remaining_first_levels)
            
            bases_per_first_level = int(bases_per_first_level + int(first_level_deficit/remaining_count))
            
            for first_level in remaining_first_levels:
            
                level_df = first_level_df[first_level_df['Group_0'] == first_level].copy()
                first_level_bases = level_df['Base_Count'].sum()
            
                if first_level_bases <= bases_per_first_level:
                    first_level_deficit-= first_level_bases
                    level_df['Allocated_Bases'] = level_df['Base_Count']
                    return_df = pd.concat([return_df,level_df])
                    first_level_df = first_level_df[first_level_df['Group_0'] != first_level].copy()
        
        # Rebalance any first levels where there are enough bases to not fully sample
        if not first_level_df.empty:
            
            remaining_first_levels = first_level_df['Group_0'].unique()
            remaining_count = len(remaining_first_levels)
            for remaining in remaining_first_levels:
                level_df = first_level_df[first_level_df['Group_0'] == remaining].copy()
                original_requested = level_df['Requested_Bases'].sum()
                new_requested = bases_per_first_level
                level_df['Adapted_Request'] = (
                    level_df['Requested_Bases'] * (new_requested / original_requested)
                ).round().astype(int)
                rebalanced_level_df = rebalance_first_level(level_df)
                return_df = pd.concat([return_df,rebalanced_level_df])

        return return_df.drop('Bases_Remaining',axis=1, errors="ignore")
                
### Main Script ###

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-b', '--base_counts', type=str, default=None, help="Base count data in CSV (Sample_ID,Read_Count,Base_Count,Forward,Reverse)")
parser.add_argument('-g', '--group_file', type=str, default=None, help="Group data in CSV (Sample_ID,Group_0,Group_1,Forward,Reverse)")
parser.add_argument('-m', '--manual_counts', type=str, default=None, help="Group data in CSV (Sample_ID,Group_0,Group_1,Base_Count,Forward,Reverse)")
parser.add_argument('-o', '--output', type=str, required=True, help="Output directory for subset reads")
parser.add_argument('-s', '--genomesize', type=int, required=True, help="Genome size estimate in bp")
parser.add_argument('-c', '--coverage', type=int, default=10, help="Desired coverage")
parser.add_argument('-p', '--out_prop', type=float, default=0.25, help="Sampling proportion for unsorted samples")
args = parser.parse_args()

# Create log file
log_dir = os.path.dirname(args.output)
log_file = f"{log_dir}/Subset_Log.txt"

if args.manual_counts:
    data_string = f"\t- Manual Count Information: {os.path.abspath(args.manual_counts)}\n"
else:
    data_string = f"\t- Base Count Information: {os.path.abspath(args.base_counts)}\n\t- Group Information: {os.path.abspath(args.group_file)}\n"

with open(log_file, 'w') as log:
    log.write("SNPRS Read Subsetter Log\n")
    log.write("-------------------------------------------------------\n\n")
    log.write(f"\t- Genome size estimate: {args.genomesize} bp\n")
    log.write(f"\t- Desired coverage: {args.coverage}X\n")
    log.write(data_string)
    log.write(f"\t- Output directory: {os.path.abspath(args.output)}\n")
    log.write("\n-------------------------------------------------------\n\n")
        
if not (args.manual_counts or (args.base_counts and args.group_file)):
    sys.exit("Must provide --manual_counts or --base_counts/--group_file")

if args.manual_counts:
    path = args.manual_counts
    with open(path) as f:
        first_line = f.readline().strip()
    expected_cols = ['Sample_ID', 'Group_0', 'Group_1', 'Base_Count', 'Forward', 'Reverse']
    if all(col in first_line for col in ['Sample_ID', 'Group_0']):
        df = pd.read_csv(path)
    else:
        sys.exit("Manual counts provided by --manual_counts must have the header: 'Sample_ID', 'Group_0', 'Group_1', 'Base_Count', 'Forward', 'Reverse'")

    base_count_df = df[['Sample_ID', 'Base_Count', 'Forward', 'Reverse']].copy()
    base_count_df['Read_Count'] = np.nan
    base_count_df = base_count_df[['Sample_ID', 'Read_Count', 'Base_Count', 'Forward', 'Reverse']]
    group_df = df[['Sample_ID', 'Group_0', 'Group_1', 'Forward', 'Reverse']].copy()

else:
    base_count_df = pd.read_csv(args.base_counts,names=['Sample_ID', 'Read_Count', 'Base_Count', 'Forward', 'Reverse'])
    group_df = pd.read_csv(args.group_file,names=['Sample_ID', 'Group_0', 'Group_1', 'Forward', 'Reverse'])


for col in ['Forward','Reverse']:
    base_count_df[col] = base_count_df[col].replace({np.nan: None})
    group_df[col] = group_df[col].replace({np.nan: None})

merged_df = pd.merge(group_df, base_count_df[['Sample_ID','Forward','Reverse','Read_Count','Base_Count']],
                     on=['Sample_ID','Forward','Reverse'], how='inner')

sample_count = merged_df.shape[0]
if len(merged_df) != len(group_df):
    raise ValueError("Mismatch after merging...")

total_bases_requested = int(args.genomesize*args.coverage)

global_outgroup_df = merged_df[merged_df['Group_0'] == "NO_GROUP"].copy()
first_level_df = merged_df[merged_df['Group_0'] != "NO_GROUP"].copy()

has_global_outgroup  = not global_outgroup_df.empty
has_first_level = not first_level_df.empty

request_df = pd.DataFrame()

if not has_first_level:
    global_outgroup_df['Requested_Bases'] = int(total_bases_requested/sample_count)
    global_outgroup_df['Final_Group'] = "Even_Sampling"
    request_df = global_outgroup_df.copy()

else:  
    if has_global_outgroup:
        global_outgroup_bases = int(total_bases_requested*args.out_prop)
        global_outgroup_count = global_outgroup_df.shape[0]
        bases_per_global_outgroup = int(global_outgroup_bases/global_outgroup_count)
        global_outgroup_df['Requested_Bases'] = bases_per_global_outgroup
        global_outgroup_df['Final_Group'] = "Global_Outgroup"
        request_df = global_outgroup_df.copy()
        
        ingroup_bases = int(total_bases_requested - global_outgroup_bases)
    else:
        ingroup_bases = total_bases_requested
    
    first_levels = first_level_df['Group_0'].unique()
    first_level_count = len(first_levels)
    bases_per_first_level = int(ingroup_bases/first_level_count)
    
    for group in first_levels:
        
        level_df = first_level_df[first_level_df['Group_0']==group].copy()
        
        outgroup_df = level_df[level_df['Group_1'] == "NO_GROUP"].copy()
        second_level_df = level_df[level_df['Group_1'] != "NO_GROUP"].copy()
        
        has_outgroup  = not outgroup_df.empty
        has_second_level = not second_level_df.empty
        
        if not has_second_level:
            level_count = outgroup_df.shape[0]
            outgroup_df['Requested_Bases'] = int(bases_per_first_level/level_count)
            outgroup_df['Final_Group'] = f"{group}_Even_Sampling"
            request_df = pd.concat([request_df,outgroup_df])
            
        else:
            second_levels = second_level_df['Group_1'].unique()
            second_level_count = len(second_levels)
            
            if has_outgroup:
                level_outgroup_bases = int(bases_per_first_level*args.out_prop)
                level_outgroup_count = outgroup_df.shape[0]
                outgroup_df['Requested_Bases'] = int(level_outgroup_bases/level_outgroup_count)
                outgroup_df['Final_Group'] = f"{group}_Outgroup"
                request_df = pd.concat([request_df,outgroup_df])
                
                level_ingroup_bases = int(bases_per_first_level - level_outgroup_bases)
                                
            else:
                level_ingroup_bases = bases_per_first_level
            
            bases_per_second_level = int(level_ingroup_bases/second_level_count)
            
            for slevel in second_levels:
                slevel_df = second_level_df[second_level_df['Group_1'] == slevel].copy()
                slevel_count = slevel_df.shape[0]
                slevel_df['Requested_Bases'] = int(bases_per_second_level/slevel_count)
                slevel_df['Final_Group'] = f"{group}_{slevel}"
                request_df = pd.concat([request_df,slevel_df])
        
request_df = request_df[["Sample_ID","Group_0","Group_1","Final_Group","Base_Count","Requested_Bases","Forward","Reverse"]].copy()
balanced_df = balance_bases(request_df)

balanced_df["Sample_Type"] = np.where(
balanced_df["Allocated_Bases"] == balanced_df["Base_Count"],
"Link",
"Sample"
)

balanced_df["Subsample_ID"] = [uuid.uuid4().hex[:12] for _ in range(len(balanced_df))]    

summary_group0 = (
    balanced_df.groupby("Group_0")[["Allocated_Bases", "Base_Count"]]
    .sum()
    .reset_index()
)
summary_group0["Level"] = "Group_0"

summary_final = (
    balanced_df.groupby(["Group_0", "Final_Group"])[["Allocated_Bases", "Base_Count"]]
    .sum()
    .reset_index()
)
summary_final["Level"] = "Final_Group"

summary_samples = (
    balanced_df.groupby(["Group_0", "Final_Group", "Sample_ID", "Forward"])[["Allocated_Bases", "Base_Count"]]
    .sum()
    .reset_index()
)
summary_samples["Level"] = "Sample"

# Make all columns align
combined_summary = pd.concat([summary_group0, summary_final, summary_samples], ignore_index=True, sort=False)

# Save to TSV
allocation_file = f"{log_dir}/Sample_Allocation.tsv"
allocation_summary_file = f"{log_dir}/Sample_Allocation_Summary.tsv"

combined_summary.to_csv(allocation_summary_file, sep="\t", index=False)
balanced_df.to_csv(allocation_file,index=False,sep="\t")
balanced_df['Output_Dir'] = os.path.abspath(args.output)
balanced_df[['Sample_ID','Subsample_ID','Forward','Reverse','Allocated_Bases','Sample_Type','Output_Dir']].to_csv(sys.stdout, index=False, header=False)