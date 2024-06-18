import os
import numpy as np

file_names = os.listdir()
file_names = list(filter(lambda x: x.endswith(".tsv"), file_names))

elements = set()

for file_name in file_names:
    
    current_element = file_name.split(".")[1]    
    elements.add(current_element)
    
for element in elements:

    # MutaBind2 read
    
    with open("MB2." + element + ".tsv", "r") as f:
        mb2_muts = f.read()
        
    mb2_muts = mb2_muts.split("\n")
    
    mut_dict = {line.split("\t")[0]: [float(line.split("\t")[1]), ] for line in mb2_muts}

    # mCSM2 read
        
    with open("mCSM." + element + ".tsv", "r") as f:
        mcsm_muts = f.read()
        
    mcsm_muts = mcsm_muts.split("\n")
    
    for mcsm_mut in mcsm_muts:
    
        current_mut = mcsm_mut.split("\t")[0]
        
        if current_mut in mut_dict:
        
            current_ddg = -float(mcsm_mut.split("\t")[1])  # NEGATIVE SIGN!! mCSM consideres ddG backwards...
            mut_dict[current_mut].append(current_ddg)

    # SAAMBE-3D read
            
    with open("SAAMBE." + element + ".tsv", "r") as f:
        saambe_muts = f.read()
        
    saambe_muts = saambe_muts.split("\n")
        
    for saambe_mut in saambe_muts:
    
        current_mut = saambe_mut.split("\t")[0]
        
        if current_mut in mut_dict:
        
            current_ddg = float(saambe_mut.split("\t")[1])
            mut_dict[current_mut].append(current_ddg)

    # Deleting mutations with incomplete lists
            
    for key in list(mut_dict.keys()):
        
        if len(mut_dict[key]) < 3:
        
            print(f"Removing {key} from {element}: {mut_dict[key]}")
            del mut_dict[key]

    # Calculating averages and standard deviations
            
    mb2_values = [mut_dict[key][0] for key in mut_dict]
    mcsm_values = [mut_dict[key][1] for key in mut_dict]
    saambe_values = [mut_dict[key][2] for key in mut_dict]
    
    mb2_avg = np.mean(mb2_values)
    mb2_std = np.std(mb2_values)
    
    mcsm_avg = np.mean(mcsm_values)
    mcsm_std = np.std(mcsm_values)
    
    saambe_avg = np.mean(saambe_values)
    saambe_std = np.std(saambe_values)

    # Calculating Z-scores and consensus Z-scores
    
    out_str = ""  # The merged file.
    for key in mut_dict:
    
        mb2_zscore = (mut_dict[key][0] - mb2_avg) / mb2_std
        mcsm_zscore = (mut_dict[key][1] - mcsm_avg) / mcsm_std
        saambe_zscore = (mut_dict[key][2] - saambe_avg) / saambe_std        
        new_value = (mb2_zscore + mcsm_zscore + saambe_zscore) / 3
        
        mut_dict[key] = new_value
        
        out_str += f"{key}\t{new_value}\n"
        
    out_str = out_str[:-1]  # Removing trailing newline.
                        
    with open("merged." + element + ".tsv", "w+") as f:
        f.write(out_str)

    # Collecting the 19 consensus Z-scores for each position
        
    avg_zscore = dict()
    for mut_id in mut_dict:
    
        if mut_id[:-1] in avg_zscore:  # position_id = mut_id[:-1]
            avg_zscore[mut_id[:-1]].append(mut_dict[mut_id])
        else:
            avg_zscore[mut_id[:-1]] = [mut_dict[mut_id], ]
            
    # Creating the .pml script
    
    out_str = "alter all, b = 1000;\n"
    
    zscore_min = float("inf")
    zscore_max = float("-inf")
    for key in avg_zscore:
    
        current_avg = np.mean(avg_zscore[key])
        
        if current_avg > zscore_max:
            zscore_max = current_avg
        if current_avg < zscore_min:
            zscore_min = current_avg
        
        out_str += f"select chain {key[1]} and resi {key[2:]}; "
        out_str += f"alter (sele), b = {current_avg};\n"
        
    out_str += f"spectrum b, blue_white_yellow, minimum={zscore_min}, maximum={zscore_max};\n"
    out_str += f"color red, b>999;"
    
    with open("script." + element + ".pml", "w+") as f:
        f.write(out_str)
    