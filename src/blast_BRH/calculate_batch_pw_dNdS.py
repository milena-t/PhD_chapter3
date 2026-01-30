"""
Make batches of dNdS to calculate at once. splits it into several jobs of orthologs
"""

import os


def make_nested_lists(dir_path, include_list = []):
    """
    make dictionary of by-species nested lists for the pairwise fasta files. if include_list is specified then only use 
    pairs that have members in that list
    """
    subdirs = [f"{dir_path}{f}/" for f in os.listdir(dir_path) if os.path.isdir(os.path.join(dir_path, f))]
    if include_list != []:
        print(f"include only species: {include_list}")

    fastas = []
    for pairdir in subdirs:
        print(f"making nested fasta lists in: {pairdir}")
        species_present = True
        if include_list != []:
            species_present = False
            for species_incl in include_list:
                print(f"\t - testing {species_incl} in {pairdir}")
                if species_incl in pairdir:
                    species_present = True
                    break
        if species_present:
            fastas.extend([f"{pairdir}{f}" for f in os.listdir(pairdir) if os.path.isfile(os.path.join(pairdir, f))])
    
    nested_dict = {}
    print(f"{len(fastas)} orthologs: {fastas[:3]} ... ")
        
    for fasta_name in fastas:
        fasta_components = fasta_name.split("/")[-1].split("_")
        species1 = f"{fasta_components[0]}_{fasta_components[1]}"
        species2 = f"{fasta_components[2]}_{fasta_components[3]}"
        if species1 in nested_dict:
            if species2 in nested_dict[species1]:
                nested_dict[species1][species2].append(fasta_name)
            else:
                nested_dict[species1][species2] = [fasta_name]
        else:
            nested_dict[species1] = {}
            nested_dict[species1][species2] = [fasta_name]

    return nested_dict


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

if __name__ == "__main__":

    # list files in dir
    datadir = f"/Users/miltr339/work/pairwise_blast_chapter_2_3/brh_tables/"
    bash_dir = f"/Users/miltr339/work/PhD_code/PhD_chapter3/bash"
    dNdS_exec = f"bash"
    pelle = False

    if not os.path.isdir(datadir):
        datadir = f"/proj/naiss2023-6-65/Milena/chapter3/dNdS_calculations/"
        bash_dir = f"/proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/bash"
        dNdS_exec = f"sbatch"
        pelle = True
    
    #############################
    ###############
    #######
    chr_type = "X"
    #######
    analysis = "LRT" ## codeml site model M1a and M2a with likelihood ratio test
    # analysis = "dNdS" ## codeml branch model
    #######
    ###############
    #############################
    ## done: A_dNdS, X_dNdS, X_LRT
    
    script_path = f"{bash_dir}/run_batch_{analysis}_{chr_type}.sh"

    ### original
    X_path = f"{datadir}brh_sequences_{chr_type}/"
    if analysis=="dNdS":
        outdir = f"{datadir}brh_results_{chr_type}_branch_model/"
    elif analysis=="LRT":
        outdir = f"{datadir}brh_results_{chr_type}/"
    else:
        raise RuntimeError(f"invalid 'analysis' option! you have {analysis} but it needs to be 'dNdS' or 'LRT' ")

    ### second run for the ancestral X-syntenic chromosome in Dcar as X
    # X_path = f"{datadir}brh_sequences_{chr_type}_Dcar_X_syntenic/"
    # outdir = f"{datadir}brh_results_{chr_type}_Dcar_X_syntenic/"

    os.chdir(outdir)
    x_paths_nested_dict = make_nested_lists(X_path, include_list=["C_magnifica","C_septempunctata","T_castaneum","T_freemani"])
    sp1_list=list(x_paths_nested_dict.keys())

    ########
    ## for testing purposes only do a few files of each pair
    ## if 0 then it takes all files, for a normal run
    num_files = 0
    ########

    for species1, subdict in x_paths_nested_dict.items():
        print(f"species1 : {species1}")
        if species1 == "C_magnifica"  or species1 == "C_septempunctata":
            print(f"run {species1}")
            pass
        elif species1 == "T_castaneum"  or species1 == "T_freemani":
            print(f"run {species1}")
            pass
        else:
            print(f"skip {species1}")
            continue

        for species2, fasta_list in subdict.items():
            if species1 == species2:
                continue
            dirname = f"{species1}_{species2}_pairwise_dNdS"
            print(f">>>> {dirname}\n")
            if not os.path.isdir(dirname):
                os.mkdir(dirname)
            os.chdir(dirname)
            wd = os.getcwd()

            if chr_type == "X" or num_files!=0:
                if num_files == 0:
                    fasta_string = " ".join([f"{fasta}" for fasta in fasta_list])
                else:
                    fasta_string = " ".join([f"{fasta}" for fasta in fasta_list[:num_files]])
                print(f"\t--> working in {wd}")
                
                if pelle:
                    jobname = f"dNdS_{chr_type}-linked_{species1}_{species2}"
                    dNdS_command = f"{dNdS_exec} -J {jobname} -o {jobname}.out {script_path} {fasta_string}"
                else:
                    dNdS_command = f"{dNdS_exec} {script_path} {fasta_string}"

                os.system(dNdS_command)
                print(f"{dNdS_command[:1000]} ...")
            
            if chr_type == "A" and num_files == 0:
                ### split into several batches because it doesn't run with too many fasta files
                fasta_overall_list = fasta_list
                i = 0
                for fasta_list in chunks(fasta_overall_list, 500):
                    fasta_string = " ".join([f"{fasta}" for fasta in fasta_list])
                    print(f"\t--> working in {wd}")
                    
                    if pelle:
                        jobname = f"{i}_dNdS_{chr_type}-linked_{species1}_{species2}"
                        dNdS_command = f"{dNdS_exec} -J {jobname} -o {jobname}.out {script_path} {fasta_string}"
                    else:
                        dNdS_command = f"{dNdS_exec} {script_path} {fasta_string}"

                    os.system(dNdS_command)
                    os.system("sleep .5") # wait a little bit after each command so that the job manager can keep up
                    print(f"{dNdS_command[:1000]} ...")
                    i +=1

            os.chdir("..")
            wd = os.getcwd()
            # print(f"\n working in {wd}")
            