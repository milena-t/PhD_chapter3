"""
Make batches of dNdS to calculate at once. splits it into several jobs of orthologs
"""

import os


def make_nested_lists(dir_path):
    """
    make dictionary of by-species nested lists for the pairwise fasta files
    """
    subdirs = [f"{dir_path}{f}/" for f in os.listdir(dir_path) if os.path.isdir(os.path.join(dir_path, f))]

    fastas = []
    for pairdir in subdirs:
        print(pairdir)
        fastas.extend([f"{pairdir}{f}" for f in os.listdir(pairdir) if os.path.isfile(os.path.join(pairdir, f))])
    
    nested_dict = {}
    print(f"{len(fastas)} orthologs: {fastas[:10]}")
        
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
    script_path = f"{bash_dir}/run_batch_dNdS.sh"
    pelle = False

    if not os.path.isdir(datadir):
        datadir = f"/proj/naiss2023-6-65/Milena/chapter3/dNdS_calculations/"
        bash_dir = f"/proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/bash"
        dNdS_exec = f"sbatch"
        script_path = f"{bash_dir}/run_batch_dNdS.sh"
        pelle = True
    
    #######
    chr_type = "X"
    #######

    ### original
    X_path = f"{datadir}brh_sequences_{chr_type}/"
    outdir = f"{datadir}brh_results_{chr_type}/"

    ### second run for the ancestral X-syntenic chromosome in Dcar as X
    # X_path = f"{datadir}brh_sequences_{chr_type}_Dcar_X_syntenic/"
    # outdir = f"{datadir}brh_results_{chr_type}_Dcar_X_syntenic/"

    os.chdir(outdir)
    x_paths_nested_dict = make_nested_lists(X_path)
    sp1_list=list(x_paths_nested_dict.keys())

    ########
    ## for testing purposes only do a few files of each pair
    ## if 0 then it takes all files, for a normal run
    num_files = 0
    ########

    for species1, subdict in x_paths_nested_dict.items():
        # if species1 != "D_carinulata":
        #     continue
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
            