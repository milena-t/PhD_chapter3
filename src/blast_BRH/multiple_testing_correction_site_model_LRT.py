"""
since the individual LRT tests for the site model are run as part of calculate_pairwise_dNdS, they do not take multiple testing error into account.
Therefore, I will re-summarize these results by extracting the p-values of all genes from the log files and assessing significance properly.
"""

from scipy.stats import false_discovery_control

# rsync -azP milenatr@pelle.uppmax.uu.se:/proj/naiss2023-6-65/Milena/chapter3/dNdS_calculations/brh_results_A/site_model_A_all_log_files.out site_model_A_all_log_files.out
# rsync -azP milenatr@pelle.uppmax.uu.se:/proj/naiss2023-6-65/Milena/chapter3/dNdS_calculations/brh_results_X/site_model_X_all_log_files.out site_model_X_all_log_files.out

def get_log_file_paths(username="miltr339"):
    out_dict = {
        "A" : f"/Users/{username}/work/chapter3/dNdS/site_model_A_all_log_files.out",
        "X" : f"/Users/{username}/work/chapter3/dNdS/site_model_X_all_log_files.out",
    }
    return out_dict


def get_pvals_dict_from_log(log_file_path:str):
    """
    get a dictionary that is {ortholog : pval}
    """
    out_dict = {}
    curr_OG_pval=False
    ortholog = ""
    with open(log_file_path, "r") as log_file:
        for log_line in log_file.readlines():
            line = log_line.strip()
            if line == "":
                continue
            if "===================" in line:
                curr_ortholog = line.split(" ")[1]
                if curr_ortholog != ortholog:
                    ortholog=curr_ortholog
                    curr_OG_pval=True
            if "p-value: " in line and curr_OG_pval:
                pval = line.split(" ")[1]
                out_dict[ortholog] = float(pval)
                curr_OG_pval=False

    return out_dict

def multiple_testing_correction(pvals_dict):
    """
    return a dictionary with corrected p-values
    """
    orthologs,pvals = zip(*pvals_dict.items())
    bp_pvals = false_discovery_control(pvals, method="bh")
    result = dict(zip(orthologs, bp_pvals))
    return result


def modify_original_file(site_classes_path, bh_pos_sel_orthologs:list, outfile):
    """
    The original file is sorted into positively or negatively selected according to the proportion of positively selected sites in the site classes table.
    this function makes a new file where p and w where positively selected site classes are deleted
    """
    summary_stats = {
        "no_pos_sel" : 0,
        "simple_pos_sel_not_after_BH" : 0,
        "BH_pos_sel" : 0
    }
    with open(site_classes_path, "r") as site_classes_uncorrected, open(outfile, "w") as outfile_bh:
        for uncorrected_line in site_classes_uncorrected.readlines():
            try:
                filepath_str, site_classes_string = uncorrected_line.strip().split(" : ")
            except:
                raise RuntimeError(f"line cannot be parsed! \n{uncorrected_line}")
            
            filepath = filepath_str.split("/")
            ortholog_name,site_class_name = filepath[-2:]
                
            if ortholog_name in bh_pos_sel_orthologs:
                # keep positive selection and write line to new outfile
                corrected_line = uncorrected_line
                summary_stats["BH_pos_sel"]+=1
            else: 
                # either already negative or no positive after correction
                line = uncorrected_line.strip()
                ortholog_path, line_list = line.split(" : ")
                
                if "no_file" in line_list:
                    continue
                
                p,w = line_list.split(";")
                p_list = p.split()
                w_list = w.split()
                assert len(p_list) == len(w_list)

                if len(p_list) == 3:
                    corrected_line = uncorrected_line # no positive selection in non-bh corrected file either, no change necessary
                    summary_stats["no_pos_sel"]+=1
                else:
                    p_fix = " ".join(p_list[:3])
                    w_fix = " ".join(w_list[:3])
                    # print(f"\t\t{ortholog_name}")
                    # print(f"\t\t{p} ---> {p_fix}")
                    # print(f"\t\t{w} ---> {w_fix}")
                    corrected_line = f"{ortholog_path} : {p_fix};{w_fix}\n"
                    summary_stats["simple_pos_sel_not_after_BH"]+=1
            outfile_bh.write(corrected_line)

    print(f"outfile written to: {outfile}.")
    print(f"\t * {summary_stats['no_pos_sel']} genes not under positive selection")
    print(f"\t * {summary_stats['BH_pos_sel']} are under positive selection after BH correction")
    print(f"\t * {summary_stats['simple_pos_sel_not_after_BH']} were positively selected according to simple p-value but are not any more after BH correction")


def write_output(pvals_dict, p_sig = 0.05, outfile_path=""):
    """ 
    write an output file with the ortholog names and TRUE/FALSE according to passed LRT and positive selection (TRUE) or failed (FALSE)
    """
    pos_sel_orthologs = []
    with open(outfile_path, "w") as outfile:
        for ortholog,pval in pvals_dict.items():
            pos_sel = "FALSE"
            if pval<p_sig:
                pos_sel="TRUE"
                pos_sel_orthologs.append(ortholog.replace("/",""))
            outfile.write(f"{ortholog}\t{pos_sel}\n")

    print(f"outfile written to: {outfile_path}")
    return pos_sel_orthologs

if __name__ == "__main__":
    username="miltr339"
    log_paths_dict = get_log_file_paths(username=username)
    p_sig = 0.05

    files = {
        "A_list" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/site_classes_A_BH_corrected.txt",
        "X_list" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/site_classes_X_BH_corrected.txt",
        "A_LRT" : f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/brh_results_A/site_classes_summary_A-linked.txt",
        "X_LRT" : f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/brh_results_X/site_classes_summary_X-linked.txt",
        "A_LRT_BH_corr" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/site_classes_summary_A_BH_corrected.txt",
        "X_LRT_BH_corr" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/site_classes_summary_X_BH_corrected.txt",
    }

    chromosomes = ["X", "A"]
    # chromosomes = ["X"]
    
    for chromosome in chromosomes:
        print(f"////////////////////////// {chromosome} //////////////////////////")
        pvals_dict = get_pvals_dict_from_log(log_file_path=log_paths_dict[chromosome])
        bh_pvals_dict = multiple_testing_correction(pvals_dict)
        pos_sel_orthologs = write_output(bh_pvals_dict, p_sig=p_sig, outfile_path=files[f"{chromosome}_list"])
        modify_original_file(site_classes_path=files[f"{chromosome}_LRT"], bh_pos_sel_orthologs=pos_sel_orthologs, outfile=files[f"{chromosome}_LRT_BH_corr"])
