"""
Plot all the branch-wise dNdS estimates from the branch model of the Bruchini phylogeny for the revision
"""

def filepaths(username="miltr339"):
    data_dir_X = f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/bruchini_seqs_revision_X"
    data_dir_A = f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/bruchini_seqs_revision_A"
    files_dict = {
        "branch_model_X" : f"{data_dir_X}run_branch_model.log",
        "site_model_X" : f"{data_dir_X}run_site_model.log",
        "site_model_beta_X" : f"{data_dir_X}run_site_model_beta.log",
        "branch_model_A" : f"{data_dir_A}run_branch_model.log",
        "site_model_A" : f"{data_dir_A}run_site_model.log",
        "site_model_beta_A" : f"{data_dir_A}run_site_model_beta.log",
    }
    return files_dict




if __name__ == "__main__":

    username = "miltr339"
    codeml_outfiles = filepaths(username=username)

    