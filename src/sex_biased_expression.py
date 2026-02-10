""" 
correlate the dNdS analysis with sex-biased genes in C. septempunctata, T. castaneum and C. maculatus
"""

import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# rsync -azP milenatr@pelle.uppmax.uu.se:/proj/naiss2023-6-65/Milena/chapter3/RNAseq/Coccinella/gene_counts/gene_counts_standard.txt /Users/miltr339/work/chapter3/DE_analysis/Csep_gene_counts.txt
# sed 's|/proj/naiss2023-6-65/Milena/chapter3/RNAseq/Coccinella/star_mapping/picard_marked_indexed/||g' Csep_gene_counts.txt | sed 's|_marked_duplicates.bam||g' > Csep_gene_counts_short_headers.txt

# rsync -azP milenatr@pelle.uppmax.uu.se:/proj/naiss2023-6-65/Milena/chapter3/RNAseq/Tribolium/gene_counts/gene_counts_standard.txt /Users/miltr339/work/chapter3/DE_analysis/Tcas_gene_counts.txt
# sed 's|/proj/naiss2023-6-65/Milena/chapter3/RNAseq/Tribolium/star_mapping/picard_marked_indexed/||g' Tcas_gene_counts.txt | sed 's|_marked_duplicates.bam||g' > Tcas_gene_counts_short_headers.txt

# rsync -azP milenatr@pelle.uppmax.uu.se:/proj/naiss2023-6-65/Milena/chapter3/RNAseq/maculatus/gene_counts/gene_counts_standard.txt /Users/miltr339/work/chapter3/DE_analysis/Cmac_gene_counts.txt
# sed 's|/proj/naiss2023-6-65/Milena/chapter3/RNAseq/maculatus/star_mapping/picard_marked_indexed/||g' Cmac_gene_counts.txt | sed 's|_marked_duplicates.bam||g' > Cmac_gene_counts_short_headers.txt

## R analysis with DESeq2 according to sebastian's github
# -Aggregate to transcript level by summing exon counts,
# -filtered on ≥3 mean counts per sample in each sex,
# -DESeq2 for DE analysis based on male vs female,
# -used vst for count normalization with variance stabilization. (Variance stabilizing transformation)

def counts_paths(username="miltr339"):
    counts_paths = {
        "Cmac" : f"/Users/{username}/work/chapter3/DE_analysis/Cmac_gene_counts_short_headers.txt",
        "Tcas" : f"/Users/{username}/work/chapter3/DE_analysis/Tcas_gene_counts_short_headers.txt",
        "Csep" : f"/Users/{username}/work/chapter3/DE_analysis/Csep_gene_counts_short_headers.txt"
    }
    VST_paths = {
        "Cmac" : f"/Users/{username}/work/chapter3/DE_analysis/Cmac_gene_counts_short_headers_vst.txt",
        "Tcas" : f"/Users/{username}/work/chapter3/DE_analysis/Tcas_gene_counts_short_headers_vst.txt",
        "Csep" : f"/Users/{username}/work/chapter3/DE_analysis/Csep_gene_counts_short_headers_vst.txt"
    }
    return counts_paths,VST_paths

def metadata_paths(username="miltr339"):
    metadata = {
        "Cmac" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/metadata/Cmac_full_SRR_list.csv",
        "Tcas" :  f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/metadata/Tcas_full_SRR_list.csv",
        "Csep" :  f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/metadata/Csep_full_SRR_list.csv"
    }
    return metadata

def lookup_tables(username="miltr339"):
    tables={
        "X" : f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/ortholog_IDs_X_transcript_IDs_association.txt",
        "A" : f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/ortholog_IDs_A_transcript_IDs_association.txt",
    }
    return tables

class OrthologTranscripts:
    def __init__(self,species1:str, species2:str, trID1:str,trID2:str, ortholog_num:str) -> None:
        self.species1=species1
        self.species2=species2
        self.trID1=trID1
        self.trID2=trID2
        self.ortholog_num=ortholog_num
        ### TODO not sure if this is the best way to do it. I essentially want a table of Cmac LFC, Cmac/Cchi dNdS and po_sel
    
    
def get_species_pairs_from_lookuptable(table_path:str):
    """
    get all species pairs that are in the orthologID_transcript lookup table
    """
    with open(table_path, "r") as table_file:
        lines = table_file.readlines()
        lines_parsed = [line.strip().split(":")[0] for line in lines]
        pairs = lines_parsed
        for i , ortholog_ID_full in enumerate(lines_parsed):
            splits = ortholog_ID_full.split("_")
            pairs[i] = "_".join(splits[:4])
        return sorted(list(set(pairs))) 
            


def read_orthologID_lookup_dict(table_path:str):
    """
    The lookup tables have the format species1_species2_orthologID:transcript1,transcript2
    read them into a dictionary with {
        species_pair : {
            orthologID : OrthologTranscripts
        }
    }
    """
    out_dict = {pair : {} for pair in get_species_pairs_from_lookuptable(table_path)}
    
    with open(table_path, "r") as table_file:
        lines = table_file.readlines()
        for line in tqdm(lines):
            line = line.strip()
            ortholog,transcrpts = line.split(":")
            t1,t2 = transcrpts.split(",")
            g1,s1,g2,s2 = ortholog.split("_")[:4]
            ol_num = ortholog.split("_")[-1]
            association = OrthologTranscripts(species1=f"{g1}_{s1}", species2=f"{g2}_{s2}", trID1=t1, trID2=t2, ortholog_num=ol_num)
            pair = f"{g1}_{s1}_{g2}_{s2}"
            out_dict[pair][]
            break
    
    return out_dict



colors_dict = {
    "sex" : {
        "male" : "#4570B0",
        "female" : "#C32B09",
    },
    "organ" : {
        "Head+thorax" : "#8D94BA", # lavender grey
        "abdomen" : "#54457F" , # dusty grape
        "body" : "#955E42", #toffee brown
        "head" : "#2E933C", # light green
        "mouthparts" : "#297045", # medium green
        "antenna" : "#204E4A", # dark green
        "legs" : "#87677B", # dusty lavender
    }
}

points_dict = {
    "sex" : {
        "male" : "v", # down triangle
        "female" : "o", # default circle
    },
    "organ" : {
        "Head+thorax" : "X", # bold X
        "abdomen" : "D" , # diamond shape
    }
}

def plot_PCA_vst_counts(counts_path:str, metadata_path:str, plot_path:str="", colors_dict=colors_dict, species_name=""):
    """
    make a PCA of the DE counts after they are vst-normalized
    """
    if plot_path=="":
        plot_path=counts_path.replace(".txt", "_PCA.png")
    
    norm_counts = pd.read_csv(counts_path, sep="\t", comment="#", index_col=0)
    # read metadata and sort/order according to counts
    metadata = pd.read_csv(metadata_path)
    metadata = metadata.set_index('Run')
    metadata = metadata.reindex(norm_counts.columns)
    print(metadata)
    # categories = list(metadata.columns)

    # create transpose
    vst_t = norm_counts.T
    vst_t.columns = vst_t.columns.astype(str) # convert all column names to string
    pca = PCA(n_components=2)
    pca_scores = pca.fit_transform(vst_t)

    pca_df = metadata.copy()
    pca_df['PC1'] = pca_scores[:,0]
    pca_df['PC2'] = pca_scores[:,1]
    # print(pca_df)

    # Explained variance
    pc1_var = pca.explained_variance_ratio_[0]
    pc2_var = pca.explained_variance_ratio_[1]

    ### plotting
    fig, ax = plt.subplots(1,1, figsize=(15, 10)) # for more than three rows

    fs = 25
    ps = fs*15 # point size

    for i, row in pca_df.iterrows():
        ax.scatter(row["PC1"], row["PC2"], marker=points_dict["sex"][row['sex']], s=ps, color=colors_dict["organ"][row['organ']])
    # ax.scatter(pca_df["PC1"], pca_df["PC2"], color=color_by_sex_vec, s=100, marker=marker_by_organ_vec)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var)", fontsize = fs)
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% var)", fontsize = fs)
    ax.tick_params(axis='x', labelsize=fs)
    ax.tick_params(axis='y', labelsize=fs)
    ylim = ax.get_ylim()
    ax.set_ylim(ylim)
    xlim = ax.get_xlim()
    ax.set_xlim(xlim)

    ## make legend
    sex_labels = { sex : points_dict['sex'][sex] for sex in pca_df['sex']}
    org_labels = { org : colors_dict['organ'][org] for org in pca_df['organ']}
    yleg = ylim[0]-1e6
    xleg = xlim[0]-1e6

    for key,value in sex_labels.items():
        ax.scatter(xleg,yleg,marker=value,s=ps,label=key, color='black')
    for key,value in org_labels.items():
        ax.scatter(xleg,yleg,color=value,s=ps,label=key, marker="s")
    ax.legend(fontsize=fs)
    plt.suptitle(f"{species_name} differential expression \nbased on {len(metadata)} samples", fontsize=fs*1.25)
    
    # Adjust layout to prevent overlap
    plt.tight_layout(rect=[0, 0.05, 1, 1])

    plt.savefig(plot_path, dpi = 300, transparent = False)
    print(f"plot saved in current working directory as: {plot_path}")


if __name__ == "__main__":

    species_list = ["Cmac","Tcas","Csep"] 
    species_list = ["Cmac"]
    species_names = {
        "Cmac" : "C. maculatus",
        "Tcas" : "T. castaneum",
        "Csep" : "C. septempunctata",
    }
    username = "miltr339"
    ortholog_ID_lookup_tables = lookup_tables(username=username)

    lookup_dict = read_orthologID_lookup_dict(ortholog_ID_lookup_tables["X"])
    print(lookup_dict)

    ### plot PCA
    if False:
        counts_paths_dict, vst_paths_dict = counts_paths(username=username)
        metadata_paths_dict = metadata_paths(username=username)
        for species in species_list:
            print(f"\n//////////////////////// {species} ////////////////////////")
            print(f"\t * {vst_paths_dict[species]}")
            print(f"\t * {metadata_paths_dict[species]}")
            plot_PCA_vst_counts(
                counts_path=vst_paths_dict[species], 
                metadata_path=metadata_paths_dict[species], 
                plot_path=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/{species}_vst_counts_PCA.png",
                species_name=species_names[species])
