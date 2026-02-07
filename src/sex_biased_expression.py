""" 
correlate the dNdS analysis with sex-biased genes in C. septempunctata, T. castaneum and C. maculatus
"""

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

# rsync -azP milenatr@pelle.uppmax.uu.se:/proj/naiss2023-6-65/Milena/chapter3/RNAseq/Coccinella/gene_counts/gene_counts_standard.txt /Users/miltr339/work/chapter3/DE_analysis/Csep_gene_counts.txt
# sed 's|/proj/naiss2023-6-65/Milena/chapter3/RNAseq/Coccinella/star_mapping/picard_marked_indexed/||g' Csep_gene_counts.txt | sed 's|_marked_duplicates.bam||g' > Csep_gene_counts_short_headers.txt

# rsync -azP milenatr@pelle.uppmax.uu.se:/proj/naiss2023-6-65/Milena/chapter3/RNAseq/Tribolium/gene_counts/gene_counts_standard.txt /Users/miltr339/work/chapter3/DE_analysis/Tcas_gene_counts.txt
# sed 's|/proj/naiss2023-6-65/Milena/chapter3/RNAseq/Tribolium/star_mapping/picard_marked_indexed/||g' Tcas_gene_counts.txt | sed 's|_marked_duplicates.bam||g' > Tcas_gene_counts_short_headers.txt

# rsync -azP milenatr@pelle.uppmax.uu.se:/proj/naiss2023-6-65/Milena/chapter3/RNAseq/maculatus/gene_counts/gene_counts_standard.txt /Users/miltr339/work/chapter3/DE_analysis/Cmac_gene_counts.txt
# sed 's|/proj/naiss2023-6-65/Milena/chapter3/RNAseq/maculatus/star_mapping/picard_marked_indexed/||g' Cmac_gene_counts.txt | sed 's|_marked_duplicates.bam||g' > Cmac_gene_counts_short_headers.txt

## R analysis with DESeq2 according to sebastian's github
# -Aggregate to transcript level by summing exon counts,
# -chose to only load and use the multimapped transcript counts,
# -filtered on ≥3 mean counts per sample in each sex,
# -DESeq2 for DE analysis based on male vs female,
# -used vst for count normalization with variance stabilization.

def counts_paths(username="miltr339"):
    counts_paths = {
        "Cmac" : f"/Users/{username}/work/chapter3/DE_analysis/Cmac_gene_counts_short_headers.txt",
        "Tcas" : f"/Users/{username}/work/chapter3/DE_analysis/Tcas_gene_counts_short_headers.txt",
        "Csep" : f"/Users/{username}/work/chapter3/DE_analysis/Csep_gene_counts_short_headers.txt"
    }
    return counts_paths

def metadata_paths(username="miltr339"):
    metadata = {
        "Cmac" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/metadata/Cmac_full_SRR_list.csv",
        "Tcas" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/metadata/Csept_full_SRR_list.csv",
        "Csep" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/metadata/Tcas_full_SRR_list.csv"
    }
    return metadata


if __name__ == "__main__":

    # species_list = ["Tcas"]
    species_list = ["Cmac","Tcas","Csep"] 
    counts_paths_dict = counts_paths()
    metadata_paths_dict = metadata_paths()

    for species in species_list:
        print(f"\n//////////////////////// {species} ////////////////////////")
        counts = pd.read_csv(counts_paths_dict[species], sep="\t", comment="#")
        # print(counts)
        metadata = pd.read_csv(metadata_paths_dict[species])
        print(metadata)