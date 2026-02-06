"""
in order to incorporate the differential expression analysis into the molecular rate analysis, I need to associate the transcript IDs from the DE analysis
with the orthogroup-pairs from the blast BRH results. all ortholog pairs are identified by their species pair members and a unique number.
"""

import os
from Bio import SeqIO


def make_lookup_dict(head_dir:str, out_path:str, skip_list:list = ["D_carinulata"]):
    """
    loop through all subdirectories (species pairs) and get all fasta files in each species pair
    save to outfile with the format ortholog_ID:transcript1,transcript2
    """
    dirs_list = [dirname for dirname in os.listdir(head_dir) if os.path.isdir(f"{head_dir}{dirname}")]
    with open(out_path, "w") as outfile:
        for pair_dir in dirs_list:
            # get all pair subdirectories
            try:
                g1,s1,g2,s2=pair_dir.split("_")
                species1 = f"{g1}_{s1}"
                species2 = f"{g2}_{s2}"
            except:
                raise RuntimeError(f"couldn't parse pair directory {pair_dir}!")
            if species1 == species2:
                continue
            elif species1 in skip_list or species2 in skip_list:
                continue
            pair_dir_full = f"{head_dir}{pair_dir}/"
            print(f"- {species1}, {species2}")

            # loop through all fastas
            fastas_list = [ortholog_fasta for ortholog_fasta in os.listdir(pair_dir_full) if os.path.isfile(f"{pair_dir_full}{ortholog_fasta}")]
            for fasta in fastas_list:
                fasta_path=f"{pair_dir_full}{fasta}"
                IDs_list = extract_ID_from_fasta(fasta_path, species1, species2)
                if len(IDs_list)==0:
                    line = f"{ortholog_ID}:not_found\n"
                elif len(IDs_list)!=2:
                    raise RuntimeError(f"fasta file {fasta_path} could not be parsed right: headers are these: \n{IDs_list}")
                elif len(IDs_list)==2:
                    ortholog_ID = fasta.replace(".fasta", "")
                    line = f"{ortholog_ID}:{IDs_list[0]},{IDs_list[1]}\n"
                outfile.write(line)    
        print(f"outfile saved to: {out_path}")


def extract_ID_from_fasta(fasta_path, species1, species2):
    """
    get a list of all fasta headers (transcript IDs) in a multifasta file 
    """
    out_list = []
    for species_record in SeqIO.parse(fasta_path, "fasta"):
        header = species_record.id.replace(">", "")
        # remove species IDs from fasta header
        header = header.split(f"_{species1}")[0]
        header = header.split(f"_{species2}")[0]
        out_list.append(header)

    return out_list
        


if __name__ == "__main__":
    
    username = "miltr339"
    outdir_X = f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/brh_sequences_X/"
    outdir_A = f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/brh_sequences_A/"

    lookup_path_X = f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/ortholog_IDs_X_transcript_IDs_association.txt"
    lookup_path_A = f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/ortholog_IDs_A_transcript_IDs_association.txt"

    # make_lookup_dict(head_dir=outdir_X, out_path=lookup_path_X)
    make_lookup_dict(head_dir=outdir_A, out_path=lookup_path_A)