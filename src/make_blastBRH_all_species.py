### run get_blast_BRH.py 

import get_blast_BRH as brh
import parse_gff as gff

def outfiles_dict(username = "miltr339"):
    blast_outdir = f"/Users/{username}/work/pairwise_blast_chapter_2_3/"
    blast_out_dict = {
        "A_obtectus" : {
            "A_obtectus" : f"{blast_outdir}A_obtectus_original_header_vs_A_obtectus_original_header.blast",
            "B_siliquastri" : f"{blast_outdir}A_obtectus_original_header_vs_B_siliquastri_original_header.blast",
            "C_chinensis" : f"{blast_outdir}A_obtectus_original_header_vs_C_chinensis_original_header.blast",
            "C_maculatus" : f"{blast_outdir}A_obtectus_original_header_vs_C_maculatus_superscaffolded_original_header.blast",
            "D_carinulata" : f"{blast_outdir}A_obtectus_original_header_vs_D_carinulata_original_header.blast",
            "D_sublineata" : f"{blast_outdir}A_obtectus_original_header_vs_D_sublienata_original_header.blast",
        },
        "B_siliquastri" : {
            "A_obtectus" : f"{blast_outdir}B_siliquastri_original_header_vs_A_obtectus_original_header.blast",
            "B_siliquastri" : f"{blast_outdir}B_siliquastri_original_header_vs_B_siliquastri_original_header.blast",
            "C_chinensis" : f"{blast_outdir}B_siliquastri_original_header_vs_C_chinensis_original_header.blast",
            "C_maculatus" : f"{blast_outdir}B_siliquastri_original_header_vs_C_maculatus_superscaffolded_original_header.blast",
            "D_carinulata" : f"{blast_outdir}B_siliquastri_original_header_vs_D_carinulata_original_header.blast",
            "D_sublineata" : f"{blast_outdir}B_siliquastri_original_header_vs_D_sublienata_original_header.blast",
        },
        "C_chinensis" : {
            "A_obtectus" : f"{blast_outdir}C_chinensis_original_header_vs_A_obtectus_original_header.blast",
            "B_siliquastri" : f"{blast_outdir}C_chinensis_original_header_vs_B_siliquastri_original_header.blast",
            "C_chinensis" : f"{blast_outdir}C_chinensis_original_header_vs_C_chinensis_original_header.blast",
            "C_maculatus" : f"{blast_outdir}C_chinensis_original_header_vs_C_maculatus_superscaffolded_original_header.blast",
            "D_carinulata" : f"{blast_outdir}C_chinensis_original_header_vs_D_carinulata_original_header.blast",
            "D_sublineata" : f"{blast_outdir}C_chinensis_original_header_vs_D_sublienata_original_header.blast",
        },
        "C_maculatus" : {
            "A_obtectus" : f"{blast_outdir}C_maculatus_superscaffolded_original_header_vs_A_obtectus_original_header.blast",
            "B_siliquastri" : f"{blast_outdir}C_maculatus_superscaffolded_original_header_vs_B_siliquastri_original_header.blast",
            "C_chinensis" : f"{blast_outdir}C_maculatus_superscaffolded_original_header_vs_C_chinensis_original_header.blast",
            "C_maculatus" : f"{blast_outdir}C_maculatus_superscaffolded_original_header_vs_C_maculatus_superscaffolded_original_header.blast",
            "D_carinulata" : f"{blast_outdir}C_maculatus_superscaffolded_original_header_vs_D_carinulata_original_header.blast",
            "D_sublineata" : f"{blast_outdir}C_maculatus_superscaffolded_original_header_vs_D_sublienata_original_header.blast",
        },
        "D_carinulata" : {
            "A_obtectus" : f"{blast_outdir}D_carinulata_original_header_vs_A_obtectus_original_header.blast",
            "B_siliquastri" : f"{blast_outdir}D_carinulata_original_header_vs_B_siliquastri_original_header.blast",
            "C_chinensis" : f"{blast_outdir}D_carinulata_original_header_vs_C_chinensis_original_header.blast",
            "C_maculatus" : f"{blast_outdir}D_carinulata_original_header_vs_C_maculatus_superscaffolded_original_header.blast",
            "D_carinulata" : f"{blast_outdir}D_carinulata_original_header_vs_D_carinulata_original_header.blast",
            "D_sublineata" : f"{blast_outdir}D_carinulata_original_header_vs_D_sublienata_original_header.blast",
        },
        "D_sublineata" : {
            "A_obtectus" : f"{blast_outdir}D_sublienata_original_header_vs_A_obtectus_original_header.blast",
            "B_siliquastri" : f"{blast_outdir}D_sublienata_original_header_vs_B_siliquastri_original_header.blast",
            "C_chinensis" : f"{blast_outdir}D_sublienata_original_header_vs_C_chinensis_original_header.blast",
            "C_maculatus" : f"{blast_outdir}D_sublienata_original_header_vs_C_maculatus_superscaffolded_original_header.blast",
            "D_carinulata" : f"{blast_outdir}D_sublienata_original_header_vs_D_carinulata_original_header.blast",
            "D_sublineata" : f"{blast_outdir}D_sublienata_original_header_vs_D_sublienata_original_header.blast",
        }
    }

    annotation_out_dir = f"/Users/{username}/work/chapter3/isoform_filtered_native_annotations/"
    annotation_dict = {
        "A_obtectus" : f"{annotation_out_dir}A_obtectus.gff",
        "B_siliquastri" : f"{annotation_out_dir}B_siliquastri.gff",
        "C_chinensis" : f"{annotation_out_dir}C_chinensis.gff",
        "C_maculatus" : f"{annotation_out_dir}C_maculatus_superscaffolded.gff",
        # "C_magnifica" : f"{annotation_out_dir}C_magnifica.gff",
        # "C_septempunctata" : f"{annotation_out_dir}C_septempunctata.gff",
        "D_carinulata" : f"{annotation_out_dir}D_carinulata.gff",
        "D_sublineata" : f"{annotation_out_dir}D_sublineata.gff",
        # "T_castaneum" : f"{annotation_out_dir}T_castaneum.gff",
        # "T_freemani" : f"{annotation_out_dir}T_freemani.gff",
    }
    return blast_out_dict, annotation_dict


if __name__ == "__main__":

    blast_out_dict, annotation_dict = outfiles_dict()
    sex_chromosomes_lists = brh.sex_chromosome_names()
    
    annotation_Class= { species : gff.parse_gff3_general(annotation_path, verbose=False) for species, annotation_path in annotation_dict.items() }

    done_sets = []

    for species1 in annotation_dict.keys():
        print(f"{species1}")
        for species2 in annotation_dict.keys():

            # don't do bidirectional species iterations. e.g. if species1 = Cmac and species2 = Bsil, skip if Bsil, Cmac already generated
            curr_set = set([species1,species2])
            if curr_set in done_sets:
                continue

            print(f" - {species2}")
            blast_infile_path1 = blast_out_dict[species1][species2]
            blast_infile_path2 = blast_out_dict[species2][species1]
            annotation1 = annotation_Class[species1]
            annotation2 = annotation_Class[species2]
            X_list1 = sex_chromosomes_lists[species1]["X"]
            X_list2 = sex_chromosomes_lists[species2]["X"]
            Y_list1 = sex_chromosomes_lists[species1]["Y"]
            Y_list2 = sex_chromosomes_lists[species2]["Y"]

            outfile = f"{blast_infile_path1}_BRH.tsv"
            besthits_infile1 = brh.read_best_hits(blast_infile_path1)
            besthits_infile2 = brh.read_best_hits(blast_infile_path2)
            
            BRH_dict = brh.get_BRHs(besthits_infile1, besthits_infile2, annotation1=annotation1, annotation2=annotation2, x_list1=X_list1, x_list2=X_list2, y_list1=Y_list1, y_list2=Y_list2, species1=species1, species2=species2,outfile_path=outfile)
            
            done_sets.append(curr_set)
            