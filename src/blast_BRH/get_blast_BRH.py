"""
Get the transcript IDs of best reciprocal hits (BRH) of two species. the input is the outfmt 6 blast results in either direction
useage: python3 get_blast_BRH.py species1_species2.blast species2_species1.blast
and the output is a tsv file with two columns, one is species1 and one is species2, and each row contains a set of gene IDs that are each others BRH
I identify the best hit by the highest bit score
"""
import sys
import parse_gff as gff
import argparse

def parse_args():
    # Create the parser
    program_description = """
script to compute best reciprocal hits and their chromosome location from a bidirectional all-vs-all proteinblast search and the corresponding species annotations
the result is a tsv file that includes the transcript IDs of the best hits and their chromosome location as X/Y/A (not actual contig IDs!)
"""
    parser = argparse.ArgumentParser(description=program_description)

    # Add the arguments
    parser.add_argument('--blast1', type=str, required=True, help='Absolute filepath to the all-vs-all blastp output where species1 is the query (outfmt 6')
    parser.add_argument('--blast2', type=str, required=True, help='Absolute filepath to the all-vs-all blastp output where species2 is the query (outfmt 6')
    parser.add_argument('-o', '--outfile', type=str, help='output filename, default: filename from blast1 + _BRH.tsv')
    parser.add_argument('--annotation1', type=str, required=True, help='Absolute filepath to the annotation that blast1 query is based on')
    parser.add_argument('--annotation2', type=str, required=True, help='Absolute filepath to the annotation that blast2 query is based on')
    parser.add_argument('--verbose', action='store_true', help="enable verbose mode")
    parser.add_argument('--X_contigs1', type=str, required=True, help='comma separated list of the names of X-linked contigs in species 1 (given directly in the command line with no spaces, not a path to a file!')
    parser.add_argument('--X_contigs2', type=str, required=True, help='comma separated list of the names of X-linked contigs in species 2 (given directly in the command line with no spaces, not a path to a file!')
    parser.add_argument('--Y_contigs1', type=str, help='comma separated list of the names of Y-linked contigs in species 1 (given directly in the command line with no spaces, not a path to a file!')
    parser.add_argument('--Y_contigs2', type=str, help='comma separated list of the names of Y-linked contigs in species 2 (given directly in the command line with no spaces, not a path to a file!')
 

    # Parse the arguments
    args = parser.parse_args()

    # set default values for non-obligatory arguments
    if not args.outfile:
        args.outfile = f"{args.blast1}_BRH.tsv"

    return args

blast_outfmt6_headers = ["qseqid", "rseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

def sex_chromosome_names():
    contig_names_dict = {
        "A_obtectus" : {
            "X" : ["CAVLJG010000002.1","CAVLJG010003236.1","CAVLJG010003544.1","CAVLJG010000099.1","CAVLJG010000155.1","CAVLJG010000244.1","CAVLJG010000377.1","CAVLJG010000488.1",],
            "Y" : ["CAVLJG010000343.1","CAVLJG010002896.1","CAVLJG010000233.1","CAVLJG010000566.1","CAVLJG010000588.1",]
        },
        "B_siliquastri" : {
            "X" : ['X'],
            "Y" : ['Y']
        },
        "C_chinensis" : {
            "X" : ["1211|quiver","1844|quiver","854|quiver","5741|quiver","2866|quiver","658|quiver","1498|quiver","1455|quiver","2404|quiver","2935|quiver","1115|quiver","370|quiver","2273|quiver","1424|quiver","1865|quiver","767|quiver","2222|quiver","1525|quiver","5023|quiver","1925|quiver","1217|quiver","2328|quiver","2475|quiver","959|quiver","537|quiver","2776|quiver","325|quiver","2576|quiver","2336|quiver","988|quiver","2252|quiver","1388|quiver","1508|quiver","1712|quiver","1260|quiver","977|quiver","2202|quiver","2223|quiver","2397|quiver","693|quiver","1092|quiver","2189|quiver","1958|quiver","1355|quiver","2241|quiver","849|quiver","703|quiver","277|quiver","518|quiver","2589|quiver","1326|quiver","2962|quiver","2341|quiver","358|quiver","462|quiver","2786|quiver","1116|quiver","525|quiver","1358|quiver","5693|quiver","1429|quiver","1253|quiver","2372|quiver","326|quiver","474|quiver","777|quiver","955|quiver","1852|quiver","718|quiver","1024|quiver","1974|quiver","2295|quiver","2356|quiver","1484|quiver","1503|quiver","3076|quiver","2091|quiver","1262|quiver","1109|quiver","1475|quiver","1695|quiver","1168|quiver","1386|quiver","2201|quiver","2320|quiver","1117|quiver","769|quiver","2050|quiver","1805|quiver","2692|quiver","411|quiver","851|quiver","5703|quiver","1585|quiver","824|quiver","1816|quiver","1370|quiver","2416|quiver","1814|quiver","1277|quiver","619|quiver","1750|quiver","2709|quiver","2664|quiver","1250|quiver","971|quiver","3020|quiver","310|quiver","1176|quiver","2510|quiver","1699|quiver","1256|quiver","1420|quiver","5727|quiver","413|quiver","1124|quiver","682|quiver","1000|quiver","1313|quiver","5708|quiver","1556|quiver","274|quiver","1787|quiver","1137|quiver","360|quiver","1469|quiver","1853|quiver","2380|quiver","1239|quiver","993|quiver","791|quiver","2540|quiver","1510|quiver","868|quiver","505|quiver","1212|quiver","376|quiver","1564|quiver","1836|quiver","1670|quiver","500|quiver","2099|quiver","353|quiver","1042|quiver","419|quiver","1314|quiver","1339|quiver","1470|quiver","1576|quiver","1717|quiver","5692|quiver","2157|quiver","700|quiver","1284|quiver","1694|quiver","2306|quiver","2712|quiver","182|quiver","1973|quiver","882|quiver","2363|quiver","2482|quiver","1640|quiver","1913|quiver","2323|quiver","1240|quiver","161|quiver","1649|quiver","1164|quiver","1054|quiver","1096|quiver","313|quiver","1815|quiver","1831|quiver","1349|quiver","151|quiver","1478|quiver","1523|quiver","1888|quiver","739|quiver","1322|quiver","2338|quiver","1798|quiver","1391|quiver","1530|quiver","1519|quiver","1651|quiver","1105|quiver","509|quiver","1308|quiver","1833|quiver","1914|quiver","1741|quiver","1080|quiver","2292|quiver","2364|quiver","643|quiver","5745|quiver","1920|quiver","1725|quiver","125|quiver","1086|quiver","2552|quiver","5749|quiver","2120|quiver","2964|quiver","5722|quiver","2045|quiver","2422|quiver","593|quiver","1496|quiver","1772|quiver","799|quiver","2690|quiver","414|quiver","1531|quiver","1443|quiver","1408|quiver","1688|quiver","1371|quiver","1501|quiver","3090|quiver","1025|quiver","5698|quiver","347|quiver","1435|quiver","476|quiver","1883|quiver","2820|quiver","5728|quiver","342|quiver","1972|quiver","1826|quiver","968|quiver","2037|quiver","1723|quiver","252|quiver","1863|quiver","2983|quiver","1947|quiver","1430|quiver","1612|quiver","1701|quiver","839|quiver","613|quiver","1979|quiver","1584|quiver","2024|quiver","1486|quiver","1097|quiver"],
            # ["1092|quiver","1124|quiver","1080|quiver","1105|quiver","1148|quiver","1339|quiver","1435|quiver","1482|quiver","1501|quiver","1565|quiver","1694|quiver","1688|quiver","1758|quiver","1618|quiver","1786|quiver","1816|quiver","1815|quiver","1817|quiver","1826|quiver","1889|quiver","1898|quiver","1908|quiver","1911|quiver","1933|quiver","2046|quiver","2054|quiver","2056|quiver","5713|quiver","2194|quiver","2226|quiver","2306|quiver","2357|quiver","2381|quiver","2392|quiver","2400|quiver","2435|quiver","2453|quiver","2513|quiver","2524|quiver","2569|quiver","2576|quiver","2580|quiver","2599|quiver","2693|quiver","2733|quiver","1210|quiver","2935|quiver","2958|quiver","2964|quiver","3034|quiver","3068|quiver","3080|quiver","3091|quiver"],
            "Y" : ["850|quiver","949|quiver","1088|quiver","1125|quiver","1159|quiver","1134|quiver","1224|quiver","1369|quiver","1410|quiver","1568|quiver","1577|quiver","1619|quiver","1634|quiver","1646|quiver","1652|quiver","1665|quiver","1681|quiver","1697|quiver","1722|quiver","1766|quiver","1783|quiver","1891|quiver","1937|quiver","1963|quiver","1790|quiver","1997|quiver","2073|quiver","2113|quiver","2163|quiver","2166|quiver","5705|quiver","2245|quiver","2259|quiver","2260|quiver","2334|quiver","2340|quiver","2382|quiver","2443|quiver","2511|quiver","2534|quiver","2573|quiver","2597|quiver","2651|quiver","2707|quiver","2766|quiver","2773|quiver","2791|quiver","2830|quiver","2875|quiver","3022|quiver","3070|quiver","3074|quiver","3075|quiver","3078|quiver"]
        },
        # "C_maculatus_Kaufmann2023" : {
        #     "X" : ['utg000057l_1','utg000114l_1','utg000139l_1','utg000191l_1','utg000326l_1','utg000359l_1','utg000532l_1','utg000602l_1'],
        #     "Y" : ['utg000322l_1','utg000312c_1','utg000610l_1','utg001235l_1']
        # },
        "C_maculatus" : {
            "X" : ['scaffold_10','scaffold_14','scaffold_23','scaffold_31','scaffold_34','scaffold_83'],
            "Y" : ['scaffold_26','scaffold_48','scaffold_103','scaffold_112','scaffold_164']
        },
        "D_carinulata" : {
            "X" : ['NC_079460.1'], ## syntenic to other X chromosomes
            # "X" : ['NC_079472.1'], ## originally identified
            "Y" : ['NC_079473.1']
        },
        "D_sublineata" : {
            "X" : ['NC_079485.1'],
            "Y" : ['NC_079486.1']
        },
        "T_castaneum" :  { 
            "X" : ['NC_087403.1'], # I think based on synteny, not identified on the NCBI
            "Y" : ['unidentified']
        },
        "T_freemani" :  { 
            "X" : ['CM039461.1'], # identified as linkage group X (LGX) on the NCBI
            "Y" : ['unidentified']
        },
        "C_septempunctata" :  { 
            "X" : ['NC_058198.1'],
            "Y" : ['unidentified']
        },
        "C_magnifica" :  { 
            "X" : ['OZ286750.1'],
            "Y" : ['unidentified']
        },
    }

    return contig_names_dict

class BestHit:
    """
    which reference ID is the best hit to the query ID
    """
    
    def __init__(self, qseqid:str, rseqid:str, bitscore:float) -> None:
        self.qseqid = qseqid
        self.rseqid = rseqid
        self.bitscore = bitscore

    def update_besthit(self, new_qseqid:str, new_rseqid:str, new_bitscore:float) -> bool:
        assert new_qseqid == self.qseqid
        if new_bitscore > self.bitscore:
            self.rseqid = new_rseqid
            self.bitscore = new_bitscore
            return True
        else:
            return False

    def __str__(self) -> str:
        return(f"  * {self.qseqid} has the best hit {self.rseqid} with a bitscore of {self.bitscore}")


def read_best_hits(blast_infile_path:str) -> dict:
    """
    make a dictionary of all the best hits instances
    { queryID : BestHit ,  ...}
    """
    best_hits_dict = {}
    with open(blast_infile_path, "r") as blast_infile:
        blast_lines = blast_infile.readlines()
        for line in blast_lines:
            line = line.strip().split("\t")
            qseqid = line[0]
            rseqid = line[1]

            if qseqid == rseqid:
                # for self-blast searches to find paralogs, skip self-hits
                continue
            try:
                bitscore = float(line[-1])
            except:
                raise RuntimeError(f"{line} \nfrom {blast_infile_path}\ncould not be parsed\n")
            if qseqid not in best_hits_dict:
                best_hits_dict[qseqid] = BestHit(qseqid=qseqid, rseqid=rseqid, bitscore=bitscore)
            else:
                best_hits_dict[qseqid].update_besthit(new_qseqid = qseqid, new_rseqid = rseqid, new_bitscore = bitscore)
    return best_hits_dict  


def get_BRHs(besthits_infile1, besthits_infile2, annotation1, annotation2, x_list1, x_list2, y_list1 = [], y_list2 = [], species1 = "", species2 = "", outfile_path:str = ""):
    """
    get a dictionary with {species1_ID : species2_ID} of all best reciprocal hits
    """
    out_dict = {}
    if type(annotation1) == str and type(annotation2)== str:
        ## use annotation filenames as species headers
        species1 = annotation1.split("/")[-1].split(".")[0]
        species2 = annotation2.split("/")[-1].split(".")[0]
        ## read annotations from filepaths
        print(f"reading annotations of species 1 and 2...")
        annotation1 = gff.parse_gff3_general(annotation1, verbose=False)
        annotation2 = gff.parse_gff3_general(annotation2, verbose=False)
    elif type(annotation1) == dict and type(annotation2) == dict:
        ## use species names from parameters and assume that the annotations are already read in
        if species1 == "" or species2 == "":
            raise RuntimeError(f"if you include parsed annotations you have to give species names in the function parameters. You have given:\n species1 = '{species1}'\n species2 = '{species2}'")
        pass
    header = f"{species1}\tchromosome\t{species2}\tchromosome\n"

    for species1_id, species1_besthit in besthits_infile1.items():
        species2_id = species1_besthit.rseqid
        try:
            species2_besthit = besthits_infile2[species2_id]
        except:
            continue
        if species1_id == species2_besthit.rseqid:
            out_dict[species1_id] = species2_id
    if outfile_path != "":
        if len(out_dict)>0:
            with open(outfile_path, "w") as outfile:
                outfile.write(header)
                for species1_id, species2_id in out_dict.items():
                    ## remove the "_1" suffix so that they can be found in the annotation
                    if species1_id[-2:] == "_1":
                        species1_id = species1_id[:-2]
                    if species2_id[-2:] == "_1":
                        species2_id = species2_id[:-2]
                    contig1 = "A"
                    contig2 = "A"
                    try:
                        ID_contig1 = annotation1[species1_id].contig
                    except:
                        raise RuntimeError(f"{species1_id} in {species1} not found in the annotation")
                    try:
                        ID_contig2 = annotation2[species2_id].contig
                    except:
                        raise RuntimeError(f"{species2_id} in {species2} not found in the annotation")
                    if ID_contig1 in x_list1:
                        contig1 = "X"
                    elif ID_contig1 in y_list1:
                        contig1 = "Y"
                    if ID_contig2 in x_list2:
                        contig2 = "X"
                    elif ID_contig2 in y_list2:
                        contig2 = "Y"
                    outfile.write(f"{species1_id}\t{contig1}\t{species2_id}\t{contig2}\n")
        print(f"outfile written to: {outfile_path}")
    return out_dict


if __name__ == "__main__":
    

    args=parse_args()
    blast_infile_path1 = args.blast1
    blast_infile_path2 = args.blast2
    outfile = args.outfile
    annotation_path1 = args.annotation1
    annotation_path2 = args.annotation2
    verbose = args.verbose
    X_list1 = args.X_contigs1.strip().split(",")
    X_list2 = args.X_contigs2.strip().split(",")
    if args.Y_contigs1 and args.Y_contigs2:
        Y_list1 = args.Y_contigs1.strip().split(",")
        Y_list2 = args.Y_contigs2.strip().split(",")
    else:
        Y_list1 = []
        Y_list2 = []

    besthits_infile1 = read_best_hits(blast_infile_path1)
    besthits_infile2 = read_best_hits(blast_infile_path2)
    # print(besthits_infile1["rna-AOBTE_LOCUS3-2_1"])

    BRH_dict = get_BRHs(besthits_infile1, besthits_infile2, annotation1=annotation_path1, annotation2=annotation_path2, x_list1=X_list1, x_list2=X_list2, y_list1=Y_list1, y_list2=Y_list2, outfile_path=outfile)

