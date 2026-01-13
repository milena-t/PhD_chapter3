# Chapter III: FastX in coleoptera

## Other code

I wrote a lot of code for this a year ago here: https://github.com/milena-t/calculate_orthogroup_dNdS

## Notes

DTOL open data release policy [here](https://www.darwintreeoflife.org/wp-content/uploads/2024/10/DToL-Open-Data-Release-Policy.docx_.pdf)

### Methods

* [Martinez-Pacheco 2020](https://academic.oup.com/gbe/article/12/11/2015/5892261) 
  * double check y linked gene presence by `blastn` against the assembly.
  * *"the best BlastN match (usually around 92–95% identity over the entire sequence) onto the annotated X chromosome of the reference genomes was considered the X gametologs"*
  * [Marques 2005](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0030357) for identifying retrogenes. blast proteins against assembly, merge nearby matches, *"query and target sequences had >50% similarity on the amino acid level and over >80% of their length* \[are\] *shared"*, verify absence of introns. with some `paml` stuff, they identified the ancestral gametolog that all retrogenes originate from (useful for FastX?)
* **[Whittle 2020](https://academic.oup.com/g3journal/article/10/3/1125/6026234): no FastX in beetles**. only compared between *T. castaneum* and *T. freemani*, X enriched for female-biased genes, no X dosage compensation in the testes. 
  * TODO read in more detail about what they say about the influence of the (lack of) dosage compensation on fastX.
  * also check methods for alignment filtering: *"It has been suggested that removal of highly divergent segments from alignments, while causing loss of some sequence regions, improves measurements of protein sequence divergence; thus, highly divergent segments were excluded using the program Gblocks v. 0.91b set at default parameters (Castresana 2000; Talavera and Castresana 2007)."*
  * Divergence time 12-47 MYA ([Angelini 2008](https://www.sciencedirect.com/science/article/pii/S1055790307002941?via%3Dihub)) which is unfortunately not super informative when trying to figure out if they are more or less far apart than *C. maculatus* and *C. chinensis*.
* [Mank 2009](https://academic.oup.com/mbe/article/27/3/661/1000994?login=true) Faster-Z in birds is mainly due to drift
  * Positive selection would be fixation of recessive male-biased mutations
* [Mank 2007](https://academic.oup.com/mbe/article/24/12/2698/978299) faster evolutionary rate of female biased genes in bird brains (not male biased)
  * Note the ZW system in birds
  * 
* [Li 2010](https://pubmed.ncbi.nlm.nih.gov/21035095/) FastZ in duplicates compared to autosomal duplicates
  * within-species comparison, make all pairwise dN/dS of all genes within a gene family (check methods specifically that they use to reduce between-sample depencence for statistical power)
  
### Expectations

* Retrotransposition of male-biased genes from X to A ([Ellegren 2011](https://www.nature.com/articles/nrg2948.pdf))
* Ampliconic regions on both X and Y, expansion of intergenic regions

### Results

* TODO read: [Vicoso & Bachtrog 2006](https://www.nature.com/articles/nrg1914) Big early paper that everyone cites about the different evolutionary forces that act differently on the X vs the autosomes
* [Navarro 2003](https://www.science.org/doi/full/10.1126/science.1080600) accelerated evolutionary rate on rearranged chromosomes between species
* (Review)[Ellegren 2011](https://www.nature.com/articles/nrg2948.pdf) review on the influence of heterogameity on sex chromosome evolution
  * *"For example, long interspersed repeat elements are enriched on both the mammalian X and the avian Z chromosome \[46,47\], whereas gene  density is lower than on autosomes in both systems as a result of intergenic expansions \[27,48\]"*
  * He also says about the selection pressure on X-linked genes in the heterogametic sex that *selection will occur more frequently* as opposed to that it is stronger, which I think doesn't change the outcome because the selection is stronger in the end compared to the autosomes due to it occuring more frequently
  * *Among the genes that generate new retrocopies, through mRNA intermediates, there is an excess of X-linked genes inserted at autosomal locations.* Might be because X linked genes are temporarily inactivated during meiosis in males (meiotic sex-chromosome inactivation **MSCI**), and genes that give selective advantage to males want to escape. The retrocopies that leave often aquire male-specific function, because, if dominant, they spread easier on the autosomes because they are temporarily inactivated on the X negating their selective advantage
  * References about Gene traffic to and from the X chromosome
    * Emerson, J. J., Kaessmann, H., Betran, E. & Long, M. Extensive gene traffic on the mammalian X chromosome. Science 303, 537–540 (2004)
    * Shiao, M. S. et al. Origins of new male germ-line functions from X-derived autosomal retrogenes in the mouse. Mol. Biol. Evol. 24, 2242–2253 (2007).
    * Vinckenbosch, N., Dupanloup, I. & Kaessmann, H. Evolutionary fate of retroposed gene copies in the human genome. Proc. Natl Acad. Sci.USA 103, 3220–3225 (2006).
    * Meisel, R. P., Han, M. V. & Hahn, M. W. A complex suite of forces drives gene traffic from Drosophila X chromosomes. Genome Biol. Evol. 1, 176–188 (2009).
    * Vibranovski, M. D., Zhang, Y. & Long, M. General gene movement off the X chromosome in the Drosophila genus. Genome Res. 19, 897–903
  (2009)
* [Yu 2026](https://www.pnas.org/doi/full/10.1073/pnas.2522417123) super conserved noncoding sex determining locus in (haplodiploid) hymenoptera

## Workflow

<details>
<summary>Flowchart for my pipeline</summary>

```mermaid
graph TD;
    species_Ass(species assemblies);
    species_Ass --> simplified_bed{{bash/make_bedfiles_for_MCScanX.sh}};
    simplified_bed -- merge all species --> mcscanx_bed(mcscanx annotation input);
    
    species_Ann(species annotations);
    species_Ann -. new annotations .-> BUSCO([annotation evaluation]);
    species_Ann --> is_filter{{bash/isoform_filter_gff.sh}};
    is_filter --> get_prot{{bash/get_fasta_from_gff.sh}};
    get_prot --> protfiles(proteinfiles of all species);
    protfiles --> blast{{all-against-all proteinblast}};
    blast -- merge all species --> mcscanx_blast(mcscanx blast input);
    

    mcscanx_bed --> run_mcscanx{{run MCScanX}};
    mcscanx_blast --> run_mcscanx
    run_mcscanx --> mcscanx_out{{collinearity file}}
    mcscanx_out --> mcscanx_plot([protein synteny plot, pairwise and riparian plot with all species])
    mcscanx_plot .-> conf_X(confirm X-chromosome identification)

    species_Ass -- NCBI and SATC --> XY_chr(identify X and Y chromosmes);
    
    blast --> XY_paralogs(identify species with XY gametologs with blast BRH);
    dNdS_groups(make orthogroups to compare dNdS);
    XY_paralogs --> dNdS_groups;

    rearrangeA(dNdS for A depending on if A is rearranged between pairs or not)
    basicFastX(whole-phylogeny 1-to-1 orthologs for basic FastX);
    gametologs(dNdS for X genes depending on if they have an ancestral Y gametolog or not);

    dNdS_groups --> gametologs;
    dNdS_groups --> basicFastX;
    dNdS_groups --> rearrangeA;
    XY_chr --> gametologs
    XY_chr --> basicFastX
    mcscanx_out --> rearrangeA

    OG_dNdS([dNdS analysis based on orthofinder output with https://github.com/milena-t/calculate_orthogroup_dNdS])
    gametologs --> OG_dNdS
    basicFastX --> OG_dNdS
    rearrangeA --> OG_dNdS

```   

</details>

### Species selection

I want to look at groups of sister-species to take into account the evolutionary distance when comparing dN/dS ratio

* Include *Tribolium freemani* as a sister species to *T. castaneum* because of Whittle2020, T. freemani assembly [here](https://www.ebi.ac.uk/ena/browser/view/GCA_022388455.1)
  * no annotation available -> but RNAseq so I am annotating it myself from scratch
* include *Coccinella magnifica* (darwin tree of life project, [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/965/644/565/GCA_965644565.1_icCocMagn1.hap1.1/)) as sister species to *Coccinella septempunctata*
  * no annotation available
  * Timetree says that they are "tne same species" so unclear what that means for their phylogenetic distance

I have annotated *T. freemani* with available RNAseq data from Whittle2020, but there was no RNAseq data for *C. septempunctata*. 

<details>
<summary>Annotation Evaluation</summary>

#### BUSCO scores

```text
C. maculatus superscaffolded (adult LOME RNAseq but no larval)
(19022 genes)

C:98.7%[S:72.2%,D:26.5%],F:0.6%,M:0.7%,n:1013	   
999	Complete BUSCOs (C)			   
731	Complete and single-copy BUSCOs (S)	   
268	Complete and duplicated BUSCOs (D)	   
6	Fragmented BUSCOs (F)			   
8	Missing BUSCOs (M)			   
1013	Total BUSCO groups searched
```

#### Single exon genes

<p float="left">
  <img src="data/annotation_evaluation/single_exon_genes_white_bg.png" width="75%" />
</p>

</details>

### sex chromosome identification

* **Bruchids**
  * *C. maculatus:* from Kaufmann et al.: 
    ```python
    { X : ['utg000057l_1','utg000114l_1','utg000139l_1','utg000191l_1','utg000326l_1','utg000359l_1','utg000532l_1','utg000602l_1'],
      Y : ['utg000322l_1','utg 000312c_1','utg 000610l_1','utg 001235l_1']}
    ## superscaffolded
    { X : ['scaffold_10','scaffold_14','scaffold_23','scaffold_31','scaffold_34','scaffold_83'],
      Y : ['scaffold_26','scaffold_48','scaffold_103','scaffold_112','scaffold_164']}
    ``` 
  
  * *C. chinensis:*
    identified by me, but the assembly is very fragmented so there is a lot of contigs.

    ```python
    { X : ["1092_quiver","1124_quiver","1080_quiver","1105_quiver","1148_quiver","1339_quiver","1435_quiver","1482_quiver","1501_quiver","1565_quiver","1694_quiver","1688_quiver","1758_quiver","1618_quiver","1786_quiver","1816_quiver","1815_quiver","1817_quiver","1826_quiver","1889_quiver","1898_quiver","1908_quiver","1911_quiver","1933_quiver","2046_quiver","2054_quiver","2056_quiver","5713_quiver","2194_quiver","2226_quiver","2306_quiver","2357_quiver","2381_quiver","2392_quiver","2400_quiver","2435_quiver","2453_quiver","2513_quiver","2524_quiver","2569_quiver","2576_quiver","2580_quiver","2599_quiver","2693_quiver","2733_quiver","1210_quiver","2935_quiver","2958_quiver","2964_quiver","3034_quiver","3068_quiver","3080_quiver","3091_quiver"],
      Y : ["850_quiver","949_quiver","1088_quiver","1125_quiver","1159_quiver","1134_quiver","1224_quiver","1369_quiver","1410_quiver","1568_quiver","1577_quiver","1619_quiver","1634_quiver","1646_quiver","1652_quiver","1665_quiver","1681_quiver","1697_quiver","1722_quiver","1766_quiver","1783_quiver","1891_quiver","1937_quiver","1963_quiver","1790_quiver","1997_quiver","2073_quiver","2113_quiver","2163_quiver","2166_quiver","5705_quiver","2245_quiver","2259_quiver","2260_quiver","2334_quiver","2340_quiver","2382_quiver","2443_quiver","2511_quiver","2534_quiver","2573_quiver","2597_quiver","2651_quiver","2707_quiver","2766_quiver","2773_quiver","2791_quiver","2830_quiver","2875_quiver","3022_quiver","3070_quiver","3074_quiver","3075_quiver","3078_quiver"]}
    ```
  * *B. siliquastri* identified by the DTOL project with the assembly
    ```python
    { X : ['X'],
      Y : ['Y']}
    ``` 
  * *A. obtectus* Identified by me and Göran's project about it
    ```python
    #### TODO double check these contig IDs
    { X : ["CAVLJG010000002.1","CAVLJG010003236.1","CAVLJG010003544.1","CAVLJG010000099.1","CAVLJG010000155.1","CAVLJG010000244.1","CAVLJG010000377.1","CAVLJG010000488.1",],
      Y : ["CAVLJG010000343.1","CAVLJG010002896.1","CAVLJG010000233.1","CAVLJG010000566.1","CAVLJG010000588.1",]}
    ``` 
    <details>
    <summary>A_obtectus chromosome name conversions</summary>
    The chromosome names in the ENA version of the assembly are weird, a whole header is like 
    ```
    >ENA|CAVLJG010000002|CAVLJG010000002.1 Acanthoscelides obtectus genome assembly, contig: chr_10
    >ENA|CAVLJG010002896|CAVLJG010002896.1 Acanthoscelides obtectus genome assembly, contig: scaffold_36
    ```
    and somehow in the mapping for the sex chromosome identification that Angela did, the contigs are identified as `chr_10` and `scaffold_36`, and not as `CAVLJG010000002.1`. This is confusing because in the original version of this assembly on uppmax after superscaffolding before it was deposited, the contigs are named like `HiC_scaffold_1`, and these scaffold numbers do *not* correspond to the ones in the ENA/NCBI headers! To be super clear I list below the complete X and Y association, which I need for further analyses because the annotation obviously uses `CAVLJG010000002.1` and not some weird tail end of the header.

    ```text
    X-contigs > 100 kb:
    * chr_10 :       CAVLJG010000002.1, 43142905 bp
    * scaffold_49 :  CAVLJG010003236.1, 333528 bp
    * scaffold_77 :  CAVLJG010003544.1, 217000 bp
    * scaffold_108 : CAVLJG010000099.1, 157500 bp
    * scaffold_113 : CAVLJG010000155.1, 147000 bp
    * scaffold_121 : CAVLJG010000244.1, 129500 bp
    * scaffold_133 : CAVLJG010000377.1, 114918 bp
    * scaffold_143 : CAVLJG010000488.1, 105500 bp

    Y-contigs > 100 kb:
    * scaffold_13 :  CAVLJG010000343.1, 1187000 bp
    * scaffold_36 :  CAVLJG010002896.1, 421500 bp
    * scaffold_120 : CAVLJG010000233.1, 131000 bp
    * scaffold_150 : CAVLJG010000566.1, 100909 bp
    * scaffold_152 : CAVLJG010000588.1, 100000 bp
    ```
    </details>

* **Diorhabada** 
  * [Diorhabda sublineata](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_026230105.1/)
  ```python
  { X : ['NC_079485.1'],
    Y : ['NC_079486.1']} 
  ``` 
  * [Diorhabda carinulata](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_026250575.1/)
  ```python
  { X : ['NC_079472.1'],
    Y : ['NC_079473.1']} 
  ``` 
  * There is also Diabrotica undecimpunctata but it does not have sex chromosomes identified, but it does have a giant genome at 1.7Gb



# Riparian plot for synteny

I use MCscanX to plot synteny, and [SynVisio](https://synvisio.github.io/) to visualize the riparian plot. *C. chinensis* is not superscaffolded and therefore not included in the plot.

<p float="left">
  <img src="data/synteny_plots/synvisio_all_riparian.png" width="100%" />
</p>

The X chromosomes are all shown at the right-most side of the plot except *D. carinulata*:
* ao1950 `CAVLJG010000002`
* bs9
* cm124 `utg000057l_1` (scaffold 10 in the superscaffolded version of Cmac)
* dc35 (!) (on the left in the plot)
* ds21

Synteny is conserved within *Bruchini* and within *Diorhabda*, but not between them.

Notably, there is one autosome missing in *A. obtectus*. Since the SynVisio website is not great to use with more species and when unplaced contigs are included, I use the setting to exclude contigs with no synteny blocks on them, which apparently excludes the one missing *A. obtectus* contig.

## No syntenic X in *Diorhabda carinulata*

Despite the broken-down synteny of the autosomes between *Bruchini* and *Diorhabda*, the X chromosomes are all syntenic except *D. carinulata* (dc 35), which is in agreement with the identified 1-to-1 orthologs below, where most pairwise comparisons have 300-400 X-linked orthologs, except any pair involving *D. carinulata*. The "old" syntenic chromosome in *D. carinulata* is `NC_079460.1` (dc13). A part of the *D. carinulata* X is syntenic to *D. sublineata* `NC_079475.1` (ds5).

# ortholog identification

I will use the complete X contig list from the beginning, *not* just the syntenic contigs shown in the riparian plot.

## orthofinder

In the original analysis with Lila we ran orthofinder and used the orthologs that are 1-to-1 across all species. With the new species selection, I tried that here as well, but it seems like the *Diorhabda* X chromosomes are not syntenic. 

<p float="left">
  <img src="data/orthofinder/orthogroups_by_contig_upset_plot.png" width="65%" />
</p>

## best reciprocal hits

All comparisons of dNdS to identify fastX are always pairwise. Therefore it doesn't matter that a gene is 1-to-1 in all the species, only in the pair where i actually do the calculation. I will use the all-vs-all proteinblast that I made for the synteny with MCScanX (see chapter 2) to identify best reciprocal hits between pairs instead.

<p float="left">
  <img src="data/fastX_ortholog_ident/BRH_A_linked_counts_heatmap.png" width="45%" />
  <img src="data/fastX_ortholog_ident/BRH_X_linked_counts_heatmap.png" width="45%" />
</p>


I also check for gametologs (genes where there is a 1-to-1 homolog on the X and Y chromosome in the same species). the within-species gametolog counts are unfortunately very low:

* *A. obtectus* : 3
* *B. siliquastri* : 4
* *C. chinensis* : 1
* *C. maculatus* : 4
* *D. carinulata* : 3
* *D. sublineata* : 0

I suspect that this is not because there are no ancestral gametologs, but because there have been more recent duplications on X and/or Y. This would make them not 1-to-1 any more and therefore not identified in this test.

## Conclusion

I will continue with best reciprocal hits, becasue the orthofinder clustering is a bit more "relaxed", meaning that a transcript is more likely to be included in an orthogorup. Since I am only using pairwise comparisons, this method may make larger orthogroups that would not qualify as 1-to-1 orthologs, and therefore reduce the sample size for our dNdS analysis.

## dNdS 

The dNdS is calculated between all pairwise comparisons for A-linked and X-linked genes identified with BRH above. I am doing unique comparisons only below, so (Cmac,Aobt) is the same as (Aobt,Cmac) which is why only one of them is shown. Also, self-comparisons make no sense here so I am excluding them as well. I show horizontal lines of the median dNdS for A and X in matching colors to the violin plot to clearly show the mean difference in dNdS. 

Most species are slowX, except the *D. carinulata* comparisons, but those have such a low sample size for X-linked genes that the mean is not reliable. 

### Statistical analysis via permutation test

Generally, all species pairs seem to show slowX, but I will do a permutation test for each species pair like this:

1. **Sample without replacement:** Merge A and X dNdS values into one list, sample new A and new X lists with original sample sizes without replacement
2. **Calculate median:** Calculate median A and median X from resampled lists, use difference A - X
3. **Normal distribution:** Compute a normal distribution from all median(A) - median(X) values, find 95% confidence intervals
4. **Plot:** Histrogram and normal distribution, compare original observed dNdS_A - dNdS_X to 95% confidence intervals of the permutations

I have chosen 10000 permutations for now, this takes less than 10 mins. 

The top right is all the permutation tests. the pink line is the measured `dNdS_A - dNdS_X`, while green is the distributions from the permutation test. Since I am plotting `dNdS_A - dNdS_X`, and I mostly observe `dNdS_A - dNdS_X > 0`, which means `dNdS_A > dNdS_X`, this indicates slowX. the bottom right is the same violin plot as above. Note the low sample sizes for *D. carinulata* X, which are "hidden" in the permutation test.

<p float="left">
  <img src="data/fastX_ortholog_ident/fastX_permutation_white_bg.png" width="100%" />
</p>

### X chromosome turnover in *D. carinulata*

As seen in the synteny plot, the X of *D. carinulata* is not syntenic with any of the other species. However, the chromosome dc13 (`NC_079460.1`) is syntenic with the remaining X chromosomes. I hypothesize that this is the result of a recent sex-chromosome turnover in *D. carinulata* where the sex determining locus has moved from the X-syntenic dc13 to dc35 (if the X chromosome is identified correctly). I will repeat the above analysis with the dc13 chromosome for *D. carinulata* instead of the dc35 one which is the actual X, and if the sex chromosome turnover is true I would expect that this still shows slowX since the turnover is much more recent than the long-term evolutionary forces that shape slowX evolution in *Chrysomelidae* (and probably *Coleoptera* in general).