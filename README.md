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


## Species selection

I want to look at groups of sister-species to take into account the evolutionary distance when comparing dN/dS ratio

* Include *Tribolium freemani* as a sister species to *T. castaneum* because of Whittle2020, T. freemani assembly [here](https://www.ebi.ac.uk/ena/browser/view/GCA_022388455.1)
  * no annotation available -> but RNAseq so I am annotating it myself from scratch
* include *Coccinella magnifica* (darwin tree of life project, [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/965/644/565/GCA_965644565.1_icCocMagn1.hap1.1/)) as sister species to *Coccinella septempunctata*
  * no annotation available
  * Timetree says that they are "tne same species" so unclear what that means for their phylogenetic distance

### sex chromosome identification

* **Bruchids**
  * *C. maculatus:* from Kaufmann et al.: 
    ```python
    { X : ['utg000057l_1','utg000114l_1','utg000139l_1','utg000191l_1','utg000326l_1','utg000359l_1','utg000532l_1','utg000602l_1'],
      Y : ['utg000322l_1','utg 000312c_1','utg 000610l_1','utg 001235l_1']}
    ``` 
  
  * *C. chinensis:*
    identified by me, but the assembly is very fragmented so there is a lot of contigs.

    ```python
    { X : ["1092_quiver","1124_quiver","1080_quiver","1105_quiver","1148_quiver","1339_quiver","1435_quiver","1482_quiver","1501_quiver","1565_quiver","1694_quiver","1688_quiver","1758_quiver","1618_quiver","1786_quiver","1816_quiver","1815_quiver","1817_quiver","1826_quiver","1889_quiver","1898_quiver","1908_quiver","1911_quiver","1933_quiver","2046_quiver","2054_quiver","2056_quiver","5713_quiver","2194_quiver","2226_quiver","2306_quiver","2357_quiver","2381_quiver","2392_quiver","2400_quiver","2435_quiver","2453_quiver","2513_quiver","2524_quiver","2569_quiver","2576_quiver","2580_quiver","2599_quiver","2693_quiver","2733_quiver","1210_quiver","2935_quiver","2958_quiver","2964_quiver","3034_quiver","3068_quiver","3080_quiver","3091_quiver"],
      Y : ["850_quiver","949_quiver","1088_quiver","1125_quiver","1159_quiver","1134_quiver","1224_quiver","1369_quiver","1410_quiver","1568_quiver","1577_quiver","1619_quiver","1634_quiver","1646_quiver","1652_quiver","1665_quiver","1681_quiver","1697_quiver","1722_quiver","1766_quiver","1783_quiver","1891_quiver","1937_quiver","1963_quiver","1790_quiver","1997_quiver","2073_quiver","2113_quiver","2163_quiver","2166_quiver","5705_quiver","2245_quiver","2259_quiver","2260_quiver","2334_quiver","2340_quiver","2382_quiver","2443_quiver","2511_quiver","2534_quiver","2573_quiver","2597_quiver","2651_quiver","2707_quiver","2766_quiver","2773_quiver","2791_quiver","2830_quiver","2875_quiver","3022_quiver","3070_quiver","3074_quiver","3075_quiver","3078_quiver"]}
    ```
  * *B. siliquastri* identified by the DTOL project with the assembly
    ```python
    { X : ['X']
      Y : ['Y']}
    ``` 
  * *A. obtectus* Identified by me and Göran's project about it
    ```python
    { X : ['CAVLJG010000002.1']# ["HiC_scaffold_10", 'scaffold_21'] # double check assembly versions for correct scaffold names, i also have 'Chr_10' as a name
      Y : ['scaffold_13', 'scaffold_86']} # some have HiC and some don't but that may just be manual renaming of the largest to chromosomes
    ``` 
* **Cochinella**
  * *C. septempunctata*
    ```python
    { X : ["NC_058198.1"]
      Y : []}
    ``` 
  * *C. magnifica*
    ```python
    { X : []
      Y : []}
    ``` 
* **Tribolium**
  * *T. castaneum*
    ```python
    { X : ["NC_087403.1"]
      Y : []}
    ``` 
  * *T. freemani*
    ```python
    { X : []
      Y : []}
    ``` 

