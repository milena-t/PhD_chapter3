# Chapter III: FastX in coleoptera

## Other code

I wrote a lot of code for this a year ago here: https://github.com/milena-t/calculate_orthogroup_dNdS

## Notes

### Methods

* [Martinez-Pacheco 2020](https://academic.oup.com/gbe/article/12/11/2015/5892261) 
  * double check y linked gene presence by `blastn` against the assembly.
  * *"the best BlastN match (usually around 92–95% identity over the entire sequence) onto the annotated X chromosome of the reference genomes was considered the X gametologs"*
  * [Marques 2005](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0030357) for identifying retrogenes. blast proteins against assembly, merge nearby matches, *"query and target sequences had >50% similarity on the amino acid level and over >80% of their length* \[are\] *shared"*, verify absence of introns. with some `paml` stuff, they identified the ancestral gametolog that all retrogenes originate from (useful for FastX?)
* **[Whittle 2020](https://academic.oup.com/g3journal/article/10/3/1125/6026234): no FastX in beetles**
* [Mank 2009](https://academic.oup.com/mbe/article/27/3/661/1000994?login=true) Faster-Z in birds is mainly due to drift
  * Positive selection would be fixation of recessive male-biased mutations
* [Li 2010](https://pubmed.ncbi.nlm.nih.gov/21035095/) FastZ in duplicates compared to autosomal duplicates
  * within-species comparison, make all pairwise dN/dS of all genes within a gene family (check methods specifically that they use to reduce between-sample depencence for statistical power)
  
### Expectations

* Retrotransposition of male-biased genes from X to A ([Ellegren 2011](https://www.nature.com/articles/nrg2948.pdf))
* Ampliconic regions on both X and Y, expansion of intergenic regions

### Results

* TODO read: [Vicoso & Bachtrog 2006](https://www.nature.com/articles/nrg1914) Big early paper that everyone cites about the different evolutionary forces that act differently on the X vs the autosomes
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
