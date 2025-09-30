#!/bin/bash -l

REPLIB=/proj/naiss2023-6-65/Milena/chapter3/repeat_libraries/Tfre
ASSEMBLY=/proj/naiss2023-6-65/Milena/chapter3/species_assemblies/Tfre_GCA_022388455.1.fasta 
# sbatch /proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/bash/run_repeatmasker.sh  $ASSEMBLY $REPLIB Tfre

ASSEMBLY=/proj/naiss2023-6-65/Milena/chapter3/species_assemblies/Cmag_GCA_965644565.1.fasta
REPLIB=/proj/naiss2023-6-65/Milena/chapter3/repeat_libraries/Cmag
sbatch /proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/bash/run_repeatmasker.sh  $ASSEMBLY $REPLIB Cmag