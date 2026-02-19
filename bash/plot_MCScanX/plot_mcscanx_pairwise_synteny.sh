#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -n 1
#SBATCH -t 10:00
#SBATCH -J plot_MCScanX
#SBATCH -o plot_MCScanX.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

# Bioinfo tools not needed on pelle
# module load java/OpenJDK_24+36
module load Java/21.0.7

cd /proj/naiss2023-6-65/Milena/chapter2/MCScanX-1_0_0/downstream_analyses

java dual_synteny_plotter \
    -g /proj/naiss2023-6-65/Milena/chapter3/MCScanX/chrysomelidae/all_species.gff \
    -s /proj/naiss2023-6-65/Milena/chapter3/MCScanX/chrysomelidae/all_species.collinearity \
    -c /proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/bash/plot_MCScanX/plot_ctl_Cmac_Bsil.ctl \
    -o /proj/naiss2023-6-65/Milena/chapter3/MCScanX/chrysomelidae/plots/plot_ctl_Cmac_Bsil.png

echo "--> done Cmac Bsil"

java dual_synteny_plotter \
    -g /proj/naiss2023-6-65/Milena/chapter3/MCScanX/chrysomelidae/all_species.gff \
    -s /proj/naiss2023-6-65/Milena/chapter3/MCScanX/chrysomelidae/all_species.collinearity \
    -c /proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/bash/plot_MCScanX/plot_ctl_Csep_Bsil.ctl \
    -o /proj/naiss2023-6-65/Milena/chapter3/MCScanX/chrysomelidae/plots/plot_ctl_Csep_Bsil.png

echo "--> done Csep Bsil"

java dual_synteny_plotter \
    -g /proj/naiss2023-6-65/Milena/chapter3/MCScanX/chrysomelidae/all_species.gff \
    -s /proj/naiss2023-6-65/Milena/chapter3/MCScanX/chrysomelidae/all_species.collinearity \
    -c /proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/bash/plot_MCScanX/plot_ctl_Aobt_Bsil.ctl \
    -o /proj/naiss2023-6-65/Milena/chapter3/MCScanX/chrysomelidae/plots/plot_ctl_Aobt_Bsil.png

echo "--> done Aobt Bsil"

java dual_synteny_plotter \
    -g /proj/naiss2023-6-65/Milena/chapter3/MCScanX/chrysomelidae/all_species.gff \
    -s /proj/naiss2023-6-65/Milena/chapter3/MCScanX/chrysomelidae/all_species.collinearity \
    -c /proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/bash/plot_MCScanX/plot_ctl_Cmac_Aobt.ctl \
    -o /proj/naiss2023-6-65/Milena/chapter3/MCScanX/chrysomelidae/plots/plot_ctl_Cmac_Aobt.png

echo "--> done Cmac Aobt"

java dual_synteny_plotter \
    -g /proj/naiss2023-6-65/Milena/chapter3/MCScanX/chrysomelidae/all_species.gff \
    -s /proj/naiss2023-6-65/Milena/chapter3/MCScanX/chrysomelidae/all_species.collinearity \
    -c /proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/bash/plot_MCScanX/plot_ctl_Dsub_Dcar.ctl \
    -o /proj/naiss2023-6-65/Milena/chapter3/MCScanX/chrysomelidae/plots/plot_ctl_Dsub_Dcar.png

echo "--> done Dsub Dcar"

java dual_synteny_plotter \
    -g /proj/naiss2023-6-65/Milena/chapter3/MCScanX/chrysomelidae/all_species.gff \
    -s /proj/naiss2023-6-65/Milena/chapter3/MCScanX/chrysomelidae/all_species.collinearity \
    -c /proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/bash/plot_MCScanX/plot_ctl_Dsub_Aobt.ctl \
    -o /proj/naiss2023-6-65/Milena/chapter3/MCScanX/chrysomelidae/plots/plot_ctl_Dsub_Aobt.png

echo "--> done Dsub Aobt"

