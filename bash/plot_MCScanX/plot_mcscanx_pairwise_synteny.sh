#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -n 1
#SBATCH -t 10:00
#SBATCH -J plot_MCScanX
#SBATCH -o plot_MCScanX.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

# Bioinfo tools not needed on pelle
module load bioinfo-tools java/OpenJDK_24+36

cd /proj/naiss2023-6-65/Milena/chapter2/MCScanX-1_0_0/downstream_analyses

java dual_synteny_plotter \
    -g /proj/naiss2023-6-65/Milena/chapter3/MCScanX/all_species_old/all_species.gff \
    -s /proj/naiss2023-6-65/Milena/chapter3/MCScanX/all_species_old/all_species.collinearity \
    -c /proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/bash/plot_MCScanX/plot_ctl_Tcas_Tfre.ctl \
    -o /proj/naiss2023-6-65/Milena/chapter3/MCScanX/plots/plot_ctl_Tcas_Tfre.png

echo "--> done Tcas Tfre"

java dual_synteny_plotter \
    -g /proj/naiss2023-6-65/Milena/chapter3/MCScanX/all_species_old/all_species.gff \
    -s /proj/naiss2023-6-65/Milena/chapter3/MCScanX/all_species_old/all_species.collinearity \
    -c /proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/bash/plot_MCScanX/plot_ctl_Csep_Cmag.ctl \
    -o /proj/naiss2023-6-65/Milena/chapter3/MCScanX/plots/plot_ctl_Csep_Cmag.png

echo "--> done Csep Cmag"


java dual_synteny_plotter \
    -g /proj/naiss2023-6-65/Milena/chapter3/MCScanX/all_species_old/all_species.gff \
    -s /proj/naiss2023-6-65/Milena/chapter3/MCScanX/all_species_old/all_species.collinearity \
    -c /proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/bash/plot_MCScanX/plot_ctl_Tcas_Bsil.ctl \
    -o /proj/naiss2023-6-65/Milena/chapter3/MCScanX/plots/plot_ctl_Tcas_Bsil.png

echo "--> done Tcas Bsil"

java dual_synteny_plotter \
    -g /proj/naiss2023-6-65/Milena/chapter3/MCScanX/all_species_old/all_species.gff \
    -s /proj/naiss2023-6-65/Milena/chapter3/MCScanX/all_species_old/all_species.collinearity \
    -c /proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/bash/plot_MCScanX/plot_ctl_Cmac_Bsil.ctl \
    -o /proj/naiss2023-6-65/Milena/chapter3/MCScanX/plots/plot_ctl_Cmac_Bsil.png

echo "--> done Cmac Bsil"

java dual_synteny_plotter \
    -g /proj/naiss2023-6-65/Milena/chapter3/MCScanX/all_species_old/all_species.gff \
    -s /proj/naiss2023-6-65/Milena/chapter3/MCScanX/all_species_old/all_species.collinearity \
    -c /proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/bash/plot_MCScanX/plot_ctl_Aobt_Bsil.ctl \
    -o /proj/naiss2023-6-65/Milena/chapter3/MCScanX/plots/plot_ctl_Aobt_Bsil.png

echo "--> done Aobt Bsil"

java dual_synteny_plotter \
    -g /proj/naiss2023-6-65/Milena/chapter3/MCScanX/all_species_old/all_species.gff \
    -s /proj/naiss2023-6-65/Milena/chapter3/MCScanX/all_species_old/all_species.collinearity \
    -c /proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/bash/plot_MCScanX/plot_ctl_Cmac_Aobt.ctl \
    -o /proj/naiss2023-6-65/Milena/chapter3/MCScanX/plots/plot_ctl_Cmac_Aobt.png

echo "--> done Cmac Aobt"