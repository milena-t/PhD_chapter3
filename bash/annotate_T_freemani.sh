# T. freemani annotation
TFRE_DIR=/proj/naiss2023-6-65/Milena/annotation_pipeline/Cmac_Lome_superscaffolded_comparison/Cmac_Lome_diverse
cd $TFRE_DIR
sbatch --job-name="T_freemani_annotation" --output="T_freemani_annotation.out" -t 5-00:00:00 \
/proj/naiss2023-6-65/Milena/annotation_pipeline/braker3_singularity_with_RNAseq_in_SNIC_TMP.sh Cmac_SI_diverse ${TFRE_DIR}/assembly_genomic.fna.masked \
$ORTHODB_ARTHROPODA ERR12383274,ERR12383248,ERR12383255,ERR12383289,ERR12383294,ERR12383314,ERR12383271,ERR12383290,ERR12383279,ERR12383266