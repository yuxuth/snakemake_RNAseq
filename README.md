# snakemake_se_RNAseq

This pepeline is inspired by [crazyhottommy](https://github.com/crazyhottommy/pyflow-RNAseq).

snakemake -j -np  99 --cluster '/Snakefile-sbatch.py'

generate workflow plot

snakemake --dag 2> /dev/null | dot -T png > ../../img/workflow_bysample.png
