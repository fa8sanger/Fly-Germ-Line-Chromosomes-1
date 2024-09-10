 #!/bin/bash
 
 #BSUB -o logs/braker3.out.%J
 #BSUB -e logs/braker3.err.%J
 #BSUB -q long
 #BSUB -n 32
 #BSUB -M40960
 #BSUB -R "select[mem>40960] rusage[mem=40960]"
 
 module load ISG/singularity/3.11.4
 export BRAKER_SIF=/software/team360/src/braker3.sif
 singularity exec -B ${PWD}:${PWD} ${BRAKER_SIF} braker.pl \
 --genome=./data/bcop_braker3_core/core.fa \
 --bam=./data/bcop_braker3_core/All_rnaseq.bam \
 --softmasking \
 --workingdir=./data/bcop_braker3_core/braker3/bcop \
 --threads 32 \
 --species=bcop_core_b3 \
 --gff3 \
 --prot_seq=./data/bcop_braker3_core/Arthropoda+contarina.simple_header.fa \
 --useexisting


