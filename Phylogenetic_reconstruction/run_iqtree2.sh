module add iqtree/2.3.4--h21ec9f0_0
cd /lustre/scratch126/tol/teams/jaron/users/fede/CDS_seqs/CleanAlignments/
mkdir LOGS
for og in OG*
 do echo $og
 bsub5000 -e LOGS/$og.err -o LOGS/$og.err iqtree2 -s $og/$og.final.aa.fa --seed 1234 -b 100 -redo
 done






