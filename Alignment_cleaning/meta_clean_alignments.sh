cd /lustre/scratch126/tol/teams/jaron/users/fede/CDS_seqs/CDS4alis
mkdir LOGS
for og in OG*fa
 do echo $og
 NAME=`echo $og | cut -d'.' -f1`
 #echo $NAME
 bsub5000 -e LOGS/$og.err2 -o LOGS/$og.err2 bash ../clean_alignments.sh $NAME
done





