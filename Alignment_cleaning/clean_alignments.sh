
og=$1
out_path=/lustre/scratch126/tol/teams/jaron/users/fede/CDS_seqs/CleanAlignments/
seqs_path=/lustre/scratch126/tol/teams/jaron/users/fede/CDS_seqs/CDS4alis/

mkdir $out_path/$og
cp -p $seqs_path/$og.fa $out_path/$og/
cd $out_path/$og
echo ""
echo "------------------------------------------------------------------------------------"
echo "Running Prequal...."
echo "------------------------------------------------------------------------------------"
/lustre/scratch126/tol/teams/jaron/users/fede/software/prequal/prequal $og.fa
echo "------------------------------------------------------------------------------------"
echo "Running MACSE...."
echo "------------------------------------------------------------------------------------"
java -Xmx3880M -jar /lustre/scratch126/tol/teams/jaron/users/fede/software/macse_v2.07.jar -prog alignSequences -seq $og.fa.dna.filtered 
echo "------------------------------------------------------------------------------------"
echo "Running BMGE...."
echo "------------------------------------------------------------------------------------"
java -Xmx3880M -jar /lustre/scratch126/tol/teams/jaron/users/fede/software/BMGE/src/BMGE.jar -i $og.fa.dna_NT.filtered -t CO -m BLOSUM80    -g 0.1 -e 0.3 -oaa $og.final.aa.fa -oco $og.final.codon.fa -onco12 $og.final.nex -b 5
echo "------------------------------------------------------------------------------------"
echo "Done!"
echo "------------------------------------------------------------------------------------"
