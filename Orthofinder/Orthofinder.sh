# Before running orthofinder I should keep only one transcript per gene (one isoform, I mean)
# I would need different scripts for each species

# scp -p one_isoform_per_gene_* farm22:/lustre/scratch125/casm/team268im/fa8/117/fede/Blast/
# cdt
# cd Orthofinder
# mkdir Proteomes_OneIsoform
# cd Proteomes
# perl one_isoform_per_gene_braker.pl     aaphi.fa       > ../Proteomes_OneIsoform/aaphi.fa
# perl one_isoform_per_gene_braker.pl     orobi.fa       > ../Proteomes_OneIsoform/orobi.fa
# perl one_isoform_per_gene_braker.pl     ling_grc.fa    > ../Proteomes_OneIsoform/ling_grc.fa
# perl one_isoform_per_gene_braker.pl     ling_core.fa   > ../Proteomes_OneIsoform/ling_core.fa
# perl one_isoform_per_gene_braker.pl     bcop_grc.fa    > ../Proteomes_OneIsoform/bcop_grc.fa
# perl one_isoform_per_gene_braker.pl     bcop_core.fa   > ../Proteomes_OneIsoform/bcop_core.fa
# perl one_isoform_per_gene_contarinia.pl contarinia.fa  > ../Proteomes_OneIsoform/contarinia.fa
# perl one_isoform_per_gene_dmel.pl       dmel.fa        > ../Proteomes_OneIsoform/dmel.fa


cdt
cd Orthofinder
mkdir Proteomes_2nd_round
cd Proteomes_2nd_round
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bcop/Braker2_grc/braker.aa bcop_grc.fa
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bcop/data/bcop_braker3_core/braker3/bcop/braker.aa bcop_core.fa
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Orobi/braker.aa orobi.fa
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Ling/Braker2_grc/braker.aa ling_grc.fa
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Ling/data/lycoriella_braker3_core/braker3/lycoriella/braker.aa ling_core.fa
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bimp/Braker2_core/braker.aa bimp_core.fa
# cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bimp/Braker2_grc/braker.aa # This failed
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bimp/data/bimp_braker3_core/braker3/bimp/braker.aa bimp_grc.fa
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Aaphi/braker.aa aaphi.fa
cp -p ../Proteomes/contarinia.fa .
cp -p ../Proteomes/dmel.fa .
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/GCA_029228625.1_FCFRP_Bhyg_1.0_translated_cds.faa phyg.fa

cdt
cd Orthofinder
mkdir Proteomes_2nd_round_OneIsoform
cd Proteomes_2nd_round
perl ../one_isoform_per_gene_braker.pl     aaphi.fa       > ../Proteomes_2nd_round_OneIsoform/aaphi.fa
perl ../one_isoform_per_gene_braker.pl     bimp_grc.fa    > ../Proteomes_2nd_round_OneIsoform/bimp_grc.fa
perl ../one_isoform_per_gene_braker.pl     bimp_core.fa   > ../Proteomes_2nd_round_OneIsoform/bimp_core.fa
perl ../one_isoform_per_gene_braker.pl     orobi.fa       > ../Proteomes_2nd_round_OneIsoform/orobi.fa
perl ../one_isoform_per_gene_braker.pl     ling_grc.fa    > ../Proteomes_2nd_round_OneIsoform/ling_grc.fa
perl ../one_isoform_per_gene_braker.pl     ling_core.fa   > ../Proteomes_2nd_round_OneIsoform/ling_core.fa
perl ../one_isoform_per_gene_braker.pl     bcop_grc.fa    > ../Proteomes_2nd_round_OneIsoform/bcop_grc.fa
perl ../one_isoform_per_gene_braker.pl     bcop_core.fa   > ../Proteomes_2nd_round_OneIsoform/bcop_core.fa
perl ../one_isoform_per_gene_contarinia.pl contarinia.fa  > ../Proteomes_2nd_round_OneIsoform/contarinia.fa
perl ../one_isoform_per_gene_dmel.pl       dmel.fa        > ../Proteomes_2nd_round_OneIsoform/dmel.fa
cp -p phyg.fa                                                         ../Proteomes_2nd_round_OneIsoform/phyg.fa

## phyg_gff=/lustre/scratch125/casm/team268im/fa8/117/fede/Braker_2nd_round/GCA_029228625.1_FCFRP_Bhyg_1.0_genomic.gff.gz
# For Phyg there is only one mRNA per gene, I think they don't have multiple isoforms
## cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker_2nd_round
## zcat GCA_029228625.1_FCFRP_Bhyg_1.0_genomic.gff.gz | cut -f3 |grep -v "^#" | sort | uniq -c | sort -n
##    6452 region
##   17881 gene
##   17881 mRNA
##   94850 CDS
##   96787 exon
## No need to do anything


##########################################################################################
# Orthofinder
##cdt
##cd Braker
##cp -p ./Orobi/braker.aa ../Orthofinder/Proteomes/orobi.fa
##cp -p ./Ling/data/lycoriella_braker3_core/braker3/lycoriella/braker.aa ../Orthofinder/Proteomes/ling_core.fa
##cp -p ./Ling/Braker2_grc/braker.aa ../Orthofinder/Proteomes/ling_grc.fa
##cp -p ./Aaphi/braker.aa ../Orthofinder/Proteomes/aaphi.fa
##cp -p ../Blast/contarinia.fa ../Orthofinder/Proteomes/contarinia.fa
##cp -p ../Blast/dmel-all-translation-r6.58.fasta ../Orthofinder/Proteomes/dmel.fa

## cdt
## cd Orthofinder
## mkdir Proteomes_OneIsoform
## cd Proteomes
## perl ../one_isoform_per_gene_braker.pl     aaphi.fa       > ../Proteomes_OneIsoform/aaphi.fa
## perl ../one_isoform_per_gene_braker.pl     orobi.fa       > ../Proteomes_OneIsoform/orobi.fa
## perl ../one_isoform_per_gene_braker.pl     ling_grc.fa    > ../Proteomes_OneIsoform/ling_grc.fa
## perl ../one_isoform_per_gene_braker.pl     ling_core.fa   > ../Proteomes_OneIsoform/ling_core.fa
## perl ../one_isoform_per_gene_braker.pl     bcop_grc.fa    > ../Proteomes_OneIsoform/bcop_grc.fa
## perl ../one_isoform_per_gene_braker.pl     bcop_core.fa   > ../Proteomes_OneIsoform/bcop_core.fa
## perl ../one_isoform_per_gene_contarinia.pl contarinia.fa  > ../Proteomes_OneIsoform/contarinia.fa
## perl ../one_isoform_per_gene_dmel.pl       dmel.fa        > ../Proteomes_OneIsoform/dmel.fa

# module add ISG/python/3.12.3
# module add mafft/7.525--h031d066_1
# cpus=16
# bsub20000 -q normal -e ortho2.err -o ortho2.err -n $cpus -R span[ptile=$cpus]  /lustre/scratch126/tol/teams/jaron/users/fede/software/OrthoFinder/orthofinder -f Proteomes_2nd_round_OneIsoform -a $cpus -t $cpus

# With MSA!
cdt
cd Orthofinder
module add ISG/python/3.12.3
module add mafft/7.525--h031d066_1
module add fasttree/2.1.11--h031d066_3
cpus=16
LSB_DEFAULT_JOBGROUP=team360
bsub60000 -q long -e ortho3.err -o ortho3.err -n $cpus -R span[ptile=$cpus]  /lustre/scratch126/tol/teams/jaron/users/fede/software/OrthoFinder/orthofinder -f Proteomes_2nd_round_OneIsoform -a $cpus -t $cpus -M msa


cut -f2,3,4,5,6,7,8,9,10,11,12 Orthogroups/Orthogroups.GeneCount.tsv | sort | uniq -c | sort -n

# En todos y 1 sola vez:
awk '$13==11' Orthogroups/Orthogroups.GeneCount.tsv |grep -v -P "\t0\t" | cut -f1 
cd /lustre/scratch126/tol/teams/jaron/users/fede/Orthofinder/Proteomes_2nd_round_OneIsoform/OrthoFinder/Results_Jul25_1
for og in OG0003564 OG0003592 OG0003593 OG0003598 OG0003603 OG0003616 OG0003624 OG0003630 OG0003631 OG0003632 OG0003658 OG0003660 OG0003661 OG0003663 OG0003666 OG0003667 OG0003670 OG0003693 OG0003694 OG0003712 OG0003713 OG0003730 OG0003744 OG0003745 OG0003748 OG0003749 OG0003750 OG0003763 OG0003784 OG0003785 OG0003859 OG0003904 OG0003921 OG0003933 OG0003954 OG0003964 OG0003965 OG0003972 OG0003978 OG0003986 OG0003987 OG0003988 OG0004001 OG0004027 OG0004039 OG0004042 OG0004058 OG0004059 OG0004072 OG0004079 OG0004086 OG0004143 OG0004144 OG0004214 OG0004247 OG0004260 OG0004269 OG0004278 OG0004279 OG0004303 OG0004305 OG0004327 OG0004332 OG0004339 OG0004345 OG0004349 OG0004355 OG0004361 OG0004368 OG0004373 OG0004377 OG0004394 OG0004404 OG0004408 OG0004409 OG0004412 OG0004414 OG0004417 OG0004419 OG0004427 OG0004434 OG0004435 OG0004436 OG0004440 OG0004445 OG0004446 OG0004447 OG0004449 OG0004450 OG0004452 OG0004453 OG0004454 OG0004465 OG0004466 OG0004467 OG0004468 OG0004469 OG0004485 OG0004488 OG0004498 OG0004505 OG0004518 OG0004519 OG0004520 OG0004552 OG0004571 OG0004648 
do echo "" `cat Gene_Trees/"$og"_tree.txt`;
done



# Orthogroups with 1 from each except for contarinia with two:
cd /lustre/scratch126/tol/teams/jaron/users/fede/Orthofinder/Proteomes_OneIsoform/OrthoFinder/Results_Jul17/Orthogroups
awk '$2==1&&$3==2&&$4==1&&$5==1&&$6==1&&$7==1' Orthogroups.GeneCount.tsv


for a in OG0003688 OG0003698 OG0003707 OG0003727 OG0003728 OG0003732 OG0003733 OG0003742 OG0003754 OG0003760 OG0003762 OG0003781 OG0003782 OG0003783 OG0003784 OG0003785 OG0003788 OG0003792 OG0003794 OG0003795 OG0003797 OG0003798 OG0003804 OG0003807 OG0003808 OG0003813 OG0003819 OG0003821 OG0003822 OG0003823 OG0003824 OG0003827 OG0003828 OG0003829 OG0003832 OG0003841 OG0003857 OG0003858 OG0003859 OG0003860 OG0003861 OG0003863 OG0003864 OG0003866 OG0003871 OG0003876 OG0003877 OG0003878 OG0003879 OG0003881 OG0003882 OG0003888 OG0003889 OG0003890 OG0003900 OG0003908 OG0003916 OG0003917 OG0003918 OG0003919 OG0003920 OG0003922 OG0003924 OG0003925 OG0003930 OG0003931 OG0003932 OG0003933 OG0003934 OG0003935 OG0003936 OG0003938 OG0003941 OG0003944 OG0003945 OG0003947 OG0003950 OG0003978 OG0003981 OG0003982 OG0004003 OG0004045 OG0004054 OG0004056 OG0004062 OG0004063 OG0004080 OG0004095 OG0004098 OG0004099 OG0004100 OG0004120 OG0004127 OG0004129 OG0004135 OG0004138 OG0004142 OG0004143 OG0004159 OG0004170 OG0004206 OG0004220 OG0004222 OG0004225 OG0004239 OG0004240 OG0004244 OG0004245 OG0004246 OG0004247 OG0004248 OG0004251 OG0004252 OG0004253 OG0004255 OG0004257 OG0004258 OG0004259 OG0004261 OG0004263 OG0004264 OG0004266 OG0004269 OG0004270 OG0004272 OG0004273 OG0004276 OG0004277 OG0004278 OG0004279 OG0004280 OG0004281 OG0004285 OG0004286 OG0004292 OG0004293 OG0004295 OG0004297 OG0004299 OG0004306 OG0004310 OG0004311 OG0004312 OG0004322 OG0004330 OG0004332 OG0004333 OG0004338 OG0004341 OG0004352 OG0004359 OG0004362 OG0004366 OG0004368 OG0004369 OG0004374 OG0004378 OG0004386 OG0004391 OG0004392 OG0004394 OG0004395 OG0004397 OG0004401 OG0004402 OG0004403 OG0004405 OG0004406 OG0004407 OG0004408 OG0004409 OG0004410 OG0004411 OG0004412 OG0004413 OG0004414 OG0004415 OG0004416 OG0004417 OG0004418 OG0004420 OG0004422 OG0004423 OG0004426 OG0004427 OG0004428 OG0004431 OG0004432 OG0004433 OG0004434 OG0004435 OG0004440 OG0004441 OG0004442 OG0004445 OG0004448 OG0004449 OG0004451 OG0004452 OG0004454 OG0004455 OG0004456 OG0004458 OG0004459 OG0004465 OG0004466 OG0004467 OG0004468 OG0004470 OG0004482 OG0004495 OG0004496 OG0004503 OG0004504 OG0004505 OG0004507 OG0004508 OG0004517 OG0004521 OG0004522 OG0004523 OG0004524 OG0004526 OG0004527 OG0004536 OG0004537 OG0004541 OG0004542 OG0004543 OG0004544 OG0004545 OG0004546 OG0004547 OG0004548 OG0004558 OG0004559 OG0004561 OG0004562 OG0004568 OG0004574 OG0004575 OG0004579 OG0004583 OG0004584 OG0004585 OG0004586 OG0004587 OG0004588 OG0004589 OG0004590 OG0004591 OG0004592 OG0004593 OG0004594 OG0004595 OG0004596 OG0004597 OG0004598 OG0004599 OG0004600 OG0004601 OG0004602 OG0004604 OG0004605 OG0004606 OG0004608 OG0004609 OG0004618 OG0004625 OG0004626 OG0004636 OG0004647 OG0004650 OG0004662 OG0004663 OG0004664 OG0004666 OG0004669 OG0004670 OG0004671 OG0004672 OG0004673 OG0004674 OG0004675 OG0004678 OG0004683 OG0004696 OG0004700 OG0004701 OG0004708 OG0004717 OG0004719 OG0004720 OG0004721 OG0004722 OG0004723 OG0004724 OG0004725 OG0004726 OG0004727 OG0004728 OG0004730 OG0004732 OG0004733 OG0004734 OG0004738 OG0004739 OG0004740 OG0004741 OG0004742 OG0004743 OG0004744 OG0004745 OG0004748 OG0004749 OG0004753 OG0004754 OG0004755 OG0004757 OG0004758 OG0004759 OG0004760 OG0004762 OG0004763 OG0004764 OG0004765 OG0004769 OG0004770 OG0004771 OG0004773 OG0004774 OG0004777 OG0004778 OG0004779 OG0004780 OG0004781 OG0004782 OG0004783 OG0004784 OG0004786 OG0004787 OG0004788 OG0004789 OG0004798 OG0004807 OG0004814 OG0004817 OG0004820 OG0004821 OG0004832 OG0004833 OG0004834 OG0004835 OG0004840 OG0004841 OG0004842 OG0004843 OG0004844 OG0004845 OG0004847 OG0004848 OG0004852 OG0004859 OG0004860 OG0004863 OG0004864 OG0004870 OG0004871 OG0004873 OG0004877 OG0004878 OG0004879 OG0004881 OG0004882 OG0004883 OG0004884 OG0004890 OG0004891 OG0004895 OG0004897 OG0004899 OG0004902 OG0004903 OG0004904 OG0004905 OG0004906 OG0004907 OG0004908 OG0004911 OG0004912 OG0004914 OG0004917 OG0004918 OG0004919 OG0004920 OG0004921 OG0004923 OG0004924 OG0004925 OG0004931 OG0004933 OG0004942 OG0004943 OG0004945 OG0004946 OG0004947 OG0004949 OG0004950 OG0004951 OG0004952 OG0004953 OG0004955 OG0004959 OG0004966 OG0004967 OG0004968 OG0004969 OG0004970 OG0004971 OG0004972 OG0004974 OG0004975 OG0004976 OG0004977 OG0004978 OG0004979 OG0004980 OG0004981 OG0004982 OG0004986 OG0004987 OG0004988 OG0004989 OG0004990 OG0004991 OG0004992 OG0004993 OG0004994 OG0004995 OG0004996 OG0004997 OG0004998 OG0005000 OG0005001 OG0005002 OG0005003 OG0005004 OG0005005 OG0005006 OG0005007 OG0005008 OG0005015 OG0005016 OG0005018 OG0005019 OG0005021 OG0005022 OG0005028 OG0005029 OG0005031 OG0005032 OG0005033 OG0005034 OG0005035 OG0005036 OG0005037 OG0005039 OG0005040 OG0005041 OG0005044 OG0005045 OG0005046 OG0005047 OG0005048 OG0005049 OG0005050 OG0005052 OG0005054 OG0005055 OG0005056 OG0005057 OG0005059 OG0005060 OG0005061 OG0005062 OG0005063 OG0005064 OG0005066 OG0005068 OG0005069 OG0005070 OG0005072 OG0005073 OG0005074 OG0005075 OG0005076 OG0005077 OG0005079 OG0005080 OG0005081 OG0005082 OG0005084 OG0005085 OG0005087 OG0005088 OG0005089 OG0005090 OG0005091 OG0005092 OG0005093 OG0005095 OG0005096 OG0005097 OG0005098 OG0005100 OG0005107 OG0005108 OG0005109 OG0005110 OG0005115 OG0005118 OG0005119 OG0005120 OG0005135 OG0005142 OG0005143 OG0005145 OG0005148 OG0005150 OG0005151 OG0005156 OG0005178 OG0005184 OG0005189 OG0005192 OG0005193 OG0005197 OG0005206 OG0005207 OG0005208 OG0005220 OG0005223 OG0005228 OG0005230 OG0005236 OG0005244 OG0005245 OG0005248 OG0005250 OG0005277 OG0005279 OG0005285 OG0005287 OG0005298 OG0005307 OG0005319 OG0005324 OG0005336 OG0005346 OG0005357 OG0005358 OG0005359 OG0005360 OG0005366 OG0005370 OG0005376 OG0005391 OG0005405 OG0005413 OG0005414 OG0005417 OG0005421 OG0005422 OG0005425
do echo $a
 echo `cat Resolved_Gene_Trees/"$a"_tree.txt` >> ALL_1per.tree
 done
