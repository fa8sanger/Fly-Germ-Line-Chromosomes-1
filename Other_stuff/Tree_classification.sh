module add ISG/python/3.12.3
python3.12 -m pip install Bio
cd /lustre/scratch126/tol/teams/jaron/users/fede/TreeClassifications
# python3.12 ../software/treefile2table_of_neighbors.py . |& less
cd /lustre/scratch126/tol/teams/jaron/users/fede/TreeClassifications
for tree in *contree
 do echo $tree
 python3.12 ../software/treefile2table_of_neighbors.py $tree >& $tree.OUT
 done






grep "families are not monophyletic" *OUT | wc -l
grep "Outgroup (dmel) absent" *OUT | wc -l
grep "No scia" *OUT | wc -l
grep "No cecidomyiidae" *OUT | wc -l

581 not monophyletic
4328 without outgroup
245 without either sciaridae or cecidomyiidae


grep -v "families are not monophyletic" *OUT | grep -v "Outgroup" | grep -v "not mono" | grep -v branch | grep -v "No ceci" | grep -v "No sci" | grep -v "NA" > ALL_CLASSIFICATIONS
 4813 genes classified, from:
 
   2488 bcop
   1007 bimp
   1306 ling

   1827 bcop	cecidomyiidae
     10 bcop	other
    651 bcop	sciaridae
    220 bimp	cecidomyiidae
    787 bimp	sciaridae
   1008 ling	cecidomyiidae
      4 ling	other
    294 ling	sciaridae

	  origin:cecid.   origin:sciarid.   other
bcop          1827               651       10
bimp           220               787        0
ling          1008               294        4 


# With the three coming from cecido:
OG0001863: FBpp0076239, https://www.uniprot.org/uniprotkb/M9PBI0/entry, CRISP-related, Cysteine-rich secretory protein-related
OG0002028: FBpp0081343, PNKP-PA, Polynucleotide kinase 3 phosphatase (PNK3P), DNA repair
OG0002412: FBpp0071971, Dcp-1, Caspase 1
OG0002438: FBpp0085885, CG5323-PA, PROTEIN FAM136A
OG0002531: FBpp0112540, Dmel\l(3)80Fg, DnaJ domain, Thioredoxin domain
OG0002754: FBpp0077874, BNIP3, BCL2/adenovirus E1B 19 kDa protein-interacting , Apoptosis-inducing protein that can overcome BCL2 suppression
OG0002755: FBpp0077149, CG3652, Yip1 domain, Golgi protein involved in vesicular transport
OG0002837: FBpp0080702, brat, brain tumour protein, translational repressor to inhibit cell proliferation. Look out for pum
OG0002865: FBpp0080689, Acn, acinus, RNA splicing, APOPTOTIC CHROMATIN CONDENSATION INDUCER IN THE NUCLEUS 1 
OG0003249: FBpp0072957, RpL28, Large ribosomal subunit protein 
OG0003273: FBpp0088438, Aprt, Phosphoribosyltransferase domain
OG0003341: 




OG0003273.final.aa.fa.contree.OUT:OG0003273.final.aa.fa.contree OG0003341.final.aa.fa.contree.OUT:OG0003341.final.aa.fa.contree 
                                                              3                                                               3 
OG0003353.final.aa.fa.contree.OUT:OG0003353.final.aa.fa.contree OG0003631.final.aa.fa.contree.OUT:OG0003631.final.aa.fa.contree 
                                                              3                                                               3 
OG0003666.final.aa.fa.contree.OUT:OG0003666.final.aa.fa.contree OG0003667.final.aa.fa.contree.OUT:OG0003667.final.aa.fa.contree 
                                                              3                                                               3 
OG0003785.final.aa.fa.contree.OUT:OG0003785.final.aa.fa.contree OG0003786.final.aa.fa.contree.OUT:OG0003786.final.aa.fa.contree 
                                                              3                                                               3 
OG0003972.final.aa.fa.contree.OUT:OG0003972.final.aa.fa.contree OG0004143.final.aa.fa.contree.OUT:OG0004143.final.aa.fa.contree 
                                                              3                                                               3 
OG0004373.final.aa.fa.contree.OUT:OG0004373.final.aa.fa.contree OG0004465.final.aa.fa.contree.OUT:OG0004465.final.aa.fa.contree 
                                                              3                                                               3 
OG0004466.final.aa.fa.contree.OUT:OG0004466.final.aa.fa.contree OG0004467.final.aa.fa.contree.OUT:OG0004467.final.aa.fa.contree 
                                                              3                                                               3 
OG0004470.final.aa.fa.contree.OUT:OG0004470.final.aa.fa.contree 
                                                              3 
