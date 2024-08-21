#!/usr/bin/env python3

import os
from collections import defaultdict
# import numpy as np
from Bio import Phylo
import sys

# species in the analysis
outgroup = 'dmel'
possible_grcs = ['bcop_grc', 'bimp_grc', 'ling_grc']
sciarid_core = ['bcop_core', 'bimp_core', 'ling_core', 'phyg']
cecidomyiidae = ['aaphi', 'contarinia', 'orobi']
possible_Sciaridae = set(['bcop', 'ling', 'phyg','bimp'])

# l_string = 'L-sciara_coprophila'
# a_string = 'A-sciara_coprophila'
# na_string = 'NA-sciara_coprophila'
# 'A-sciara_coprophila'
# sciaridae = set(['phytosciara_flavipes', 'trichosia_splendens'])
# cecidomyiidae = set(['mayetiola_destructor', 'porricondyla_nigripennis', 'catotricha_subobsoleta', 'lestremia_cinerea'])

def is_monophyletic_sciaridae(clade):
    return(all([tip.name.split('_')[0] in possible_Sciaridae for tip in clade]))

def is_monophyletic_cecidomyiidae(clade):
    return(all([tip.name.split('_')[0] in cecidomyiidae or 'grc' in tip.name for tip in clade]))

def tip2sp_name(tip):
    return(tip.name.split('_')[0])

def tip2gene(tip):
    return("_".join(tip.name.split('_')[-1:]).rstrip("'"))

def tree2assigments(input_newick):
    tree_name = input_newick.split('/')[-1:][0].split('_')[0]
    tree = Phylo.read(input_newick, "newick")

    sp2node_name = defaultdict(list)
    for tip in tree.get_terminals():
        sp = tip.name.split('_')[0]
        sp2node_name[sp].append(tip)

    # This commented out bit is testing for monophyly of all sciarid and all cecidomyiidae branches
    # print('parsing: ' + input_newick)
    # sp2node_name = defaultdict(list)
    # present_Sciaridae = set()
    # sciaridae_tips = []
    # cecidomyiidae_tips = []
    # for tip in tree.get_terminals():
    #     sp = tip.name.split('_')[0]
    #     is_grc = 'grc' in tip.name
    #     # sp = "_".join(tip.name.split('_')[:2]).rstrip("'")
    #     sp2node_name[sp].append(tip)
    #     if sp in cecidomyiidae:
    #         cecidomyiidae_tips.append(tip)

    #     if sp in possible_Sciaridae and not is_grc:
    #         # defining which sciaridae are present in the treefile
    #         present_Sciaridae.add(sp)
    #         sciaridae_tips.append(tip)

    # if not all([tip.name.split('_')[0] in possible_Sciaridae for tip in tree.common_ancestor(sciaridae_tips).get_terminals()]):
    #     print("Non-monophyletic sciaridae") # this tests if all the non-GRC Sciarid genes are monophyletic, while disregarding GRCs for the purpose of the test (which is the reason why we had to use this oneline instead of is_monophyletic)

    # if not all([tip.name.split('_')[0] in cecidomyiidae or 'grc' in tip.name for tip in tree.common_ancestor(cecidomyiidae_tips).get_terminals()]):
    #     print("Non-monophyletic cecidomyiidae") # the same test for cecidomyiidae, but explicitly filtering grc genes (because they are not cecidomyiidae);

    if not sp2node_name['dmel']:
        sys.stderr.write("Outgroup (Dmel) absent\n")
        sys.stdout.write("\t".join([tree_name, 'NA', 'NA', 'other', 'NA']) + '\n')
        return(0)


    tree.root_with_outgroup(sp2node_name['dmel'][0])
    
    if tree.root.clades[0] != sp2node_name['dmel'][0]:
        gnat_subtree = tree.root.clades[0]
    else:
        gnat_subtree = tree.root.clades[1]

    gnats1 = gnat_subtree.clades[0].get_terminals()
    gnats2 = gnat_subtree.clades[1].get_terminals()

    sciaridae_subtree_identified = False
    cecidomyiidae_subtree_identified = False
    if is_monophyletic_sciaridae(gnats1):
        sciaridae_subtree = gnats1
        sciaridae_subtree_identified = True
    elif is_monophyletic_cecidomyiidae(gnats1):
        cecidomyiidae_subtree = gnats1
        cecidomyiidae_subtree_identified = True
    else:
        sys.stderr.write('one gnat subclade is not monophyletic\n')
        sys.stdout.write("\t".join([tree_name, 'NA', 'NA', 'other', 'NA']) + '\n')
        return(0)

    if is_monophyletic_sciaridae(gnats2) and not sciaridae_subtree_identified:
        sciaridae_subtree = gnats2
        sciaridae_subtree_identified = True
    elif is_monophyletic_cecidomyiidae(gnats2) and not cecidomyiidae_subtree_identified:
        cecidomyiidae_subtree = gnats2
        cecidomyiidae_subtree_identified = True
    else:
        sys.stderr.write('one gnat subclade is not monophyletic\n')
        sys.stdout.write("\t".join([tree_name, 'NA', 'NA', 'other', 'NA']) + '\n')
        return(0)
    
    if not (cecidomyiidae_subtree_identified and sciaridae_subtree_identified):
        sys.stderr.write('sciaridae or cecidomyiidae are not monophyletic\n')
        classification = 'other'
        sys.stdout.write("\t".join([tree_name, 'NA', 'NA', 'Other', 'NA']) + '\n')
        # for tip in tree.get_terminals():
        #     print

    for sciarid_tip in sciaridae_subtree:
        if 'grc' in sciarid_tip.name:
            sp = tip2sp_name(sciarid_tip)
            classification = 'sciaridae'
            gene = tip2gene(sciarid_tip)
            branch_length = str(sciarid_tip.branch_length)
            sys.stdout.write("\t".join([tree_name, sp, gene, classification, branch_length]) + '\n')

    for cecido_tip in cecidomyiidae_subtree:
        if 'grc' in cecido_tip.name:
            sp = tip2sp_name(cecido_tip)
            classification = 'cecidomyiidae'
            gene = tip2gene(cecido_tip)
            branch_length = str(cecido_tip.branch_length)
            sys.stdout.write("\t".join([tree_name, sp, gene, classification, branch_length]) + '\n')

    return(0)


# input_dir = sys.argv[1]
input_dir = 'data/testing_trees'
# tree_files = [i for i in os.listdir(input_dir) if i.endswith('treefile')]
tree_files = os.listdir('data/testing_trees/')

# with open('tables/L-busco-phylogenies-summary.tsv', 'w') as tab:
sys.stdout.write('BUSCO_id\ttype\tscfs\tgene_tree_location\tbootstraps\tbranch_lengths\n')

for file in tree_files:
    # file = 'OG0003006_tree.txt'
    gene = file.split('.')[0] # this is wrong
    input_newick = input_dir + '/' + file
    tree2assigments(input_newick)
    # sys.stdout.write(gene + '\t' + tree2assigments(input_newick) + '\n')


