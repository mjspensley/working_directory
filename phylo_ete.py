# import required packages
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from ete4 import Tree
from ete4.treeview import TreeStyle
from ete4.treeview import NodeStyle
import re
import matplotlib as plt

# read alignment file
with open("phylo.aln-fasta") as aln: 
    alignment = AlignIO.read(aln, "fasta")

# calculate distance matrix
calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(alignment)

# build phylogenetic tree
constructor = DistanceTreeConstructor(calculator)
nj_tree = constructor.nj(distance_matrix)

# save tree in newick format
Phylo.write(nj_tree, "phylo_tree.nwk", "newick")

# clean newick file ready for visualisation
with open("phylo_tree.nwk", "r") as f:
    newick_content = f.read()
cleaned_newick = re.sub(r'\)Inner\d+:', '):', newick_content)   # removes internal node names

# save cleaned newick file
with open("phylo_tree_cleaned.nwk", "w") as f:
    f.write(cleaned_newick)

# load newick file with ete toolkit
ete_tree = Tree("phylo_tree_cleaned.nwk")

# remove blue dots, thicker lines
ns = NodeStyle()
ns["size"] = 0              # remove node circles (set size to 0)
ns["hz_line_width"] = 1     # horizontal line thickness
ns["vt_line_width"] = 1     # vertical line thickness

# apply style to all nodes
for node in ete_tree.traverse():
    node.set_style(ns)

# visualise tree circular mode
ts = TreeStyle()
ts.mode = "c"               # circular mode
ts.show_leaf_name = True    # show species names

# Render the tree with higher resolution
ete_tree.render("phylo_tree.png", 
                w=1000,
                h=1000,
                dpi=300,
                tree_style=ts)
print("phylo_tree.png")