import ahocorasick
import numpy as np
from rev_script import *
from coverage_map import *

# Open the psm file and generate a list of psm peptides
with open("psm.tsv") as fin:
    # Ignores the first line because it is the header
    psm_list = fin.read().split("\n")[1:]
    # The column is seperated by tabs and the third one is the peptide we need
    # The last line is an empty line, thus len(psm_list)-1 is used
    peptide_list = [psm_list[peptide].split("\t")[2] for peptide in range(0,len(psm_list)-1)]

# Open the Uniport fasta file for proteins
with open("UniProt_Human.fasta") as fin:
    # Split the protein entries with ">" in the beginning of the description
    split_list = [protein for protein in fin.read()[1:].split('\n>')]
    # Uniport list
    uniport_list = [split_list[uniport].split("|")[1] for uniport in range(0, len(split_list))]
    # Protein sequence with join for multiple lines
    protein_list = ["".join(split_list[uniport].split("\n")[1:]) for uniport in range(0, len(split_list))]

# Long string for aho-corasick
aho_seq = "|".join(protein_list)
# a list of zeroes for aho-corasick protein identify
zero_list = np.zeros(len(aho_seq), dtype= int).tolist()
# initialize the protein index
ind = int(0)
# iterate through the aho-string, if it meets "|", the index should add 1
for i in range(0,len(zero_list)):
    if aho_seq[i] != "|":
        zero_list[i] += ind
    else:
        ind += 1

# Create Automaton
automaton = ahocorasick.Automaton()

for idx, key in enumerate(peptide_list):
    automaton.add_word(key, (idx, key))

# Convert the trie to an Aho-Corasick automaton to enable Aho-Corasick search
automaton.make_automaton()

found_dict = {}
# Do the search of the aho-corasick
for end_index, (insert_order, original_value) in automaton.iter(aho_seq):
    start_index = end_index - len(original_value) + 1
    assert aho_seq[start_index:start_index + len(original_value)] == original_value
    uniport_id = uniport_list[zero_list[end_index]]
    if uniport_id in found_dict:
        found_dict[uniport_id].append(original_value)
    else:
        found_dict[uniport_id] = []
        found_dict[uniport_id].append(original_value)

if __name__ == "__main__":
    print(found_dict)
    print("")
    print(split_list[0])
    print(rev_string(split_list[0], 1))
    print(rev_string(split_list[0], 2))
    print(rev_string(split_list[0], 3))
    protein_coverage_map(uniport_list[0], found_dict, protein_list, uniport_list)