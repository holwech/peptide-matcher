# %%
from typing import Dict, List, Optional, Tuple
from tqdm import tqdm

# %%

with open('../data/SARS_Cov_2_peptide_seqs_9mers.txt') as file:
    peptides_file: List[str] = file.readlines()

with open('../data/UP000005640_9606.fasta') as file:
    proteins_file: List[str] = file.readlines()

# %%
def fasta_to_dict(file: List[str]) -> Dict[str, str]:
    mapped_proteins: Dict[str, str] = {}

    key: str = "";
    for row in file:
        if (row[0] == ">"):
            key = row.split("|")[1]
            mapped_proteins[key] = ""
        else:
            mapped_proteins[key] += row.strip("\n")
    return mapped_proteins
            
# %%
# Convert to bytes might speed up?
proteins = fasta_to_dict(proteins_file)
print(len(proteins.items()))

# %%
peptides = [line.split("\t")[-1].strip("\n") for line in peptides_file]
print(len(peptides))

# %%
"TPAAPAPAP" in proteins["P04637"]

# %%
# O(m*n*k) -> O(n^3)
def slow_search_all(peptides: List[str], proteins: Dict[str, str]) -> List[Tuple[str, str]]:
    matches: List[Tuple[str, str]] = []
    for peptide in tqdm(peptides):
        for key, protein in proteins.items():
            # Uses a mix of Boyer-More and Horspool
            # O(n)?
            if peptide in protein:
                matches.append((key, peptide))
    return matches

# %%
%timeit -n 1 -r 1 slow_search_all(peptides, proteins)
slow_search_results = slow_search_all(peptides, proteins)

# %%
from ahocorapy.keywordtree import KeywordTree

# O(mk)
def aho_corasick_search(peptides: List[str], proteins: Dict[str, str]) -> List[Tuple[str, str]]:
    matches: List[Tuple[str, str]] = []

    kwtree = KeywordTree(case_insensitive=True)
    for peptide in peptides:
        kwtree.add(peptide)
    kwtree.finalize()
    for key, protein in tqdm(proteins.items()):
        match = kwtree.search(protein)
        if match != None:
            matches.append((key, match[0])) 
    return matches

# %%
%timeit -n 1 -r 1 aho_corasick_search(peptides, proteins)
aho_corasick_search_results = aho_corasick_search(peptides, proteins)

# %%
# O(m+n) -> O(n)?
# Building hash set is O(m), looping through all proteins is O(n)
# If string is UTF-8, string slicing is O(k), making this O(m+n*k) -> O(n^2)
# If strings are short, the cost of slicing is neglectable
def joachim_search(peptides: List[str], proteins: Dict[str, str], peptide_len: int = 9) -> List[Tuple[str, str]]:
    matches: List[Tuple[str, str]] = []

    peptides_set = set(peptides)
    for key, protein in tqdm(proteins.items()):
        for i in range(0, len(protein) - peptide_len + 1):
            sequence = protein[i:i + peptide_len]
            if sequence in peptides_set:
                return True
                matches.append((key, sequence))
    return matches

# %%
%timeit -n 1 -r 1 joachim_search(peptides, proteins)
joachim_search_results = joachim_search(peptides, proteins)

# %%
def sort_list_of_lists(results: List[List]) -> None:
    for result in results:
        result.sort()

sort_list_of_lists([
    slow_search_results,
    aho_corasick_search_results,
    joachim_search_results
])

# %%
len(slow_search_results)

# %%
len(aho_corasick_search_results)

# %%
len(joachim_search_results)
# %%
slow_search_results
# %%
aho_corasick_search_results
# %%
joachim_search_results

# %%
