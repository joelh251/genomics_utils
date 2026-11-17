"""
For a fasta file with NCBI gene ID as the only component in the header, produces a fasta file with the longest transcript per gene.
"""

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input fasta file with NCBI gene ID headers")
parser.add_argument("-o", "--output", help="Path to save primary transcripts")
args = parser.parse_args()

if not args.input or not args.output:
    exit('Must supply "--input" and "--output"')

try:
    inFile = open(args.input, "r")
except:
    exit('Could not open "--input"')
try:
    outFile = open(args.output, "w")
except:
    exit('Could not open "--output"')

inSeqs = []
pair = ["", ""]
for line in inFile:
    if line.startswith(">"):
        if "" not in pair:
            inSeqs.append(pair)
        pair = ["", ""]
        pair[0] = line.strip()
    else:
        pair[1] += line.strip()
if "" not in pair:
    inSeqs.append(pair)
inFile.close()

genes = list(set([pair[0] for pair in inSeqs]))

best = {}
for pair in inSeqs:
    gene, seq = pair[0], pair[1]
    if gene not in best or len(seq) > len(best[gene][1]):
        best[gene] = pair

outSeqs = [best[gene] for gene in genes if gene in best]

for outSeq in outSeqs:
    outFile.write(outSeq[0] + "\n")
    outFile.write(outSeq[1] + "\n")
