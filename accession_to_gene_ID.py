"""
Translates protein fasta file headers to NCBI gene IDs
"""

import argparse
import re
import csv

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Fasta file to translate headers of")
parser.add_argument("-c", "--cds", help="CDS file to build dictionary from")
parser.add_argument("-d", "--dict", help="Pre-generated dictionary file")
parser.add_argument("-o", "--output", help="Path to save translated fasta file")
parser.add_argument("-s", "--save", help="Path to save generated dictionary")
args = parser.parse_args()

if args.dict and args.cds:
    exit('Mutually exclusive arguments "--cds" and "--dict" supplied, only 1 may be supplied')

refSeqAccessionPattern = r'(NW_\d+(?:\.\d+)?)'
genePattern = r'\[gene=((?:[^\[\]]|\[[^\[\]]*\])*)\]'
geneIDPattern = r'\[db_xref=GeneID:((?:[^\[\]]|\[[^\[\]]*\])*)\]'
protPattern = r'\[protein=((?:[^\[\]]|\[[^\[\]]*\])*)\]'
protIDPattern = r'\[protein_id=((?:[^\[\]]|\[[^\[\]]*\])*)\]'

if args.input and not args.output:
    if args.input:
        exit('No "--output" specified')
    if args.output:
        exit('No "--input" specified')

if args.cds:
    try:
        file = open(args.cds, "r")
    except:
        exit("CDS file could not be opened")
    
    cdsHeaders = []
    cdsDict = {}
    for line in file:
        if line.startswith(">"):
            cdsHead = [None, None, None, None, None]
            refSeqMatch = re.search(refSeqAccessionPattern, line)
            if refSeqMatch:
                cdsHead[0] = refSeqMatch.group(1)

            geneMatch = re.search(genePattern, line)
            if geneMatch:
                cdsHead[1] = geneMatch.group(1)

            geneIDMatch = re.search(geneIDPattern, line)
            if geneIDMatch:
                cdsHead[2] = geneIDMatch.group(1)

            protMatch = re.search(protPattern, line)
            if protMatch:
                cdsHead[3] = protMatch.group(1)

            protIDMatch = re.search(protIDPattern, line)
            if protIDMatch:
                cdsHead[4] = protIDMatch.group(1)

            cdsHeaders.append(cdsHead)
            cdsDict[cdsHead[4]] = cdsHead[0:3]
            # Unavoidable issue here. There are duplicated protein accessions.
            # Dicts don't allow duplicated keys, so these will get overwritten
            # Not much I can really do about it

    file.close()

if args.save:
    with open(args.save, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["refseq_accession", "gene", "NCBI_gene_id", "protein", "protein_accession"])
        writer.writerows(cdsHeaders)

if args.dict:
    try:
        file = open(args.dict, "r", newline="")
    except:
        exit("Dictionary file could not be opened")

    reader = csv.reader(file)
    next(reader)
    rows = [row for row in reader]
    cdsDict = {}
    for row in rows:
        cdsDict[row[4]] = row[0:3]

    file.close()

protAccessionPattern = r'(XP_\d+(?:\.\d+)?)'

if args.input and args.output:
    try:
        inFile = open(args.input, "r")
    except:
        exit('Could not open "--input"')
    try:
        outFile = open(args.output, "w")
    except:
        exit('Could not open "--output"')

    pair = ["", ""]
    for line in inFile:
        if line.startswith(">"):
            if "" not in pair:
                pair[1] += "\n"
                outFile.write(pair[0])
                outFile.write(pair[1])
            pair = ["", ""]
            protIDMatch = re.search(protAccessionPattern, line)
            if protIDMatch:
                geneID = cdsDict[protIDMatch.group(1)][2]
                pair[0] = ">" + geneID + "\n"
        else:
            pair[1] += line.strip()
    if "" not in pair:
                    pair[1] += "\n"
                    outFile.write(pair[0])
                    outFile.write(pair[1])
    inFile.close()
    outFile.close()