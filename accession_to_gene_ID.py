"""
Translates protein fasta file headers to NCBI gene IDs
"""

import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Fasta file to translate headers of")
parser.add_argument("-c", "--cds", help="CDS file to build dictionary from")
parser.add_argument("-o", "--output", help="Path to save translated fasta file")
args = parser.parse_args()


refSeqAccessionPattern = r'(NW_\d+(?:\.\d+)?)'
genePattern = r'\[gene=((?:[^\[\]]|\[[^\[\]]*\])*)\]'
geneIDPattern = r'\[db_xref=GeneID:((?:[^\[\]]|\[[^\[\]]*\])*)\]'
protPattern = r'\[protein=((?:[^\[\]]|\[[^\[\]]*\])*)\]'
protIDPattern = r'\[protein_id=(XP_(?:[^\[\]]|\[[^\[\]]*\])*)\]'


if args.cds:
    try:
        file = open(args.cds, "r")
    except:
        exit("CDS file could not be opened")
    
    cdsHeaders = []
    cdsDict = {}
    duplicatesDict = {}
    for line in file:
        if line.startswith(">"):
            cdsHead = [None, None, None, None, None]

            protIDMatch = re.search(protIDPattern, line)
            if not protIDMatch:
                next
            cdsHead[4] = protIDMatch.group(1)
            
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

            cdsHeaders.append(cdsHead)

            if cdsHead[4] in cdsDict:
                if cdsHead[4] in duplicatesDict:
                    duplicatesDict[cdsHead[4]].append(cdsHead[0:3])
                else:
                    duplicatesDict[cdsHead[4]] = [cdsHead[0:3]]
            else:
                cdsDict[cdsHead[4]] = cdsHead[0:3]

    file.close()

for key in duplicatesDict:
    print(str(key) + " : " + str(duplicatesDict[key]))
    print(key)

print(len(duplicatesDict))

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

    duplicateMode = False
    pair = ["", ""]
    for line in inFile:
        if line.startswith(">"):
            if "" not in pair and not duplicateMode:
                pair[1] += "\n"
                outFile.write(pair[0])
                outFile.write(pair[1])

            if "" not in pair and duplicateMode:
                for duplicate in duplicatesDict[pair[0]]:
                    head = ">" + duplicate[2] + "\n"
                    outFile.write(head)
                    outFile.write(pair[1] + "\n")
                    duplicateMode = False

            pair = ["", ""]
            protIDMatch = re.search(protAccessionPattern, line)
            if protIDMatch:
                protRef = protIDMatch.group(1)
                if protRef in duplicatesDict:
                    pair[0] = protRef
                    duplicateMode = True
                else:
                    geneID = cdsDict[protRef][2]
                    pair[0] = ">" + geneID + "\n"
        else:
            pair[1] += line.strip()

    if "" not in pair and not duplicateMode:
        pair[1] += "\n"
        outFile.write(pair[0])
        outFile.write(pair[1])

    if "" not in pair and duplicateMode:
        for duplicate in duplicatesDict[pair[0]]:
            head = ">" + duplicate[2] + "\n"
            outFile.write(head)
            outFile.write(pair[1] + "\n")
            duplicateMode = False

    inFile.close()
    outFile.close()

