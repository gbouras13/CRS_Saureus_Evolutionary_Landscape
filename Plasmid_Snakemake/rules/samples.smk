"""
Function for parsing the 'Plasmids' config and identifying samples and read files
"""

from itertools import chain

def samplesFromDirectory(dir):
    """Parse samples from a directory"""
    outDict = {}
    # https://stackoverflow.com/questions/11860476/how-to-unnest-a-nested-list
    samples= glob_wildcards(os.path.join(dir,'{sample}.fasta'))
    samples2 = chain(*samples)
    for sample in samples2:
        outDict[sample] = {}
        fasta = os.path.join(dir,f'{sample}.fasta')
        if os.path.isfile(fasta):
            outDict[sample]['fasta'] = fasta
        else:
            sys.stderr.write("\n"
                             "    FATAL: Error globbing files."
                             f"    {fasta} \n"
                             "    does not exist. Ensure consistent formatting and file extensions."
                             "\n")
            sys.exit(1)
    return outDict

# def samplesFromCsv(csvFile):
#     """Read samples and files from a TSV"""
#     outDict = {}
#     with open(csvFile,'r') as tsv:
#         for line in tsv:
#             l = line.strip().split(',')
#             if len(l) == 4:
#                 outDict[l[0]] = {}
#                 if os.path.isfile(l[1]) and os.path.isfile(l[2]) and os.path.isfile(l[3]):
#                     outDict[l[0]]['LR'] = l[1]
#                     outDict[l[0]]['R1'] = l[2]
#                     outDict[l[0]]['R2'] = l[3]
#                 else:
#                     sys.stderr.write("\n"
#                                      f"    FATAL: Error parsing {csvFile}. One of \n"
#                                      f"    {l[1]} or \n"
#                                      f"    {l[2]} or \n"
#                                      f"    {l[3]} \n"
#                                      "    does not exist. Check formatting, and that \n" 
#                                      "    file names and file paths are correct.\n"
#                                      "\n")
#                     sys.exit(1)
#     return outDict

def parseSamples(readFileDir):
    # for reading from directory
    if os.path.isdir(readFileDir):
      sampleDict = samplesFromDirectory(readFileDir)
    # if os.path.isfile(csvfile):
    #     sampleDict = samplesFromCsv(csvfile)
    # else:
    #     sys.stderr.write("\n"
    #                      f"    FATAL: {csvfile} is neither a file nor directory.\n"
    #                      "\n")
    #     sys.exit(1)
    if len(sampleDict.keys()) == 0:
        sys.stderr.write("\n"
                         "    FATAL: We could not detect any samples at all.\n"
                         "\n")
        sys.exit(1)
    return sampleDict


