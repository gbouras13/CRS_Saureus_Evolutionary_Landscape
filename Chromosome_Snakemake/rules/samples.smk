"""
Function for parsing the metadata csv file config and identifying samples and read files
"""

from itertools import chain



def gessFromCsv(csvFile):
    """Read gess samples and files from a CSV"""
    outDict = {}
    with open(csvFile,'r') as csv:
        for line in csv:
            l = line.strip().split(',')
            if len(l) == 12:
                outDict[l[0]] = {}
                if  os.path.isfile(l[3])  and os.path.isfile(l[4]) and os.path.isfile(l[5]) and os.path.isfile(l[6]) and os.path.isfile(l[7]) and os.path.isfile(l[8]) and os.path.isfile(l[9]) and os.path.isfile(l[10]) and os.path.isfile(l[11])   :
                    outDict[l[0]]['T0'] = l[1]
                    outDict[l[0]]['T1'] = l[2]
                    outDict[l[0]]['T0_fasta'] = l[3]
                    outDict[l[0]]['T1_fasta'] = l[4]
                    outDict[l[0]]['T0_gbk'] = l[5]
                    outDict[l[0]]['T1_R1'] = l[6]
                    outDict[l[0]]['T1_R2'] = l[7]
                    outDict[l[0]]['T0_gff'] = l[8]
                    outDict[l[0]]['T1_gff'] = l[9]
                    outDict[l[0]]['T1_gbk'] = l[10]
                    outDict[l[0]]['T1_long'] = l[11]
                else:
                    sys.stderr.write("\n"
                                     f"    FATAL: Error parsing {csvFile}. The  file \n"
                                     f"    {l[3]} \n"
                                     "    does not exist. Check formatting, and that \n" 
                                     "    file names and file paths are correct.\n"
                                     "\n")
                    sys.exit(1)
    return outDict



def parseGess(csvfile):
    # for reading from directory
    #if os.path.isdir(readFileDir):
    #   sampleDict = samplesFromDirectory(readFileDir)
    if os.path.isfile(csvfile):
        sampleDict = gessFromCsv(csvfile)
    else:
        sys.stderr.write("\n"
                         f"    FATAL: {csvfile} is neither a file nor directory.\n"
                         "\n")
        sys.exit(1)
    if len(sampleDict.keys()) == 0:
        sys.stderr.write("\n"
                         "    FATAL: We could not detect any samples at all.\n"
                         "\n")
        sys.exit(1)
    return sampleDict
