import sys
from random import random, randint

"""
Read a VCF file and outputs a dictionary mapping positions to tuples of (ref, alt).

Input:
    Path to VCF file

Output:
    Dictionary mapping position to reference and alteration, i.e. { pos:( ref, alt ) ... }
"""
def readVCF(vcf):
    readVCF = open(vcf, "r").readlines()
    vcfList = [ line.split("\t") for line in readVCF if line[0] != '#' ]
    
    return { int(ele[1]):(ele[3], ele[4]) for ele in vcfList }


"""
Generate an altered version (one with mutations) of the chromosome given as
input based on the VCF file given as input.

Input:
    fasta: .FA file as string
    vcf: vcf dictionary generated with readVCF()

Output:
    chrm2 as a string generated based on input files
"""
def generateChrm(fasta, vcf):
    chrm2 = list(fasta)

    for pos in vcf:
        ref = vcf[pos][0]
        alt = vcf[pos][1]
        if chrm2[pos] == ref:
            chrm2[pos] = alt

    return ''.join(chrm2)


"""
Generate reads based on the two chromosomes given as input. This is done
as follows:
    - Choose chrm1 or chrm2 at random
    - Choose a start position at random
    - Go over all letters in the read and change each letter with given prob.
    - Repeat for the number of reads

Input:
    chrm1: A string consisting of the reference chromosome
    chrm2: A string consisting of the mutated chromosome
    numberReads: The number of reads we want to generate
    readLength: The length of the reads
    errorProb: The error probability

Output:
    A list containing all reads generated as strings
"""
def generateReads(chrm1, chrm2, numberReads=100, readLength=1000, errorProb=0.1):

    reads = []
    nucleotides = ['A', 'T', 'C', 'G']

    for i in range(numberReads):
        currRead = []
        r = random()
        index = randint(0, len(chrm1) - readLength - 1)

        if r < 0.5:
            for a in range(readLength):
                ind = index + a
                if random() < errorProb:
                    n = [ nuc for nuc in nucleotides if nuc != chrm1[ind] ]
                    currRead.append(n[randint(0,2)])
                else:
                    currRead.append(chrm1[ind])
            reads.append( (index, ''.join(currRead)) )
        else:
            for b in range(readLength):
                ind = index + b
                if random() < errorProb:
                    n = [ nuc for nuc in nucleotides if nuc != chrm2[ind] ]
                    currRead.append(n[randint(0,2)])
                else:
                    currRead.append(chrm2[ind])
            reads.append( (index, ''.join(currRead)) )

    return reads


"""
Takes a list of reads and generates a SAM file and saves it to the given path.

Input:
    reads: A list of reads
    pathOut: The output path to which the SAM file is saved
"""
def generateSAM(reads, chrm, readLength=1000):
    ## Header for .sam-file
    sam = ['@HD\tVN:1.4\tGO:none\tSO:coordinate']     
    sam.append('@SQ\tSN:chr1\tLN:%d' % len(chrm))
    sam.append('@RG\tID:1\tPL:Platform\tLB:library\tSM:sample')

    ## Initialize static values
    RNAME = 'chr1'
    RNEXT = '='
    PNEXT = '0'
    MAPQ = '255' # TODO: Should this be calculated? '255' indicates that the value is not available
    TLEN = '%d' % len(chrm) # len of crhm from which I generate reads
    FLAG = '0' # TODO: Is this right?
    QUAL = 'J'*readLength 
    CIGAR = '%dM' % readLength

    ## Generate lines for .sam-file
    for read in reads:
        POS = '%d' % read[0]
	QNAME = 'r00%d' % read[0]
        SEQ = read[1]

        ## Append line to sam list
        sam.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL) )

    return sam


"""
## Input:
    argv[0]: .FA file
    argv[1]: .VCF file
    argv[2]: Output path (i.e. path/name of .sam-file)

Output:
    .SAM file
"""
def main(argv):
    ## Save FASTA-file as string
    fastaList = open(argv[0], "r").readlines()
    chrm1 = ''.join(fastaList[1:]).replace("\n", "").replace("N", "").upper()

    ## Generate altered version based on VCF
    vcf_dict = readVCF(argv[1])
    chrm2 = generateChrm(chrm1, vcf_dict)

    ## Generate reads based on algorithm
    # default: generateReads(chrm1, chrm2, numberReads=1000000, readLength=1000, errorProb=0.1)
    reads = generateReads(chrm1, chrm2)

    ## Output as .SAM
    sam = generateSAM(reads, chrm1)

    ## Write output to .sam-file
    fileOut = open(argv[2], "w")
    for line in sam:
        fileOut.write(line)
        fileOut.write("\n")
    fileOut.close()

if __name__ == "__main__":
    main(sys.argv[1:])
