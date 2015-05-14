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
        ref, alt = vcf[pos][0], vcf[pos][1]
        if chrm2[pos] == vcf[pos][0]:
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
def generateAllReads(chrm1, chrm2, numberReads, readLength, errorProb):

    reads = []
    nucleotides = ['A', 'T', 'C', 'G']

    for i in range(numberReads):
        r = random()
        index = randint(0, len(chrm1) - readLength - 1)

        if r < 0.5:
            reads.append( (index, ''.join(generateRead(chrm1, readLength, index, errorProb)) ) )
        else:
            reads.append( (index, ''.join(generateRead(chrm2, readLength, index, errorProb) ) ) )

    return reads


"""
Generate a single read. Used in generateAllReads().

Input:
    chrm: The chromosome to use
    readLength: The read length
    index: The index of the chromosome from which the read should start
    errorProb: The sequencing error probability

Output:
    A list containing the read generated
"""
def generateRead(chrm, readLength, index, errorProb):

    nucleotides = ['A', 'T', 'C', 'G']
    currRead = []

    for i in range(readLength):
        ind = index + i
        if random() < errorProb:
            n = [ nuc for nuc in nucleotides if nuc != chrm[ind] ]
            currRead.append(n[randint(0,2)])
        else:
            currRead.append(chrm[ind])

    return currRead


"""
Takes a list of reads and generates a SAM file and saves it to the given path.

Input:
    reads: A list of reads
    chrmLength: The length of the chromosome
"""
def generateSAM(reads, chrmLength, readLength=1000):

    ## Header for .sam-file
    sam = ['@HD\tVN:1.4\tGO:none\tSO:coordinate']     
    sam.append('@SQ\tSN:chr1\tLN:%d' % chrmLength)
    sam.append('@RG\tID:1\tPL:Platform\tLB:library\tSM:sample')

    ## Initialize static values
    FLAG = '0'
    RNAME = 'chr1'
    MAPQ = '255' # '255' indicates that the value is not available
    CIGAR = '%dM' % readLength
    RNEXT = '='
    PNEXT = '0'
    TLEN = '%d' % chrmLength # len of crhm from which the reads are generated
    QUAL = 'J'*readLength 

    ## Generate lines for .sam-file
    for read in reads:
        POS = '%d' % read[0]
        QNAME = 'r00%d' % read[0]
        SEQ = read[1]

        ## Append line to sam list
        sam.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL) )

    return sam


"""
Simulate DNA reads from the given parameters and saves the simulated reads
as a .sam-file.
## Input:
    argv[0]: .FA file
    argv[1]: .VCF file
    argv[2]: Number of reads
    argv[3]: Length of reads
    argv[4]: Sequencing error probability
    argv[5]: Output path (i.e. path/name of .sam-file)

Output:
    .sam-file: Saves the simulated reads to a Sequence Alignment/Map-file, following the
    conventions of the SAM-format
"""
def main(argv):
    numberReads = int(argv[2])
    readLength = int(argv[3])
    errorProb = float(argv[4])

    ## Save FASTA-file as string
    fastaList = open(argv[0], "r").readlines()
    chrm1 = ''.join(fastaList[1:]).replace("\n", "").replace("N", "").upper()

    ## Generate altered version based on VCF
    vcf_dict = readVCF(argv[1])
    chrm2 = generateChrm(chrm1, vcf_dict)

    ## Generate reads based on algorithm
    reads = generateAllReads(chrm1, chrm2, numberReads, readLength, errorProb)

    ## Output as .sam
    sam = generateSAM(reads, len(chrm1), readLength)

    ## Write output to .sam-file
    fileOut = open(argv[5], "w")
    for line in sam:
        fileOut.write(line)
        fileOut.write("\n")
    fileOut.close()

if __name__ == "__main__":
    main(sys.argv[1:])
