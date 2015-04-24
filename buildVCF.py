import sys
from random import randint, random


"""
Reads a FASTA file and generates a VCF file with mutations at approx. every 400
nucleotide in the FASTA file.

Input:
    Path to a FASTA file with chromosome data (.fa file)
    Path of the output file

Output:
    a VCF file with mutation information

"""
## Takes two parameters: input filename and output filename
def main(argv):
    fileIn = argv[0]
    fileOut = argv[1]
    fastaList = open(fileIn, "r").readlines()
    chrom = fastaList[0]
    sample = '0/1'
    filter = 'PASS'
    info = '.'
    format = 'GT'
    fasta = ''.join(fastaList[1:]).replace("\n", "").replace("N", "").upper()

    nucleotides = ['A', 'T', 'C', 'G']
    header = ['##fileformat=VCFv4.1', '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']
    listOut = [["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Venter"]]

    rand = randint(0, 100)

    for i in range(len(fasta)):
        #prob = 0.0025
        #if random() < prob:
        if i % (400+rand) == 0:
            rand = randint(0, 100)
            n = [ nuc for nuc in nucleotides if nuc != fasta[i] ]
            alt = n[randint(0,2)]
            listOut.append([chrom, i, id, fasta[i], alt, '.', filter, info, format, sample])

    fo = open(fileOut, "w")

    for h in header:
        fo.write(h+"\n")

    for line in listOut:
        for ele in line:
            fo.write(ele+"\t")
        fo.write("\n")

    fo.close()

if __name__ == "__main__":
    main(sys.argv[1:])

