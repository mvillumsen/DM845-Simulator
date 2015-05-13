import sys

"""
...
"""

def main(origVCF, phasedVCF):
	original = open(origVCF, 'r').readlines()
	phased = open(phasedVCF, 'r').readlines()

	res_original = [ line for line in original if line[0] != '#' and len(line.split("\t")) > 1 ]
	res_phased = [ line for line in phased if line[0] != '#' and len(line.split("\t")) > 1 ]
	
	match, diff = 0, 0

	if len(res_original) != len(res_phased):
		print "VCF-files have different lenghts"
		print phasedVCF, len(res_phased)
		return
	
	for i in range(len(res_original)):
		orig = res_original[i].replace("\n", "").split("\t")
		phas = res_phased[i].replace("\n", "").split("\t")
		
		if orig[9] == phas[9]:
			match += 1
		else:
			diff += 1

	fileOut = 'output/results/' + phasedVCF.split("/")[-1].split(".")[0] + '.txt'
	f = open(fileOut, "w")
	f.write('matches\tdiffs\n')
	f.write('%d\t%d\n' % (match, diff))
	f.close()


if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2])
