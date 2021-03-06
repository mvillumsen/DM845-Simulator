# DNA Read Simulator
##### Made for Advanced Algorithms for Computational Biology (DM845) at SDU

This read simulator is made with the purpose of generating long reads which will be used with whatshap, a software for phasing genomic variants using DNA sequencing reads, also called haplotype assembly.

* The source code can be found in src/. In order to generate a VCF-file describing a mutated version of the reference chromosome we can use the following script:

```shell
python buildVCF.py chr1.fa chr1.vcf
```

where chr1.fa is the reference chromosome and chr1.vcf is the file that will be generated by the script.

In order to simulate DNA reads we can use the following script:

```shell
python sim.py chr1.fa chr1.vcf numberReads readLength errorProb chr1.sam
```

where numberReads and readLength are integers stating the number of reads to simulate and the length of the reads respectively. errorProb is a decimal num- ber that states the sequencing error probability and chr1.sam is the sequence alignment/map file that will be generated which describes the simulated reads.

In order to evaluate the results and count the number of correctly phased SNPs we can use the following script:

```shell
python evaluate.py chr1.vcf phased.vcf
```


* doc/ contains the final report of the project.

* output/ contains the output of WhatsHap (output/100k/ and output/500k/), plots generated with src/results.R (output/plots/), results of evaluate.py (output/results/) and the number of unphasable positions (output/unphasable100k.txt and output/unphasable500k.txt)

* The source code for WhatsHap can be found at: https://bitbucket.org/whatshap/whatshap
* The documentation for WhatsHap can be found at: https://whatshap.readthedocs.org/
