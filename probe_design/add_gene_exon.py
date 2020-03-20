import sys

genes_file = open(sys.argv[1],'r')
pos_file = open(sys.argv[2],'r')

genes = {}

for line1 in genes_file:
	
	column = line1.rstrip().split()

	chrom = column[0]
	pos = column[1]
	gene = column[2]
	exon = column[3]
	strand = column[4]
        d0_exp = column[5]
        d28_exp = column[6]

        genes[chrom+":"+pos] = gene + ":" + exon + ":" + strand + ":" + d0_exp + ":" + d28_exp

for line2 in pos_file:

	column = line2.rstrip().split()

	chrom = column[0]
	pos = column[1]

	info = genes[chrom+":"+pos].split(":")

	print chrom + "\t" + pos + "\t" + info[0] + "\t" + info[1] + "\t" + info[2] + "\t" + info[3] + "\t" + info[4]
