import sys

#gene = "ACTB"

gene_file = open(sys.argv[1],'r')
exon_file = open(sys.argv[2],'r')

genes = {} 

for line in gene_file:

    #gene = line.rstrip()
    #genes[gene] = 1
    column = line.rstrip().split()
    genes[column[0]] = column[1]+":"+column[2]
    
for line in exon_file:

    column = line.rstrip().split()
    curr_gene = column[3]

    if curr_gene in genes:
            chrom = column[0]
            start = int(column[1])
            end = int(column[2])
            exon = column[4]
            strand = column[5]

            d0_exp = genes[curr_gene].split(":")[0]
            d28_exp = genes[curr_gene].split(":")[1]

            for i in range(start,end+1):
                    print chrom + "\t" + str(i) + "\t" + curr_gene + "\t" + exon + "\t" + strand + "\t" + d0_exp + "\t" + d28_exp


