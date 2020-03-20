import sys

gene_exp_file = open(sys.argv[1],'r')
thresh = int(sys.argv[2])

for line in gene_exp_file:

    if "external_gene_name" in line:
        continue

    column = line.rstrip().split(",")
    gene = column[1][1:-1]
    day0_avg = (int(column[7]) + int(column[8]) + int(column[9]))/3
    day28_avg = (int(column[31]) + int(column[32]) + int(column[33]))/3

    if day0_avg > thresh and day28_avg > thresh:
        #print gene, day0_avg, day28_avg
        print gene + "\t" + str(day0_avg) + "\t" + str(day28_avg)
