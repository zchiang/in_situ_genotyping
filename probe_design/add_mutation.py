import sys

snp_file = open(sys.argv[1],'r')
vcf_file = open(sys.argv[2],'r')

snps = {}

for line in snp_file:

    column = line.rstrip().split()

    chrom = column[0]
    pos = column[1]

    snps[chrom+":"+pos] = 1


for line in vcf_file:

    if line[0] == "#":
        continue

    column = line.rstrip().split()

    chrom = column[0]
    pos = column[1]

    if chrom+":"+pos in snps:
        
        ref = column[3]
        alt = column[4]

        snps[chrom+":"+pos] = ref+":"+alt

snp_file.seek(0,0)

for line in snp_file:

    column = line.rstrip().split()

    chrom = column[0]
    pos = column[1]

    ref = snps[chrom+":"+pos].split(":")[0]
    alt = snps[chrom+":"+pos].split(":")[1]

    column.append(ref)
    column.append(alt)

    print "\t".join(column)

