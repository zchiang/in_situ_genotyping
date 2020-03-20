import sys

probe_file = open(sys.argv[1],'r')
exon_file = open(sys.argv[2],'r')

region_size = 15
comp = {"A":"T","C":"G","G":"C","T":"A","#":"#","+":"+"}

def rev_comp(seq):

    rev_comp_seq = ""

    for i in reversed(seq):
        rev_comp_seq += comp[i.upper()]
    return rev_comp_seq

for line in probe_file:

    #print line.rstrip()

    column = line.rstrip().split()
    chrom = column[0][1:-1]
    pos = int(column[1])
    gene = column[2][1:-1]
    strand = column[4][1:-1]

    exon_file.seek(0,0)
    exons = []

    for line2 in exon_file:
        
        column2 = line2.rstrip().split()
        exon_gene = column2[3]
        if gene == exon_gene:
            exons.append(line2.rstrip())

    exon_num = 0
    snp_exon = ""
    snp_exon_num = -1

    for exon in exons:

        exon_column = exon.split()
        exon_start = int(exon_column[1])
        exon_end = int(exon_column[2])

        if exon_start <= pos and pos <= exon_end:
            snp_exon = exon_column
            snp_exon_num = exon_num
            snp_pos = pos-exon_start
            #print exon_column[0:5]

        exon_num += 1

    sel = int(column[13])
    dibase_ref = ""
    dibase_alt = ""
    
    if (sel == 5 or sel == 8) and strand == "+":
        dibase_ref = column[9][1:3]
        dibase_alt = column[10][1:3]
        dibase_left = pos-1
        dibase_right = pos
        snp_pos_left = snp_pos-1
        snp_pos_right = snp_pos
    elif (sel == 3) and strand == "+": 
        dibase_ref = column[9][2:4]
        dibase_alt = column[10][2:4]
        dibase_left = pos
        dibase_right = pos+1
        snp_pos_left = snp_pos
        snp_pos_right = snp_pos+1
    if (sel == 5 or sel == 8) and strand == "-":
        dibase_ref = column[9][1:3]
        dibase_alt = column[10][1:3]
        dibase_left = pos
        dibase_right = pos+1
        snp_pos_left = snp_pos
        snp_pos_right = snp_pos+1
    elif (sel == 3) and strand == "-": 
        dibase_ref = column[9][2:4]
        dibase_alt = column[10][2:4]
        dibase_left = pos-1
        dibase_right = pos
        snp_pos_left = snp_pos-1
        snp_pos_right = snp_pos

    region_left = ""
    region_right = ""

    if dibase_left - region_size >= int(snp_exon[1]):
        region_left = snp_exon[6][snp_pos_left-15-1:snp_pos_left-1]
    else:
        curr_exon_seq = snp_exon[6][0:snp_pos_left-1]
        if snp_exon_num-1 > -1:
            prev_exon = exons[snp_exon_num-1].split()
            prev_exon_seq = prev_exon[6][len(curr_exon_seq)-region_size:]
        else:
            #print "NO PREVIOUS EXONS"
            prev_exon_seq = "".join(["#"]*(region_size-len(curr_exon_seq)))
        region_left = prev_exon_seq + "+" + curr_exon_seq


    if dibase_right + region_size < int(snp_exon[2]):
        region_right = snp_exon[6][snp_pos_right:snp_pos_right+15]
    else:
        curr_exon_seq = snp_exon[6][snp_pos_right:]
        if snp_exon_num+1 <= len(exons):
            next_exon = exons[snp_exon_num+1].split()
            next_exon_seq = next_exon[6][1:region_size-len(curr_exon_seq)]
        else:
            #print "NO FOLLOWING EXONS"
            next_exon_seq = "".join(["#"]*(region_size-len(curr_exon_seq)))
        region_right = curr_exon_seq + "+" + next_exon_seq

    if strand == "-":
        tmp = region_left
        region_left = rev_comp(region_right)
        region_right = rev_comp(tmp)
        #region_left = rev_comp(region_left)
        #region_right = rev_comp(region_right)

    column.append(region_left.upper())
    column.append(dibase_ref + "/" + dibase_alt)
    column.append(region_right.upper())
    #print region_left.upper() + "\t" + dibase_ref + "/" + dibase_alt + "\t" + region_right.upper()
    print "\t".join(column)

