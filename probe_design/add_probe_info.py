import sys

comp = {"A":"T","C":"G","G":"C","T":"A"}
eff = {"dA/pdA":1,"dA/pdC":0,"dA/pdG":0,"dA/pdT":1,"dC/pdA":1,"dC/pdC":0,"dC/pdG":0,"dC/pdT":1,"dG/pdA":1,"dG/pdC":0,"dG/pdG":0,"dG/pdT":1,"dT/pdA":1,"dT/pdC":1,"dT/pdG":0,"dT/pdT":1}

for line in sys.stdin:

        column = line.rstrip().split()

        # remove indels and *s

        if len(column[7]) > 1 or len(column[8]) > 1 or column[7] == "*" or column[8] == "*":
            column.extend(["NA","NA","NA","NA"])
            print "\t".join(column)
            continue

        strand = column[4]

        # change ref seq to uppercase

        column[9] = column[9].upper()

        if strand == "+":

            ref_seq = column[9]
            alt_seq = column[9][0] + column[8] + column[9][2]

        elif strand == "-":

            ref_seq = comp[column[9][2]] + comp[column[9][1]] + comp[column[9][0]]
            alt_seq = comp[column[9][2]] + comp[column[8]] + comp[column[9][0]]

        probe5_ref = "d" + comp[ref_seq[1]] + "/pd" + comp[ref_seq[0]]
        probe5_alt = "d" + comp[alt_seq[1]] + "/pd" + comp[alt_seq[0]]        
        probe3_ref = "d" + comp[ref_seq[2]] + "/pd" + comp[ref_seq[1]]
        probe3_alt = "d" + comp[alt_seq[2]] + "/pd" + comp[alt_seq[1]]

        probe5 = ref_seq[0:2] + "->" + probe5_ref + "(" + str(eff[probe5_ref]) + "):" + alt_seq[0:2] + "->" + probe5_alt + "(" + str(eff[probe5_alt]) + ")"
        probe3 = ref_seq[1:3] + "->" + probe3_ref + "(" + str(eff[probe3_ref]) + "):" + alt_seq[1:3] + "->" + probe3_alt + "(" + str(eff[probe3_alt]) + ")"

        if (eff[probe5_ref] and eff[probe5_alt]) and (eff[probe3_ref] and eff[probe3_alt]):
            sel = 8
        elif eff[probe5_ref] and eff[probe5_alt]:
            sel = 5
        elif eff[probe3_ref] and eff[probe3_alt]:
            sel = 3
        else:
            sel = 0

        column[9] = ref_seq
        column.extend([alt_seq,probe5,probe3,str(sel)])
        print "\t".join(column)


