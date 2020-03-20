import sys
import re

repeats = ["AAAAAA","CCCCCC","GGGGGG","TTTTTT"]

for line in sys.stdin:

    column = line.rstrip().split()

    left_region = column[14]
    dibase = column[15]
    right_region = column[16]

    left_region = re.sub('\+', '', left_region)
    right_region = re.sub('\+', '', right_region)
    
    has_repeat = 0
    left_gc = -1
    right_gc = -1

    if "#" in left_region or "#" in right_region:
        column.extend([str(has_repeat),str(left_gc),str(right_gc)])
        print "\t".join(column)
        continue

    for repeat in repeats:
        if repeat in left_region or repeat in right_region:
            has_repeat = 1
            continue

    left_gc = left_region.count("C") + dibase[0].count("C") + left_region.count("G") + dibase[0].count("G")
    right_gc = right_region.count("C") + dibase[1].count("C") + right_region.count("G") + dibase[1].count("G")

    column.extend([str(has_repeat),str(left_gc),str(right_gc)])
    print "\t".join(column)
