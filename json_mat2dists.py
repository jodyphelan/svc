import sys
import json
from collections import defaultdict

data = json.load(open(sys.argv[1]))
print data.keys()

set_variants = defaultdict(set)
for chrom,pos in data["loci"]:
	for s in data["vars"][chrom][str(pos)]:
		if data["vars"][chrom][str(pos)][s]!="-" and data["vars"][chrom][str(pos)][s]!="N":
			set_variants[s].add((chrom,pos))

	
for s1 in data["samples"]:
	for s2 in data["samples"]:
		dist = len(set_variants[s1]-set_variants[s2]) + len(set_variants[s2]-set_variants[s1])
		print "%s\t%s\t%s" % (s1,s2,dist)
