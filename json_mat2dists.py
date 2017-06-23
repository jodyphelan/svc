import sys
import json
from collections import defaultdict

data = json.load(open(sys.argv[1]))

set_variants = defaultdict(set)
for chrom,pos in data["loci"]:
	for s in data["vars"][chrom][str(pos)]:
		if data["vars"][chrom][str(pos)][s]!="-" and data["vars"][chrom][str(pos)][s]!="N":
			set_variants[s].add((chrom,pos))

print "\t%s" % ("\t".join(data["samples"]))	
for s1 in data["samples"]:
	dists = []
	for s2 in data["samples"]:
		dist = len(set_variants[s1]-set_variants[s2]) + len(set_variants[s2]-set_variants[s1])
		dists.append(dist)
	print "%s\t%s" % (s1,"\t".join([str(x) for x in dists]))
