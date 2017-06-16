import json
import sys
from collections import defaultdict

json_file = sys.argv[1]
outfile = sys.argv[2]
jvar = json.load(open(json_file))
seqs = defaultdict(list)
for chrom,pos in jvar["loci"]:
	temp = jvar["vars"][chrom][str(pos)]
	temp_vars = "".join([x for x in temp.values() if x!="N" and x!="-"])
	if "+" in temp_vars or "-" in temp_vars: continue
	for s in jvar["samples"]:
		if s in temp:
			if temp[s]=="-":
				seqs[s].append("N")
			else:
				seqs[s].append(temp[s])
		else:
			seqs[s].append(jvar["ref"][chrom][str(pos)])

with open(outfile,"w") as o:
	for s in jvar["samples"]:
		o.write(">%s\n%s\n" % (s,"".join(seqs[s])))
