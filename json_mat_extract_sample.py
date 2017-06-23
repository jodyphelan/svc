import sys
import json
from collections import defaultdict

sample = sys.argv[2]
data = json.load(open(sys.argv[1]))

sample_vars = defaultdict(dict)
for chrom,pos in data["loci"]:
	if sample in data["vars"][chrom][pos]:
		sample_vars[chrom][pos] = data["vars"][chrom][pos][sample]

print sample_vars

