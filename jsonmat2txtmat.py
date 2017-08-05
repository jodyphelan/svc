import sys
import json

infile = sys.argv[1]
outfile = sys.argv[2]

mat = json.load(open(infile))

F = open(outfile,"w")
F.write("#chr\tpos\tref\tinfo\ttype\t%s\n" % ("\t".join(mat["samples"])))
for chrom,pos in mat["loci"]:
	F.write("%s\t%s\t%s\t.\t.\t%s\n" % (chrom,pos,mat["ref"][chrom][pos],"\t".join([mat["vars"][chrom][pos][s] if s in mat["vars"][chrom][pos] else mat["ref"][chrom][pos] for s in mat["samples"]])))
