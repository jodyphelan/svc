#! /usr/bin/python
import sys
import subprocess
from collections import defaultdict,Counter
import json

def nucmer_align(fafile,ref,prefix):
	subprocess.call("nucmer %s %s -p %s.raw" % (ref,fafile,prefix),shell=True)
	subprocess.call("delta-filter -q %s.raw.delta > %s.delta" % (prefix,prefix),shell=True)
	subprocess.call("show-snps -CTHlr %s.delta > %s.vars" % (prefix,prefix),shell=True)


def get_var_pos(s):
	varpos = defaultdict(set)
	calls = defaultdict(lambda:defaultdict(str))
	ref = defaultdict(lambda:defaultdict(str))
	rawindels = defaultdict(lambda: defaultdict(lambda: defaultdict(str)))
	temp_alt_nuc = ""
	temp_ref_pos = ""
	for l in open("%s.vars" % s):
		arr = l.rstrip().split()	
		if temp_alt_nuc != arr[2]: temp_ref_pos = arr[0]
		temp_alt_nuc = arr[2]
		if arr[1] != "." and arr[2] != ".":
			varpos[arr[10]].add(int(arr[0]))
			calls[arr[10]][arr[0]] = arr[2]
			ref[arr[10]][arr[0]] = arr[1]
		else:
			if arr[2]==".":
				rawindels[arr[10]][temp_ref_pos]["ref"]+=arr[1]
				rawindels[arr[10]][temp_ref_pos]["alt"]+=arr[2]
			else:
				rawindels[arr[10]][arr[0]]["ref"]+=arr[1]
				rawindels[arr[10]][arr[0]]["alt"]+=arr[2]
	for chrom in rawindels:
		for pos in rawindels[chrom]:
			rseq = rawindels[chrom][pos]["ref"]
			aseq = rawindels[chrom][pos]["alt"]
			indel_len = len(rseq)
			indel_type = "+" if rseq[0]=="." else "-"
			indel_seq = aseq if rseq[0]=="." else rseq
			indel_pos = int(pos)-1
			ref_nt = ref_fa_dict[chrom][indel_pos-1:indel_pos]
			indel_str = "%s%s%s%s" % (ref_nt,indel_type,indel_len,indel_seq)
			varpos[chrom].add(int(indel_pos))
			calls[chrom][str(indel_pos)] = indel_str				
			ref[chrom][str(indel_pos)] = ref_nt

	loci = []
	for chrom in varpos:
		for i in [str(x) for x in sorted(list(varpos[chrom]))]:
			loci.append((chrom,i))
	return loci,calls,ref

def delta2cov(prefix,ref_fa_dict):
	cov = defaultdict(lambda:defaultdict(int))
	results = defaultdict(lambda:defaultdict(int))
	for l in subprocess.Popen("show-coords -THr %s.delta" % prefix,shell=True,stdout=subprocess.PIPE).stdout:
		arr = l.rstrip().split()
		for i in range(int(arr[0]),int(arr[1])):
			cov[arr[7]][i]+=1
	for chrom in ref_fa_dict:
		for pos in range(1,len(ref_fa_dict[chrom])):
			if cov[chrom][pos]!=1:
				results[chrom][str(pos)] = cov[pos]
	return results

def fa2dict(filename):
        print "Loading sequences"
        fa_dict = {}
        seq_name = ""
        for l in open(filename):
                line = l.rstrip()
                if line[0] == ">":
                        seq_name = line[1:].split()[0]
                        fa_dict[seq_name] = []
                else:
                        fa_dict[seq_name].append(line)
        result = {}
        for seq in fa_dict:
                result[seq] = "".join(fa_dict[seq])
        return result



qry_fa = sys.argv[1]
ref_fa = sys.argv[2]
prefix = sys.argv[3]

ref_fa_dict = fa2dict(ref_fa)
outfile = "%s.var.json" % (prefix)
nucmer_align(qry_fa,ref_fa,prefix)
cov_dict = delta2cov(prefix,ref_fa_dict)

var_pos,calls,ref_calls = get_var_pos(prefix)
for chrom in cov_dict:
	for pos in cov_dict[chrom]:
		calls[chrom][pos] = "N"
subprocess.call("rm %s.vars %s.delta %s.raw.delta" % (prefix,prefix,prefix),shell=True)

for s in calls:
	print "Variants:%s\nMissing:%s" % (len([x for x in calls[s].values() if x!="N" and x!="-"]),len([x for x in calls[s].values() if x=="N" or x=="-"]))

 
json.dump(calls,open(outfile,"w"))
