import sys
import subprocess
from collections import defaultdict,Counter
from tqdm import tqdm
import json
import os.path
import gzip 

scriptDir = os.path.dirname(os.path.realpath(__file__))

vcfutils = scriptDir+"/bin/vcfutils.pl"
samtools = scriptdir+"/bin/samtools"
bcftools = scriptdir+"/bin/bcftools"
sort_alt = scriptdir+"/bin/sort-alt"

def check_files(programs):
	for p in programs:
		if not os.path.isfile(programs[p]):
			print "Can't find %s at %s" % (p,programs[p])
			quit()

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

def recode_indel(ref,alt):
	if len(ref)<len(alt):
		ldiff = len(alt)-len(ref)
		return "%s+%s%s" % (ref[:1],ldiff,ref[1:1+ldiff])
	else:
		ldiff = len(ref)-len(alt)
		return "%s-%s%s" % (ref[:1],ldiff,ref[1:1+ldiff])
	return indel


def call_variants(ref,bam,prefix,threads):
	fai = "%s.fai" % ref
	bai = "%s.bai" % bam
	vcf_file = "%s.vcf.gz" % prefix
	if not os.path.isfile(fai): subprocess.call("%s faidx %s" % (samtools,ref),shell=True)
	if not os.path.isfile(bai): subprocess.call("%s index %s" % (samtools,bam),shell=True)

	subprocess.call("%s splitchr -l 200000 %s | xargs -P%s -i sh -c \"%s mpileup -ugf %s -r {} %s | %s call --ploidy 1 -vmO z -o %s.part.{}.vcf.gz\"" % (vcfutils,fai,threads,samtools,ref,bam,bcftools,prefix),shell=True)
	subprocess.call("bcftools concat `ls %s.part*.vcf.gz | %s -N` | gzip -c > %s" % (vcf_concat,prefix,sort_alt,vcf_file),shell=True)
	subprocess.call("rm %s.part*.vcf.gz" % prefix,shell=True)

	base_calls = defaultdict(lambda : defaultdict(dict))
	for l in subprocess.Popen("%s depth -aa --reference %s -d 10 -Q 60 %s"  % (samtools,ref,bam),shell=True,stdout=subprocess.PIPE).stdout:
		arr = l.rstrip().split()
		if int(arr[2])<10:
			base_calls[arr[0]][int(arr[1])] = "-"
			
	
	for l in gzip.open(vcf_file,"rb"):
		if l[0]=="#": continue
		arr = l.rstrip().split()
		if arr[5]<200: continue
		if len(arr[3])>1 or len(arr[4])>1: 
			recode_indel(arr[3],arr[4])
		else:
			base_calls[arr[0]][arr[1]] = arr[4]
	return base_calls

	
		
				

if len(sys.argv)!=5:
	print "bam2varjson.py <ref> <bam> <prefix> <threads>"
	quit()

ref_file = sys.argv[1]
bam = sys.argv[2]
prefix = sys.argv[3]
threads = sys.argv[4]

outfile = "%s.var.json" % (prefix)

final_calls = call_variants(ref_file,bam,prefix,threads)
for s in final_calls:
	print "Variants:%s\nMissing:%s" % (len([x for x in final_calls[s].values() if x!="N" and x!="-"]),len([x for x in final_calls[s].values() if x=="N" or x=="-"]))
json.dump(final_calls,open(outfile,"w"))
