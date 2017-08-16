#! /usr/bin/python
import sys
import os.path
import subprocess

def index_ref(ref):
	if not os.path.isfile(ref+".pac"):
		subprocess.call("bwa index %s" % ref,shell=True)

def bwa_mapping(ref,prefix,threads):
	tp1 = prefix+"_1_trimmed_paired.fq"
	tp2 = prefix+"_2_trimmed_paired.fq"
	tu1 = prefix+"_1_trimmed_unpaired.fq"
	tu2 = prefix+"_2_trimmed_unpaired.fq"
	pbam = prefix+".paired.bam"
	ubam1 = prefix+".single_1.bam"
	ubam2 = prefix+".single_2.bam"

	print "Running BWA for $sample"
	index_ref(ref)
	bwa_prefix = "bwa mem -t %s -c 100 -R '@RG\\tID:$sample\\tSM:$sample\\tPL:Illumina' -M -T 50" % threads
	print "Mapping paired reads"
	subprocess.call("%s %s %s %s 2>> err | sambamba view -t %s -S -f bam /dev/stdin 2>> err | sambamba sort -o %s -t %s /dev/stdin > log 2>> err" % (bwa_prefix,ref,tp1,tp2,threads,pbam,threads),shell=True)
	print "Mapping unpaired reads"
	subprocess.call("%s %s %s 2>> err | sambamba view -t %s -S -f bam /dev/stdin 2>> err | sambamba sort -o %s -t %s /dev/stdin > log 2>> err" % (bwa_prefix,ref,tu1,threads,ubam1,threads),shell=True)
	subprocess.call("%s %s %s 2>> err | sambamba view -t %s -S -f bam /dev/stdin 2>> err | sambamba sort -o %s -t %s /dev/stdin > log 2>> err" % (bwa_prefix,ref,tu2,threads,ubam2,threads),shell=True)
	print "Merging bams"
	subprocess.call("sambamba merge -t %s %s.unsorted.bam %s %s %s 2>> log" % (threads,prefix,pbam,ubam1,ubam2),shell=True)
	print "Sorting bams"
	subprocess.call("sambamba sort -t %s -o %s.bam %s.unsorted.bam 2>> log" % (threads,prefix,prefix),shell=True)
		


def trim(r1,r2,prefix,threads):
	print "Trimming reads"
	subprocess.call("java -jar /opt/storage/pathogenseq/trimmomatic.jar PE -threads %s -phred33 %s %s %s_1_trimmed_paired.fq %s_1_trimmed_unpaired.fq %s_2_trimmed_paired.fq %s_2_trimmed_unpaired.fq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 2>> log" % (threads,r1,r2,prefix,prefix,prefix,prefix),shell=True)



if len(sys.argv)<6:
	print "fastq2bam.py <ref> <r1> <r2> <prefix> <threads>"
	quit()



script,ref,r1,r2,prefix,threads = sys.argv[:6]

print "ref=%s\nread1=%s\nread2=%s\nprefix=%s\nthreads=%s" % (ref,r1,r2,prefix,threads)
trim(r1,r2,prefix,threads)
bwa_mapping(ref,prefix,threads)
if "--noclean" not in sys.argv:
	subprocess.call("rm %s_*.fq  %s.unsorted* %s.paired.*  %s.single* log err" % (prefix,prefix,prefix,prefix),shell=True) 
