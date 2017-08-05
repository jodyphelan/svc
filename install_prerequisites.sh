mkdir bin

#htslib
cd htslib
autoheader
autoconf
./configure --disable-lzma;make
mv tabix bgzip ../bin/
cd ../

#bcftools
cd bcftools
make
mv bcftools ../bin
cd ../

#samtools
cd samtools
make
mv samtools ../bin
cd ../

mv vcfutils.pl bin/


#sort 
wget http://lh3lh3.users.sourceforge.net/download/sort-20101217.tar.bz2
tar -xvf sort-20101217.tar.bz2 
cd sort-20101217
make
mv sort ../bin/sort-alt
cd ../


