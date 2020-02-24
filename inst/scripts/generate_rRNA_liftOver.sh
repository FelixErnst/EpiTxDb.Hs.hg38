#!/bin/bash

cd ../extdata/

# download tools
mkdir tools
cd ./tools
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/twoBitInfo ./twoBitInfo
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faToTwoBit  ./faToTwoBit
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/blat/blat ./blat
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/axtChain ./axtChain
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/chainMergeSort ./chainMergeSort
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/chainSplit ./chainSplit
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/chainSort ./chainSort
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/chainNet ./chainNet
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/netChainSubset ./netChainSubset
cd ..

# generate 2bit files and chrominfo files
./tools/faToTwoBit rRNA_hg19.fasta rRNA_hg19.2bit
./tools/faToTwoBit rRNA_hg38.fasta rRNA_hg38.2bit
./tools/twoBitInfo rRNA_hg19.2bit rRNA_hg19.chromInfo
./tools/twoBitInfo rRNA_hg38.2bit rRNA_hg38.chromInfo

# generate rRNA.hg19Tohg38.liftOver
mkdir hg19Tohg38
cd ./hg19Tohg38
../tools/blat ../rRNA_hg19.2bit ../rRNA_hg38.2bit -minIdentity=98 rRNA.hg19Tohg38.psl -noHead -minScore=100
../tools/axtChain -linearGap=medium -psl rRNA.hg19Tohg38.psl ../rRNA_hg19.2bit ../rRNA_hg38.2bit rRNA.hg19Tohg38.chain

mkdir chainMerge
../tools/chainMergeSort rRNA.hg19Tohg38.chain | ../tools/chainSplit chainMerge stdin -lump=50

cat chainMerge/*.chain > all.chain
../tools/chainSort all.chain all.sorted.chain

mkdir net
../tools/chainNet all.sorted.chain ../rRNA_hg19.chromInfo ../rRNA_hg38.chromInfo net/all.net /dev/null

../tools/netChainSubset net/all.net all.chain rRNA.hg19Tohg38.liftOver

cp rRNA.hg19Tohg38.liftOver ../rRNA.hg19Tohg38.liftOver
cd ..

# generate rRNA.hg19Tohg38.liftOver
mkdir hg38Tohg19
cd ./hg38Tohg19
../tools/blat ../rRNA_hg38.2bit ../rRNA_hg19.2bit -minIdentity=98 rRNA.hg38Tohg19.psl -noHead -minScore=100
../tools/axtChain -linearGap=medium -psl rRNA.hg38Tohg19.psl ../rRNA_hg38.2bit ../rRNA_hg19.2bit rRNA.hg38Tohg19.chain

mkdir chainMerge
../tools/chainMergeSort rRNA.hg38Tohg19.chain | ../tools/chainSplit chainMerge stdin -lump=50

cat chainMerge/*.chain > all.chain
../tools/chainSort all.chain all.sorted.chain

mkdir net
../tools/chainNet all.sorted.chain ../rRNA_hg38.chromInfo ../rRNA_hg19.chromInfo net/all.net /dev/null

../tools/netChainSubset net/all.net all.chain rRNA.hg38Tohg19.liftOver

cp rRNA.hg38Tohg19.liftOver ../rRNA.hg38Tohg19.liftOver
cd ..

# cleanup
rm -R ./hg19Tohg38
rm -R ./hg38Tohg19
rm *.2bit
rm *.chromInfo
rm -R ./tools
