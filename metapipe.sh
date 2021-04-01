#!/bin/bash

date=$(date +%Y-%m-%d)
contigs=$1
threads=10
sample=$2

#software
blastn=/home/bioinf/software/ncbi-blast-2.8.1+/bin/blastn
ntdb=/data/DBs/NT/nt #if not then: --->update_blastdb.pl nt
#from: https://github.com/sdwfrost/Simple-LCA changed the line 20 as >>>> a = list(map(str.strip, taxid.split("|")))
# wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip and unzipping it on the same directory
addtaxonomydir=/home/bioinf/software/Simple-LCA
blasthit=/home/bioinf/.egjb/Results_Minimus2_Prodigal/Megahit_IDBA/blastn/blasthitcounter.pl

#seting the environment
mkdir ${sample}_metapipe_$date
cd ${sample}_metapipe_$date

#make a blastdb from the Metagneome sequece
mkdir ${sample}_DB
echo "making database......"
makeblastdb -in ../$contigs -dbtype nucl -parse_seqids -out ./${sample}_DB/${sample}_database

#1blast
echo "Blasting........... "
blastn -query ../$contigs -out ${sample}_blastvsNT.txt -db $ntdb -outfmt '6 qseqid stitle sacc staxids pident qcovs evalue bitscore' -num_threads $threads -task blastn &
PID=$!
i=3
sp="/-\|"
echo -n ' '
while [ -d /proc/$PID ]
do
  printf "\b${sp:i++%${#sp}:1}"
done
##optional with max hits and max hsp =1
#blastn -query ${sample} -out ${sample}_blastvsNT_maxT_100_maxH1.txt -db $ntdb -outfmt '6 qseqid stitle sacc staxids pident qcovs evalue bitscore' -num_threads $threads -task blastn -max_target_seqs 100 -max_hsps 1

#2. add taxonomy
grep -v "[0-9];[0-9]" ${sample}_blastvsNT.txt  > ${sample}_blastvsNT_sing.txt
python $addtaxonomydir/add_taxonomy.py -i ${sample}_blastvsNT_sing.txt -t $addtaxonomydir/rankedlineage.dmp -m $addtaxonomydir/merged.dmp -o ${sample}_blastvsNT_sing_addtax.out

#3.blasthitcounter
evalue=20
perl $blasthit ${sample}_blastvsNT_sing_addtax.out $evalue > ${sample}_blastvsNT_sing_addtax_allTax.txt
cat ${sample}_*_max.txt | awk '{print $1"\t"$2"\t"$3"\t"$4}' > ${sample}_majorityruleassignation.txt
echo 'Extracting fasta information ..........'
for i in $(ls *_max.txt);do awk '{print $1}' $i > list2extract; blastdbcmd -entry_batch list2extract -dbtype nucl -db ${sample}_DB/${sample}_database >  $i.fasta;done
rm list2extract
echo "Done"
