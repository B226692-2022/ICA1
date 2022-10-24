#!/bin/bash

destination="./"
cp -v /localdisk/data/BPSM/ICA1/fastq/*.fq.gz $destination
cp -v /localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed $destination
cp -v /localdisk/data/BPSM/ICA1/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz $destination
cp -v /localdisk/data/BPSM/ICA1/fastq/*.fqfiles $destination
#Print new file with old name and new
awk 'BEGIN{FS="\t"}{print $6, $9 = $1"-"$2"-"$4"h-"$5"_1.fq"}' Tco.fqfiles > Tco.fqfiles1.txt
awk 'BEGIN{FS="\t"}{print $7, $9 = $1"-"$2"-"$4"h-"$5"_2.fq"}' Tco.fqfiles > Tco.fqfiles2.txt

#Rename files based on new tables (original name in column 1 - New in column 2)
awk '{f=$1;sub("-.*\\.","."$2".",$1);system("mv "f" "$1)}' Tco.fqfiles1.txt
awk '{f=$1;sub("-.*\\.","."$2".",$1);system("mv "f" "$1)}' Tco.fqfiles2.txt

#Fastqc check -t specifies use of more cores, --extract unzips fastqc summary

fastqc -t 16 --extract *.fq.gz

#Summarise files 

cat ./Tco*/summary.txt > summarymega.txt

total=$(grep "PASS\|FAIL\|WARN" summarymega.txt | wc -l)
passes=$(grep "PASS" summarymega.txt | wc -l)
warn=$(grep "WARN" summarymega.txt | wc -l)
fail=$(grep "FAIL" summarymega.txt | wc -l)


echo "Overall Summary"
echo "${passes} Passes"
echo "${warn} Warnings"
echo "${fail} Fails"

# Per base sequence quality

awk 'BEGIN{FS="\t";}{if($2 == "Per base sequence quality"){print $0, $1, $2;}}' summarymega.txt > 1.txt
total=$(grep "PASS\|FAIL\|WARN" 1.txt | wc -l)
passes=$(grep "PASS" 1.txt | wc -l)
warn=$(grep "WARN" 1.txt | wc -l)
fail=$(grep "FAIL" 1.txt | wc -l)

echo "Overall Summary for: Per base sequence quality
        ${passes} Passes
        ${warn} Warnings
        ${fail} Fails" > pbsq.txt

# Per base sequence content
awk 'BEGIN{FS="\t";}{if($2 == "Per base sequence content"){print $0, $1, $2;}}' summarymega.txt > 2.txt
total=$(grep "PASS\|FAIL\|WARN" 2.txt | wc -l)
passes=$(grep "PASS" 2.txt | wc -l)
warn=$(grep "WARN" 2.txt | wc -l)
fail=$(grep "FAIL" 2.txt | wc -l)

echo "Overall Summary for: Per base sequence content
        ${passes} Passes
        ${warn} Warnings
        ${fail} Fails" > pbsc.txt
#duplicates

awk 'BEGIN{FS="\t";}{if($2 == "Sequence Duplication Levels"){print $0, $1, $2;}}' summarymega.txt > 3.txt
total=$(grep "PASS\|FAIL\|WARN" 3.txt | wc -l)
passes=$(grep "PASS" 3.txt | wc -l)
warn=$(grep "WARN" 3.txt | wc -l)
fail=$(grep "FAIL" 3.txt | wc -l)

echo "Overall Summary for: Sequence Duplication Levels
        ${passes} Passes
        ${warn} Warnings
        ${fail} Fails" > sdl.txt

#Overrepresented

awk 'BEGIN{FS="\t";}{if($2 == "Overrepresented sequences"){print $0, $1, $2;}}' summarymega.txt > 4.txt
total=$(grep "PASS\|FAIL\|WARN" 4.txt | wc -l)
passes=$(grep "PASS" 4.txt | wc -l)
warn=$(grep "WARN" 4.txt | wc -l)
fail=$(grep "FAIL" 4.txt | wc -l)

echo "Overall Summary for: Overrepresented sequences
        ${passes} Passes
        ${warn} Warnings
        ${fail} Fails" > os.txt


## options


PS3='Select option for more detail: '
qcchecks=("Per base sequence quality" "Per base sequence content" "Sequence Duplication Levels" "Overrepresented sequences" "Continue")
select opt in "${qcchecks[@]}"; do
    case $opt in
        "Per base sequence quality")
            echo "$opt"
            echo "$(<pbsq.txt)"
            ;;
        "Per base sequence content")
            echo "$opt"
            echo "$(<pbsc.txt)"
            ;;
        "Sequence Duplication Levels")
            echo "$opt"
            echo "$(<sdl.txt)"
            ;;
        "Overrepresented sequences")
            echo "$opt"
            echo "$(<os.txt)"
            ;;
        "Continue")
            echo "Thanks for using FastQC quality checker - Goodbye Al"
            break
            ;;
        *) echo "$REPLY is not an option Al!";;
    esac
done


# Index full sequence with bowtie

echo "Indexing reference genome"
bowtie2-build TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz congogenome

# align read pairs with ref genome
# --no-unal means reads that dont align aren't written to output. -p = no. of cores -x is ref genome 
#This needs to be a loop so...

 for f in *1.fq.gz
    do 
    base=$(basename $f "_1.fq.gz")
    echo "Aligning ${base} to reference genome"
    bowtie2 -p 16 -x congogenome -1 ${base}_1.fq.gz -2 ${base}_2.fq.gz | samtools view -b -o ${base}.bam -
    done

#Sorting and indexing files - nt neccesary for intersect so omitted from program

#Sort bam file 

#samtools sort output.bam -o output.sorted.bam

#for b in *.bam
#    do
#    base=$(basename $b ".bam")
#    samtools sort -@ 16 ${base}.bam -o ${base}sorted.bam 
#    done


#Index bam file 

#for s in *sorted.bam
 #   do 
 #   base=$(basename $s "sorted.bam")
 #   samtools index -@ 16 ${base}sorted.bam
 #   done


#Use bed tools intersect and output results to text file
echo "Counting intersects in group: Clone1 - 0h..."
bedtools intersect -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b *Clone1-0h*.bam -c > clone1-0hresults.txt

echo "Calculating mean in group: Clone1 - 0h..." 
cc10f=$(ls *Clone1-0h*.bam | wc -l)
awk -v cc10f="$cc10f" 'BEGIN{FS="\t"}{print $4, $5, $6, $7 = $6/cc10f}' clone1-0hresults.txt > clone1-0hmean.txt

###

echo "Counting intersects in group: Clone1 - 24h Uninduced..."
bedtools intersect -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b *Clone1-24h-Uninduced*.bam -c > clone1-24h-Uninduced-results.txt


echo "Calculating mean in group: Clone1 - 24h Uninduced..." 
cc10f=$(ls *Clone1-24h*.bam | wc -l)
awk -v cc10f="$cc10f" 'BEGIN{FS="\t"}{print $4, $5, $6, $7 = $6/cc10f}' clone1-24h-Uninduced-results.txt > clone1-24h-Uninduced-mean.txt

##

echo "Counting intersects in group: Clone1 - 24h Induced..."
bedtools intersect -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b *Clone1-24h-Induced*.bam -c > clone1-24h-Induced-results.txt


echo "Calculating mean in group: Clone1 - 24h Induced..." 
cc10f=$(ls *Clone1-24h*.bam | wc -l)
awk -v cc10f="$cc10f" 'BEGIN{FS="\t"}{print $4, $5, $6, $7 = $6/cc10f}' clone1-24h-Induced-results.txt > clone1-24h-Induced-mean.txt

###

echo "Counting intersects in group: Clone1 - 48h Uninduced..."
bedtools intersect -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b *Clone1-48h-Uninduced*.bam -c > clone1-48h-Uninduced-results.txt

echo "Calculating mean in group: Clone1 - 48h Uninduced..." 
cc10f=$(ls *Clone1-48h*.bam | wc -l)
awk -v cc10f="$cc10f" 'BEGIN{FS="\t"}{print $4, $5, $6, $7 = $6/cc10f}' clone1-48h-Uninduced-results.txt > clone1-48h-Uninduced-mean.txt

##

echo "Counting intersects in group: Clone1 - 48h Induced..."
bedtools intersect -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b *Clone1-48h-Induced*.bam -c > clone1-48h-Induced-results.txt

echo "Calculating mean in group: Clone1 - 48h Induced..." 
cc10f=$(ls *Clone1-48h*.bam | wc -l)
awk -v cc10f="$cc10f" 'BEGIN{FS="\t"}{print $4, $5, $6, $7 = $6/cc10f}' clone1-48h-Induced-results.txt > clone1-48h-Induced-mean.txt

###

echo "Counting intersects in group: Clone2 - 0h..."
bedtools intersect -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b *Clone2-0h*.bam -c > clone2-0hresults.txt

echo "Calculating mean in group: Clone - 0h..." 
cc10f=$(ls *Clone2-0h*.bam | wc -l)
awk -v cc10f="$cc10f" 'BEGIN{FS="\t"}{print $4, $5, $6, $7 = $6/cc10f}' clone2-0hresults.txt > clone2-0hmean.txt

###

echo "Counting intersects in group: Clone2 - 24h Uninduced..."
bedtools intersect -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b *Clone2-24h-Uninduced*.bam -c > clone2-24h-Uninduced-results.txt


echo "Calculating mean in group: Clone2 - 24h Uninduced..." 
cc10f=$(ls *Clone2-24h*.bam | wc -l)
awk -v cc10f="$cc10f" 'BEGIN{FS="\t"}{print $4, $5, $6, $7 = $6/cc10f}' clone2-24h-Uninduced-results.txt > clone2-24h-Uninduced-mean.txt

##

echo "Counting intersects in group: Clone2 - 24h Induced..."
bedtools intersect -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b *Clone2-24h-Induced*.bam -c > clone2-24h-Induced-results.txt


echo "Calculating mean in group: Clone2 - 24h Induced..." 
cc10f=$(ls *Clone2-24h*.bam | wc -l)
awk -v cc10f="$cc10f" 'BEGIN{FS="\t"}{print $4, $5, $6, $7 = $6/cc10f}' clone2-24h-Induced-results.txt > clone2-24h-Induced-mean.txt

###

echo "Counting intersects in group: Clone2 - 48h Uninduced..."
bedtools intersect -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b *Clone2-48h-Uninduced*.bam -c > clone2-48h-Uninduced-results.txt

echo "Calculating mean in group: Clone2 - 48h Uninduced..." 
cc10f=$(ls *Clone2-48h*.bam | wc -l)
awk -v cc10f="$cc10f" 'BEGIN{FS="\t"}{print $4, $5, $6, $7 = $6/cc10f}' clone2-48h-Uninduced-results.txt > clone2-48h-Uninduced-mean.txt

##

echo "Counting intersects in group: Clone2 - 48h Induced..."
bedtools intersect -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b *Clone2-48h-Induced*.bam -c > clone2-48h-Induced-results.txt

echo "Calculating mean in group: Clone2 - 48h Induced..." 
cc10f=$(ls *Clone2-48h*.bam | wc -l)
awk -v cc10f="$cc10f" 'BEGIN{FS="\t"}{print $4, $5, $6, $7 = $6/cc10f}' clone2-48h-Induced-results.txt > clone2-48h-Induced-mean.txt


#####

echo "Counting intersects in group: WT - 0h..."
bedtools intersect -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b *WT-0h*.bam -c > WT-0hresults.txt

echo "Calculating mean in group: WT - 0h..." 
cc10f=$(ls *WT-0h*.bam | wc -l)
awk -v cc10f="$cc10f" 'BEGIN{FS="\t"}{print $4, $5, $6, $7 = $6/cc10f}' WT-0hresults.txt > WT-0hmean.txt

###

echo "Counting intersects in group: WT - 24h Uninduced..."
bedtools intersect -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b *WT-24h-Uninduced*.bam -c > WT-24h-Uninduced-results.txt


echo "Calculating mean in group: WT - 24h Uninduced..." 
cc10f=$(ls *WT-24h*.bam | wc -l)
awk -v cc10f="$cc10f" 'BEGIN{FS="\t"}{print $4, $5, $6, $7 = $6/cc10f}' WT-24h-Uninduced-results.txt > WT-24h-Uninduced-mean.txt

##

echo "Counting intersects in group: WT - 24h Induced..."
bedtools intersect -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b *WT-24h-Induced*.bam -c > WT-24h-Induced-results.txt


echo "Calculating mean in group: WT - 24h Induced..." 
cc10f=$(ls *WT-24h*.bam | wc -l)
awk -v cc10f="$cc10f" 'BEGIN{FS="\t"}{print $4, $5, $6, $7 = $6/cc10f}' WT-24h-Induced-results.txt > WT-24h-Induced-mean.txt

###

echo "Counting intersects in group: WT - 48h Uninduced..."
bedtools intersect -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b *WT-48h-Uninduced*.bam -c > WT-48h-Uninduced-results.txt

echo "Calculating mean in group: WT - 48h Uninduced..." 
cc10f=$(ls *WT-48h*.bam | wc -l)
awk -v cc10f="$cc10f" 'BEGIN{FS="\t"}{print $4, $5, $6, $7 = $6/cc10f}' WT-48h-Uninduced-results.txt > WT-48h-Uninduced-mean.txt

##

echo "Counting intersects in group: WT - 48h Induced..."
bedtools intersect -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b *WT-48h-Induced*.bam -c > WT-48h-Induced-results.txt

echo "Calculating mean in group: WT - 48h Induced..." 
cc10f=$(ls *WT-48h*.bam | wc -l)
awk -v cc10f="$cc10f" 'BEGIN{FS="\t"}{print $4, $5, $6, $7 = $6/cc10f}' WT-48h-Induced-results.txt > WT-48h-Induced-mean.txt

 
#Unable to find method to divide by no. of replicates, therefore final mean result is incorect

#Below: Attempting to merge files and cxalculated fold change using awk: Unable to find workaround for "Fatal division by zero" 

#join WT-0hmean.txt clone1-0hmean.txt > newf1.txt
#join newf1.txt clone2-0hmean.txt > wtc1c20.txt
#awk 'BEGIN{FS="\t"}{print $1, $2 = $2/$2, $3 = $4/$2, $4 = $7/$2, $5 = $10/$2}' wtc1c20.txt > groupwisecomp0hunin.txt

mkdir results
cp *mean.txt /results
rm -f *fq.gz
rm -f *.bam
