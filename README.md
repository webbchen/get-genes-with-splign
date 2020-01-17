# get-genes-with-splign

This pipeline describes recovery of CDS of genes by aligning CDS with [NCBI splign](https://www.ncbi.nlm.nih.gov/sutils/splign/splign.cgi) to genomic sequence.
CDS recovery with splign is more accurate than using coordinates obtained with [NCBI blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi) for example which might lose the ends of sequences, making accurate translation impossible.


### rename contigs and assemblies in assembly to make downstream processing easier

The asssemblies are renamed and the identifiers simplified
ZtST99CH3D7_genomic.fna
ZtST99CH3D1_genomic.fna
ZtST99CH1E4_genomic.fna
ZtST99CH1A5_genomic.fna

The search target is a file with CDSs of genes within a QTL:
QTL_CDSs.fasta

the contigs or chromosomes are renamed as follows:
```
for assembly in /projects/Z_tritici/assemblies_splign/Zt*_genomic.fna
do
ID=$(echo $(basename ${assembly}) | sed -E 's/(Zt[A-Z0-9]+)(_genomic.fna)/\1/')
echo $ID
python /projects/Z_tritici/final_assemblies/replace_fasta_ID2.py -i $assembly -I chromosome -o /projects/Z_tritici/assemblies_splign/${ID}_renamed.fasta
done
```
## Splign
### databases
splign requires a blast database to be made for the genomic sequence as well as for the search target CDSs

```
for assembly in /projects/Z_tritici/assemblies_splign/Zt*CH*_renamed.fasta
do
echo $assembly
makeblastdb -dbtype nucl -parse_seqids -in $assembly
done

makeblastdb -dbtype nucl -parse_seqids -in /projects/Z_tritici/QTL_CDSs.fasta
```
### compart
Compart produces pre-alignmnets to be fed into splign.
This and the following steps should not be run on servers where low read-write speeds are an issue as these steps produce lots of intermediary files.

```
mkdir /tmp/annew/compart-files
for assembly in /projects/Z_tritici/assemblies_splign/Zt*CH*_renamed.fasta
do
ID=$(basename $assembly)
echo $ID
/scratch/software/ncbi_cxx--22_0_0/GCC630-DebugMT64/bin/compart -qdb /projects/Z_tritici/QTL_CDSs.fasta -sdb $assembly > /tmp/annew/compart-files/${ID}.compartments
done
```
### make directories with individual fasta sequences and indices for each isolate


```
# catenate assembly with CDSs to be searched
# make a directory for the split files for each isolate

for assembly in /projects/Z_tritici/assemblies_splign/Zt*CH*_renamed.fasta
do
ID=$(echo $(basename ${assembly}) | sed -E 's/(Zt[A-Z0-9]+)(_renamed.fasta)/\1/')
echo $ID
mkdir /tmp/annew/${ID}QTL_CDS_fastas
cat /projects/Z_tritici/QTL_CDSs.fasta $assembly > /tmp/annew/${ID}_splign_fastas.fasta
done

# split catenated files into individual fastas

for file in /tmp/annew/Zt*_splign_fastas.fasta
do
ID=$(echo $(basename ${file}) | sed -E 's/(Zt[A-Z0-9]+)(_splign_fastas.fasta)/\1/')
echo $ID
python2 /projects/Z_tritici/split_fastas.py $file -p /tmp/annew/${ID}QTL_CDS_fastas/${ID}_
done

# index directories with single files

for directory in /tmp/annew/Zt*QTL_CDS_fastas
do
/scratch/software/ncbi_cxx--22_0_0/GCC630-DebugMT64/bin/splign -mklds $directory
done
```
### run splign itself

Outputs are a text file similar to tabular blast output and a file showing the alignments of the CDSs\' exons including translations where an open reading frame is found
Alignments are shown for both strands in 3-5 and 5-3 orientation.

```
mkdir /tmp/annew/splign_results
for directory in /tmp/annew/Zt*CH*QTL_CDS_fastas
do
ID=$(echo $(basename ${directory}) | sed -E 's/(Zt[A-Z0-9]+)(QTL_CDS_fastas)/\1/')
echo $ID
/scratch/software/ncbi_cxx--22_0_0/GCC630-DebugMT64/bin/splign -ldsdir $directory -comps /tmp/annew/compart-files/${ID}_renamed.fasta.compartments -aln /tmp/annew/splign_results/${ID}.aln > /tmp/annew/splign_results/${ID}splign_output.txt 
done
```

reformat the splign output table to .bed and .faidx formats for downstream processing.

```
for file in /tmp/annew/splign_results/*splign_output.txt
do
ID=$(echo $(basename ${file}) | sed -E 's/(Zt[A-Z0-9]+)(splign_output.txt)/\1/')
echo $ID
python /projects/Z_tritici/splign_to_faidx.py $file -o /tmp/annew/splign_results/${ID}splign_reformat                                   
done
```
### recover CDS sequences, concatenate and translate them

Use bedtools to get sequences for cdss out:

```
for assembly in /projects/Z_tritici/assemblies_splign/Zt*CH*_renamed.fasta
do
ID=$(echo $(basename ${assembly}) | sed -E 's/(Zt[A-Z0-9]+)(_renamed.fasta)/\1/')
echo $ID
/scratch/software/bedtools2/bin/fastaFromBed -s -name -fi $assembly -bed /tmp/annew/splign_results/${ID}splign_reformat.bed | fold -w 60 > /tmp/annew/QTL_CDS/${ID}splign_CDS.fasta
done
```
These CDSs are in either direction, depending on which strand they were found. they need to be catenated and , where applicable, the reverse complement found for translation in the fiorst frame

```
mkdir /tmp/annew/test_mRNA
for file in /tmp/annew/QTL_CDS/Zt*CH*splign_CDS.fasta
do
ID=$(echo $(basename ${file}) | sed -E 's/(Zt[A-Z0-9]+)(splign_CDS.fasta)/\1/')
echo $ID
python /projects/Z_tritici/splign-fastaFromBed-CDSs_to_mRNA.py $file -o /tmp/annew/test_mRNA/${ID}_splign_mRNA.fasta
done
```

translate and prefix isolate name to sequence identifiers as follows:
```
>Zt1_QTL_7_1
```

```
mkdir /tmp/annew/QTL_translated_test
for file in test_mRNA/Zt*CH*_splign_mRNA.fasta
do
ID=$(echo $(basename ${file}) | sed -E 's/(Zt[A-Z0-9]+)(_splign_mRNA.fasta)/\1/')
echo $ID
/scratch/software/emboss-latest/EMBOSS-6.6.0/bin/transeq -trim -sequence $file -outseq QTL_translated_test/${ID}_QTL_protein.fasta
done

for file in QTL_translated_test/Zt*CH*_QTL_protein.fasta
do
ID=$(echo $(basename ${file}) | sed -E 's/(Zt[A-Z0-9]+)(_QTL_protein.fasta)/\1/')
echo $ID
python fasta_IDv1.py -i $file -I $ID -o QTL_translated_renamed/${ID}_QTL_protein_renamed.fasta
done
```

### Catenate the renamed protein files and sort like proteins into separate files

```
cat QTL_translated_renamed/*QTL_protein_renamed.fasta > all_QTL_proteins.fasta
```

Make a list of all gene identifiers present in the file and run the script whcih re-sorts the files.


```
sed -E 's/(^>Zt[A-Z0-9]+_)([A-Za-z0-9]+[_]*)/\2/' all_QTL_proteins.fasta | grep "_" | cut -d' ' -f1 | sort | uniq > genelist.txt

mkdir /tmp/annew/QTL_parcels_ext
python make_parcels_for_alignment_large_v2.py all_QTL_proteins.fasta genelist.txt -a protein -o QTL_parcels_ext
rm fasta_parcels.idx
```
### clean files up

The sfilenames as well as the sequence identifiers themselves contain braces, - and + which ought to be removed

rename files:
```
rename 's/\(//' QTL_parcels_ext/*.fasta
rename 's/\)//' QTL_parcels_ext/*.fasta
rename 's/\+/_plus/' QTL_parcels_ext/*.fasta
rename 's/\-/_minus/' QTL_parcels_ext/*.fasta
```
combine the plus and minus orientations for each protein into one file:
```
mkdir QTL_parcels_combined_ext
for file in /tmp/annew/QTL_parcels_ext/*plus_1.fasta
do
ID=$(echo $(basename ${file}) | sed -E 's/([A-Za-z0-9_]*)(plus_1.fasta)/\1/')
echo $ID
cat $file $(ls $(dirname $file)/${ID}minus_1.fasta) > /tmp/annew/QTL_parcels_combined_ext/${ID}combined.fasta
done
```
Remove any empty sequences which only have headers and an \'unknown description\':
Finally, remove troublesome characterss from headers:

```
mkdir /tmp/annew/QTL_parcels_combined_cleaned_ext
for file in /tmp/annew/QTL_parcels_combined_ext/*.fasta
do
ID=$(echo $(basename ${file}) | sed -E 's/(.*)(.fasta)/\1/')
echo $ID
sed '/unknown description/d' $file > /tmp/annew/QTL_parcels_combined_cleaned_ext/${ID}_cleaned.fasta
done

# remove protected chars from headers:

for file in /tmp/annew/QTL_parcels_combined_cleaned_ext/*.fasta
do
sed -E -i 's/\(\-\)/_minus/g' $file
sed -E -i 's/\(\+\)/_plus/g' $file
done
```
The resulting files should contain the protein sequences for each of the query proteins from each of the isolates. There can be two copies if the CDS was found on both strands.
