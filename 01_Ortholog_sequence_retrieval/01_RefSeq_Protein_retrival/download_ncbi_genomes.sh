#!/bin/sh

set -u #exit on onbound variable
#TAXLIST=("Daphnia pulex" "Drosophila melanogaster" "Anopheles gambiae" "Pediculus humanus"
#"Ixodes scapularis" "Apis mellifera"  "Bombyx mori")
#TAXLIST=("Strigamia maritima")
TAXLIST=$@ # provide multiple taxa on the cmd-line
for TAX in "${TAXLIST[@]}" ; do
echo getting genome for: $TAX
   mkdir "$TAX"

   GENOME=$(esearch -db genome -query txid${TAX}[Organism:exp] |
       efetch -format docsum | tee "${TAX}.genome.esearch.docsum")
   ACC=`echo $GENOME | xtract -pattern DocumentSummary  -element Assembly_Accession`
   NAME=`echo $GENOME |  xtract -pattern DocumentSummary -element Assembly_Name`
   echo authoritative genome: $ACC $NAME
   RESULT=$(esearch -db assembly -query "$ACC" |
       efetch -format docsum | tee "${TAX}.assembly.esearch.docsum")
   FTPP=`echo $RESULT | xtract -pattern DocumentSummary  -element FtpPath_GenBank`
   TAXID=`echo $RESULT | xtract -pattern DocumentSummary  -element Taxid`
   echo FtpPath: $FTPP
   BASENAME=`basename $FTPP`
   FTPPATHG=$FTPP/$BASENAME'_genomic.fna.gz'
   FTPPATHP=$FTPP/$BASENAME'_protein.faa.gz'
   echo Downloading $FTPPATHG ...

   ## get genome data
   wget $FTPPATHG
   BASENAME=`basename $FTPPATHG`
   gunzip -f $BASENAME
   BASENAME=`echo $BASENAME | sed s/.gz//`
   mv $BASENAME "./$TAX/$BASENAME"

   echo Downloading $FTPPATHP ...
   ## get protein data
   wget $FTPPATHP  # this may throw an error if there is no proteome file
   BASENAME=`basename $FTPPATHP`
   gunzip -f $BASENAME
   BASENAME=`echo $BASENAME | sed s/.gz//`
   mv $BASENAME "./$TAX/$BASENAME"

   mv *.docsum ./DocumentSummary/.

done
