#!/usr/bin/env bash

export METADATA_FILE=$1

# Download fasta sequences from a list of sequence accessions (skip first line)
echo "> Running download script..."
mkdir -p seq
mkdir -p output
{
  read
    while IFS=$'\t' read -r accession country collection_date host strain isolate
    do
      curl "https://www.ebi.ac.uk/ena/browser/api/fasta/${accession}?download=true" --output seq/${accession}.fasta
    done
} < ${METADATA_FILE}
echo "> Running download script... [DONE]"

echo "> Combining all sequences to a single fasta file and copying to hMPXV..."
cat seq/* > sequences.fasta
cp sequences.fasta data/hMPXV
echo "> Combining all sequences to a single fasta file and copying to hMPXV... [DONE]"

echo "> Running nextclade assignment on sequences..."
nextclade run --input-dataset data/hMPXV --output-all=output sequences.fasta
echo "> Running nextclade assignment on sequences... [DONE]"