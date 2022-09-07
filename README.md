# ena-mpx-prep

This repository details preparatory steps required to host an INSDC-focused Nextstrain Monkeypox instance.

###Usage

1. Ensure that you have installed requirements.
2. Run `nextclade dataset get --name 'hMPXV' --output-dir 'data/hMPXV'`. This downloads the pre-requisites in order to assign clades for your sequences.
3. Run `python metadata_processing.py`. This will download metadata and sequences and then prepare the metadata for ingestion by the nextstrain workflow.

###Requirements

- Python3+
- Bioconda
- Nextstrain, including nextclade