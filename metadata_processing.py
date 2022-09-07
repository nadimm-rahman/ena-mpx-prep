#!/usr/bin/env python

__author__ = 'Nadim Rahman'

import datetime, io, re, requests, subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from requests.auth import HTTPBasicAuth


ena_searches = {
    'sequence': {'search_fields': ['country', 'collection_date', 'host', 'strain', 'isolate', 'first_public', 'collected_by'], 'query': 'tax_tree(10244) AND country="*" AND collection_date="*"', 'result_type': 'sequence', 'data_portal': 'ena'},
}
format = 'fasta'

class retrieve_data:
    def __init__(self, ena_search, username="", password=""):
        self.ena_search = ena_search        # A dictionary of that includes: query and result type (to search), data portal (to search in) and search fields to return
        self.username = username
        self.password = password
        self.BASE_PORTAL_API_SEARCH_URL = 'https://www.ebi.ac.uk/ena/portal/api/search'
        self.ena_headers = {
            'accept': '*/*',
        }       # Define headers for the requests

    def req(URL, headers, params, user="", password=""):
        """
        Run a request and retrieve the output
        :param URL: URL used for the search
        :param headers: Headers to be used in the search
        :param params: Parameters for the request
        :param user: Username (only if authentication required)
        :param password: Username password (only if authentication required)
        :return: A response object with the results of the search
        """
        if user == "":
            # If a username was not provided
            response = requests.get(URL, headers=headers, params=params)  # No authentication required in the query
        else:
            # If a username was provided
            response = requests.get(URL, headers=headers, params=params,
                                    auth=HTTPBasicAuth(user, password))  # Authentication required in the query
            print(response.url)
        return response

    def build_request_params(**kwargs):
        """
        Build parameters for the request search
        :param query_url: query URL to use in the search
        :param start: starting point for results to be retrieved
        :return: A tuple of tuples of the parameters to be used in the request search
        """
        search_params = {}
        for key, value in kwargs.items():
            if type(value) is list:
                search_params[key] = ",".join(value)
            else:
                search_params[key] = value
        return search_params

    def download_fastas(self, date):
        """
        Download FASTA sequences from list of accessions. Invoke shell script.
        :param date: Date of analysis download
        :return:
        """
        ena_metadata = 'ENA_Search_{}_{}.tsv'.format(self.ena_search['result_type'], date)
        process = subprocess.Popen(['./download_fasta.sh', ena_metadata],
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        print(str(stdout))
        print("[ERRORS]: "+str(stderr))

    def coordinate_retrieval(self):
        """
        Run the retrieval of ENA data
        :return: Data frame
        """
        today = datetime.date.today()
        date = today.strftime('%d%m%Y')
        print('> Running data request... [{}]'.format(datetime.datetime.now()))
        if 'authentication' in self.ena_search.keys():
            # If there is an authentication flag for the result type to search for
            self.ena_search_params = retrieve_data.build_request_params(dataPortal=self.ena_search['data_portal'],
                                                          fields=self.ena_search['search_fields'],
                                                          result=self.ena_search['result_type'],
                                                          dccDataOnly=True,
                                                          limit=0)  # Create the parameter tuple
            print(self.ena_search_params)
            self.ena_search_result = retrieve_data.req(self.BASE_PORTAL_API_SEARCH_URL, self.ena_headers, self.ena_search_params, self.username, self.password)
        else:
            self.ena_search_params = retrieve_data.build_request_params(dataPortal=self.ena_search['data_portal'],
                                                          fields=self.ena_search['search_fields'],
                                                          query=self.ena_search['query'],
                                                          result=self.ena_search['result_type'],
                                                          excludeAccessions='LC722946',
                                                          limit=0)
            print(self.ena_search_params)
            self.ena_search_result = retrieve_data.req(self.BASE_PORTAL_API_SEARCH_URL, self.ena_headers, self.ena_search_params)      # Search the query
        self.ena_results = pd.read_csv(io.StringIO(self.ena_search_result.content.decode('UTF-8')), sep="\t")      # Save results in a dataframe
        output_file = 'ENA_Search_{}_{}.tsv'.format(self.ena_search['result_type'], date)
        self.ena_results.to_csv(output_file, sep="\t", index=False)      # Save search results to a dataframe
        print('> Running data request... [DONE] [{}]'.format(datetime.datetime.now()))

        print('> Downloading FASTA sequences... [{}]'.format(datetime.datetime.now()))
        self.download_fastas(date)
        print('> Downloading FASTA sequences... [DONE] [{}]'.format(datetime.datetime.now()))
        return self.ena_results


class process_fastas:
    def __init__(self, sequence_metadata, sequence_file, format):
        self.sequence_metadata = sequence_metadata
        self.sequence_file = sequence_file
        self.format = format

    def get_id(self, record):
        result = re.search('\|(.*)\|', record.id)
        return result.group(1)

    def get_metadata(self):
        accession_info = self.sequence_metadata.loc[self.sequence_metadata['accession'] == self.id]

        # Handle fields that must contain information by default
        country = accession_info.iloc[0]['country']
        collection_date = accession_info.iloc[0]['collection_date']
        host = accession_info.iloc[0]['host']

        # Obtain appropriate strain information
        if accession_info['strain'].isnull().values.any():
            if accession_info['isolate'].isnull().values.any():
                strain = id
            else:
                strain = accession_info.iloc[0]['isolate']
        else:
            strain = accession_info.iloc[0]['strain']
        return country, collection_date, host, strain

    def create_record(self, seq):
        new_record = SeqRecord(
            Seq(seq),
            id=str(self.id),
            description=""
        )
        return new_record

    def process(self):
        records = []            # List for all sequence records
        for record in SeqIO.parse(self.sequence_file, self.format):
            self.id = self.get_id(record)
            seq = record.seq

            print('*' * 50)
            print(self.id)
            print('*' * 50)

            # country, collection_date, host, strain = self.get_metadata()
            new_record = self.create_record(seq)
            records.append(new_record)
        SeqIO.write(records, "insdc_cleaned_sequences.fasta", "fasta")            # Save all sequence records in FASTA file


if __name__ == '__main__':
    for key, value in ena_searches.items():
        data_retrieval = retrieve_data(ena_searches[key])     # Instantiate class with information
        ena_results = data_retrieval.coordinate_retrieval()

    # ena_results = pd.read_csv('ENA_Search_sequence_07092022.tsv', sep="\t")       # Temp: for bug fixing

    nextclade_metadata = pd.read_csv('output/nextclade.tsv', sep="\t")
    nextclade_metadata['accession'] = nextclade_metadata['seqName'].str.extract('ENA\|([A-Z0-9]+)\|.*$').fillna('')
    nextclade_metadata['genbank_accession_rev'] = nextclade_metadata['seqName'].str.extract('([A-Z0-9]+\.[0-9]).+$').fillna('')

    ena_results[['country', 'region']] = ena_results['country'].str.split(':', expand=True)
    ena_results = ena_results.rename(columns={'collection_date': 'date', 'collected_by': 'institution', 'first_public': 'date_submitted'})

    # Process fasta sequence files
    processor = process_fastas(ena_results, 'sequences.fasta', format)            # Instantiate class with information
    process = processor.process()

    merged_meta = ena_results.merge(nextclade_metadata, how='inner', on='accession')
    merged_meta.to_csv('metadata.tsv', sep="\t", index=False)