#!/usr/bin/python

import re, sys, math, os, json, argparse, ConfigParser
from collections import defaultdict, Counter

def process_maf_file(MAF, transcipt_dict, sample_whitelist, gene_whitelist):
    # 1. Identify file indices
    #
    indice_list = None

    with open(MAF) as maf_file:
        for line in maf_file:
            if line.startswith('#'):
                continue

            line = line.rstrip().split('\t')

            # Indentify indices of relevant columns
            if not indice_list:
                indice_list = define_indices(line)

            [gene, sample, class_type, mutation_type, valid_stat,
                mut_stat, location, transcript_id, codon, aachange] = [line[i] for i in indice_list]

            # take only first three segments of name if sample is from TCGA
            if sample[:4] == 'TCGA':
                sample = '-'.join((sample.split('-'))[:3])

            # Javascript can't have "." in gene names
            gene = gene.replace(".", "-") 

            return 0


def define_indices(header_line):
    '''
    Using the column names, find the index of each relevant piece of data.
    '''
    # Dict for easy matching of column name to index, e.g. 'hugo':0
    # keys must be set to lowercase to match column names
    indice_dict = {key.lower():index for index, key in enumerate(header_line)}

    # List of possible column names for data we care about
    hugo = ["hugo_symbol"]  # GENE NAME
    sample = ["tumor_sample_barcode"]	# SAMPLE NAME
    class_type = ["variant_classification"]	# SILENT OR NOT
    valid_stat = ["validation_status"] # FOR WILDTYPE CHECKING
    mut_stat = ["mutation_status"] # FOR GERMLINE CHECKING
    mutation_type = ["variant_type"] # MUTATION TYPE
    location = ["start_position"] # FOR POSITION OF MUTATION
    codon = ["codon_change", "c_position", "c_position_wu", "chromchange", "amino_acids"] # CODON
    aachange = ["protein_change", "amino_acid_change", "aachange", "amino_acid_change_wu", "hgvsp_short"] # Protein change
    transcript_id = ["refseq_mrna_id", "transcript_name", "transcript_name_wu", "transcriptid", "transcript_id"]

    header_name_list = [hugo, sample, class_type, mutation_type, valid_stat, mut_stat,
                        location, transcript_id, codon, aachange]

    indice_list = []
    for header_options in header_name_list:
        existing_headers = set(header_options) & set(indice_dict.keys())
        if len(existing_headers) == 0 :
            raise IndexError("Names %s not found" % str(header_options))
        else:
            indice_list.append(indice_dict[existing_headers.pop()])
    assert len(indice_list) == 10
    return indice_list

def get_amino_acid_change():
    pass

def get_parser():
    '''
    Parse arguments.
    '''
    transcriptFile='transcript-lengths.json'

    parser = argparse.ArgumentParser(description='Parse MAF files')
    parser.add_argument('-mf', '--maf_file', required=True, help='MAF file.')
    parser.add_argument('-sf', '--sample_file', required=False, help="Sample whitelist.")
    parser.add_argument('-tf', '--transcript_file', default=transcriptFile, help="JSON file of transcript lengths.")
    parser.add_argument('-ia', '--inactive_types', nargs="*", type=str, help="Inactivating mutation types.",
        default=["frame_shift_ins", "nonstop_mutation", "nonsense_mutation", "splice_site", "frame_shift_del"])
    parser.add_argument('-o', '--output_prefix', default=None, help='Output prefix.')

    return parser

def get_config():
    configPars = ConfigParser.SafeConfigParser()
    found = configPars.read('config.cfg')
    if not found:
        raise IOError("Config file not found!")
    raw_config = configPars._sections

    return raw_config

def run(args,config):
    sample_whitelist = config['Whitelists']['sample']
    gene_whitelist = config['Whitelists']['gene']

    transcript_dict = json.load(open(args.transcript_file))
    process_maf_file(args.maf_file, transcript_dict, sample_whitelist,gene_whitelist)


if __name__ == '__main__': run(get_parser().parse_args(sys.argv[1:]),get_config())
