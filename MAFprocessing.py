#!/usr/bin/python

import re, sys, math, os, json, argparse, ConfigParser
from collections import defaultdict, Counter

def process_maf_file(MAF, transcipt_dict, sample_whitelist, gene_whitelist, config):
    # 1. Identify file indices
    #
    indice_list = None

    gene_to_sample = defaultdict(lambda: defaultdict())
    sample_to_gene = defaultdict(lambda: defaultdict())

    exclude_classes, exclude_mutations = get_mutation_exclusions(config)

    with open(MAF) as maf_file:
        for line in maf_file:
            if line.startswith('#'):
                continue

            line = line.rstrip().split('\t')

            # Indentify indices of relevant columns
            if not indice_list:
                indice_list = define_indices(line)
                continue

            [gene, sample, variant_class_type, mutation_type, valid_stat,
                mut_stat, location, transcript_id, codon, aa_change] = [line[i] for i in indice_list]

            # If a whitelist is provided, skip any genes/samples not in the list
            if gene_whitelist and gene not in gene_whitelist:
                continue
            if sample_whitelist and sample not in sample_whitelist:
                continue

            # take only first three segments of name if sample is from TCGA
            if sample[:4] == 'TCGA':
                sample = '-'.join((sample.split('-'))[:3])

            # Javascript can't have "." in gene names
            gene = gene.replace(".", "-")



            if (variant_class_type not in exclude_classes 
             and mutation_type not in exclude_mutations):
                original_amino_acid, new_amino_acid, amino_acid_location = get_amino_acid_change(
                                                aa_change, mutation_type, variant_class_type, codon)
                




    return 0

def get_mutation_exclusions(config):
    '''
    Return the two exclusion sets, mutation types and mutation classes.
    If none are provided in the config, defaults are used.
    '''

    if not config.get('maf', 'mutation_types'):
        exclude_mutations = set(["Silent", "Intron", "3'UTR", "5'UTR", "IGR", "Intron", "lincRNA"])
    else:
        exclude_mutations = set(','.split(config.get('maf', 'mutation_types')))

    if not config.get('maf', 'mutation_classes'):
        exclude_classes = set(["Germline"])
    else:
        exclude_classes = set(','.split(config.get('maf', 'mutation_classes')))

    return exclude_classes, exclude_mutations

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
    variant_class_type = ["variant_classification"]	# SILENT OR NOT
    valid_stat = ["validation_status"] # FOR WILDTYPE CHECKING
    mut_stat = ["mutation_status"] # FOR GERMLINE CHECKING
    mutation_type = ["variant_type"] # MUTATION TYPE
    location = ["start_position"] # FOR POSITION OF MUTATION
    codon = ["codon_change", "c_position", "c_position_wu", "chromchange", "amino_acids"] # CODON
    aachange = ["protein_change", "amino_acid_change", "aachange", "amino_acid_change_wu", "hgvsp_short"] # Protein change
    transcript_id = ["refseq_mrna_id", "transcript_name", "transcript_name_wu", "transcriptid", "transcript_id"]

    header_name_list = [hugo, sample, variant_class_type, mutation_type, valid_stat, 
                        mut_stat, location, transcript_id, codon, aachange]

    indice_list = []
    for header_options in header_name_list:
        existing_headers = set(header_options) & set(indice_dict.keys())
        if len(existing_headers) == 0 :
            raise IndexError("Names %s not found" % str(header_options))
        else:
            indice_list.append(indice_dict[existing_headers.pop()])
    assert len(indice_list) == 10
    return indice_list

def get_amino_acid_change(aa_change, mutation_type, variant_class_type, codon):
    '''
    Attempt to parse amino acid change and change location.
    '''

    return 1, 2, 3

def get_parser():
    '''
    Parse arguments.
    '''
    transcript_file='transcript-lengths.json'

    parser = argparse.ArgumentParser(description='Parse MAF files')
    parser.add_argument('-ia', '--inactive_types', nargs="*", type=str, help="Inactivating mutation types.",
        default=["frame_shift_ins", "nonstop_mutation", "nonsense_mutation", "splice_site", "frame_shift_del"])
    parser.add_argument('-o', '--output_prefix', default=None, help='Output prefix.')

    return parser

def get_config():
    config = ConfigParser.SafeConfigParser()
    found = config.read('config.cfg')
    if not found:
        raise IOError("Error: Config file not found!")

    if not os.path.isfile(config.get('maf', 'file')):
        raise IOError('Error: MAF file not found. Please check location in configuration file')

    if not os.path.isfile(config.get('transcript', 'database')):
        raise IOError('Error: Transcript database file not found. Please check location in configuration file')

    return config

def run(args, config):
    '''
    Main function.
    '''

    sample_whitelist = config.get('whitelists', 'sample')
    gene_whitelist = config.get('whitelists', 'gene')
    maf_file = config.get('maf', 'file')

    transcript_dict = json.load(open(config.get('transcript', 'database')))
    process_maf_file(maf_file, transcript_dict, sample_whitelist, gene_whitelist, config)

    _ = args.output_prefix


if __name__ == '__main__':
    run(get_parser().parse_args(sys.argv[1:]), get_config())
