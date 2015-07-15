#!/usr/bin/python

import re, sys, math, os, json, argparse, ConfigParser
from collections import defaultdict

def process_maf_file(maf_path, transcipt_dict, sample_whitelist, gene_whitelist, config):

    indice_list = None
    transcript_db = None
    unparseable_aa = []
    missing_transcripts = set()

    # Used for MAGI output
    gene_to_sample = defaultdict(lambda: defaultdict())

    # Used for CoMEt/HotNet2 output
    sample_to_gene = defaultdict(set)

    exclude_mutations = set(','.split(config.get('maf', 'mutation_types_blacklist')))
    exclude_status = set(','.split(config.get('maf', 'mutation_status_blacklist')))
    exclude_validation = set(','.split(config.get('maf', 'validation_status_blacklist')))

    with open(maf_path) as maf_file:
        for line in maf_file:

            # Ignore comments
            if line.startswith('#'):
                continue

            line = line.rstrip().split('\t')

            # Indentify indices of relevant columns
            if not indice_list:
                indice_list = define_indices(line)
                continue

            [gene, sample, variant_class_type, variant_type, valid_status,
                mutation_status, transcript_id, codon, aa_change] = [line[i] for i in indice_list]

            # If a whitelist is provided, skip any genes/samples not in the list
            if gene_whitelist and gene not in gene_whitelist:
                continue
            if sample_whitelist and sample not in sample_whitelist:
                continue

            # Identify which database the transcript is from. If none can be found,
            # add transcript id to missing transcript set and move to next line
            if not transcript_db:
                for database_name, transcript_list in transcipt_dict.iteritems():
                    if transcript_id in transcript_list:
                        transcript_db = database_name
                        break 
                if not transcript_db:
                    missing_transcripts.add(transcript_id)
                    continue


            # take only first three segments of name if sample is from TCGA
            if sample[:4] == 'TCGA':
                sample = '-'.join((sample.split('-'))[:3])

            # Javascript can't have "." in gene names
            gene = gene.replace(".", "-")

            if (mutation_status not in exclude_status 
                    and variant_class_type not in exclude_mutations
                    and valid_status not in exclude_validation):
                try:
                    original_amino_acid, new_amino_acid, amino_acid_location = get_amino_acid_change(
                                                aa_change, variant_type, variant_class_type, codon)
                except ValueError:
                    unparseable_aa.append(aa_change)
                    continue

                if original_amino_acid and new_amino_acid and amino_acid_location:
                    # check if transcript is in database, if so add to final
                    pass



    return 0

def define_indices(header_line):
    '''
    Using the column names, find the index of each relevant piece of data.
    '''
    # Dict for easy matching of column name to index, e.g. 'hugo':0
    # keys must be set to lowercase to match column names
    indice_dict = {key.lower():index for index, key in enumerate(header_line)}

    # List of possible column names for data we care about
    gene_name = ["hugo_symbol"]  # GENE NAME
    sample = ["tumor_sample_barcode"]	# SAMPLE NAME
    variant_class_type = ["variant_classification"]	# SILENT OR NOT
    valid_status = ["validation_status"] # FOR WILDTYPE/INVALID CHECKING
    mut_stat = ["mutation_status"] # FOR GERMLINE CHECKING
    variant_type = ["variant_type"] # MUTATION TYPE
    codon = ["codon_change", "c_position", "c_position_wu", "chromchange", "amino_acids"] # CODON
    aachange = ["protein_change", "amino_acid_change", "aachange", "amino_acid_change_wu", "hgvsp_short"] # Protein change
    transcript_id = ["refseq_mrna_id", "transcript_name", "transcript_name_wu", "transcriptid", "transcript_id"]

    header_name_list = [gene_name, sample, variant_class_type, variant_type, valid_status, 
                        mut_stat, transcript_id, codon, aachange]

    indice_list = []
    for header_options in header_name_list:
        existing_headers = set(header_options) & set(indice_dict.keys())
        if len(existing_headers) == 0:
            raise IndexError("Names %s not found" % str(header_options))
        else:
            indice_list.append(indice_dict[existing_headers.pop()])
            
    assert len(indice_list) == 9
    return indice_list

def get_amino_acid_change(aa_change, variant_type, variant_class_type, codon):
    '''
    Attempt to parse amino acid change and change location.
    '''
    aa_new, aa_original, aa_location = None, None, None

    if aa_change and aa_change not in ("N/A", ".", "NULL"):
        if variant_type in ("SNP", "DNP", "TNP", "ONP"):
            if variant_class_type in ("Missense_Mutation", "Nonsense_Mutation",
                                      "RNA", "Translation_Start_Site", "Start_Codon_SNP"):
                aa_original, aa_new, aa_location = snp_mutation(aa_change)
            elif variant_class_type == "Splice_Site":
                aa_original, aa_new, aa_location = splice_site_mutation(aa_change, codon)
            else:
                sys.stderr.write("New mutation effect can't be parsed: %s\n" % variant_class_type)
                exit(1)

        elif variant_type in ("DEL", "INS"):
            aa_original, aa_new, aa_location = ins_del_mutation(aa_change, codon)
        else:
            sys.stderr.write("New mutation type can't be parsed: %s\n" % variant_type)
            exit(1)

    return aa_original, aa_new, aa_location

def get_parser():
    '''
    Parse arguments.
    '''

    parser = argparse.ArgumentParser(description='Parse MAF files')

    parser.add_argument('-ia', '--inactive_types', nargs="*",
                        type=str, help="Inactivating mutation types.",
                        default=["frame_shift_ins", "nonstop_mutation", "nonsense_mutation",
                                 "splice_site", "frame_shift_del"])

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

    with open(config.get('transcript', 'database')) as t_file:
        transcript_dict = json.load(t_file)

    process_maf_file(maf_file, transcript_dict, sample_whitelist, gene_whitelist, config)

    _ = args.output_prefix


if __name__ == '__main__':
    run(get_parser().parse_args(sys.argv[1:]), get_config())


## BELOW IS OLD CODE, included just to get this working. Need to survey
## MAF files to get an idea of how to improve amino acid change info extraction
################################################################################
# Parse the different mutation formats to get the amino acid change
def snp_mutation( aa_change ):
    # parsing SNP, DNP and TNP
    aao, aan, aaloc = "", "", ""
    
    # capture formats like p.A222Q, p.*222Q, p.222Q, p.A222*
    if re.search(r'p.([a-zA-z*_]+)?(\d+)([a-zA-Z*_]+)', aa_change): 
        aao = re.search(r'p.([a-zA-Z*_]+)?(\d+)([a-zA-Z*_]+)', aa_change).group(1)
        aan = re.search(r'p.([a-zA-Z*]+)?(\d+)([a-zA-Z*_]+)', aa_change).group(3)
        aaloc = re.search(r'p.([a-zA-Z*_]+)?(\d+)([a-zA-Z*_]+)', aa_change).group(2)  
        if not aao:
            aao = ""
    else:
        raise ValueError("Error format can't be parsed: " + aa_change + "\n")

    return aao, aan, aaloc

def splice_site_mutation(aa_change, codon):
    aao, aan, aaloc = "", "", ""

    # capture formats like p.A222_splice or e20-1 or p.321_splice
    if re.search(r'^p\.[A-Z]\d+_splice', aa_change): #p.A222_splice
        aao = aa_change[2]
        aan = 'splice'
        aaloc = re.search(r'^p\.[A-Z](\d+)_splice', aa_change).group(1)
    elif re.search(r'^p\.\d+_splice', aa_change): # p.333_splice
        aao = 'splice'
        aan = 'splice'
        aaloc = re.search(r'^p\.(\d+)_splice', aa_change).group(1)
    
    elif re.search(r'^e', aa_change) and re.search(r'c\.(\d+)(\+|\-)(\d+)', codon):
        aao = 'splice'
        aan = 'splice'
        aaloc = int(math.ceil( float(re.search(r'^c\.(\d+)(\+|\-)(\d+)', codon).group(1)) / 3))
    else:
        raise ValueError("Error format can't be parsed: " + aa_change + "\n")

    return aao, aan, aaloc

def ins_del_mutation(aa_change, codon):
    aao, aan, aaloc = "", "", ""

    # capture formats like p.A222Q, p.*222Q, p.222Q, p.A222*, inframe_shift, or splice format
    if re.search(r'p.([a-zA-Z*_]+)?(\d+)([a-zA-Z*_]+)', aa_change):
        aao = re.search(r'p.([a-zA-Z*_]+)?(\d+)([a-zA-Z*_]+)', aa_change).group(1)
        aan = re.search(r'p.([a-zA-Z*_]+)?(\d+)([a-zA-Z*_]+)', aa_change).group(3)
        aaloc = re.search(r'p.([a-zA-Z*_]+)?(\d+)([a-zA-Z*_]+)', aa_change).group(2)
        if not aao:
            aao = ""
    elif re.search(r'p.(\d+)_(\d+)(\w+)', aa_change):
        aao = re.search(r'p.(\d+)_(\d+)(\w+)', aa_change).group(3)
        aan = re.search(r'p.(\d+)_(\d+)(\w+)', aa_change).group(3)
        aaloc = re.search(r'p.(\d+)_(\d+)(\w+)', aa_change).group(1)                                                
        aaloc2 = re.search(r'p.(\d+)_(\d+)(\w+)', aa_change).group(2)
        if aaloc == aaloc2:
            aaloc2 = -1
    elif re.search(r'p.-(\d+)fs', aa_change):
        aaloc = re.search(r'p.-(\d+)fs', aa_change).group(1)
        return None, aaloc, None
    elif re.search(r'^e', aa_change) and re.search(r'^c\.(\d+)(\+|\-)(\d+)', codon): 
        aao = 'splice'
        aan = 'splice'
        aaloc =  math.ceil(float(re.search(r'^c\.(\d+)(\+|\-)(\d+)', codon).group(1))/3)                
    else:
        raise ValueError("Error format can't be parsed: " + aa_change + "\n")
    
    return aao, aan, aaloc