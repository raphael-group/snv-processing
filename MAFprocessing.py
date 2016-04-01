#!/usr/bin/python

import re, sys, math, os, json, argparse, ConfigParser
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

def process_maf_file(maf_path, transcript_dict, sample_whitelist_file, gene_whitelist_file, config):


    if gene_whitelist_file:
        gene_whitelist = parse_whitelist(gene_whitelist_file)
    if sample_whitelist_file:
        sample_whitelist = parse_whitelist(sample_whitelist_file)

    required_indices = None
    magi_indices = None
    stats = {'total_mutations':0, 'processed_mutations':0, 'missing_transcripts':set(), 
             'samples':defaultdict(lambda: 0), 'genes':defaultdict(lambda: 0), 'mutation_types':defaultdict(lambda: 0),
             'unknown_mutations':set()}

    # Used for MAGI output
    gene_to_sample = defaultdict(lambda: defaultdict(list))

    # Used for CoMEt/HotNet2 output
    sample_to_gene = defaultdict(set)
    exclude_mutations = set(config.get('options', 'mutation_types_blacklist').lower().split(' '))
    exclude_status = set(config.get('options', 'mutation_status_blacklist').lower().split(' '))
    exclude_validation = set(config.get('options', 'validation_status_blacklist').lower().split(' '))

    with open(maf_path) as maf_file:
        for line in maf_file:

            # Ignore comments
            if line.startswith('#'):
                continue

            stats['total_mutations'] += 1

            line = line.rstrip().split('\t')

            # Indentify indices of relevant columns
            if not required_indices:
                try:
                    required_indices, magi_indices = define_indices(line)
                except IndexError:
                    print "\nCould not parse all indices of " + str(config.get('options','file'))+"\n"
                    print "Fatal error occured, required headers missing. Exiting."
                    exit(1)
                continue

            # populate data from line
            try:
                [gene, sample, variant_class_type, variant_type, valid_status, 
                    mutation_status] = [line[i] for i in required_indices]
                if magi_indices:
                    [transcript_id, codon, aa_change] = [line[i] for i in magi_indices]
                else:
                    gene_to_sample = None
            except IndexError:
                print "Index out of range in file " + str(config.get('options','file'))
                print "Index list:"
                print str(required_indices)
                print str(magi_indices)
                print "Line: "
                print line
                exit(1)

            # take only first three segments of name if sample is from TCGA
            if sample[:4] == 'TCGA':
                sample = '-'.join((sample.split('-'))[:3])

            # If a whitelist is provided, skip any genes/samples not in the list
            if gene_whitelist_file and gene not in gene_whitelist:
                continue
            if sample_whitelist_file and sample not in sample_whitelist:
                continue

            # Javascript can't have "." in gene names
            gene = gene.replace(".", "-")

            if (mutation_status.lower() not in exclude_status 
                    and variant_class_type.lower() not in exclude_mutations
                    and valid_status.lower() not in exclude_validation):

                sample_to_gene[sample].add(gene)

                stats['samples'][sample] += 1
                stats['genes'][gene] += 1
                stats['processed_mutations'] += 1
                stats['mutation_types'][variant_class_type] += 1

                # If magi_indices is None, this means that the data required to produce
                # the MAGI output could not be found in the MAF file. This is either because
                # the column names were not recognized, or the data was simply not present.

                # *************************************************************************** #

                if magi_indices == None:
                    continue



                if '.' in transcript_id:
                    transcript_id = transcript_id.split('.')[0]

                # Search all transcript databases for transcript id and
                # return corresponding length if found. If none is found,
                # add to missing transcripts set and ignore this mutation
                length = None
                for _, database in transcript_dict.items():
                    if transcript_id in database:
                        length = database[transcript_id]
                if not length:
                    stats['missing_transcripts'].add(transcript_id)
                    length = 'transcript_not_found'

                if not transcript_id:
                    transcript_id = 'no_transcript_id'

                try:
                    original_amino_acid, new_amino_acid, amino_acid_location = get_amino_acid_change(
                                                aa_change, variant_type, variant_class_type, codon)
                except ValueError:
                    stats['unknown_mutations'].add(' '.join([aa_change, variant_type, variant_class_type]))
                    continue

                if original_amino_acid and new_amino_acid and amino_acid_location:


                    gene_to_sample[gene][sample].append({'transcript':transcript_id,
                         'length':length, 'locus':amino_acid_location, 'mutation_type': variant_class_type,
                         'o_amino_acid':original_amino_acid, 'n_amino_acid':new_amino_acid})



    return gene_to_sample, sample_to_gene, stats

def parse_whitelist(whitelist_file):

    white_set = set()
    with open(whitelist_file, "r") as white_file:
        for line in white_file:
            white_set.add(line.rstrip())
    return white_set


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

    # header_name_list = [gene_name, sample, variant_class_type, variant_type, valid_status, 
                        # mut_stat, transcript_id, codon, aachange]
    required_headers = [gene_name, sample, variant_class_type, variant_type, valid_status, mut_stat]
    magi_required_headers = [transcript_id, codon, aachange]

    required_indices = []
    magi_indices = []

    # Find indices needed to process HotNet2 and CoMEt outputs. Without these, it
    # is not possible to process the MAF at all.
    for header_options in required_headers:
        existing_headers = set(header_options) & set(indice_dict.keys())
        if len(existing_headers) == 0:
            print "Header %s not found, no output possible. Aborting." % str(header_options)
            raise IndexError()
        else:
            required_indices.append(indice_dict[existing_headers.pop()])
    
    # Find indices needed to process MAGI output. 
    for header_options in magi_required_headers:
        if magi_indices == None:
            continue
        existing_headers = set(header_options) & set(indice_dict.keys())
        if len(existing_headers) == 0:
            print "Names %s not found, MAGI ouput not possible, attempting to output other formats if selected" % str(header_options)
            magi_indices = None
        else:
            magi_indices.append(indice_dict[existing_headers.pop()])
            
    return required_indices, magi_indices

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
            elif variant_class_type in ('5\'Flank', 'Read-through', "Nonstop_Mutation"):
                return None, None, None
            else:
                sys.stderr.write("New mutation effect can't be parsed: %s\n" % variant_class_type)
                raise ValueError

        elif variant_type in ("DEL", "INS"):
            aa_original, aa_new, aa_location = ins_del_mutation(aa_change, codon)
        else:
            print aa_change + ' ' + variant_type + ' ' + variant_class_type + ' ' + codon
            sys.stderr.write("New mutation type can't be parsed: %s\n" % variant_type)
            raise ValueError

    return aa_original, aa_new, aa_location

def write_magi(gene_to_sample, config):
    out_dir = config.get('options','output_dir')
    out_name = config.get('options', 'prefix')+"_magi.tsv"

    with open(os.path.join(out_dir, out_name), "w") as outfile:
        if gene_to_sample == None:
            outfile.write("THERE WAS AN ERROR PROCESSING THIS MAF FILE.\n")
            outfile.write("Data required for MAGI processing was not found.")
            return
        header = "#Gene\tSample\tTranscript\tTranscript_Length\tLocus\t"\
                     "Mutation_Type\tOriginal_Amino_Acid\tNew_Amino_Acid\n"
        outfile.write(header)
        for gene, samplelist in sorted(gene_to_sample.items()):
            for sample in samplelist:
                for mut in samplelist[sample]:
                    outfile.write('\t'.join([gene,sample, mut['transcript'], str(mut['length']),
                        str(mut['locus']), mut['mutation_type'], mut['o_amino_acid'], mut['n_amino_acid']])+'\n')

def write_other(sample_to_gene, config, out_type):
    """
    Outputs a mutation file compatible with both HotNet2 and CoMEt.
    """
    out_dir = config.get('options','output_dir')
    out_name = config.get('options', 'prefix')+"_"+out_type+".tsv"

    with open(os.path.join(out_dir,out_name), "w") as outfile:
        for sample in sorted(sample_to_gene):
            if len(sample_to_gene[sample]) > 0:
                outfile.write(sample + '\t' + '\t'.join(sorted(list(sample_to_gene[sample])))+'\n')

def output_stats(stats, config):
# stats = {'total_mutations':0, 'processed_mutations:':0, 'missing_transcripts':set(), 
#      'samples':set(), 'genes':set(), 'mutation_types':defaultdict(lambda: 0),
#      'unknown_mutations':set()}

    output_list = []
    output_list.append("*** Summary statistics ***")
    output_list.append("*  Total mutations in file:\t" + str(stats['total_mutations']))
    output_list.append("*  Mutations successfully processed:\t" + str(stats['processed_mutations']))
    output_list.append("*  Success ratio:\t" + str(float(stats['processed_mutations']) / stats['total_mutations']))
    output_list.append("*  Unique samples:\t" + str(len(stats['samples'])))
    output_list.append("*  Unique genes:\t" + str(len(stats['genes'])))
    output_list.append("*  Mutation types and totals:")
    for mut_type, quantity in stats['mutation_types'].items():
        output_list.append(mut_type + '\t' + str(quantity))
    if len(stats['missing_transcripts']) > 0 :
        output_list.append("*  Missing transcripts:")
        for transcript in stats['missing_transcripts']:
            output_list.append("**\t" + transcript)
    if len(stats['unknown_mutations']) > 0 :
        output_list.append("*  Unknown mutations:")
        for mutation in stats['unknown_mutations']:
            output_list.append("**\t" + mutation)

    out_dir = config.get('options','output_dir')
    out_name = config.get('options', 'prefix')+'_summary.txt'


    with open(os.path.join(out_dir, out_name), 'w') as write_file:
        for line in output_list:
            write_file.write(line+'\n')

def visualize_data(stats, gene_to_sample, config):

    # ADD AMINO ACID CHANGE
     #####################################################
    #           Number of each kind of mutation           # 
     #####################################################

    plot_items = []
    plot_value = []

    for mutation, quantity in stats['mutation_types'].items():
        plot_items.append(mutation)
        plot_value.append(quantity)

    y_pos = np.arange(len(plot_items))
    plt.barh(y_pos, np.array(plot_value), align='center', alpha=0.8)
    plt.yticks(y_pos, plot_items)
    plt.xlabel('Total mutations')

    plt.tight_layout()
    plt.savefig(config.get('options','prefix')+'_mut_count.svg')
    plt.close()

     #####################################################
    #             Top 40 Amino Acid Mutations             # 
     #####################################################

    change_count = defaultdict(int)

    for gene in gene_to_sample:
        for sample in gene_to_sample[gene]:
            for mut in gene_to_sample[gene][sample]:
                key = mut['o_amino_acid'] + ' to ' + mut['n_amino_acid']
                change_count[key] += 1


    top_names = sorted(change_count, key=lambda k: change_count[k], reverse=True)[:40]

    top_values = []
    for mut in top_names:
        top_values.append(change_count[mut])

    plt.figure(figsize=(12,6))

    x_pos = np.arange(len(top_names))
    plt.xticks(x_pos, top_names,rotation=45)
    plt.xlabel('Mutation')
    plt.ylabel('Count')
    plt.title('Top 40 Amino Acid Changes')

    plt.bar(x_pos, np.array(top_values), align='center', alpha=0.8)
    plt.tight_layout()
    plt.savefig(config.get('options','prefix')+'_type_count.svg')
    plt.close()

def get_parser():
    '''
    Parse arguments.
    '''

    parser = argparse.ArgumentParser(description='Parse MAF files')

    parser.add_argument('-f', '--file', help='Path to MAF file to be processed')
    parser.add_argument('-d', '--database', help='Path to transcripts database file')
    parser.add_argument('-s', '--statistics', action='store_true')
    parser.add_argument('-v', '--visualization', action='store_true')
    parser.add_argument('-o', '--output_dir',help="Folder to output processed MAF data")
    parser.add_argument('-t', '--type', help="Output format options: magi,hotnet2,comet", type=str.lower,
                         choices=['magi', 'hotnet2', 'comet'],nargs='+')
    parser.add_argument('-p', '--prefix', help='String to prefix output files.')


    return parser

def get_config(args):

    config = ConfigParser.SafeConfigParser()
    found = config.read('maf.cfg')
    if not found:
        raise IOError("Error: Config file not found!")

    # Integrate command line arguments into configuration
    for arg in vars(args):
        attribute = getattr(args, arg)
        if attribute:
            config.set('options', arg, str(attribute))


    if not os.path.isfile(config.get('options', 'file')):
        raise IOError('Error: MAF file not found. Please check location in configuration file')

    if not os.path.isfile(config.get('options', 'database')):
        raise IOError('Error: Transcript database file not found. Please check location in configuration file')

    if not config.has_option('options','prefix') or config.get('options','prefix') == '':
        config.set('options','prefix',os.path.basename(config.get('options', 'file')).split('.')[0][:10])

    return config

def run(config):
    '''
    Main function.
    '''

    sample_whitelist = config.get('options', 'sample')
    gene_whitelist = config.get('options', 'gene')
    maf_file = config.get('options', 'file')

    with open(config.get('options', 'database')) as t_file:
        transcript_dict = json.load(t_file)

    gene_to_sample, sample_to_gene, stats = process_maf_file(maf_file, transcript_dict, sample_whitelist, gene_whitelist, config)

    if not os.path.exists(config.get('options', 'output_dir')):
        os.makedirs(config.get('options', 'output_dir'))
    if 'magi' in config.get('options', 'type'):
        write_magi(gene_to_sample, config)
    if 'comet' in config.get('options', 'type'):
        write_other(sample_to_gene, config, 'comet')
    if 'hotnet2' in config.get('options', 'type'):
        write_other(sample_to_gene, config, 'hotnet2')

    if config.getboolean('options','statistics'):
        output_stats(stats, config)

    if config.getboolean('options','visualization'):
        visualize_data(stats, gene_to_sample, config)


################################################################################
# Parse the different mutation formats to get the amino acid change

def snp_mutation(aa_change):
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




if __name__ == '__main__':
    run(get_config(get_parser().parse_args(sys.argv[1:])))