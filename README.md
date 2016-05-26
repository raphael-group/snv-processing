# SNV Processing

Purpose 
---
Proces SNV data (.maf files) for use in CoMEt, MAGI, and HotNet2.

Requires
---
* Python 2.7 (tested with 2.7.9)
* Matplotlib (tested with 1.4.2)
* Numpy (tested with 1.8.2)
* Transcript database file

Quick Start
---

`python MAFprocessing.py -f <path to MAF file> -d <path to transcripts database> `

Output
---

#### HotNet2
A tab separated file, each row is a sample and associated gene mutations. The first column is the sample name, each entry after is a gene. The name is {prefix}_hotnet2.tsv. The prefix is a optional argument. If no prefix is supplied, the default is 'output'.

#### CoMEt
The same as HotNet2, named {prefix}_comet.tsv.

#### MAGI
A eight column, tab separated file with a header. Each row (in order) consists of a gene, a sample, the transcript ID, the transcript length, the locus of the mutation, the mutation type, the original amino acid, and the new amino acid.


Usage
---
All options can be set by command line or via a configuration file. The 
command line will override configuration file settings, and the names and usage for
command line and configuration file options are the same. 

### Required ###

The following options must be specified by the user either at the command line or in the config file, and cannot be left to defaults:

Config/long argument | Short argument   | Input type | Description 
:-------------------------------------| :----- |:----- |:-----
database  | -d | Path to transcript database |Path to the transcript database, a json file with a key/value pair of { 'database' : {'transcript' : length (integer)} }. See example in this repository (transcript-lenghts.json).
file     | -f | Path to MAF file | Path to the maf file to be processed.

