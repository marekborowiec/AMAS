#! /usr/bin/env python3

#   Program to calculate various statistics on a multiple sequence alignment
#   and allow efficient manipulation of phylogenomic data sets

#   Copyright (C) 2015 Marek Borowiec

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
  
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
This stand-alone program allows manipulations of multiple sequence
alignments. It supports sequential FASTA, PHYLIP, NEXUS, and interleaved PHYLIP 
and NEXUS formats for DNA and aino acid sequences. It can print summary statistics,
convert among formats, and concatenate alignments.

Current statistics include the number of taxa, alignment length, total number
of matrix cells, overall number of undetermined characters, percent of missing 
data, AT and GC contents (for DNA alignments), number and proportion of 
variable sites, number and proportion of parsimony informative sites,
and counts of all characters present in the relevant (nucleotide or amino acid) alphabet.
"""


import argparse, multiprocessing as mp, re, sys
from random import sample
from os import path, remove
from collections import defaultdict, Counter
from itertools import compress

def proportion(x):
    # needed to prevent input of invalid floats in trim mode
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x

class ParsedArgs:

    def __init__(self):
        parser = argparse.ArgumentParser(
            usage='''AMAS <command> [<args>]

The AMAS commands are:
  concat      Concatenate input alignments
  convert     Convert to other file format
  replicate   Create replicate data sets for phylogenetic jackknife
  split       Split alignment according to a partitions file
  summary     Write alignment summary
  remove      Remove taxa from alignment
  translate   Translate DNA alignment into protein alignment
  trim        Remove columns from alignment

Use AMAS <command> -h for help with arguments of the command of interest
'''
        )

        parser.add_argument(
            "command", 
            help="Subcommand to run"
        )

        # parse_args defaults to [1:] for args, but you need to
        # exclude the rest of the args too, or validation will fail
        self.args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, self.args.command):
            print('Unrecognized command')
            parser.print_help()
            exit(1)
        # use dispatch pattern to invoke method with same name
        getattr(self, self.args.command)()

    def add_common_args(self, parser):
        # define required arguments for every command
        requiredNamed = parser.add_argument_group('required arguments')
        parser.add_argument(
            "-e",
            "--check-align",
            dest = "check_align",
            action = "store_true",
            default = False,
            help = "Check if input sequences are aligned. Default: no check"
        )
        parser.add_argument(
        # parallelization is used for file parsing and calculating summary stats
            "-c",
            "--cores",
            dest = "cores",
            default = 1,
            help = "Number of cores used. Default: 1"
        )
        
        requiredNamed.add_argument(
            "-i",
            "--in-files",
            nargs = "+",
            type = str,
            dest = "in_files",
            required = True,
            help = """Alignment files to be taken as input.
            You can specify multiple files using wildcards (e.g. --in-files *fasta)"""
        )
        requiredNamed.add_argument(
            "-f",
            "--in-format",
            dest = "in_format",
            required = True,
            choices = ["fasta", "phylip", "nexus", "phylip-int", "nexus-int"],
            help = "The format of input alignment"
        )
        requiredNamed.add_argument(
            "-d",
            "--data-type",
            dest = "data_type",
            required = True,
            choices = ["aa", "dna"],
            help = "Type of data"
        )

    def trim(self):
        # convert command

        parser = argparse.ArgumentParser(
            description="Trim alignment by occupancy. Optionally removes sites that are not parsimony informative. \n CAUTION: when running on amino acids stop codons marked with * will be treated as missing data!",
        )
        parser.add_argument(
            "-u",
            "--out-format",
            dest = "out_format",
            choices = ["fasta", "phylip", "nexus", "phylip-int", "nexus-int"],
            default = "fasta",
            help = "File format for the output alignment. Default: fasta"
        )
        parser.add_argument(
            "-o",
            "--trim-out",
            dest = "trim_out",
            help = "File name for the trimmed alignment when providing a single file as input."
        )
        parser.add_argument(
            "-t",
            "--trim-fraction",
            type = proportion, 
            dest = "trim_fraction",
            default = 0.6,
            help = "Columns in the alignments with occupancy lower than this value will be removed. Default: 0.6"
        )
        parser.add_argument(
            "-p",
            "--retain-only-parsimony-sites",
            dest = "parsimony_check",
            action = "store_true",
            default = False,
            help = "Only write parsimony informative columns in trimmed alignment Default: write all columns"
        )
        # add shared arguments
        self.add_common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        return args

    def summary(self):
       # summary command
        parser = argparse.ArgumentParser(
            description="Write alignment summary",
        )
        parser.add_argument(
            "-o",
            "--summary-out",
            dest = "summary_out",
            default = "summary.txt",
            help = "File name for the alignment summary. Default: 'summary.txt'"
        )
        parser.add_argument(
            "-s",
            "--by-taxon",
            dest = "by_taxon_summary",
            action = "store_true",
            default = False,
            help = "In addition to alignment summary, write by sequence/taxon summaries. Default: Don't write"
        )
        # add shared arguments
        self.add_common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        return args

    def concat(self):
        # concat command
        parser = argparse.ArgumentParser(
            description="Concatenate input alignments"
        )
        parser.add_argument(
            "-p",
            "--concat-part",
            dest = "concat_part",
            default = "partitions.txt",
            help = "File name for the concatenated alignment partitions. Default: 'partitions.txt'"
        ) 
        parser.add_argument(
            "-t",
            "--concat-out",
            dest = "concat_out",
            default = "concatenated.out",
            help = "File name for the concatenated alignment. Default: 'concatenated.out'"
        )
        parser.add_argument(
            "-u",
            "--out-format",
            dest = "out_format",
            choices = ["fasta", "phylip", "nexus", "phylip-int", "nexus-int"],
            default = "fasta",
            help = "File format for the output alignment. Default: fasta"
        )
        parser.add_argument(
            "-y",
            "--part-format",
            dest = "part_format",
            choices = ["nexus", "raxml", "unspecified"],
            default = "unspecified",
            help = "Format of the partitions file. Default: 'unspecified'"
        ) 
        # add shared arguments
        self.add_common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        return args

    def convert(self):
        # convert command
        parser = argparse.ArgumentParser(
            description="Convert to other file format",
        )
        parser.add_argument(
            "-u",
            "--out-format",
            dest = "out_format",
            choices = ["fasta", "phylip", "nexus", "phylip-int", "nexus-int"],
            default = "fasta",
            help = "File format for the output alignment. Default: fasta"
        )
        # add shared arguments
        self.add_common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        return args

    def replicate(self):
        # replicate command
        parser = argparse.ArgumentParser(
            description="Create replicate datasets for phylogenetic jackknife",
        )
        parser.add_argument(
            "-r",
            "--rep-aln",
            nargs = 2,
            type = int,
            dest = "replicate_args",
            help = "Create replicate data sets for phylogenetic jackknife [replicates, no alignments for each replicate]",
            required = True
        ) 
        parser.add_argument(
            "-u",
            "--out-format",
            dest = "out_format",
            choices = ["fasta", "phylip", "nexus", "phylip-int", "nexus-int"],
            default = "fasta",
            help = "File format for the output alignment. Default: fasta"
        ) 
        # add shared arguments
        self.add_common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        return args

    def split(self):
        # split command
        parser = argparse.ArgumentParser(
            description="Split alignment according to a partitions file",
        )
        parser.add_argument(
            "-l",
            "--split-by",
            dest = "split_by",
            help = "File name for partitions to be used for alignment splitting.",
            required = True
        )
        parser.add_argument(
            "-j",
            "--remove-empty",
            dest = "remove_empty",
            action = "store_true",
            default = False,
            help = "Remove taxa with sequences composed of only undetermined characters? Default: Don't remove"
        )
        parser.add_argument(
            "-u",
            "--out-format",
            dest = "out_format",
            choices = ["fasta", "phylip", "nexus", "phylip-int", "nexus-int"],
            default = "fasta",
            help = "File format for the output alignment. Default: fasta"
        )
        # add shared arguments
        self.add_common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        return args

    def translate(self):
        # translate command
        parser = argparse.ArgumentParser(
            description="Translate a protein-coding DNA alignment into amino acids",
        )
        parser.add_argument(
            "-b",
            "--code",
            type = int,
            dest = "genetic_code",
            choices = [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 16, 21, 22, 23, 24, 25, 26],
            default = 1,
            help = "NCBI genetic code to use: 1. The Standard Code, 2. The Vertebrate Mitochondrial Code, \
3. The Yeast Mitochondrial Code, 4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code, \
5. The Invertebrate Mitochondrial Code, 6. The Ciliate, Dasycladacean and Hexamita Nuclear Code, \
9. The Echinoderm and Flatworm Mitochondrial Code, 10. The Euplotid Nuclear Code, 11. The Bacterial, Archaeal and Plant Plastid Code, \
12. The Alternative Yeast Nuclear Code, 13. The Ascidian Mitochondrial Code, 14. The Alternative Flatworm Mitochondrial Code, \
16. Chlorophycean Mitochondrial Code, 21. Trematode Mitochondrial Code, 22. Scenedesmus obliquus Mitochondrial Code, \
23. Thraustochytrium Mitochondrial Code, 24. Pterobranchia Mitochondrial Code, 25. Candidate Division SR1 and Gracilibacteria Code, \
26. Pachysolen tannophilus Nuclear Code. Default: 1."
        )
        parser.add_argument(
            "-k",
            "--reading-frame",
            type = int,
            dest = "reading_frame",
            choices = [1, 2, 3],
            default = 1,
            help = "Number specifying reading frame; i.e. '2' means codons start at the second character of the alignment. Default: 1",
        )
        parser.add_argument(
            "-u",
            "--out-format",
            dest = "out_format",
            choices = ["fasta", "phylip", "nexus", "phylip-int", "nexus-int"],
            default = "fasta",
            help = "File format for the output alignment. Default: fasta"
        )
        # add shared arguments
        self.add_common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        return args

    def remove(self):
        # remove taxa command
        parser = argparse.ArgumentParser(
            description="Remove taxa from alignment",
        )
        parser.add_argument(
            "-x",
            "--taxa-to-remove",
            nargs = "+",
            type = str,
            dest = "taxa_to_remove",
            help = "Taxon/sequence names to be removed.",
            required = True
        )
        parser.add_argument(
            "-u",
            "--out-format",
            dest = "out_format",
            choices = ["fasta", "phylip", "nexus", "phylip-int", "nexus-int"],
            default = "fasta",
            help = "File format for the output alignment. Default: fasta"
        )
        parser.add_argument(
            "-g",
            "--out-prefix",
            dest = "out_prefix",
            default = "reduced_",
            help = "File name prefix for the concatenated alignment. Default: 'reduced_'"
        )
        # add shared arguments
        self.add_common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        return args


    def get_args_dict(self):
    # store arguments in a dictionary
        command = self.args.__dict__
        arguments = getattr(self, self.args.command)().__dict__
        argument_dictionary = command.copy()
        argument_dictionary.update(arguments)
        
        return argument_dictionary
    

class FileHandler:
    """Define file handle that closes when out of scope"""

    def __init__(self, file_name):
        self.file_name = file_name

    def __enter__(self):
        try:
            self.in_file = open(self.file_name, "r")
        except FileNotFoundError:
            print("ERROR: File '" + self.file_name + "' not found.")
            sys.exit()
        return self.in_file

    def __exit__(self, *args):
        self.in_file.close()
    
    def get_file_name(self):
        return self.file_name
        
class FileParser:
    """Parse file contents and return sequences and sequence names"""
    
    def __init__(self, in_file):
        self.in_file = in_file
        with FileHandler(in_file) as handle:
            self.in_file_lines = handle.read().rstrip("\r\n")    

    def fasta_parse(self):
        # use regex to parse names and sequences in sequential fasta files
        matches = re.finditer(
            r"^>(.*[^$])([^>]*)",
            self.in_file_lines, re.MULTILINE
        )
        records = {}
 
        for match in matches:
            name_match = match.group(1).replace("\n","")
            seq_match = match.group(2).replace("\n","").upper()
            seq_match = self.translate_ambiguous(seq_match)
            records[name_match] = seq_match

        return records
    
    def phylip_parse(self):
        # use regex to parse names and sequences in sequential phylip files
        matches = re.finditer(
            r"^(\s+)?(\S+)\s+([A-Za-z*?.{}-]+)",
            self.in_file_lines, re.MULTILINE
        )

        records = {}

        for match in matches:
            name_match = match.group(2).replace("\n","")
            seq_match = match.group(3).replace("\n","").upper()
            seq_match = self.translate_ambiguous(seq_match)
            records[name_match] = seq_match

        return records    

    def phylip_interleaved_parse(self):
    # use regex to parse names and sequences in interleaved phylip files
        tax_chars_matches = re.finditer(
            r"^(\s+)?([0-9]+)[ \t]+([0-9]+)",
            self.in_file_lines, re.MULTILINE
        )
        name_matches = re.finditer(
            r"^(\s+)?(\S+)[ \t]+[A-Za-z*?.{}-]+",
            self.in_file_lines, re.MULTILINE
        )
        seq_matches = re.finditer(
            r"(^(\s+)?\S+[ \t]+|^)([A-Za-z*?.{}-]+)$",
            self.in_file_lines, re.MULTILINE
        )
        # get number of taxa and chars
        for match in tax_chars_matches:
            tax_match = match.group(2)
            chars_match = match.group(3)

        # initiate lists for taxa names and sequence strings on separate lines
        taxa = []
        sequences = []
        # initiate a dictionary for the name:sequence records
        records = {}
        # initiate a counter to keep track of sequences strung together
        # from separate lines
        counter = 0
        
        for match in name_matches:
            name_match = match.group(2).replace("\n","")
            taxa.append(name_match)

        for match in seq_matches:
            seq_match = match.group(3).replace("\n","").upper()
            seq_match = self.translate_ambiguous(seq_match)
            sequences.append(seq_match)
        # try parsing PHYLUCE-style interleaved phylip
        if len(taxa) != int(tax_match):
            taxa = []
            sequences = []
            matches = re.finditer(
                r"(^(\s+)?(\S+)( ){2,}|^\s+)([ A-Za-z*?.{}-]+)",
                self.in_file_lines, re.MULTILINE
            )
            
            for match in matches:
                try:
                    name_match = match.group(3).replace("\n","")
                    taxa.append(name_match)
                except AttributeError:
                    pass
                seq_match = match.group(5).replace("\n","").upper()
                seq_match = "".join(seq_match.split())
                seq_match = self.translate_ambiguous(seq_match)
                sequences.append(seq_match)

        for taxon_no in range(len(taxa)):
            sequence = ""
            for index in range(counter,len(sequences),len(taxa)):
                sequence += sequences[index] 
           
            records[taxa[taxon_no]] = sequence
            counter += 1 
            
        return records
        
    def nexus_parse(self):
        # use regex to parse names and sequences in sequential nexus files
        # find the matrix block
        matches = re.finditer(
            r"(\s+)?(MATRIX\n|matrix\n|MATRIX\r\n|matrix\r\n)(.*?;)",
            self.in_file_lines, re.DOTALL
        )
        
        records = {}
        # get names and sequences from the matrix block

        for match in matches:
            matrix_match = match.group(3)
            seq_matches = re.finditer(
                 r"^(\s+)?[']?(\S+\s\S+|\S+)[']?\s+([A-Za-z*?.{}-]+)($|\s+\[[0-9]+\]$)",
                 matrix_match, re.MULTILINE
             )

            for match in seq_matches:
                name_match = match.group(2).replace("\n","")
                seq_match = match.group(3).replace("\n","").upper()
                seq_match = self.translate_ambiguous(seq_match)
                records[name_match] = seq_match

        return records
        
    def nexus_interleaved_parse(self):
    # use regex to parse names and sequences in sequential nexus files
    # find the matrix block
        matches = re.finditer(
            r"(\s+)?(MATRIX\n|matrix\n|MATRIX\r\n|matrix\r\n)(.*?;)",
            self.in_file_lines, re.DOTALL
        )
        # initiate lists for taxa names and sequence strings on separate lines
        taxa = []
        sequences = []
        # initiate a dictionary for the name:sequence records
        records = {}

        for match in matches:
            matrix_match = match.group(3)
            # get names and sequences from the matrix block
            seq_matches = re.finditer(
                r"^(\s+)?[']?(\S+\s\S+|\S+)[']?\s+([A-Za-z*?.{}-]+)($|\s+\[[0-9]+\]$)",
                matrix_match, re.MULTILINE
            )

            for match in seq_matches:
                name_match = match.group(2)
                if name_match not in taxa:
                    taxa.append(name_match)
                seq_match = match.group(3)
                
                sequences.append(seq_match)

        # initiate a counter to keep track of sequences strung together
        # from separate lines
        counter = 0

        for taxon_no in range(len(taxa)):

            full_length_sequence = "".join([sequences[index] for index in range(counter,len(sequences),len(taxa))])    
            records[taxa[taxon_no]] = self.translate_ambiguous(full_length_sequence).replace("\n","").upper()
            counter += 1 

        return records

    def translate_ambiguous(self, seq):
        # translate ambiguous characters from curly bracket format
        # to single letter format 
        # also remove spaces from sequences
        seq = seq.replace("{GT}","K")
        seq = seq.replace("{AC}","M")
        seq = seq.replace("{AG}","R")
        seq = seq.replace("{CT}","Y")
        seq = seq.replace("{CG}","S")
        seq = seq.replace("{AT}","W")
        seq = seq.replace("{CGT}","B")
        seq = seq.replace("{ACG}","V")
        seq = seq.replace("{ACT}","H")
        seq = seq.replace("{AGT}","D")
        seq = seq.replace("{GATC}","N")
        seq = seq.replace(" ","")

        return seq

    def partitions_parse(self):
        # parse partitions file using regex
        matches = re.finditer(r"^(\s+)?([^ =]+)[ =]+([\\0-9, -]+)", self.in_file_lines, re.MULTILINE)
        
        # initiate list to store dictionaries with lists
        # of slice positions as values
        partitions = []
        add_to_partitions = partitions.append
        for match in matches:
            # initiate dictionary of partition name as key
            dict_of_dicts = {}
            # and list of dictionaries with slice positions
            list_of_dicts = []
            add_to_list_of_dicts = list_of_dicts.append
            # get parition name and numbers from parsed partition strings
            partition_name = match.group(2)
            numbers = match.group(3)
            # find all numbers that will be used to parse positions
            positions = re.findall(r"([^ ,]+)", numbers)
            
            for position in positions:
                # create dictionary for slicing input sequence
                # conditioning on whether positions are represented
                # by range, range with stride, or single number
                pos_dict = {}
        
                if "-" in position:
                    m = re.search(r"([0-9]+)-([0-9]+)", position)
                    pos_dict["start"] = int(m.group(1)) - 1
                    pos_dict["stop"] = int(m.group(2))
                else:
                    pos_dict["start"] = int(position) - 1
                    pos_dict["stop"] = int(position)
        
                if "\\" in position:
                    pos_dict["stride"] = 3
                elif "\\" not in position:
                    pos_dict["stride"] = 1
        
                add_to_list_of_dicts(pos_dict)
                
            dict_of_dicts[partition_name] = list_of_dicts 
            add_to_partitions(dict_of_dicts)

        return partitions

 
class Alignment:
    """Gets in parsed sequences as input and summarizes their stats"""
    
    def __init__(self, in_file, in_format, data_type):    
        # initialize alignment class with parsed records and alignment name as arguments,
        # create empty lists for list of sequences, sites without
        # ambiguous or missing characters, and initialize variable for the number
        # of parsimony informative sites
        self.in_file = in_file
        self.in_format = in_format
        self.data_type = data_type

        self.parsed_aln = self.get_parsed_aln()
        
    def __str__(self):
        return self.get_name

    def get_aln_input(self):
        # open and parse input file
        aln_input = FileParser(self.in_file)
        return aln_input

    def get_parsed_aln(self):
        # parse according to the given format
        aln_input = self.get_aln_input()
        if self.in_format == "fasta":
            parsed_aln = aln_input.fasta_parse()
        elif self.in_format == "phylip":
            parsed_aln = aln_input.phylip_parse()
        elif self.in_format == "phylip-int":
            parsed_aln = aln_input.phylip_interleaved_parse()
        elif self.in_format == "nexus":
            parsed_aln = aln_input.nexus_parse()
        elif self.in_format == "nexus-int":
            parsed_aln = aln_input.nexus_interleaved_parse()

        return parsed_aln
        
    def summarize_alignment(self):
        # call methods to create sequences list, matrix, sites without ambiguous or
        # missing characters; get and summarize alignment statistics
        summary = []
        self.length = str(self.get_alignment_length())
        self.matrix = self.matrix_creator()
        self.no_missing_ambiguous = self.get_sites_no_missing_ambiguous()
        self.variable_sites = self.get_variable()
        self.prop_variable = self.get_prop_variable()
        self.parsimony_informative = self.get_parsimony_informative()
        self.prop_parsimony = self.get_prop_parsimony()
        self.missing_records = self.get_missing_from_parsed()
        name = str(self.get_name())
        taxa_no = str(self.get_taxa_no())
        cells = str(self.get_matrix_cells())
        missing = str(self.get_missing())
        missing_percent = str(self.get_missing_percent())
        self.check_data_type()
        summary = [name, taxa_no, self.length, cells, missing, missing_percent, \
         str(self.variable_sites), str(self.prop_variable), str(self.parsimony_informative), str(self.prop_parsimony)]
        return summary

    def summarize_alignment_by_taxa(self):
        # get summary for all taxa/sequences in alignment
        per_taxon_summary = []
        taxa_no = self.get_taxa_no()
        self.missing_records = self.get_missing_from_parsed()
        self.length = self.get_alignment_length()
        lengths = (self.length for i in range(taxa_no))
        name = self.get_name()
        names = (name for i in range(taxa_no))
        taxa_names = (taxon.replace(" ","_").replace(".","_").replace("'","") \
         for taxon, missing_count, missing_percent in self.missing_records)
        missing = (missing_count for taxon, missing_count, missing_percent in self.missing_records)
        missing_percent = (missing_percent for taxon, missing_count, missing_percent in self.missing_records)
        self.check_data_type()
        per_taxon_summary = (names, taxa_names, lengths, missing, missing_percent)
        zipped = list(zip(*per_taxon_summary))
        return zipped

    def get_char_summary(self):
        # get summary of frequencies for all characters
        characters = []
        counts = []
        add_to_chars = characters.append
        add_to_counts = counts.append
        char_count_dicts = self.get_counts()
        for char in self.alphabet:
            add_to_chars(char)
            if char in char_count_dicts.keys():
                add_to_counts(str(char_count_dicts[char]))
            else:
                add_to_counts("0")
        return characters, counts

    def get_taxon_char_summary(self):
        # get summary of frequencies for all characters
        records = (self.append_count(char_dict) for taxon, char_dict in self.get_counts_from_parsed())
        return records

    def append_count(self, char_dict):
        count_list = []
        for char in self.alphabet:
            if char in char_dict.keys():
                count_list.append(char_dict[char])
            else:
                count_list.append(0)
        return count_list
     
    def matrix_creator(self):
        # decompose character matrix into a two-dimensional list
        matrix = [list(sequence) for sequence in self.parsed_aln.values()]
        return matrix

    def get_column(self, i):
        # get site from the character matrix
        return [row[i] for row in self.matrix]
        
    def all_same(self, site):
        # check if all elements of a site are the same
         return not site or site.count(site[0]) == len(site)

    def get_sites_no_missing_ambiguous(self):
        # get each site without missing or ambiguous characters
         no_missing_ambiguous_sites = [self.get_site_no_missing_ambiguous(column) for column in range(self.get_alignment_length())]
         return no_missing_ambiguous_sites

    def get_site_no_missing_ambiguous(self, column):
        site = self.get_column(column)
        return [char for char in site if char not in self.missing_ambiguous_chars]

    def replace_missing(self, column):
        return ["-" if x in self.missing_chars else x for x in self.get_column(column)]

    def get_trim_selection(self, trim_fraction, parsimony_check):
        # this checks each column of alignment for minimum occupancy
        self.matrix = self.matrix_creator()
        trim_vector = []
        for column in range(self.get_alignment_length()):
            site = self.replace_missing(column)
            occ = (len(site) - site.count("-")) / len(site)
            if parsimony_check:
                unique_chars = set(site)
                try:
                    unique_chars.remove("-")
                except KeyError: 
                    pass # this occurs if we have no missing data
                pattern = [base for base in unique_chars if site.count(base) >= 2]
                trim_vector.append(len(pattern) >= 2 and occ >= trim_fraction)
            else:
                trim_vector.append(occ >= trim_fraction)
        return trim_vector
   
    def get_variable(self):
        # if all elements of a site without missing or ambiguous characters 
        # are not the same, consider it variable
        variable = len([site for site in self.no_missing_ambiguous \
         if not self.all_same(site)])
        return variable
    
    def get_parsimony_informative(self):
        # if the count for a unique character in a site is at least two, 
        # and there are at least two such characters in a site without missing
        # or ambiguous characters, consider it parsimony informative
        parsimony_informative = 0
        for site in self.no_missing_ambiguous:
            unique_chars = set(site)
            pattern = [base for base in unique_chars if site.count(base) >= 2]
            no_patterns = len(pattern)
            
            if no_patterns >= 2:
                parsimony_informative += 1
        return parsimony_informative
    
    def get_prop_variable(self):
        # get proportion of variable sites to all sites
        prop_variable = self.variable_sites / int(self.length)
        return round(prop_variable, 3)
        
    def get_prop_parsimony(self):
        # get proportion of parsimony informative sites to all sites
        prop_parsimony = self.parsimony_informative / int(self.length)
        return round(prop_parsimony, 3)

    def get_name(self):
        # get input file name
        in_filename = path.basename(self.in_file)
        return in_filename
        
    def get_taxa_no(self):
        # get number of taxa
        return len(self.parsed_aln.values())
    
    def get_alignment_length(self):
        # get alignment length by just checking the first seq length
        # this assumes that all sequences are of equal length
        return len(next(iter(self.parsed_aln.values())))

    def get_matrix_cells(self):
    # count all matrix cells
        self.all_matrix_cells = len(self.parsed_aln.values()) \
         * int(self.length)
        return self.all_matrix_cells

    def get_missing(self):
        # count missing characters from the list of missing for all sequences
        self.missing = sum(count for taxon, count, percent in self.missing_records)
        return self.missing

    def get_missing_percent(self):
        # get missing percent
        missing_percent = round((self.missing / self.all_matrix_cells * 100), 3)
        return missing_percent
        
    def get_missing_from_parsed(self):
        # get missing count and percent from parsed alignment
        # return a list of tuples with taxon name, count, and percent missing
        self.missing_records = sorted([(taxon, self.get_missing_from_seq(seq), self.get_missing_percent_from_seq(seq)) \
         for taxon, seq in self.parsed_aln.items()])
        return self.missing_records
    
    def get_missing_from_seq(self, seq):
        # count missing characters for individual sequence
        missing_count = sum(seq.count(char) for char in self.missing_chars)
        return missing_count

    def get_missing_percent_from_seq(self, seq):
        # get missing percent from individual sequence
        missing_seq_percent = round((self.get_missing_from_seq(seq) / self.get_alignment_length() * 100), 3)
        return missing_seq_percent

    def get_counts(self):
        # get counts of each character in the used alphabet for all sequences
        counters = [Counter(chars) for taxon, chars in self.get_counts_from_parsed()]
        all_counts = sum(counters, Counter())
        counts_dict = dict(all_counts)
        return counts_dict

    def get_counts_from_parsed(self):
        # get counts of all characters from parsed alignment
        # return a list of tuples with taxon name and counts
        return sorted([(taxon, self.get_counts_from_seq(seq)) \
         for taxon, seq in self.parsed_aln.items()])

    def get_counts_from_seq(self, seq):
        # get all alphabet chars count for individual sequence
        char_counts = {char : seq.count(char) for char in self.alphabet}
        return char_counts

    def check_data_type(self):
        # check if the data type is correct; only one seq to save on computation
        seq = next(iter(self.parsed_aln.values()))
        self.check = any(char in self.non_alphabet for char in seq)
        if self.check is True:
            print("WARNING: found non-" + self.data_type + " characters. "\
             "Are you sure you specified the right data type?")


class AminoAcidAlignment(Alignment):
    """Alphabets specific to amino acid alignments"""

    alphabet = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R", \
     "S","T","V","W","Y","B","J","Z","X",".","*","-","?"]
    missing_ambiguous_chars = ["B","J","Z","X",".","*","-","?"]
    missing_chars = ["X",".","*","-","?"]
    non_alphabet = ["O"]

    def get_summary(self):
        # get alignment summary specific to amino acids
        data = self.summarize_alignment()
        new_data = data + list(self.get_char_summary()[1])
        return new_data

    def get_taxa_summary(self):
        # get per-taxon/sequence alignment summary specific to amino acids
        data = self.summarize_alignment_by_taxa()
        aa_summary = (data, self.get_taxon_char_summary())
        zipped_list = list(zip(*aa_summary))
        new_data = [list(data_tupl) + chars for data_tupl, chars in zipped_list]
        return new_data
           
class DNAAlignment(Alignment):
    """Alphabets specific to DNA alignments"""
    
    alphabet = ["A","C","G","T","K","M","R","Y","S","W","B","V","H","D","X", \
     "N", "O", "-","?"]
    missing_ambiguous_chars = ["K","M","R","Y","S","W","B","V","H","D","X", \
     "N", "O", "-","?"] 
    missing_chars = ["X","N","O","-","?"]
    non_alphabet = ["E", "F", "I", "L", "P", "Q", "J", "Z", ".", "*"]    

    def get_summary(self):
        # get alignment summarry specific to nucleotide
        data = self.summarize_alignment()
        new_data = data + self.get_atgc_content() \
         + list(self.get_char_summary()[1])
        return new_data
        
    def get_taxa_summary(self):
        # get per-taxon/sequence alignment summary specific to nucleotides
        data = self.summarize_alignment_by_taxa()
        dna_summary = (data, self.get_list_from_atgc(), self.get_taxon_char_summary())
        zipped_list = list(zip(*dna_summary))
        new_data = [list(data_tupl) + list(atgc) + chars for data_tupl, atgc, chars in zipped_list]
        return new_data

    def get_atgc_content(self):
        # get AC and GC contents for all sequences
        # AT content is the first element of AT, GC content tuple
        # returned by get_atgc_from_seq()
        atgc_records = self.get_atgc_from_parsed()
        at_content = round(sum(atgc[0] for taxon, atgc in atgc_records) \
         / self.get_taxa_no(), 3)
        gc_content = round(1 - float(at_content), 3)
        
        atgc_content = [str(at_content), str(gc_content)]
        return atgc_content

    def get_list_from_atgc(self):
        records = (atgc for taxon, atgc in self.get_atgc_from_parsed())
        return records

    def get_atgc_from_parsed(self):
        # get AT and GC contents from parsed alignment dictionary
        # return a list of tuples with taxon name, AT content, and GC content
        return sorted([(taxon, self.get_atgc_from_seq(seq)) \
         for taxon, seq in self.parsed_aln.items()])
        
    def get_atgc_from_seq(self, seq):
        # get AT and GC contents from individual sequences

        at_count = seq.count("A") + seq.count("T") + seq.count("W")
        gc_count = seq.count("G") + seq.count("C") + seq.count("S")
        
        try:
            at_content = round(at_count / (at_count + gc_count), 3)
            gc_content = round(1 - float(at_content), 3)

        except ZeroDivisionError:
            at_content = 0
            gc_content = 0

        return at_content, gc_content

class MetaAlignment():
    """Class of multiple sequence alignments"""
 
    def __init__(self, **kwargs):
        # set defaults and get values from kwargs
        self.in_files = kwargs.get("in_files")
        self.in_format = kwargs.get("in_format")
        self.data_type = kwargs.get("data_type")
        self.command = kwargs.get("command")
        self.concat_out = kwargs.get("concat_out", "concatenated.out")
        self.check_align = kwargs.get("check_align", False)
        self.cores = kwargs.get("cores")
        self.by_taxon_summary = kwargs.get("by_taxon_summary")
     
        if self.command == "replicate":
            self.no_replicates = kwargs.get("replicate_args")[0]
            self.no_loci = kwargs.get("replicate_args")[1]

        if self.command == "split":
            self.split = kwargs.get("split_by")
            self.remove_empty = kwargs.get("remove_empty", False)

        if self.command == "remove":
            self.species_to_remove = kwargs.get("taxa_to_remove")
            self.reduced_file_prefix = kwargs.get("out_prefix")
            self.check_taxa = kwargs.get("check_taxa", False)

        if self.command == "translate":
            self.reading_frame = kwargs.get("reading_frame")
            self.genetic_code = kwargs.get("genetic_code")
        if self.command == "trim":
            self.trim_fraction = kwargs.get("trim_fraction")
            self.trim_out = kwargs.get("trim_out")
            self.parsimony_check = kwargs.get("parsimony_check", False)

        self.alignment_objects = self.get_alignment_objects()
        self.parsed_alignments = self.get_parsed_alignments()

        # The code list:
        self.codes_list = """
        1. The Standard Code
        2. The Vertebrate Mitochondrial Code
        3. The Yeast Mitochondrial Code
        4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
        5. The Invertebrate Mitochondrial Code
        6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
        9. The Echinoderm and Flatworm Mitochondrial Code
        10. The Euplotid Nuclear Code
        11. The Bacterial, Archaeal and Plant Plastid Code
        12. The Alternative Yeast Nuclear Code
        13. The Ascidian Mitochondrial Code
        14. The Alternative Flatworm Mitochondrial Code
        16. Chlorophycean Mitochondrial Code
        21. Trematode Mitochondrial Code
        22. Scenedesmus obliquus Mitochondrial Code
        23. Thraustochytrium Mitochondrial Code
        24. Pterobranchia Mitochondrial Code
        25. Candidate Division SR1 and Gracilibacteria Code
        26. Pachysolen tannophilus Nuclear Code
        """
        
        # 1: The Standard Code
        self.gencode_NCBI_1 = {
        "TTT" : "F", # Phe
        "TCT" : "S", # Ser
        "TAT" : "Y", # Tyr
        "TGT" : "C", # Cys
        "TTC" : "F", # Phe
        "TCC" : "S", # Ser
        "TAC" : "Y", # Tyr
        "TGC" : "C", # Cys
        "TTA" : "L", # Leu
        "TCA" : "S", # Ser
        "TAA" : "*", # Ter
        "TGA" : "*", # Ter
        "TTG" : "L", # Leu i
        "TCG" : "S", # Ser
        "TAG" : "*", # Ter
        "TGG" : "W", # Trp
        "CTT" : "L", # Leu
        "CCT" : "P", # Pro
        "CAT" : "H", # His
        "CGT" : "R", # Arg
        "CTC" : "L", # Leu
        "CCC" : "P", # Pro
        "CAC" : "H", # His
        "CGC" : "R", # Arg
        "CTA" : "L", # Leu
        "CCA" : "P", # Pro
        "CAA" : "Q", # Gln
        "CGA" : "R", # Arg
        "CTG" : "L", # Leu i
        "CCG" : "P", # Pro
        "CAG" : "Q", # Gln
        "CGG" : "R", # Arg
        "ATT" : "I", # Ile
        "ACT" : "T", # Thr
        "AAT" : "N", # Asn
        "AGT" : "S", # Ser
        "ATC" : "I", # Ile
        "ACC" : "T", # Thr
        "AAC" : "N", # Asn
        "AGC" : "S", # Ser
        "ATA" : "I", # Ile
        "ACA" : "T", # Thr
        "AAA" : "K", # Lys
        "AGA" : "R", # Arg
        "ATG" : "M", # Met i
        "ACG" : "T", # Thr
        "AAG" : "K", # Lys
        "AGG" : "R", # Arg
        "GTT" : "V", # Val
        "GCT" : "A", # Ala
        "GAT" : "D", # Asp
        "GGT" : "G", # Gly
        "GTC" : "V", # Val
        "GCC" : "A", # Ala
        "GAC" : "D", # Asp
        "GGC" : "G", # Gly
        "GTA" : "V", # Val
        "GCA" : "A", # Ala
        "GAA" : "E", # Glu
        "GGA" : "G", # Gly
        "GTG" : "V", # Val
        "GCG" : "A", # Ala
        "GAG" : "E", # Glu
        "GGG" : "G", # Gly
        "---" : "-", # Gap
        "???" : "?", # Unk
        "NNN" : "X", # Unk
        }
        
        # 2: The Vertebrate Mitochondrial Code
        self.gencode_NCBI_2 = self.gencode_NCBI_1.copy()
        self.gencode_NCBI_2["AGA"] = "*" # Ter
        self.gencode_NCBI_2["AGG"] = "*" # Ter
        self.gencode_NCBI_2["ATA"] = "M" # Met
        self.gencode_NCBI_2["TGA"] = "W" # Trp
        
        # 3: The Yeast Mitochondrial Code
        self.gencode_NCBI_3 = self.gencode_NCBI_1.copy()
        self.gencode_NCBI_3["ATA"] = "M" # Met
        self.gencode_NCBI_3["CTT"] = "T" # Thr
        self.gencode_NCBI_3["CTC"] = "T" # Thr
        self.gencode_NCBI_3["CTA"] = "T" # Thr
        self.gencode_NCBI_3["CTG"] = "T" # Thr
        self.gencode_NCBI_3["TGA"] = "W" # Trp
        
        del self.gencode_NCBI_3["CGA"]
        del self.gencode_NCBI_3["CGC"]
        
        # 4: The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
        self.gencode_NCBI_4 = self.gencode_NCBI_1.copy()
        self.gencode_NCBI_4["TGA"] = "W" # Trp
        
        # 5: The Invertebrate Mitochondrial Code
        self.gencode_NCBI_5 = self.gencode_NCBI_1.copy()
        self.gencode_NCBI_5["AGA"] = "S" # Ser
        self.gencode_NCBI_5["AGG"] = "S" # Ser
        self.gencode_NCBI_5["ATA"] = "M" # Met
        self.gencode_NCBI_5["TGA"] = "W" # Trp
        
        # 6: The Ciliate, Dasycladacean and Hexamita Nuclear Code
        self.gencode_NCBI_6 = self.gencode_NCBI_1.copy()
        self.gencode_NCBI_6["TAA"] = "Q" # Gln
        self.gencode_NCBI_6["TAG"] = "Q" # Gln
        
        # 9: The Echinoderm and Flatworm Mitochondrial Code
        self.gencode_NCBI_9 = self.gencode_NCBI_1.copy()
        self.gencode_NCBI_9["AAA"] = "N" # Asn
        self.gencode_NCBI_9["AGA"] = "S" # Ser
        self.gencode_NCBI_9["AGG"] = "S" # Ser
        self.gencode_NCBI_9["TGA"] = "W" # Trp
        
        # 10: The Euplotid Nuclear Code
        self.gencode_NCBI_10 = self.gencode_NCBI_1.copy()
        self.gencode_NCBI_10["TGA"] = "C" # Cys
        
        # 11: The Bacterial, Archaeal and Plant Plastid Code
        self.gencode_NCBI_11 = self.gencode_NCBI_1.copy()
        
        # 12: The Alternative Yeast Nuclear Code
        self.gencode_NCBI_12 = self.gencode_NCBI_1.copy()
        self.gencode_NCBI_12["CTG"] = "S" # Ser
        
        # 13: The Ascidian Mitochondrial Code 
        self.gencode_NCBI_13 = self.gencode_NCBI_1.copy()
        self.gencode_NCBI_13["AGA"] = "G" # Gly
        self.gencode_NCBI_13["AGG"] = "G" # Gly
        self.gencode_NCBI_13["ATA"] = "M" # Met
        self.gencode_NCBI_13["TGA"] = "W" # Trp
        
        # 14: The Alternative Flatworm Mitochondrial Code
        self.gencode_NCBI_14 = self.gencode_NCBI_1.copy()
        self.gencode_NCBI_14["AAA"] = "N" # Asn
        self.gencode_NCBI_14["AGA"] = "S" # Ser
        self.gencode_NCBI_14["AGG"] = "S" # Ser
        self.gencode_NCBI_14["TAA"] = "Y" # Tyr
        self.gencode_NCBI_14["TGA"] = "W" # Trp
        
        # 16: Chlorophycean Mitochondrial Code
        self.gencode_NCBI_16 = self.gencode_NCBI_1.copy()
        self.gencode_NCBI_16["TAG"] = "L" # Leu
        
        # 21: Trematode Mitochondrial Code
        self.gencode_NCBI_21 = self.gencode_NCBI_1.copy()
        self.gencode_NCBI_21["TGA"] = "W" # Trp
        self.gencode_NCBI_21["ATA"] = "M" # Met
        self.gencode_NCBI_21["AGA"] = "S" # Ser
        self.gencode_NCBI_21["AGG"] = "S" # Ser
        self.gencode_NCBI_21["AAA"] = "N" # Asn
        
        # 22: Scenedesmus obliquus Mitochondrial Code
        self.gencode_NCBI_22 = self.gencode_NCBI_1.copy()
        self.gencode_NCBI_22["TCA"] = "*" # Ter
        self.gencode_NCBI_22["TAG"] = "L" # Leu
        
        # 23: Thraustochytrium Mitochondrial Code
        self.gencode_NCBI_23 = self.gencode_NCBI_1.copy()
        self.gencode_NCBI_23["TTA"] = "*" # Ter
        
        # 24: Pterobranchia Mitochondrial Code
        self.gencode_NCBI_24 = self.gencode_NCBI_1.copy()
        self.gencode_NCBI_24["AGA"] = "S" # Ser
        self.gencode_NCBI_24["AGG"] = "K" # Lys
        self.gencode_NCBI_24["TGA"] = "W" # Trp
        
        # 25: Candidate Division SR1 and Gracilibacteria Code
        self.gencode_NCBI_25 = self.gencode_NCBI_1.copy()
        self.gencode_NCBI_25["TGA"] = "G" # Gly
        
        # 26: Pachysolen tannophilus Nuclear Code
        self.gencode_NCBI_26 = self.gencode_NCBI_1.copy()
        self.gencode_NCBI_26["CTG"] = "A" # Ala
        
        self.codes = {
        1 : self.gencode_NCBI_1,
        2 : self.gencode_NCBI_2,
        3 : self.gencode_NCBI_3,
        4 : self.gencode_NCBI_4,
        5 : self.gencode_NCBI_5,
        6 : self.gencode_NCBI_6,
        9 : self.gencode_NCBI_9,
        10 : self.gencode_NCBI_10,
        11 : self.gencode_NCBI_11,
        12 : self.gencode_NCBI_12,
        13 : self.gencode_NCBI_13,
        14 : self.gencode_NCBI_14,
        16 : self.gencode_NCBI_16,
        21 : self.gencode_NCBI_21,
        22 : self.gencode_NCBI_22,
        23 : self.gencode_NCBI_23,
        24 : self.gencode_NCBI_24,
        25 : self.gencode_NCBI_25,
        26 : self.gencode_NCBI_26
        }

    def translate_dna_to_aa(self, seq, translation_table, frame):
    # translate DNA string into amino acids
        # where the last codon starts
        last_codon_start = len(seq) - 2
        # where the first codon starts
        if frame == 1:
            first = 0
        elif frame == 2:
            first = 1
        elif frame == 3:
            first = 2
        # create protein sequence by growing list
        protein = []
        add_to_protein = protein.append
        for start in range(first, last_codon_start, 3):
            codon = seq[start : start + 3]
            aa = translation_table.get(codon.upper(), 'X')
            add_to_protein(aa)           

        return "".join(protein)

    def translate_dict(self, source_dict):
        translation_table = self.codes.get(self.genetic_code)
        translated_dict = {}
        for taxon, seq in sorted(source_dict.items()):
            translated_seq = self.translate_dna_to_aa(seq, translation_table, self.reading_frame)
            if "*" in translated_seq:
                print("WARNING: stop codon(s), indicated as *, found in {} sequence".format(taxon)) 
            translated_dict[taxon] = translated_seq    

        return translated_dict

    def get_translated(self, translation_table, reading_frame):
        if int(self.cores) == 1:
            translated_alignments = [self.translate_dict(alignment) for alignment in self.parsed_alignments]            
        elif int(self.cores) > 1:
            pool = mp.Pool(int(self.cores))
            translated_alignments = pool.map(self.translate_dict, self.parsed_alignments)

        return translated_alignments

    def trim_dict(self, alignment):
        trim_vector = alignment.get_trim_selection(self.trim_fraction, self.parsimony_check)
        aln_dict = alignment.parsed_aln
        for key in aln_dict:
            aln_dict[key] = ''.join(list(compress(aln_dict[key], trim_vector)))

        return aln_dict

    def get_trimmed(self, trim_fraction, parsimony_check):
        if int(self.cores) == 1:
            trimmed_alignments = [self.trim_dict(alignment) for alignment in self.alignment_objects]            
        elif int(self.cores) > 1:
            pool = mp.Pool(int(self.cores))
            trimmed_alignments = pool.map(self.trim_dict, self.alignment_objects)

        return trimmed_alignments

        
    def remove_unknown_chars(self, seq):
        # remove unknown characters from sequence
        new_seq = seq.replace("?","").replace("-","")
        
        return new_seq

    def remove_empty_sequences(self, split_alignment):
        # remove taxa from alignment if they are composed of only empty sequences
        new_alignment = {taxon : seq for taxon, seq in split_alignment.items() if self.remove_unknown_chars(seq)}

        return new_alignment

    def get_partitions(self, partitions_file):
        # parse and get partitions from partitions file 
        partitions = FileParser(partitions_file)
        parsed_partitions = partitions.partitions_parse()
        
        return parsed_partitions

    def get_alignment_object(self, alignment):
        # parse according to the given alphabet
        if self.data_type == "aa":
            aln = AminoAcidAlignment(alignment, self.in_format, self.data_type)
        elif self.data_type == "dna":
            aln = DNAAlignment(alignment, self.in_format, self.data_type)
        return aln

    def get_alignment_objects(self):
        # get alignment objects on which statistics can be computed
        # use multiprocessing if more than one core specified
        if int(self.cores) == 1:
            alignments = [self.get_alignment_object(alignment) for alignment in self.in_files]            
        elif int(self.cores) > 1:
            pool = mp.Pool(int(self.cores))
            alignments = pool.map(self.get_alignment_object, self.in_files)
        return alignments

    def get_parsed_alignments(self):
        # get parsed dictionaries with taxa and sequences
        parsed_alignments = []
        add_to_parsed_alignments = parsed_alignments.append
        for alignment in self.alignment_objects:
            parsed = alignment.parsed_aln
            add_to_parsed_alignments(parsed)
            # checking if every seq has the same length or if parsed is not empty; exit if false
            if self.check_align == True:
                equal = all(x == [len(list(parsed.values())[i]) for i in range(0,len(list(parsed.values())))][0] 
                 for x in [len(list(parsed.values())[i]) for i in range(0,len(list(parsed.values())))])
                if equal is False:
                    print("ERROR: Sequences in input are of varying lengths. Be sure to align them first.")
                    sys.exit()

            if not parsed.keys() or not any(parsed.values()):
                print("ERROR: Parsed sequences of " + alignment.in_file + " are empty. "\
                 "Are you sure you specified the right input format and/or that input is a valid alignment?")
                sys.exit()
 
        return parsed_alignments

    def get_partitioned(self, partitions_file):
        # partition alignment according to a partitions file
        partitions = self.get_partitions(partitions_file)
        alignment = self.parsed_alignments[0]

        # initiate list of newly partitioned alignments 
        list_of_parts = []
        add_to_list_of_parts = list_of_parts.append
        for partition in partitions:
            # loop over all parsed partitions, adding taxa and sliced sequences
            for name, elements in partition.items():    
                new_dict = {}
         
                for taxon, seq in alignment.items():
                    new_seq = ""
         
                    for dictionary in elements:
                        new_seq = new_seq + seq[dictionary["start"]:dictionary["stop"]:dictionary["stride"]]
                        new_dict[taxon] = new_seq

            if self.remove_empty:
            # check if remove empty sequences
                no_empty_dict = self.remove_empty_sequences(new_dict)  
                add_to_list_of_parts({name : no_empty_dict})
            else:
            # add partition name : dict of taxa and sequences to the list
                add_to_list_of_parts({name : new_dict})
    
        return list_of_parts

    def get_summaries(self):
        # get summaries for all alignment objects

        # define different headers for aa and dna alignments
        aa_header = [
            "Alignment_name",
            "No_of_taxa",
            "Alignment_length",
            "Total_matrix_cells",
            "Undetermined_characters",
            "Missing_percent",
            "No_variable_sites",
            "Proportion_variable_sites",
            "Parsimony_informative_sites",
            "Proportion_parsimony_informative"
        ]

        dna_header = [
            "Alignment_name",
            "No_of_taxa",
            "Alignment_length",
            "Total_matrix_cells",
            "Undetermined_characters",
            "Missing_percent",
            "No_variable_sites",
            "Proportion_variable_sites",
            "Parsimony_informative_sites",
            "Proportion_parsimony_informative",
            "AT_content",
            "GC_content"
        ]

        alignments = self.alignment_objects
        parsed_alignments = self.parsed_alignments
        freq_header = [char for char in alignments[0].alphabet]
        
        if self.data_type == "aa":
            header = aa_header + freq_header
        elif self.data_type == "dna":
            header = dna_header + freq_header

        # use multiprocessing if more than one core specified
        if int(self.cores) == 1:
            summaries = [alignment.get_summary() for alignment in alignments]            
        elif int(self.cores) > 1:
            pool = mp.Pool(int(self.cores))
            summaries = pool.map(self.summarize_alignments, alignments)
        return header, summaries

    def summarize_alignments(self, alignment):
        # helper function to summarize alignments
        summary = alignment.get_summary()
        return summary

    def get_taxon_summaries(self):
        # get per-sequence summaries for all alignment objects

        # define different headers for aa and dna alignments
        aa_header = [
            "Alignment_name",
            "Taxon_name",
            "Sequence_length",
            "Undetermined_characters",
            "Missing_percent"
        ]

        dna_header = [
            "Alignment_name",
            "Taxon_name",
            "Sequence_length",
            "Undetermined_characters",
            "Missing_percent",
            "AT_content",
            "GC_content"
        ]

        alignments = self.alignment_objects
        parsed_alignments = self.parsed_alignments
        freq_header = alignments[0].alphabet
        
        if self.data_type == "aa":
            header = aa_header + freq_header
        elif self.data_type == "dna":
            header = dna_header + freq_header

        # use multiprocessing if more than one core specified
        if int(self.cores) == 1:
            summaries = [alignment.get_taxa_summary() for alignment in alignments]            
        elif int(self.cores) > 1:
            pool = mp.Pool(int(self.cores))
            summaries = pool.map(self.summarize_alignments_taxa, alignments)
           
        return header, summaries

    def summarize_alignments_taxa(self, alignment):
        # helper function to summarize alignments by taxon
        summary = alignment.get_taxa_summary()
        return summary

    def write_summaries(self, file_name):
        # write summaries to file

        self.file_overwrite_error(file_name)        

        summary_file = open(file_name, "w")
        summary_out = self.get_summaries()
        header = '\t'.join(summary_out[0])
        new_summ = ['\t'.join(summary) for summary in summary_out[1]]
        summary_file.write(header + '\n')
        summary_file.write('\n'.join(new_summ))
        summary_file.close()
        print("Wrote summaries to file '" + file_name + "'")

    def write_taxa_summaries(self):
        # write by-taxon summaries to file
        for index, in_file_name in enumerate(self.in_files):
            out_file_name = in_file_name + "-seq-summary.txt"
            self.file_overwrite_error(out_file_name)
            summary_file = open(out_file_name, "w")
            summary_out = self.get_taxon_summaries()
            header = '\t'.join(summary_out[0])
            summ = [[str(col) for col in element] for element in summary_out[1][index]]            
            new_summ = ['\t'.join(row) for row in summ]
            summary_file.write(header + '\n')
            summary_file.write('\n'.join(new_summ))
            summary_file.close()
       
    def get_replicate(self, no_replicates, no_loci):
        # construct replicate data sets for phylogenetic jackknife
        replicates = []
        add_to_replicates = replicates.append
        counter = 1
        for replicate in range(no_replicates):
            
            try:
                random_alignments = sample(self.parsed_alignments, no_loci)
            except ValueError:
                print("ERROR: You specified more loci per replicate than there are in your input.")
                sys.exit()

            random_alignments = sample(self.parsed_alignments, no_loci)
            concat_replicate = self.get_concatenated(random_alignments)[0]
            add_to_replicates(concat_replicate)
            counter += 1
        
        return replicates 

    def get_concatenated(self, alignments):
        # concatenate muntiple input alignments
        # create empty dictionary of lists
        concatenated = defaultdict(list)

        # first create list of taxa in all alignments
        # you need this to insert empty seqs in
        # the concatenated alignment
        all_taxa = []
        for alignment in alignments:
            for taxon in alignment.keys():
                if taxon not in all_taxa:
                    all_taxa.append(taxon)

        # start counters to keep track of partitions
        partition_counter = 1
        position_counter = 1
        # get dict for alignment name and partition
        partitions = {}

        for alignment in alignments:        
            
            # get alignment length from a random taxon
            partition_length = len(alignment[list(alignment.keys())[0]])
            # get base name of each alignment for use when writing partitions file
            # NOTE: the base name here is whatever comes before fist perion in the file name
            alignment_name = self.alignment_objects[partition_counter - 1].get_name().split('.')[0]
            # add a prefix to the partition names
            partition_name = "p" + str(partition_counter) + "_" + alignment_name
            
            start = position_counter
            position_counter += partition_length
            end = position_counter - 1
            partitions[partition_name] = str(start) + "-" + str(end)
            partition_counter += 1
            
            # get empty sequence if there is missing taxon
            # getting length from first element of list of keys
            # created from the original dict for this alignment
            empty_seq = '?' * partition_length

            for taxon in all_taxa:

                if taxon not in alignment.keys():
                    concatenated[taxon].append(empty_seq)
                else:
                    concatenated[taxon].append(alignment[taxon])
 
        concatenated = {taxon:''.join(seqs) for taxon, seqs in concatenated.items()}
        
        return concatenated, partitions

    def remove_from_alignment(self, alignment, species_to_remove, index):
        # remove taxa from alignment
        aln_name = self.get_alignment_name_no_ext(index)
        for taxon in species_to_remove:
            if taxon not in alignment.keys():
                print("WARNING: Taxon '" + taxon + "' not found in '" + aln_name + "'.\nIf you expected it to be there, make sure to replace all taxon name spaces with underscores and that you are not using quotes.")

            new_alignment = {species: seq for species, seq in alignment.items() if species not in species_to_remove}

            aln_tuple = (aln_name, new_alignment)

        return aln_tuple

    def remove_taxa(self, species_to_remove):
        new_alns = {}
        for index, alignment in enumerate(self.parsed_alignments):
            aln_name, aln_dict = self.remove_from_alignment(alignment, species_to_remove, index)
            # check if alignment is not empty:
            if aln_dict:
                new_alns[aln_name] = aln_dict
            else:
                print("ERROR: You asked to remove all taxa from the alignment " + aln_name + ". No output file will be written.")

        return new_alns

    def print_fasta(self, source_dict):
        # print fasta-formatted string from a dictionary        
        fasta_string = ""
        # each sequence line will have 80 characters 
        n = 80
        
        for taxon, seq in sorted(source_dict.items()):
            # split dictionary values to a list of string, each n chars long
            seq = [seq[i:i+n] for i in range(0, len(seq), n)]
            # in case there are unwanted spaces in taxon names
            taxon = taxon.replace(" ","_").strip("'")
            fasta_string += ">" + taxon + "\n"
            for element in seq:
                fasta_string += element + "\n"

        return fasta_string

    def print_phylip(self, source_dict):
        # print phylip-formatted string from a dictionary
        taxa_list = list(source_dict.keys())
        no_taxa = len(taxa_list)
        # figure out the max length of a taxon for nice padding of sequences
        pad_longest_name = len(max(taxa_list, key=len)) + 3
        # get sequence length from a random value
        seq_length = len(next(iter(source_dict.values())))
        header = str(len(source_dict)) + " " + str(seq_length)
        phylip_string = header + "\n"
        for taxon, seq in sorted(source_dict.items()):
            taxon = taxon.replace(" ","_").strip("'")
            # left-justify taxon names relative to sequences
            phylip_string += taxon.ljust(pad_longest_name, ' ') + seq + "\n"
 
        return phylip_string

    def print_phylip_int(self, source_dict):
        # print phylip interleaved-formatted string from a dictionary        
        taxa_list = list(source_dict.keys())
        no_taxa = len(taxa_list)
        pad_longest_name = len(max(taxa_list, key=len)) + 3
        seq_length = len(next(iter(source_dict.values())))
        header = str(len(source_dict)) + " " + str(seq_length)
        phylip_int_string = header + "\n\n"
        # this will be a list of tuples to hold taxa names and sequences
        seq_matrix = []
        
        # each sequence line will have 500 characters
        n = 500
        
        # recreate sequence matrix
        add_to_matrix = seq_matrix.append
        for taxon, seq in sorted(source_dict.items()):
            add_to_matrix((taxon, [seq[i:i+n] for i in range(0, len(seq), n)]))

        first_seq = seq_matrix[0][1]
        for index, item in enumerate(first_seq):
            for taxon, sequence in seq_matrix:
                if index == 0:
                    phylip_int_string += taxon.ljust(pad_longest_name, ' ') + sequence[index] + "\n"
                else:
                    phylip_int_string += sequence[index] + "\n"
            phylip_int_string += "\n"
         
        return phylip_int_string

    def print_nexus(self, source_dict):
        # print nexus-formatted string from a dictionary    
        if self.data_type == "aa" or self.command == "translate":
            data_type = "PROTEIN"
        elif self.data_type == "dna":
            data_type = "DNA"
        
        taxa_list = list(source_dict.keys())
        no_taxa = len(taxa_list)
        pad_longest_name = len(max(taxa_list, key=len)) + 3
        seq_length = len(next(iter(source_dict.values())))
        header = str(len(source_dict)) + " " + str(seq_length)
        nexus_string = "#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS  NTAX=" + str(no_taxa) +\
         " NCHAR=" + str(seq_length) + ";\n\tFORMAT DATATYPE=" + data_type +\
          "  GAP = - MISSING = ?;\n\tMATRIX\n"

        for taxon, seq in sorted(source_dict.items()):
            taxon = taxon.replace(" ","_").strip("'")
            nexus_string += "\t" + taxon.ljust(pad_longest_name, ' ') + seq + "\n"
        nexus_string += "\n;\n\nEND;"
        
        return nexus_string

    def print_nexus_int(self, source_dict):
        # print nexus interleaved-formatted string from a dictionary

        if self.data_type == "aa":
            data_type = "PROTEIN"
        elif self.data_type == "dna":
            data_type = "DNA"
        
        taxa_list = list(source_dict.keys())
        no_taxa = len(taxa_list)
        pad_longest_name = len(max(taxa_list, key=len)) + 3
        seq_length = len(next(iter(source_dict.values())))
        header = str(len(source_dict)) + " " + str(seq_length)
        # this will be a list of tuples to hold taxa names and sequences
        seq_matrix = []
        
        nexus_int_string = "#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS  NTAX=" +\
         str(no_taxa) + " NCHAR=" + str(seq_length) + ";\n\tFORMAT   INTERLEAVE" +\
          "   DATATYPE=" + data_type + "  GAP = - MISSING = ?;\n\tMATRIX\n"

        n = 500
        
        # recreate sequence matrix
        add_to_matrix = seq_matrix.append
        for taxon, seq in sorted(source_dict.items()):
            add_to_matrix((taxon, [seq[i:i+n] for i in range(0, len(seq), n)]))

        first_seq = seq_matrix[0][1]
        for index, item in enumerate(first_seq):
            for taxon, sequence in seq_matrix:
                if index == 0:
                    nexus_int_string += taxon.ljust(pad_longest_name, ' ') + sequence[index] + "\n"
                else:
                    nexus_int_string += sequence[index] + "\n"
            nexus_int_string += "\n"

        nexus_int_string += "\n;\n\nEND;"
        
        return nexus_int_string

    def natural_sort(self, a_list):
        # create a function that does 'human sort' on a list
        convert = lambda text: int(text) if text.isdigit() else text.lower() 
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)] 
        return sorted(a_list, key = alphanum_key)

    def print_unspecified_partitions(self):
        # print partitions for concatenated alignment
        part_string = ""
        part_dict = self.get_concatenated(self.parsed_alignments)[1]
        part_list = self.natural_sort(part_dict.keys())
        for key in part_list:
            part_string += key + " = " + str(part_dict[key]) + "\n"
        return part_string

    def print_nexus_partitions(self):
        # print partitions for concatenated alignment
        part_string = ""
        part_dict = self.get_concatenated(self.parsed_alignments)[1]
        part_list = self.natural_sort(part_dict.keys())
        # write beginning of nexus sets
        part_string += "#NEXUS\n\n"
        part_string += "Begin sets;\n"
        for key in part_list:
            part_string += "\tcharset " + key + " = " + str(part_dict[key]) + ";\n"
        part_string += "End;"
        return part_string

    def print_raxml_partitions(self, data_type):
        # print partitions for concatenated alignment
        part_string = ""
        part_dict = self.get_concatenated(self.parsed_alignments)[1]
        part_list = self.natural_sort(part_dict.keys())
        if data_type == "dna":
            for key in part_list:
                part_string += "DNA, " + key + " = " + str(part_dict[key]) + "\n"
        if data_type == "aa":
            for key in part_list:
                part_string += "WAG, " + key + " = " + str(part_dict[key]) + "\n"
        return part_string

    def write_partitions(self, file_name, part_format):
        # write partitions file for concatenated alignment
        self.file_overwrite_error(file_name)        
        part_file = open(file_name, "w")
        if part_format == "nexus":
            part_file.write(self.print_nexus_partitions())
        if part_format == "raxml":
            part_file.write(self.print_raxml_partitions(self.data_type))
        if part_format == "unspecified":
            part_file.write(self.print_unspecified_partitions())
        print("Wrote partitions for the concatenated file to '" + file_name + "'")

    def get_extension(self, file_format):
        # get proper extension string
        if file_format == "phylip":
            extension = "-out.phy"
        elif file_format == "phylip-int":
            extension = "-out.int-phy"
        elif file_format == "fasta":
            extension = "-out.fas"
        elif file_format == "nexus":
            extension = "-out.nex"
        elif file_format == "nexus-int":
            extension = "-out.int-nex"
        
        return extension

    def file_overwrite_error(self, file_name):
        # print warning when overwriting a file
        if path.exists(file_name):
            print("WARNING: You are overwriting '" + file_name + "'")

    def write_formatted_file(self, file_format, file_name, alignment):
        # write the correct format string into a file        
        out_file = open(file_name, "w")
        if file_format == "phylip":
            out_file.write(self.print_phylip(alignment))
        elif file_format == "fasta":
            out_file.write(self.print_fasta(alignment))
        elif file_format == "phylip-int":
            out_file.write(self.print_phylip_int(alignment))
        elif file_format == "nexus":
            out_file.write(self.print_nexus(alignment))
        elif file_format == "nexus-int":
            out_file.write(self.print_nexus_int(alignment))
        out_file.close()

    def get_alignment_name(self, i, extension):
        # get file name
        file_name = self.alignment_objects[i].get_name() + extension

        return file_name

    def get_alignment_name_no_ext(self, i):
        # get file name without extension
        file_name = self.alignment_objects[i].get_name()

        return file_name

    def write_concat(self, file_format):
        # write concatenated alignment into a file
        concatenated_alignment = self.get_concatenated(self.parsed_alignments)[0]
        file_name = self.concat_out
        self.file_overwrite_error(file_name)
        self.write_formatted_file(file_format, file_name, concatenated_alignment)
        
        print("Wrote concatenated sequences to " + file_format + " file '" + file_name + "'")

    def write_convert(self, index, alignment, file_format, extension):
        # write converted alignment into a file
        file_name = self.get_alignment_name(index, extension)
        self.file_overwrite_error(file_name)        
        self.write_formatted_file(file_format, file_name, alignment)

    def write_replicate(self, index, alignment, file_format, extension):
        # write replicate alignment into a file
        file_name = "replicate" + str(index + 1) + "_" + str(self.no_loci) + "-loci" + extension
        self.file_overwrite_error(file_name)                        
        self.write_formatted_file(file_format, file_name, alignment)

    def write_split(self, index, item, file_format, extension):
        # write split alignments from partitions file
        # bad practice with the dicts; figure out better solution
        try:       
            file_name = str(self.in_files[0].split('.')[0]) + "_" + list(item.keys())[0] + extension
            alignment = list(item.values())[0]
            self.file_overwrite_error(file_name)
            self.write_formatted_file(file_format, file_name, alignment)
        except ValueError:
            print("WARNING: There was no data to write for file '" + file_name + "'. Perhaps a partition composed of missing data only?")
            remove(file_name)
            raise ValueError

    def write_reduced(self, file_format, extension):
        # write alignment with taxa removed into a file
        prefix =  self.reduced_file_prefix
        alns = self.remove_taxa(self.species_to_remove)
        for file_name, aln_dict in alns.items():
            out_file_name = prefix + file_name + extension
            self.file_overwrite_error(out_file_name)          
            self.write_formatted_file(file_format, out_file_name, aln_dict)
        return len(alns)

    def write_translated(self, index, alignment, file_format, extension):
        # write alignments translated into amino acids
        prefix = "translated_"
        file_name = self.get_alignment_name(index, extension)
        out_file_name = prefix + file_name + extension
        self.file_overwrite_error(out_file_name)   
        self.write_formatted_file(file_format, out_file_name, alignment)

    def write_trimmed(self, index, alignment, file_format, extension):
        # write trimmed alignments
        if self.trim_out: 
            out_file_name = self.trim_out
        else:
            prefix = "trimmed_"
            file_name = self.get_alignment_name(index, extension)
            out_file_name = prefix + file_name
        self.file_overwrite_error(out_file_name)   
        self.write_formatted_file(file_format, out_file_name, alignment)

    def write_out(self, action, file_format):
        # write other output files depending on command (action) 
        extension = self.get_extension(file_format)

        if action == "concat":        
            self.write_concat(file_format)    

        elif action == "convert":        
            length = len(self.alignment_objects)
            [self.write_convert(i, alignment, file_format, extension) \
             for i, alignment in enumerate(self.parsed_alignments)]
            print("Converted " + str(length) + " files from " + self.in_format + " to " + file_format)

        elif action == "replicate":
            [self.write_replicate(i, alignment, file_format, extension) \
             for i, alignment in enumerate(self.get_replicate(self.no_replicates, self.no_loci))]

            print("Constructed " + str(self.no_replicates) + " replicate data sets, each from " \
             + str(self.no_loci) + " alignments")

        elif action == "split":
            list_of_alignments = self.get_partitioned(self.split)
            length = len(list_of_alignments)
            err_indx = 0
            for i, item in enumerate(list_of_alignments):
                try:
                    self.write_split(i, item, file_format, extension)
                except ValueError:
                    err_indx += 1
                    pass
            print("Wrote " + str(length - err_indx) + " " + str(file_format) + " files from partitions provided")

        elif action == "remove":
            aln_no = self.write_reduced(file_format, extension)
            if aln_no:
                print("Wrote " + str(aln_no) + " " + str(file_format) + " files with reduced taxon set")

        elif action == "translate":
            if self.data_type == "aa":
                print("ERROR: cannot translate; you said your alignment already contains amino acids")
                sys.exit()
            translated_alignment_dicts = self.get_translated(self.genetic_code, self.reading_frame)
            length = len(self.alignment_objects)
            [self.write_translated(i, alignment, file_format, extension) \
             for i, alignment in enumerate(translated_alignment_dicts)]
            print("Translated " + str(length) + " files to amino acid sequences")

        elif action == "trim": # self.trim_fraction, self.parsimony_check
            trimmed_alignment_dicts = self.get_trimmed(self.trim_fraction, self.parsimony_check)
            length = len(self.alignment_objects)
            [self.write_trimmed(i, alignment, file_format, extension) \
             for i, alignment in enumerate(trimmed_alignment_dicts)]
            print("Trimmed", str(length), "file(s) to have", self.trim_fraction, "minimum occupancy per alignment column") 


def main():
    
    # initialize parsed arguments and meta alignment objects
    kwargs = run()
    meta_aln = MetaAlignment(**kwargs)
       
    if meta_aln.command == "summary":
        meta_aln.write_summaries(kwargs["summary_out"])
    if meta_aln.by_taxon_summary:
        print("Printing taxon summaries")
        meta_aln.write_taxa_summaries()
    if meta_aln.command == "convert":
        meta_aln.write_out("convert", kwargs["out_format"])
    if meta_aln.command == "concat":
        meta_aln.write_out("concat", kwargs["out_format"])
        meta_aln.write_partitions(kwargs["concat_part"], kwargs["part_format"])
    if meta_aln.command == "replicate":
        meta_aln.write_out("replicate", kwargs["out_format"])
    if meta_aln.command == "split":
        meta_aln.write_out("split", kwargs["out_format"])
    if meta_aln.command == "remove":
        meta_aln.write_out("remove", kwargs["out_format"])
    if meta_aln.command == "translate":
        meta_aln.write_out("translate", kwargs["out_format"])
    if meta_aln.command == "trim":
        meta_aln.write_out("trim", kwargs["out_format"])

        # meta_aln.write_out("translate", kwargs["out_format"])
  
def run():

    # initialize parsed arguments
    config = ParsedArgs()
    # get arguments
    config_dict = config.get_args_dict()
    return config_dict
    
if __name__ == '__main__':
        
        main()
