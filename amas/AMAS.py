#! /usr/bin/env python3

#   Program to calculate various statistics on a multiple sequence alignment

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
This stand-alone program calculates various statistics on a multiple sequence
alignment. It supports sequential FASTA, PHYLIP, NEXUS, and interleaved PHYLIP 
and NEXUS formats for DNA and aino acid sequences.

Current statistics include the number of taxa, alignment length, total number
of matrix cells, overall number of undetermined characters, percent of missing 
data, AT and GC contents (for DNA alignments), number and proportion of 
variable sites, number and proportion of parsimony informative sites,
and proportions of all characters relative to matrix size.
"""

import sys
from sys import argv
import re

Usage = """
Usage: AMAS.py <input_file> <format> <alphabet>

Supported formats: "fasta", "phylip", "nexus", "phylip-int", "nexus-int"
Supported alphabets: "aa", "dna"
"""


class FileHandler:
    """Define file handle that closes when out of scope"""

    def __init__(self, file_name):
        self.file_name = file_name

    def __enter__(self):
        self.in_file = open(self.file_name, "r")
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
        matches = re.finditer(r"^>(.*[^$])([^>]*)", self.in_file_lines, re.MULTILINE)
        records = {}
 
        for match in matches:
            name_match = match.group(1).replace("\n","")
            seq_match = match.group(2).replace("\n","").upper()
            seq_match = self.translate_ambiguous(seq_match)
            records[name_match] = seq_match

        return records
    
    def phylip_parse(self):
    # use regex to parse names and sequences in sequential phylip files
        matches = re.finditer(r"^(\s+)?(\S+)\s+([A-Za-z*?.{}-]+)", \
         self.in_file_lines, re.MULTILINE)
        records = {}

        for match in matches:
            name_match = match.group(2).replace("\n","")
            seq_match = match.group(3).replace("\n","").upper()
            seq_match = self.translate_ambiguous(seq_match)
            records[name_match] = seq_match
        return records    

    def phylip_interleaved_parse(self):
    # use regex to parse names and sequences in interleaved phylip files
        name_matches = re.finditer(r"^(\s+)?(\S+)[ \t]+[A-Za-z*?.{}-]+", \
          self.in_file_lines, re.MULTILINE)
        seq_matches = re.finditer(r"(^(\s+)?\S+[ \t]+|^)([A-Za-z*?.{}-]+)$", \
         self.in_file_lines, re.MULTILINE)
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
        matches = re.finditer(r"(\s+)?(MATRIX\n|matrix\n|MATRIX\r\n|matrix\r\n)(.*?;)", \
         self.in_file_lines, re.DOTALL)
        records = {}
        # get names and sequences from the matrix block

        for match in matches:
            matrix_match = match.group(3)
            seq_matches = \
             re.finditer(r"^(\s+)?[']?(\S+\s\S+|\S+)[']?\s+([A-Za-z*?.{}-]+)($|\s+\[[0-9]+\]$)", \
              matrix_match, re.MULTILINE)

            for match in seq_matches:
                name_match = match.group(2).replace("\n","")
                seq_match = match.group(3).replace("\n","").upper()
                seq_match = self.translate_ambiguous(seq_match)
                records[name_match] = seq_match

        return records
        
    def nexus_interleaved_parse(self):
    # use regex to parse names and sequences in sequential nexus files
    # find the matrix block
        matches = re.finditer(r"(\s+)?(MATRIX\n|matrix\n|MATRIX\r\n|matrix\r\n)(.*?;)", \
         self.in_file_lines, re.DOTALL)
        # initiate lists for taxa names and sequence strings on separate lines
        taxa = []
        sequences = []
        # initiate a dictionary for the name:sequence records
        records = {}
        # initiate a counter to keep track of sequences strung together
        # from separate lines
        counter = 0

        for match in matches:
            matrix_match = match.group(3)
            # get names and sequences from the matrix block
            seq_matches = \
             re.finditer(r"^(\s+)?[']?(\S+\s\S+|\S+)[']?\s+([A-Za-z*?.{}-]+)($|\s+\[[0-9]+\]$)", \
              matrix_match, re.MULTILINE)

            for match in seq_matches:
                name_match = match.group(2).replace("\n","")
                if name_match not in taxa:
                    taxa.append(name_match)
                seq_match = match.group(3).replace("\n","").upper()
                seq_match = self.translate_ambiguous(seq_match)
                sequences.append(seq_match)

        for taxon_no in range(len(taxa)):

            full_length_sequence = ""
            for index in range(counter,len(sequences),len(taxa)):
                full_length_sequence += sequences[index]
            
            counter += 1 
            records[taxa[taxon_no]] = full_length_sequence

        return records

    def translate_ambiguous(self, seq):
    # translate ambiguous characters from curly bracket format
    # to single letter format 
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
        return seq

 
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
        self.list_of_seqs = []
        self.no_missing_ambiguous_sites = []
        self.parsimony_informative = 0
        
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
        else:
            print(Usage)

        return parsed_aln

    def summarize_alignment(self):
    # call methods to create sequences list, matrix, sites without ambiguous or
    # missing characters; get and summarize alignment statistics
        summary = []
        self.seq_grabber()
        self.matrix_creator()
        self.get_sites_no_missing_ambiguous()
        variable_sites = str(self.get_variable())
        prop_variable = str(self.get_prop_variable())
        parsimony_informative = str(self.get_parsimony_informative())
        prop_parsimony = str(self.get_prop_parsimony())
        name = str(self.get_name())
        taxa_no = str(self.get_taxa_no())
        length = str(self.get_alignment_length())
        cells = str(self.get_matrix_cells())
        missing = str(self.get_missing())
        missing_percent = str(self.get_missing_percent())
        summary = [name, taxa_no, length, cells, missing, missing_percent, \
         variable_sites, prop_variable, parsimony_informative, prop_parsimony]
        return summary

    def get_freq_summary(self):
    # get summary of frequencies for all characters
        characters = []
        frequencies = []
        
        for item in self.get_frequencies():
            for char, freq in item.items():
                characters.append(str(char))
                frequencies.append(str(freq))
        return characters, frequencies           
        
    def seq_grabber(self):
    # create a list of sequences from parsed dictionary of names and seqs 
        parsed_aln = self.get_parsed_aln()
        self.list_of_seqs = [seq for name, seq in parsed_aln.items()]
        return self.list_of_seqs
               
    def matrix_creator(self):
    # decompose character matrix into a two-dimensional list
        self.matrix = [[character for character in sequence] \
         for sequence in self.list_of_seqs]
        return self.matrix

    def get_column(self, i):
    # get site from the character matrix
        return [row[i] for row in self.matrix]
        
    def all_same(self, site):
    # check if all elements of a site are the same
        return all(base == site[0] for base in site)
        
    def get_sites_no_missing_ambiguous(self):
    # get each site without missing or ambiguous characters  
        for column in range(self.get_alignment_length()):
            site = self.get_column(column)
            site = [char for char in site if char not in self.missing_ambiguous_chars]
            self.no_missing_ambiguous_sites.append(site)
        return self.no_missing_ambiguous_sites
        
    def get_variable(self):
    # if all elements of a site without missing or ambiguous characters 
    # are not the same, consider it variable
        self.variable = len([site for site in self.no_missing_ambiguous_sites \
         if not self.all_same(site)])      
        return self.variable
    
    def get_parsimony_informative(self):
    # if the count for a unique character in a site is at least two, 
    # and there are at least two such characters in a site without missing
    # or ambiguous characters, consider it parsimony informative 
        for site in self.no_missing_ambiguous_sites:
            unique_chars = set(site)
            
            pattern = [base for base in unique_chars if site.count(base) >= 2]
            no_patterns = len(pattern)
            
            if no_patterns >= 2:
                self.parsimony_informative += 1
        return self.parsimony_informative
    
    def get_prop_variable(self):
    # get proportion of variable sites to all sites
        prop_variable = self.variable / len(self.list_of_seqs[0])
        return round(prop_variable, 3)
        
    def get_prop_parsimony(self):
    # get proportion of parsimony informative sites to all sites
        prop_parsimony = self.parsimony_informative / len(self.list_of_seqs[0])
        return round(prop_parsimony, 3)

    def get_name(self):
        return self.in_file
        
    def get_taxa_no(self):
        return len(self.list_of_seqs)
    
    def get_alignment_length(self):
        return len(self.list_of_seqs[0])

    def get_matrix_cells(self):
        self.all_matrix_cells = len(self.list_of_seqs) \
         * len(self.list_of_seqs[0])
        return self.all_matrix_cells

    def get_missing_percent(self):
        missing_percent = round((self.missing / self.all_matrix_cells * 100), 3)
        return missing_percent
        
    def get_missing(self):
        self.missing = sum(sum(seq.count(char) for seq in self.list_of_seqs) \
         for char in self.missing_chars)
        return self.missing
    
    def get_frequencies(self):
    # get frequencies of each character in the used alphabet
        frequencies = []
        
        for char in self.alphabet:
            count = sum(seq.count(char) for seq in self.list_of_seqs) \
             / self.all_matrix_cells
            frequencies.append({char : round(count, 3)})
        return frequencies


class AminoAcidAlignment(Alignment):
    """Summary specific to aa alignments"""

    alphabet = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R", \
     "S","T","V","W","Y","B","J","Z","X",".","*","-","?"]
    missing_ambiguous_chars = ["B","J","Z","X",".","*","-","?"]
    missing_chars = ["X",".","*","-","?"]

    def get_summary(self):
        data = self.summarize_alignment()
        new_data = data + list(self.get_freq_summary()[1])
        
        header = ["Alignment_name", "No_of_taxa", "Alignment_length", \
         "Total_matrix_cells", "Undetermined_characters", "Missing_percent", \
          "No_variable_sites", "Proportion_variable_sites", \
           "Parsimony_informative_sites", "Proportion_parsimony_informative"]
        
        freq_header = list(self.get_freq_summary()[0])
        new_header = header + freq_header
        print("\t".join(new_header))
        print("\t".join(new_data))

           
class DNAAlignment(Alignment):
    """Summary specific to DNA alignments"""
    
    alphabet = ["A","C","G","T","K","M","R","Y","S","W","B","V","H","D","X", \
     "N", "O", "-","?"]
    missing_ambiguous_chars = ["K","M","R","Y","S","W","B","V","H","D","X", \
     "N", "O", "-","?"] 
    missing_chars = ["X","N","O","-","?"]

    def get_summary(self):
        data = self.summarize_alignment()
        
        new_data = data + self.get_atgc_content() \
         + list(self.get_freq_summary()[1])
        header = ["Alignment_name", "No_of_taxa", "Alignment_length", \
         "Total_matrix_cells", "Undetermined_characters", "Missing_percent", \
          "No_variable_sites", "Proportion_variable_sites", "Parsimony_informative_sites", \
            "Proportion_parsimony_informative", "AT_content", "GC_content"]
        
        freq_header = list(self.get_freq_summary()[0])
        new_header = header + freq_header
        
        print("\t".join(new_header))
        print("\t".join(new_data))
        
    def get_atgc_content(self):
    # get AC and GC contents
        atgc_content = []
        
        at_count = sum((seq.count("A") + seq.count("T") + seq.count("W")) \
         for seq in self.list_of_seqs)
        gc_count = sum((seq.count("G") + seq.count("C") + seq.count("S")) \
         for seq in self.list_of_seqs)
        
        at_content = str(round(at_count / (at_count + gc_count), 3))
        gc_content = str(round(1 - float(at_content), 3))
        
        atgc_content.extend((at_content, gc_content))
        return atgc_content


def main():
    # print usage instructions if number of arguments given is incorrect
    if len(argv) is not 4:
        print(Usage)
        sys.exit()

    # define variables from arguments given
    script, in_file, in_format, data_type = argv

    # parse according to the given alphabet
    if data_type == "aa":
        aln = AminoAcidAlignment(in_file, in_format, data_type)
    elif data_type == "dna":
        aln = DNAAlignment(in_file, in_format, data_type)
    else:
        print(Usage)

    # get alignment summary
    aln.get_summary()


if __name__ == '__main__':
    main()
