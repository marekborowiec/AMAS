#! /usr/bin/env python3

"""
AMAS_v0.2.py by Marek Borowiec

This stand-alone program calculates various statistics on a multiple sequence
alignment. It supports sequential FASTA, PHYLIP, and NEXUS formats for DNA and 
amino acid sequences.

Current statistics include the number of taxa, alignment length, total number
of matrix cells, overall number of undetermined characters, percent of missing 
data, AT and GC contents (for DNA alignments), number and proportion of 
variable sites, number and proportion of parsimony informative sites,
and relative proportions of all characters relative to matrix size.
"""

Usage = """
Usage: AMAS_v0.2.py <input_file> <format> <alphabet>

Supported formats: "fasta", "phylip", "nexus"
Supported alphabets: "aa", "dna"
"""

from sys import argv
import re


class FileParser:
    """Parse file contents and return sequences and sequence names"""

    def __init__(self, file_name):
    # initialize file parser with file name and read in file contents
        self.file_name = file_name
        self.in_file = open(file_name, "r")
        self.in_file_lines = self.in_file.read().rstrip("\r\n")
    
    def get_file_name(self):
        return self.file_name

    def fasta_parse(self):
    # use regex to parse names and sequences in sequential fasta files
        matches = re.finditer(r"^>(.*[^$])([^>]*)", self.in_file_lines, re.MULTILINE)
        records = {}
        
        for match in matches:
            name_match = match.group(1).replace("\n","")
            seq_match = match.group(2).replace("\n","").upper()
            seq_match = self.translate_ambiguous(seq_match)
            records[name_match] = seq_match
        #print(records)
        return records
    
    def phylip_parse(self):
    # use regex to parse names and sequences in sequential phylip files
        matches = re.finditer(r"^(\s+)?(\S+)\s+([A-Za-z*?{}-]+)", self.in_file_lines, re.MULTILINE)
        records = {}
        
        for match in matches:
            name_match = match.group(2).replace("\n","")
            seq_match = match.group(3).replace("\n","").upper()
            seq_match = self.translate_ambiguous(seq_match)
            records[name_match] = seq_match
        #print(records)
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
            #print(matrix_match)
            seq_matches = \
             re.finditer(r"^(\s+)?[']?(\S+\s\S+|\S+)[']?\s+([A-Za-z*?{}-]+)($|\s+\[[0-9]+\]$)", \
             matrix_match, re.MULTILINE)
            
            for match in seq_matches:
                name_match = match.group(2).replace("\n","")
                seq_match = match.group(3).replace("\n","").upper()
                seq_match = self.translate_ambiguous(seq_match)
                records[name_match] = seq_match
        #print(records.keys())
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
    
    def __init__(self, parsed_records, aln_name):
    
    # initialize alignment class with parsed records and alignment name as arguments,
    # create empty lists for list of sequences, character matrix, sites without
    # ambiguous or missing characters, and initialize variables for the number
    # of variable sites and parsimony informative sites 
    
        self.parsed_records = parsed_records
        self.aln_name = aln_name
       
        self.list_of_seqs = []
        self.matrix = []
        self.no_missing_ambiguous_sites = []
        self.variable = 0
        self.parsimony_informative = 0
        
    def __str__(self):
        return self.get_name
            
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
        for name, seq in self.parsed_records.items():
            self.list_of_seqs.append(seq)
        return self.list_of_seqs
               
    def matrix_creator(self):
    # decompose character matrix into a two-dimensional list
        new = []
       
        for sequence in self.list_of_seqs:
            for character in sequence:
                new.append(character)
            self.matrix.append(new)
            new = []
        return self.matrix

    def get_column(self, i):
    # get site from the character matrix
        return [row[i] for row in self.matrix]
        
    def all_same(self, items):
    # check if all elements of a list are the same
        return all(x == items[0] for x in items)
        
    def not_missing_ambiguous(self, character):
    # check if character is not missing or ambiguous
        if character not in self.missing_ambiguous_chars:
            return True
            
    def get_sites_no_missing_ambiguous(self):
    # get each site without missing or ambiguous characters  
        for column in range(self.get_alignment_length()):
            site = self.get_column(column)
            site = list(char for char in site if self.not_missing_ambiguous(char))
            self.no_missing_ambiguous_sites.append(site)
        return self.no_missing_ambiguous_sites
        
    def get_variable(self):
    # if all elements of a site without missing or ambiguous characters 
    # are not the same, consider it variable
        for site in self.no_missing_ambiguous_sites:
            if not self.all_same(site):
                #print(site)
                self.variable += 1
        
        return self.variable
    
    def get_parsimony_informative(self):
    # if the count for a unique character in a site is at least two, 
    # and there are at least two such characters in a site without missing
    # or ambiguous characters, consider it parsimony informative 
        for site in self.no_missing_ambiguous_sites:
            pattern = []
            unique_chars = set(site)
            
            for base in unique_chars:
                freq = site.count(base)
                if freq >= 2:
                    pattern.append(base)

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
        return self.aln_name
        
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
        
# print usage instructions if number of agruments given is incorrect
if len(argv) is not 4:
    print(Usage)

# define variables from arguments given
script, in_file, in_format, data_type = argv

# open and parse input file
aln_input = FileParser(in_file)

# parse according to the given format
if in_format == "fasta":
    parsed_aln = aln_input.fasta_parse()
elif in_format == "phylip":
    parsed_aln = aln_input.phylip_parse()
elif in_format == "nexus":
    parsed_aln = aln_input.nexus_parse()
else:
    print(Usage)

# parse according to the given alphabet
if data_type == "aa":
    aln = AminoAcidAlignment(parsed_aln, in_file)
elif data_type == "dna":
    aln = DNAAlignment(parsed_aln, in_file)
else:
    print(Usage)

# get alignment summary
aln.get_summary()
