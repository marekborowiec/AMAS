#! /usr/bin/env python3

"""
AMAS_v0.1.py by Marek Borowiec

This stand-alone program calculates various statistics on a multiple sequence
alignment. It supports sequential FASTA, PHYLIP, and NEXUS formats for DNA and 
amino acid sequences.

Current statistics include the number of taxa, alignment length, total number
of matrix cells, overall number of undetermined characters, percent of missing 
data, AT and GC contents (for DNA alignments), and relative proportions of all
characters relative to matrix size.
"""

Usage = """
Usage: AMAS_v0.1.py <input_file> <format> <alphabet>

Supported formats: "fasta", "phylip", "nexus"
Supported alphabets: "aa", "dna"
"""

from sys import argv
import re


class FileParser:
    """Parse file contents and return sequences and sequence names"""

    def __init__(self, file_name):
        self.file_name = file_name
        self.in_file = open(file_name, "r")
        self.in_file_lines = self.in_file.read().rstrip("\r\n")
    
    def get_file_name(self):
        return self.file_name

    def fasta_parse(self):
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
        matches = re.finditer(r"(\s+)?(MATRIX\n|matrix\n|MATRIX\r\n|matrix\r\n)(.*?;)", \
        self.in_file_lines, re.DOTALL)
        records = {}
        
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
       self.parsed_records = parsed_records
       self.aln_name = aln_name
       
       self.list_of_seqs = []
        
    def __str__(self):
        return self.get_name
            
    def summarize_alignment(self):
        summary = []
        self.seq_grabber()
        name = str(self.get_name())
        taxa_no = str(self.get_taxa_no())
        length = str(self.get_alignment_length())
        cells = str(self.get_matrix_cells())
        missing = str(self.get_missing())
        missing_percent = str(self.get_missing_percent())
        summary = [name, taxa_no, length, cells, missing, missing_percent]
        return summary

    def get_freq_summary(self):
        characters = []
        frequencies = []
        
        for item in self.get_frequencies():
            for char, freq in item.items():
                characters.append(str(char))
                frequencies.append(str(freq))
        return characters, frequencies           
        
    def matrix_creator(self):	           
       self.matrix = []
       self.new = []
       
       for sequence in self.list_of_seqs:
           for character in sequence:
               self.new.append(character)
           self.matrix.append(self.new)
       self.new = []
       return self.matrix

    def seq_grabber(self):
       for name, seq in self.parsed_records.items():
          self.list_of_seqs.append(seq)
       return self.list_of_seqs
       
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
        frequencies = []
        
        for char in self.alphabet:
            count = sum(seq.count(char) for seq in self.list_of_seqs) \
             / self.all_matrix_cells
            frequencies.append({char : round(count, 3)})
        return frequencies
            

class AminoAcidAlignment(Alignment):
    """Summary specific to aa alignments"""

    alphabet = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R", \
     "S","T","V","W","Y","B","Z","X","*","-","?"]
    missing_chars = ["X","*","-","?"]

    def get_summary(self):
        data = self.summarize_alignment()
        new_data = data + list(self.get_freq_summary()[1])
        
        header = ["Alignment_name", "No_of_taxa", "Alignment_length", \
        "Total_matrix_cells", "Undetermined_characters", "Missing_percent"]
        
        freq_header = list(self.get_freq_summary()[0])
        new_header = header + freq_header
        print("\t".join(new_header))
        print("\t".join(new_data))

           
class DNAAlignment(Alignment):
    """Summary specific to DNA alignments"""
    
    alphabet = ["A","C","G","T","K","M","R","Y","S","W","B","V","H","D","X", \
     "N", "O", "-","?"]
    missing_chars = ["X","N","O","-","?"]

    def get_summary(self):
        data = self.summarize_alignment()
        
        new_data = data + self.get_atgc_content() \
         + list(self.get_freq_summary()[1])
        header = ["Alignment_name", "No_of_taxa", "Alignment_length", \
         "Total_matrix_cells", "Undetermined_characters", "Missing_percent", \
          "AT_content", "GC_content"]
        
        freq_header = list(self.get_freq_summary()[0])
        new_header = header + freq_header
        
        print("\t".join(new_header))
        print("\t".join(new_data))
        
    def get_atgc_content(self):
        atgc_content = []
        
        at_count = sum((seq.count("A") + seq.count("T") + seq.count("W")) \
         for seq in self.list_of_seqs)
        gc_count = sum((seq.count("G") + seq.count("C") + seq.count("S")) \
         for seq in self.list_of_seqs)
        
        at_content = str(round(at_count / (at_count + gc_count), 3))
        gc_content = str(round(1 - float(at_content), 3))
        
        atgc_content.extend((at_content, gc_content))
        return atgc_content
        

if len(argv) is not 4:
    print(Usage)

script, in_file, in_format, data_type = argv

aln_input = FileParser(in_file)

if in_format == "fasta":
    parsed_aln = aln_input.fasta_parse()
elif in_format == "phylip":
    parsed_aln = aln_input.phylip_parse()
elif in_format == "nexus":
    parsed_aln = aln_input.nexus_parse()
else:
    print(Usage)

if data_type == "aa":
    aln = AminoAcidAlignment(parsed_aln, in_file)
elif data_type == "dna":
    aln = DNAAlignment(parsed_aln, in_file)
else:
    print(Usage)
   
aln.get_summary()
