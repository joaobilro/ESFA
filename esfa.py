# !/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# esfa.py
#
# Version beta1.0
#
# This Python 3 script allows the user to extract codon sites from alignments. This was
# created mainly to be coupled with CodeML, as a way to analyse site-by-site results in
# a user-friendly manner. In order for it to properly work, it needs to be given a list
# as an input, with all of the genes and sites that the user wants to analyse. The list
# should follow the same structure as the following:
# 
# list.txt
# 
# gene_a: 1, 2, 3, 4, 5, 6
# gene_b: 20, 21, 22, 23, 24
# gene_c: 123, 3, 1245, 235, 5
# ...
# 
# The codon sites should be separated by a comma and a blank space, and they do not need
# to be in order. Make sure that the gene name (for example, gene_a) corresponds with the
# name of the alignment file. The extension should not be included, but check if it is 
# supported by the program. Also make sure that LF line breaks were used to separate the
# lines in the list.
#           _
#         ><_> 
#     
#        
# MIT License
# 
# Copyright (c) 2024 João Bilro (joaobilro)
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import argparse
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os


esfa = argparse.ArgumentParser(description="This Python 3 script allows the user to extract codon"
                               "sites from alignments. The files should be in the .FASTA or .PHYLIP format."
                               "For more help, please refer to the description of the arguments.")

esfa.add_argument("--list", "-l", dest="input_gene_list", required=True, type=str, help="The file path to the list, in .txt format")

esfa.add_argument("--dir", "-d", dest="genes_dir", required=True, type=str, help="The path to the directory containing the alignment files, in .fasta or .phy format.")

args = esfa.parse_args()

class Extraction:
    """Contains the necessary functions to extract data from each of the desired alignments.""" 

    def __init__(self, gene_list, main_directory): 
        self.gene_list = gene_list
        self.main_directory = main_directory
        self.gene_sites = self.parse_gene_list()
        self.all_species = self.get_species()

    def parse_gene_list(self):
        """Parses the gene list .txt file to find which of the alignments (and codon sites, correspondingly) it is supposed to parse through."""
        gene_sites = {}

        with open(self.gene_list, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue   ### Skips empty lines or comments
            
                ### Get alignment names and codon sites
                gene_name, sites = line.split(":")
                gene_name = gene_name.strip()
                sites = [int(s.strip()) for s in sites.split(",")]
                gene_sites[gene_name] = sites
        
        return gene_sites
    
    def get_alignments(self):
        """Extracts the path to the desired alignments."""
        alignment_files = []

        os.chdir(self.main_directory)
        for path, subdirs, files in os.walk(self.main_directory):
            for file in files:
                ### List only the FASTA and PHYLIP files
                if file.endswith(".fasta") or file.endswith(".phy"):
                    alignment_files.append(os.path.join(path, file))
        
        genes_of_interest = []
        for file_path in alignment_files:
            for gene_name in self.gene_sites.keys():
                if gene_name in os.path.basename(file_path):
                    genes_of_interest.append(file_path)
        
        return genes_of_interest
    
    def get_species(self):
        """Gets all of the species found in the desired alignments, so that columns can be correctly populated."""
        all_species = set()

        for file_path in self.get_alignments():
            alignment = AlignIO.read(file_path, "fasta" if file_path.endswith(".fasta") else "phylip")
            species_names = {record.id for record in alignment}
            all_species.update(species_names)

        return all_species
    
    def get_sites(self, output_fasta = "SitesPerGene.fasta", concatenated_fasta = "ConcatenatedSites.fasta"):
        """Extracts the corresponding codon sites from the desired alignments."""

        species_sequences = {species: [] for species in self.all_species}

        ### Get gene info
        gene_data = {} 
        for file_path in self.get_alignments():
            gene_name = next(gene for gene in self.gene_sites if gene in os.path.basename(file_path))
            sites = self.gene_sites[gene_name]

            print(f"Processing sites for gene {gene_name}...")

            alignment = AlignIO.read(file_path, "fasta" if file_path.endswith(".fasta") else "phylip")
            alignment_dict = {record.id: str(record.seq) for record in alignment}

            ### For each species, make sure there is a sequence
            for species in self.all_species:
                if species not in alignment_dict:
                    alignment_dict[species] = "-" * len(alignment[0].seq)   ### Fill missing taxa columns with gaps
            
            ### Store site info
            gene_data[gene_name] = alignment_dict

        for species in self.all_species:
            concatenated_seq = []
            for gene_name, alignment_dict in gene_data.items():
                seq = alignment_dict.get(species, "-" * len(alignment_dict[next(iter(alignment_dict))]))    ### Fill missing taxa columns with gaps
                extracted_seq = ["-"] * (len(self.gene_sites[gene_name]) * 3)
                for idx, site in enumerate(self.gene_sites[gene_name]):
                    codon_start = (site - 1) * 3
                    codon = seq[codon_start:codon_start + 3]
                    extracted_seq[idx * 3: (idx + 1) * 3] = str(codon)
                concatenated_seq.append("".join(extracted_seq)) 

            ### Write concat FASTA
            with open(concatenated_fasta, "a") as concat_fasta:
                concat_record = SeqRecord(Seq("".join(concatenated_seq)),
                                          id = species,
                                          description = "Concatenated extracted sites")
                SeqIO.write([concat_record], concat_fasta, "fasta")
        
        ### Write intermediate FASTA
        with open(output_fasta, "w") as out_fasta:
            for gene_name, alignment_dict in gene_data.items():
                for species, seq in alignment_dict.items():
                    extracted_seq = ["-"] * (len(self.gene_sites[gene_name]) * 3)
                    for idx, site in enumerate(self.gene_sites[gene_name]):
                        codon_start = (site - 1) * 3
                        codon = seq[codon_start:codon_start + 3]
                        extracted_seq[idx * 3: (idx + 1) * 3] = str(codon)

                        ### Create new SeqRecord
                        new_record = SeqRecord(Seq("".join(extracted_seq)),
                                               id = f"{species}_{gene_name}",
                                               description = f"Extracted sites for {gene_name} from {os.path.basename(file_path)}")
                        SeqIO.write([new_record], out_fasta, "fasta")

if __name__ == "__main__":
    ### Extraction
    extraction = Extraction(args.input_gene_list, args.genes_dir)
    extraction.get_sites()
    print(f"Extraction complete. Data saved to SitesPerGene.fasta and ConcatenatedSites.fasta")