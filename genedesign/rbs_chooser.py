import re
from typing import List, Dict, Tuple, Set
import pandas as pd
# !pip install biopython
from Bio import SeqIO
from collections import defaultdict
from seq_utils.Translate import Translate

def extract_genes_info(genbank_file):
    gene_dict = defaultdict(dict)  # Dictionary to store gene info
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "gene":
                locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
                gene_name = feature.qualifiers.get("gene", [None])[0]

                # CDS information
                cds_feature = None
                for cds in record.features:
                    if cds.type == "CDS" and cds.qualifiers.get("locus_tag") == [locus_tag]:
                        cds_feature = cds
                        break

                if cds_feature:
                    start, end = cds_feature.location.start, cds_feature.location.end
                    strand = cds_feature.location.strand
                    if strand == 1:  # Forward strand
                        utr_start = max(0, start - 50)
                        utr_seq = record.seq[utr_start:start]
                    else:  # Reverse strand, we need to reverse complement
                        utr_start = end
                        utr_seq = record.seq[utr_start:utr_start + 50].reverse_complement()

                    cds_seq = cds_feature.extract(record.seq)
                    # Save the gene information in the dictionary
                    gene_dict[locus_tag] = {
                        "gene": gene_name,
                        "UTR": utr_seq,
                        "CDS": cds_seq
                    }
    return gene_dict

def prune_proteomics_data(file_path):
    # Read the proteomics data file into a DataFrame
    df = pd.read_csv(file_path, sep='\t', header=None, names=['Locus_Tag', 'Abundance'])

    # Remove any unwanted prefixes from 'Locus_Tag'
    df['Locus_Tag'] = df['Locus_Tag'].str.replace('511145.', '', regex=False)

    # Convert the 'Abundance' column to numeric values, coercing any errors
    df['Abundance'] = pd.to_numeric(df['Abundance'], errors='coerce')

    # Drop rows where 'Abundance' is NaN
    df = df.dropna(subset=['Abundance'])

    # Sort the DataFrame by 'Abundance' in descending order
    df = df.sort_values(by='Abundance', ascending=False)

    # Calculate the number of top 5% entries
    top_5_percent_count = max(1, int(len(df) * 0.05))

    # Get the top 5% entries
    top_5_percent_df = df.head(top_5_percent_count)

    # Convert to list of locus_tag:abundance pairs
    pruned_list = top_5_percent_df.set_index('Locus_Tag')['Abundance'].to_dict()

    return pruned_list

def merge_rbs_cds_with_proteomics(file_path, gb_file_path):
    # Prune the proteomics data to get top 5% entries
    pruned_proteomics = prune_proteomics_data(file_path)

    # Extract gene information from the GenBank file
    genes_info = extract_genes_info(gb_file_path)

    # Create a new dictionary for the merged RBS/CDS data
    merged_data = {}

    # Iterate over pruned proteomics data
    for locus_tag in pruned_proteomics.keys():
        if locus_tag in genes_info:
            # Include the RBS/CDS data for matching locus tags
            merged_data[locus_tag] = genes_info[locus_tag]

    return merged_data

def hairpin_counter(sequence, min_stem=3, min_loop=4, max_loop=9):
    """
    Counts the number of potential hairpin structures in a DNA sequence. Hairpins are common secondary structures
    in nucleic acids where a sequence of nucleotides can fold back on itself to form a double-stranded stem with a single-stranded loop.

    The algorithm searches for regions within the sequence where a segment can base-pair with its reverse complement separated by a loop.
    This function scans for such occurrences by examining every possible substring as a potential stem and ensuring the intervening
    sequence, which would form the loop, meets the specified length requirements.

    Parameters:
        sequence (str): The DNA sequence to analyze.
        min_stem (int): Minimum number of bases in the stem for stable hairpin. The stem must be able to form at least this many
                        complementary base pairs to be considered.
        min_loop (int): Minimum number of bases in the loop. This prevents the formation of overly tight hairpins which may not be biologically relevant.
        max_loop (int): Maximum number of bases in the loop. This constrains the loop to a realistic size, preventing unlikely structures.

    Returns:
        int: The count of potential hairpin structures detected, indicating regions where secondary structure formation might inhibit biological processes like transcription or translation.

    This method does not account for the thermodynamic stability of the predicted hairpins, focusing solely on their potential for formation based on sequence complementarity and specified geometrical constraints.
    """
    count = 0
    seq_len = len(sequence)

    # Iterate through each base to consider it as a start of a stem
    for i in range(seq_len):
        # Only consider end points that would fit within the loop constraints
        for j in range(i + min_stem + min_loop, min(i + min_stem + max_loop + 1, seq_len)):
            stem1 = sequence[i:i+min_stem]
            stem2 = sequence[j:j+min_stem]

            # Check if the stems are complementary (reverse complement match)
            stem2_rc = ''.join(['ATCG'['TAGC'.index(n)] for n in stem2[::-1]])

            if stem1 == stem2_rc:
                count += 1

    return count

def calculate_edit_distance(s1, s2):
    """
    Compute the edit distance between two strings using a dynamic programming approach based on the Smith-Waterman algorithm for local alignment.

    Parameters:
        s1 (str): The first string to compare.
        s2 (str): The second string to compare.

    Returns:
        int: The edit distance between the two strings, defined as the minimum number of edits (insertions, deletions, or substitutions) required to transform one string into the other.
    """
    s1_len = len(s1)
    s2_len = len(s2)
    dist = [[0] * (s2_len + 1) for _ in range(s1_len + 1)]

    # Initialize distances for transformations involving empty strings
    for i in range(s1_len + 1):
        dist[i][0] = i
    for j in range(s2_len + 1):
        dist[0][j] = j

    # Compute distances
    for i in range(1, s1_len + 1):
        for j in range(1, s2_len + 1):
            if s1[i - 1] == s2[j - 1]:
                dist[i][j] = dist[i - 1][j - 1]
            else:
                dist[i][j] = 1 + min(dist[i - 1][j], dist[i][j - 1], dist[i - 1][j - 1])

    return dist[s1_len][s2_len]




##########################################################################################################################################


from genedesign.models.rbs_option import RBSOption
import re
from typing import List, Dict, Tuple, Set

class RBSOption:
    """
    A class representing an RBS option, including its UTR, gene name, CDS, and the first six amino acids.

    Attributes:
        utr (str): The 5' UTR sequence.
        gene_name (str): The name of the gene associated with this RBS option.
        cds (str): The coding sequence.
        first_six_aas (str): The first six amino acids of the associated protein.
    """
    def __init__(self, utr: str, gene_name: str, cds: str, first_six_aas: str):
        self.utr = utr
        self.gene_name = gene_name
        self.cds = cds
        self.first_six_aas = first_six_aas

    def __repr__(self):
        return f"RBSOption(utr='{self.utr}', gene_name='{self.gene_name}', cds='{self.cds}', first_six_aas='{self.first_six_aas}')"

class RBSChooser:
    """
    A class to choose the best RBS for a given CDS sequence.
    """

    def __init__(self):
        self.rbs_options: List[RBSOption] = []  # Class variable to hold RBS options

    def shine_dalgarno_rbs_finder(self, cds: str) -> str:
        """
        Finds the RBS sequence based on the Shine-Dalgarno sequence.

        Parameters:
            cds (str): The coding sequence from which to find the RBS.

        Returns:
            str: The RBS sequence, or an empty string if no RBS is found.
        """
        # Define common Shine-Dalgarno sequences (can be modified as needed)
        sd_sequences = ["AGGAGG", "GGA", "GGAG", "AGG"]  # Include common variants

        # Search for the Shine-Dalgarno sequence in the CDS
        for sd_seq in sd_sequences:
            # Look for the RBS in the CDS, allowing for some variability
            pattern = re.compile(f"(?={sd_seq})")  # Look for the position of the RBS
            match = pattern.search(cds)
            if match:
                # Return the found RBS sequence
                return cds[match.start():match.start() + len(sd_seq)]
        return ""  # No RBS found

    def initiate(self, file_path: str, gb_file_path: str) -> None:
        """
        Initialization method for RBSChooser.
        This method constructs RBSOption instances and populates the rbs_options list.

        Parameters:
            file_path (str): The path to the proteomics data file.
            gb_file_path (str): The path to the GenBank file.
        """
        # Merge RBS/CDS data with proteomics data
        merged_data = merge_rbs_cds_with_proteomics(file_path, gb_file_path)

        # Construct RBS options from merged data
        for locus_tag, rbs_cds in merged_data.items():
            rbs_sequence = self.shine_dalgarno_rbs_finder(str(rbs_cds['CDS']))  # Extract the CDS from rbs_cds
            if rbs_sequence:  # Only proceed if an RBS was found
                first_six_aas = self.get_first_six_aas(rbs_cds['CDS'])  # Assuming you have a function for this
                rbs_option = RBSOption(utr=rbs_sequence, gene_name=locus_tag, cds=rbs_cds['CDS'], first_six_aas=first_six_aas)
                self.rbs_options.append(rbs_option)

        # Print out the number of RBS options created
        print(f"Number of RBS options created: {len(self.rbs_options)}")
        # for option in self.rbs_options[:4]:  # Print the first 4 RBS options as examples
        #     print(option)

    def get_first_six_aas(self, cds: str) -> str:
        """
        Translates the first 18 bases of the CDS to get the first 6 amino acids.

        Parameters:
            cds (str): The coding sequence.

        Returns:
            str: The first 6 amino acids.
        """
        translator = Translate()  # Assuming you have a Translate class defined elsewhere
        translator.initiate()
        return translator.run(cds[:18])  # Translate first 18 bases

    def run(self, cds: str, ignores: Set[RBSOption]) -> RBSOption:
        """
        Executes the RBS selection process for the given CDS.

        Parameters:
        - cds (str): The coding sequence to pair with an RBS.
        - ignores (Set[RBSOption]): A set of RBSOption instances to ignore during selection.

        Returns:
        - RBSOption: The selected RBSOption that best pairs with the given CDS.
        """
        # Step 1: Exclude RBS from the 'ignores' list
        valid_options = [option for option in self.rbs_options if option not in ignores]

        # Step 2: Compare peptides
        input_peptide = self.get_first_six_aas(cds)
        peptide_similarity_scores = {}

        for option in valid_options:
            edit_distance = calculate_edit_distance(input_peptide, option.first_six_aas)
            peptide_similarity_scores[option] = edit_distance

        # Print out the peptide similarity scores
        # for option, score in peptide_similarity_scores.items():
        #     print(f"RBS Option: {option}, Edit Distance: {score}")

        # Step 3: Evaluate potential hairpin structures
        hairpin_scores = {}

        for option in valid_options:
            combined_sequence = option.utr + cds  # Combine RBS and CDS sequences
            hairpin_count = hairpin_counter(combined_sequence)  # Count potential hairpin structures
            hairpin_scores[option] = hairpin_count  # Store the hairpin score for each option

        # Print out hairpin scores
        # for option, count in hairpin_scores.items():
        #     print(f"RBS Option: {option}, Hairpin Count: {count}")

        # Step 4: Choose the best RBS option based on criteria
        best_option = None
        best_score = float('inf')  # Initialize best score as infinite for comparison

        for option in valid_options:
            hairpin_count = hairpin_scores[option]
            peptide_score = peptide_similarity_scores[option]

            # Define a combined score based on your criteria
            combined_score = hairpin_count + peptide_score  # Adjust this formula as needed

            if combined_score < best_score:
                best_score = combined_score
                best_option = option

        return best_option  # Return the best RBS option found
    



# # Instantiate RBSChooser and initiate
# chooser = RBSChooser()
# file_path = '/data/511145-WHOLE_ORGANISM-integrated.txt'  # Path to the proteomics data file
# gb_file_path = '/data/sequence.gb'  # Path to the GenBank file
# chooser.initiate(file_path, gb_file_path)

# # Test CDS
# cds = "ATGGCTAGCAAATACGATTTTACAATATAA"

# # Initialize empty ignores set
# ignores = set()

# # First run, no ignored RBS options
# selected_rbs_1 = chooser.run(cds, ignores)
# print(f"Selected RBS: {selected_rbs_1}")

# # Repeat to confirm determinacy
# selected_rbs_2 = chooser.run(cds, ignores)
# print(f"Selected RBS: {selected_rbs_2}")

# # Add the returned RBSOption to ignores and run again, confirm a different selected RBS
# ignores.add(selected_rbs_1)
# selected_rbs_3 = chooser.run(cds, ignores)
# print(f"Selected RBS after ignoring {selected_rbs_1}: {selected_rbs_3}")



# from dataclasses import dataclass

# @dataclass
# class Translate:
#     """
#     Translates a DNA sequence into a protein sequence using the standard genetic code, halting at the first stop codon encountered and throwing an error for invalid codons.

#     Attributes:
#         codon_table (dict): Maps each DNA codon to its corresponding single-letter amino acid code.
#     """
#     codon_table: dict = None

#     def initiate(self) -> None:
#         """
#         Initializes the codon table with the genetic code for translating nucleotide triplets into amino acids.
#         """
#         self.codon_table = {
#             "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
#             "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
#             "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
#             "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
#             "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
#             "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
#             "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
#             "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
#             "TAT": "Y", "TAC": "Y", "TAA": "Stop", "TAG": "Stop",
#             "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
#             "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
#             "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
#             "TGT": "C", "TGC": "C", "TGA": "Stop", "TGG": "W",
#             "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
#             "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
#             "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
#         }

#     def run(self, dna_sequence: str) -> str:
#         """
#         Translates a DNA sequence into a protein sequence using the codon table.

#         Parameters:
#             dna_sequence (str): The DNA sequence to translate.

#         Returns:
#             str: The corresponding amino acid sequence.

#         Raises:
#             ValueError: If the DNA sequence length is not a multiple of three, contains untranslated sequence after a stop codon, or contains invalid codons.
#         """
#         if len(dna_sequence) % 3 != 0:
#             raise ValueError("The DNA sequence length must be a multiple of 3.")

#         protein = []
#         for i in range(0, len(dna_sequence), 3):
#             codon = dna_sequence[i:i+3]
#             if codon not in self.codon_table:
#                 raise ValueError(f"Invalid codon '{codon}' encountered in DNA sequence.")
#             amino_acid = self.codon_table[codon]
#             if amino_acid == "Stop":
#                 if i + 3 != len(dna_sequence):
#                     raise ValueError("Untranslated sequence after stop codon.")
#                 break
#             protein.append(amino_acid)

#         return ''.join(protein)