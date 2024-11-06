from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
import random
import pandas as pd
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.gc_content_checker import GCContentChecker

proteomics_file_path = "/Users/josephjoo/Downloads/bioe134/bioe134-234-transcriptdesigner-project-3-josephhjoo/genedesign/data/511145-WHOLE_ORGANISM-integrated (1).txt"
genome_file_path = "/Users/josephjoo/Downloads/bioe134/bioe134-234-transcriptdesigner-project-3-josephhjoo/genedesign/data/sequence (1).gb"
codon_usage_file_path = "/Users/josephjoo/Downloads/bioe134/bioe134-234-transcriptdesigner-project-3-josephhjoo/genedesign/data/codon_usage.txt"
codon_usage_file = pd.read_csv(codon_usage_file_path, sep='\t', names=["Codon", "Amino_Acid", "Fraction", "Frequency"])

class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence and chooses an RBS using optimized codons.
    """

    def __init__(self):
        self.aminoAcidToCodon = {}
        self.rbsChooser = None
        self.codon_usage = {}
        self.codon_to_aa = {}

    def initiate(self) -> None:
        """
        Initializes the codon table and the RBS chooser.
        """
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate(proteomics_file_path, genome_file_path)

        # Codon table with highest CAI codon for each amino acid (for E. coli)
        self.aminoAcidToCodon = {
            'A': "GCG", 'C': "TGC", 'D': "GAT", 'E': "GAA", 'F': "TTC",
            'G': "GGT", 'H': "CAC", 'I': "ATC", 'K': "AAA", 'L': "CTG",
            'M': "ATG", 'N': "AAC", 'P': "CCG", 'Q': "CAG", 'R': "CGT",
            'S': "TCT", 'T': "ACC", 'V': "GTT", 'W': "TGG", 'Y': "TAC"
        }

        self.codon_to_aa = {codon: aa for aa, codon in self.aminoAcidToCodon.items()}

        # Populating codon_usage dict
        for _, row in codon_usage_file.iterrows():
            codon = row.iloc[0]
            amino_acid = row.iloc[1]
            fraction = row.iloc[2]

            if amino_acid not in self.codon_usage:
                self.codon_usage[amino_acid] = {}
            self.codon_usage[amino_acid][codon] = fraction
        
        # print(self.codon_usage)

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Translates the peptide sequence to DNA with Monte Carlo-based codon optimization and selects an RBS.
        
        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.
        
        Returns:
            Transcript: The transcript object with the selected RBS and optimized codon sequence.
        """
        # Translate peptide to codons
        codons = [self.aminoAcidToCodon[aa] for aa in peptide]
        codons_string = "".join(codons)
        # print("Initial codons:", codons_string)

        # Define maximum attempts and initialize attempt counter
        max_attempts = 100
        best_codons = codons_string
        max_checks_passed = 0

        # Optimize the first three codons separately (EDGE CASE)
        attempt_count = 0
        while attempt_count < max_attempts:
            optimized_first = self.monte_carlo_codon_choice(codons_string[:9], 0)
            # print(f"Attempt {attempt_count + 1}: Optimized first codons: {optimized_first}")

            if self.apply_checkers(optimized_first):
                # print(f"Checks passed for first codons: {self.apply_checkers(optimized_first)}")
                codons_string = optimized_first + codons_string[9:]
                break
            attempt_count += 1

        # Main sliding window optimization for the middle codons
        window_size = 12
        for start_index in range(3, len(codons_string) - window_size - 6, 3):
            attempt_count = 0
            window = codons_string[start_index:start_index + window_size]

            while attempt_count < max_attempts:
                optimized_window = self.monte_carlo_codon_choice(window, start_index)
                check_results = self.apply_checkers(optimized_window)

                if all(check_results):
                    # print(f"Middle window optimized: {optimized_window}")
                    codons_string = (
                        codons_string[:start_index] +
                        optimized_window +
                        codons_string[start_index + window_size:]
                    )
                    break
                else:
                    num_checks_passed = sum(check_results)
                    if num_checks_passed > max_checks_passed:
                        best_codons = (
                            codons_string[:start_index] +
                            optimized_window +
                            codons_string[start_index + window_size:]
                        )
                        max_checks_passed = num_checks_passed

                attempt_count += 1

        # Optimize the last six codons separately (EDGE CASE)
        start_index = len(codons_string) - 18
        attempt_count = 0
        while attempt_count < max_attempts:
            optimized_last = self.monte_carlo_codon_choice(codons_string[start_index:], start_index)
            if self.apply_checkers(optimized_last):
                codons_string = codons_string[:start_index] + optimized_last
                break
            attempt_count += 1

        # If no valid codon string passes all checks, use the best found so far
        optimized_codons = best_codons if attempt_count >= max_attempts else codons_string

        # Append the stop codon
        optimized_codons += "TAA"

        # Choose an RBS
        selectedRBS = self.rbsChooser.run(optimized_codons, ignores)

        # Return the Transcript object
        return Transcript(selectedRBS, peptide, optimized_codons)


    def monte_carlo_codon_choice(self, codons_string, start_index):
        """
        Selects optimized codons for positions 3-6 in the given codon sequence window
        based on usage frequency using Monte Carlo sampling.

        Parameters:
            codons_string (str): The original codon sequence string.
            start_index (int): The starting index of the current window in codons_string.

        Returns:
            str: Updated codons_string with the optimized codons in positions 3-6 of the window.
        """
        new_codons_list = []
        for i in range(2, 6):  # Corresponds to 3rd to 6th codon (index 2 to 5)
            codon = codons_string[start_index + i * 3:start_index + (i + 1) * 3]
            # print("Original codon:", codon)

            # If the codon is not recognized, skip to the next one
            if codon not in self.codon_to_aa:
                new_codons_list.append(codon)
                continue
            
            amino_acid = self.codon_to_amino_acid(codon)
            relevant_dict = self.codon_usage.get(amino_acid, {})
            
            # Monte Carlo sampling for the new codon choice
            pick_list = [key for key, value in relevant_dict.items() for _ in range(int(value * 100))]
            selected_codon = random.choice(pick_list) if pick_list else codon  # fallback if pick_list is empty
            
            new_codons_list.append(selected_codon)

        # Combine the original sequence with the new optimized codons
        optimized_window = (
            codons_string[:start_index + 6] +  # Include codons up to the second position being optimized
            "".join(new_codons_list) +         # Optimized codons (3rd to 6th)
            codons_string[start_index + 18:]    # Remainder of the original codons after the optimized section
        )

        return optimized_window

    def apply_checkers(self, codons: str) -> tuple:
        forbidden_sequence = ForbiddenSequenceChecker()
        forb_seq_bool = forbidden_sequence.run(codons)[0]
        hair_check_bool = hairpin_checker(codons)[0]
        internal_promoter = PromoterChecker()
        internal_promoter.initiate()
        int_prom, _ = internal_promoter.run(codons)

        gc_checker = GCContentChecker() # Checks if the GC content of a codon sequence is within a specified range
        gc_check_bool = gc_checker.run(codons)

        return forb_seq_bool, hair_check_bool, int_prom, gc_check_bool

    def codon_to_amino_acid(self, codon: str) -> str:
        """
        Maps a given codon string to its corresponding amino acid.

        Parameters:
            codon (str): The codon string to convert.

        Returns:
            str: The corresponding amino acid, or an empty string if not found.
        """
        return self.codon_to_aa.get(codon, "")
    
    def get_codons_from_peptide(peptide):
        # Example codon table
        codon_table = {
            'A': "GCG", 'C': "TGC", 'D': "GAT", 'E': "GAA", 'F': "TTC",
            'G': "GGT", 'H': "CAC", 'I': "ATC", 'K': "AAA", 'L': "CTG",
            'M': "ATG", 'N': "AAC", 'P': "CCG", 'Q': "CAG", 'R': "CGT",
            'S': "TCT", 'T': "ACC", 'V': "GTT", 'W': "TGG", 'Y': "TAC"
        }
        
        codons = []
        for aa in peptide:
            if aa in codon_table:
                # Choose a random codon for the amino acid
                selected_codon = random.choice(codon_table[aa])
                codons.append(selected_codon)
                # print(f"Original codon: {selected_codon}")  # Debugging line
            else:
                None
                # print(f"Amino acid {aa} not found in codon table.")  # Handle unknown AA

        return codons


    
    # Additional checkers (forbidden_sequence_checker, internal_promoter_checker, secondary_structure_checker)
    # would be defined here as separate methods.

    class GCContentChecker:
        """
        Checks if the GC content of a codon sequence is within a specified range.
        """

        def __init__(self, min_gc=0.4, max_gc=0.6):
            self.min_gc = min_gc
            self.max_gc = max_gc

        def run(self, codons: str) -> bool:
            """
            Calculates the GC content and checks if it is within the acceptable range.

            Parameters:
                codons (str): The codon sequence as a string.

            Returns:
                bool: True if GC content is within the range, False otherwise.
            """
            gc_count = codons.count("G") + codons.count("C")
            gc_content = gc_count / len(codons)
            return self.min_gc <= gc_content <= self.max_gc


if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    
    peptide = "MYPFIRTARMTV"
    peptide2 = "MSKGEEKLMNPKKKTRVSASASTVWTWVASDFGASDFGPQRVWCALLMAGS"
    peptide3 = "MKTWQRPLVFHYCGANDSIR"
    
    designer = TranscriptDesigner()
    designer.initiate()  # Pass the file paths here

    ignores = set()
    transcript = designer.run(peptide2, ignores)
    
    # Print out the transcript information
    print(transcript)