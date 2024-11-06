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