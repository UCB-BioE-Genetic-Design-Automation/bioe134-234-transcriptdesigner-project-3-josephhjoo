import pytest
from genedesign.checkers.gc_content_checker import GCContentChecker  # Adjust the import based on your file structure

@pytest.fixture
def gc_checker():
    return GCContentChecker(min_gc=0.4, max_gc=0.6)

def test_gc_content_below_range(gc_checker):
    # Sequence with GC content below 40%
    low_gc_sequence = "ATATATAT"
    assert not gc_checker.run(low_gc_sequence), "Expected GC content to be below range"

def test_gc_content_within_range(gc_checker):
    # Sequence with GC content within 40-60%
    within_gc_sequence = "GCGCATAT"
    assert gc_checker.run(within_gc_sequence), "Expected GC content to be within range"

def test_gc_content_above_range(gc_checker):
    # Sequence with GC content above 60%
    high_gc_sequence = "GCGCGCGC"
    assert not gc_checker.run(high_gc_sequence), "Expected GC content to be above range"

def test_gc_content_exact_min(gc_checker):
    # Sequence with exact minimum GC content (40%)
    exact_min_gc_sequence = "GCGTATAT"
    assert gc_checker.run(exact_min_gc_sequence), "Expected GC content to be exactly at min range"

def test_gc_content_exact_max(gc_checker):
    # Sequence with exact maximum GC content (60%)
    exact_max_gc_sequence = "GCGCGTAT"
    assert gc_checker.run(exact_max_gc_sequence), "Expected GC content to be exactly at max range"

