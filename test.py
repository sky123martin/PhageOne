from application.utility import *
import pandas as pd
import unittest
import requests
import random

class BRED(unittest.TestCase):
    def test_download_phage_fasta(self):
        pass

    def test_fasta_to_DNA(self):
        pass

    def test_collect_gene_info(self):
        pass

    def test_find_primers(self):
        pass

    def find_amplicon(self):
        amplicon = test_find_amplicon("ATGGAGGGGCGATG", "TGG", "ATC")
        self.assertEqual(amplicon, "AGGGGC")
        
        # test exception throwing when unable to find primers
        with self.assertRaises(PrimerNotFound):
            test_find_amplicon("AAAAAAAAAA", "TGG", "ATC")
        
    def test_find_editing_substrate(self):
        # test insertion
        substrate, edited_DNA = find_editing_substrate("abcdefghi", 4, 4, new_DNA = "xyz", homo_length = 2)
        self.assertEqual(substrate, "cdxyzef")
        self.assertEqual(edited_DNA, "abcdxyzefghi")

        # test swap
        substrate, edited_DNA = find_editing_substrate("abcdefghi", 4, 6, new_DNA = "xyz", homo_length = 2)
        self.assertEqual(substrate, "bcxyzgh")
        self.assertEqual(edited_DNA, "abcxyzghi")

        # test deletion
        substrate, edited_DNA = find_editing_substrate("abcdefghi", 4, 6, new_DNA = "", homo_length = 2)
        self.assertEqual(substrate, "bcgh")
        self.assertEqual(edited_DNA, "abcghi")


if __name__ == '__main__':
    unittest.main()