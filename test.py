from application.utility import *
import pandas as pd
import unittest
import requests
import random
import subprocess
import string
from random import choice

class BRED(unittest.TestCase):
    def setUp(self):
        #proc = subprocess.check_call("rm -R data")
        proc = subprocess.check_call("mkdir -p data/genes_by_phage", shell=True)
        proc = subprocess.check_call("mkdir -p data/fasta_files", shell=True)

    def tearDown(self):
        #proc = subprocess.check_call("rm -R data")
        pass

    def test_fasta_to_DNA(self):
        DNA = fasta_to_DNA("D29")

    def test_download_phage_fasta(self):
        download_phage_fasta("D29")

        out = download_phage_fasta("not a real phage name")
        self.assertEqual(out, "unable to find phage")

    def test_download_phage_genes(self):
        df = download_phage_genes("D29")
        self.assertEqual(df["gene number"].count(), 77)

    def test_find_amplicon(self):
        amplicon = find_amplicon("ATGGAGGGGCGATG", "TGG", "ATC")
        self.assertEqual(amplicon, "AGGGGC")
        
        # test exception throwing when unable to find primers
        with self.assertRaises(PrimerNotFound):
            find_amplicon("AAAAAAAAAA", "TGG", "ATC")
        
    def test_find_editing_substrate(self):
        # test insertion
        substrate, edited_DNA = find_editing_substrate("abcdefghi", 4, 4, template_DNA = "xyz", homo_length = 2)
        self.assertEqual(substrate, "cdxyzef")
        self.assertEqual(edited_DNA, "abcdxyzefghi")

        # test swap
        substrate, edited_DNA = find_editing_substrate("abcdefghi", 4, 6, template_DNA = "xyz", homo_length = 2)
        self.assertEqual(substrate, "bcxyzgh")
        self.assertEqual(edited_DNA, "abcxyzghi")

        # test deletion
        substrate, edited_DNA = find_editing_substrate("abcdefghi", 4, 6, template_DNA = "", homo_length = 2)
        self.assertEqual(substrate, "bcgh")
        self.assertEqual(edited_DNA, "abcghi")

    def test_all_possible_subsets(self):
        seq = "abcde"
        primer_length = 3
        out = all_possible_subsets(seq, primer_length, "c")
        expected = pd.DataFrame(np.array([["abc","c"],["bcd","c"],["cde","c"]]), columns=['primer_seq', 'primer'])
        self.assertEqual(True, expected.equals(out))
    
    def test_primer3_calculate_tm(self):
        seq = 'GTAAAACGACGGCCAGT'
        tm = primer3_calculate_tm(seq)
        self.assertEqual(tm, 49.16808228911765)

    def test_primer3_calculate_hairpin(self):
        # test initial calc
        seq = 'CCCCCATCCGATCAGGGGG'
        hp = primer3_calculate_hairpin(seq)
        self.assertEqual(hp, -2770.524493437486)

        # test generated hairpin
        seq = 'CCCCCATCCGATCAGGGGG'
        hairpin_seq = seq + reverse(complement(seq))
        hp = primer3_calculate_hairpin(hairpin_seq)
        self.assertEqual(hp, -17242.790467175044)

        # test generated hairpin
        seq = "AGTAGAGTCTCGTGGACGGT" #"".join([choice("CGTA") for _ in range(20)])
        hp = primer3_calculate_hairpin(seq+seq)
        self.assertEqual(hp, -1210.7754890635188)

    def test_primer3_calculate_heterodimer(self):
        seq1 = "TTATACCCCCGGAATCGGCTACGGGCCATACAGG"
        seq2 = "CCCGGTATCTGTCAAGGTGATCTACGC"
        dh_diff = primer3_calculate_heterodimer(seq1, seq2)

        seq1 = "TTATACCCCCGGAATCGGCTACGGGCCATACAGG"
        seq2 = reverse(complement(seq1[10:20]))
        dh_same_subset = primer3_calculate_heterodimer(seq1, seq2)

        seq1 = "TTATACCCCCGGAATCGGCTACGGGCCATACAGG"
        seq2 = reverse(complement(seq1))
        dh_same = primer3_calculate_heterodimer(seq1, seq2)

        self.assertEqual(dh_diff, -5038.009493440593)
        self.assertEqual(dh_same_subset, -9512.098480306246)
        self.assertEqual(dh_same, -35428.989427781344)
        self.assertEqual(dh_diff>dh_same, True)
        self.assertEqual(dh_same<dh_same_subset, True)

    def test_find_possible_primers(self):
        # test replacement
        DNA = "abcdefghij"+"klmnop"+"qstuvwxyz" # 10 bp on each side
        edited_DNA = "abcdefghij"+"1234567"+"qstuvwxyz"
        bp_position_start = 11
        bp_position_stop = 16
        edit_type = "replacement"
        template_DNA = "1234567"
        primer_length = 3
        search_size = 6
        out = find_possible_primers(DNA, edited_DNA, bp_position_start, bp_position_stop, edit_type, template_DNA, primer_length, search_size, True)
        f = all_possible_subsets("bcdefg", primer_length, "f")
        r = all_possible_subsets("tuvwxy", primer_length, "r")
        c = all_possible_subsets("jklmno", primer_length, "c")
        e = all_possible_subsets("j12345", primer_length, "e")
        expected = pd.concat([f, r, c, e])
        self.assertEqual(True, expected.equals(out))

        # test insertion
        DNA = "abcdefghij"+"qstuvwxyz" # 10 bp on each side
        edited_DNA = "abcdefghij"+"1234567"+"qstuvwxyz"
        bp_position = 11
        edit_type = "insertion"
        template_DNA = "1234567"
        primer_length = 4
        search_size = 6
        out = find_possible_primers(DNA, edited_DNA, bp_position, bp_position, edit_type, template_DNA, primer_length, search_size, True)
        f = all_possible_subsets("abcdef", primer_length, "f")
        r = all_possible_subsets("vwxyz", primer_length, "r")
        c = all_possible_subsets("ijqs", primer_length, "c")
        e = all_possible_subsets("123456", primer_length, "e")
        expected = pd.concat([f, r, c, e])
        self.assertEqual(True, expected.equals(out))

        # test deletion
        DNA = "abcdefghij"+"klmnop"+"qstuvwxyz" # 10 bp on each side
        edited_DNA = "abcdefghij"+"qstuvwxyz"
        bp_position_start = 11
        bp_position_stop = 16
        edit_type = "deletion"
        template_DNA = ""
        primer_length = 4
        search_size = 6
        out = find_possible_primers(DNA, edited_DNA, bp_position_start, bp_position_stop, edit_type, template_DNA, primer_length, search_size, True)
        f = all_possible_subsets("abcdef", primer_length, "f")
        r = all_possible_subsets("uvwxyz", primer_length, "r")
        c = all_possible_subsets("ijklmn", primer_length, "c")
        e = all_possible_subsets("ijqs", primer_length, "e")
        expected = pd.concat([f, r, c, e])
        self.assertEqual(True, expected.equals(out))
    
    def test_find_primers(self):
        DNA = fasta_to_DNA("D29")
        bp_position_start = 500
        bp_position_stop = 800
        edit_type = "replacement"
        template_DNA = "GTCTCCGAGCGATCTATCCACGACCAATTTGACATGGGTGCGCCGTTTGTAAAGGCCGTGGACAAAG"
        substrate, edited_DNA = find_editing_substrate(DNA, bp_position_start, bp_position_stop, template_DNA)
        primer_sets = find_primers(DNA, edited_DNA, bp_position_start, bp_position_stop, edit_type, template_DNA)
        
        # check that control primers do not bind with edited DNA
        self.assertEqual([0], primer_sets["c edited start"].unique())

        # check that expirement primers do not bind with unedited DNA
        self.assertEqual([0], primer_sets["e start"].unique())

        DNA = fasta_to_DNA("D29")
        bp_position_start = 15000
        bp_position_stop = 16000
        edit_type = "deletion"
        substrate, edited_DNA = find_editing_substrate(DNA, bp_position_start, bp_position_stop)
        primer_sets = find_primers(DNA, edited_DNA, bp_position_start, bp_position_stop, edit_type)
        
        # check that control primers do not bind with edited DNA
        self.assertEqual([0], primer_sets["c edited start"].unique())

        # check that expirement primers do not bind with unedited DNA
        self.assertEqual([0], primer_sets["e start"].unique())

        DNA = fasta_to_DNA("D29")
        bp_position_start = 10000
        bp_position_stop = 10000
        template_DNA = "GTCTCCGAGCGATCTATCCACGACCAATTTGACATGGGTGCGCCGTTTGTAAAGGCCGTGGACAAAGCAGAACCCCCGGCACCGAGGGGGGCCGGGGGCCTGTTGGCGTGACGCCGTTGCTGCGCTCTTGGGGGTAGACGCGTCTCTCAGGCTTCGCGGTCGGCGAGGTCGGCGCCGAATATGCGATCGGCGATGGCTTCGGTGAGCGCAGACATGCTGAGCACCGCGTTGACTTGCACGAGGTCGTCGCCCGCCAGCAGCGCGCGGGTGTCGACGTACCACTTGCTGTTGCGGGTCTGCTTCTGGAGTTCTTCGGTGAGCGCGGCGCGCACTCGTAGCCGGATGTTGGATTCGACGGTCACTGGTTCTCCATTTCGGCGTCGGACAGCTGCCAGATGTGC"
        edit_type = "insertion"
        substrate, edited_DNA = find_editing_substrate(DNA, bp_position_start, bp_position_stop,template_DNA)
        primer_sets = find_primers(DNA, edited_DNA, bp_position_start, bp_position_stop, edit_type, template_DNA)
        # check that control primers do not bind with edited DNA
        self.assertEqual([0], primer_sets["c edited start"].unique())

        # check that expirement primers do not bind with unedited DNA
        self.assertEqual([0], primer_sets["e start"].unique())

if __name__ == '__main__':
    unittest.main()