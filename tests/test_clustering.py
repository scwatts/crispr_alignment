import unittest
import unittest.mock


import crispr_spacer_alignment as csa


class TestClustering(unittest.TestCase):

    def setUp(self):
        crisprs = unittest.mock.Mock()
        crisprs.spacers_seqs = ['aaaaaaaaaaa',
                                'aaaaaaaaaat',
                                'aaaaaaaaaac',
                                'aaaaaaaaaag',
                                'ggggggggggg']
    def test_native_absolute(self):
        args = unittest.mock.Mock()
        args.identity = 0.90
        args.native_orientation = True
        args.absolute_identity = False

        csa.cluster_spacer_sequences(self.crisprs, args)
