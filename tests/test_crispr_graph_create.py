import copy
import unittest


import crispr_spacer_alignment as csa


import igraph


class TestCrisprGraphCreate(unittest.TestCase):

    def setUp(self):
        self.crisprs = list()
        for i in range(4):
            crispr = csa.Crispr(i, None, None, None)
            crispr.spacers = ['a', 'b', 'c', 'd']
            self.crisprs.append(crispr)

        self.crisprs_misorder = copy.deepcopy(self.crisprs)
        crispr_misorder = self.crisprs_misorder[3]
        crispr_misorder.spacers[2] = 'd'
        crispr_misorder.spacers[3] = 'c'

    def test_ordered_graph(self):
        graph = csa.create_spacer_graph(self.crisprs)
        nodes = [n['name'] for n in graph.vs]
        self.assertEqual(nodes, ['a', 'b', 'c', 'd'])

    def test_misordered_graph(self):
        graph = csa.create_spacer_graph(self.crisprs_misorder)
        nodes = [n['name'] for n in graph.vs]
        self.assertEqual(nodes, ['a', 'b', 'c', 'd'])
