import copy
import unittest


import crispr_spacer_alignment as csa


import igraph


class TestTopoOrder(unittest.TestCase):

    def setUp(self):
        self.order_one = ['a', 'b', 'c', 'd']
        self.order_many = [copy.copy(self.order_one) for i in range(4)]

        self.order_many_misorder = copy.deepcopy(self.order_many)
        misorder_list = self.order_many_misorder[3]
        misorder_list[2], misorder_list[3] = misorder_list[3], misorder_list[2]

        self.graph_ordered = self.create_graph(self.order_many)
        self.graph_misordered = self.create_graph(self.order_many_misorder)


    @staticmethod
    def create_graph(symbols_lists):
        graph = igraph.Graph(directed=True)
        for vertex in sorted({s for so in symbols_lists for s in so}):
            graph.add_vertex(vertex)

        for i, symbols in enumerate(symbols_lists):
            for i in range(len(symbols) - 1):
                graph.add_edge(symbols[i], symbols[i+1], name=i)
        return graph


    def test_ordered(self):
        strong_order, misorders = csa.get_spacer_order(self.graph_ordered)
        self.assertEqual(strong_order, [0, 1, 2, 3])
        self.assertEqual(misorders, None)


    def test_misordered(self):
        strong_order, misorders = csa.get_spacer_order(self.graph_misordered)
        self.assertEqual(strong_order, [0, 1, 2, 3])
        self.assertEqual(misorders, {2: [['d', 'c']]})
