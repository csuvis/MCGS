import csv

import networkx as nx


def loadGraph(graph_name):
    """Load the graph file and return the graph object.

    Parameters
    ----------
    graph_name : string
        The name of the graph.

    Returns
    -------
    G : Graph object
        The graph data object with nodes, edges and other graph attributes.

    Notes
    -----
    The parameter graph_name is the name of graph, and it does not contain any suffix or extra
    information. For example, we can load the graph "Cpan" by using loadGraph("Cpan")
    instead of loadGraph("Cpan.csv").

    """
    with open('./dataset/{}_node.csv'.format(graph_name), 'r',
              encoding='utf-8') as fp:
        reader = csv.reader(fp)
        nodes = list(int(_[0]) for _ in reader)
    with open('./dataset/{}_edge.csv'.format(graph_name), 'r',
              encoding='utf-8') as fp:
        reader = csv.reader(fp)
        edges = list((int(_[0]), int(_[1])) for _ in reader if _[0] != _[1])
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    return G
