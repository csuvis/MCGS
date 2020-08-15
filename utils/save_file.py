import re


def saveAsSVG(dataset, G):
    """Save the sampled graph data as SVG type file.

    Parameters
    ----------
    dataset : string
        The name of the graph.
    G : Graph object
        The graph data object with nodes, edges and other graph attributes.

    Returns
    -------
    None

    Notes
    -----
    This method is based on the SVG file of the original graph in the directory
    'dataset'. Please store the basic SVG file in the directory 'dataset' in advance.
    And the result file will be stored in the directory 'output'.

    """
    # load the original SVG file information
    with open("./dataset/{}.svg".format(dataset), 'r', encoding='utf-8')as fp:
        svg_text = fp.readline()
        svg_header = ''
        while ('version="1.1"' not in svg_text):
            svg_header += svg_text
            svg_text = fp.readline()
        svg_header += svg_text
        svg_text = fp.read()

    # information records to be written
    svg_nodes = []
    svg_edges = []

    # edges and nodes information in the original SVG file
    nodes = re.findall('<circle .*\n.*\n.*/>', svg_text)
    edges = re.findall('<path .*\n.*\n.*/>', svg_text)

    # edges and nodes in the sampled graph
    G_nodes = G.nodes()
    G_edges = G.edges()

    # traverse the original graph data to obtain the SVG information of the sampled graph
    for _ in edges:
        edge_num = re.findall('class="id_([0-9]*) id_([0-9]*)"', _)[0]
        edge_id = (int(edge_num[0]), int(edge_num[1]))
        if (edge_id[0], edge_id[1]) in G_edges or (edge_id[1], edge_id[0]) in G_edges:
            svg_edges.append(_)
    for _ in nodes:
        node_id = int(re.findall('class="id_([0-9]*)"', _)[0])
        if node_id in G_nodes:
            svg_nodes.append(_)

    # save the SVG information of the sampled graph
    with open("./output/{}_MCGS.svg".format(dataset), 'w', encoding='utf-8')as fp:
        fp.write(svg_header)
        fp.write('    <g>\n        {}\n</g>\n'.format('\n        '.join(svg_edges)))
        fp.write('    <g>\n        {}\n</g>\n'.format('\n        '.join(svg_nodes)))
        fp.write("</svg>")
