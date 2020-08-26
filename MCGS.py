import math
import random

import networkx as nx
import numpy as np


class MCGS(object):
    """Create an mino-centric graph sampling (MCGS) algorithm model.

    Parameters
    ----------
    None

    Returns
    -------
    m : MCGS object
        An MCGS object represents the MCGS sampling algorithm instance model without
        properties. It has only one open interface method named "run_sampling".
        This method can be applied to the graph data and complete the sampling process.
        Besides, the run_sampling method can receive specific parameters, such as
        sampling rate, alpha, beta and loss weight, for multi-dimensional sampling of
        the graph data.

    Examples
    --------
    Create an MCGS object.

    >>> MCGS_model = MCGS()

    Run the sampling algorithm on graph data using the run_sampling() method
    with a sampling rate of 0.5 and other parameters by default.

    >>> G = nx.Graph()
    >>> G.add_nodes_from([1, 2, 3, 4, 5, 6, 7])
    >>> G.add_edges_from([(1, 2), (1, 3), (1, 4), (1, 5), (5, 6), (5, 7), (6, 7)])
    >>> list(G.nodes())
    [1, 2, 3, 4, 5, 6, 7]
    >>> Gs = MCGS_model.run_sampling(G, 0.5)  # sample the graph with a sampling rate of 0.5
    >>> list(Gs.nodes())
    [1, 2, 4, 5]

    More information about the run_sampling method will be introduced later.

    """

    def run_sampling(self, G, rate, alpha=1, beta=2, loss_weight=None):
        """Run the sampling algorithm on graph data.

        Parameters
        ----------
        G : Graph object
            The original graph data object with nodes, edges and other graph attributes.
        rate : float
            Sampling rate, namely, the proportion of the nodes preserved in the sample.
        alpha : int (optional, default=1)
            Controlling parameter of minority structure preservation, which controls
            the preservation of minority structures with the ratio of rate / alpha.
        beta : int (optional, default=2)
            Controlling parameter of neighborhood structures preservation, which controls
            the preservation of one-step neighbors of important minority structures with
            the ratio of rate / beta.
        loss_weight : list (optional, default=None)
            Weight coefficient list for importance evaluation of majority structures.
            It's a list of three numerical values measuring the weight of DD, NCC,
            and JI index respectively in the loss funtion of majority structure sampling.
            If None, then the weight coefficients will be set as 1:0:0 by default.

        Returns
        -------
        Gs : Graph object
            Sampling graph, an induced subgraph view of the original graph. The graph
            structure cannot be changed, but node/edge attributes can and are shared
            with the original graph.

        Examples
        --------
        Initialize the MCGS model and the original graph.
        >>> MCGS_model = MCGS()
        >>> G = nx.Graph()
        >>> G.add_nodes_from([1, 2, 3, 4, 5, 6, 7])
        >>> G.add_edges_from([(1, 2), (1, 3), (1, 4), (1, 5), (5, 6), (5, 7), (6, 7)])
        >>> list(G.nodes())
        [1, 2, 3, 4, 5, 6, 7]

        Run the MCGS sampling on graph G with a sampling rate of 0.5 and other parameters by default.
        >>> Gs = MCGS_model.run_sampling(G, 0.5)
        >>> list(Gs.nodes())
        [1, 2, 4, 5]
        >>> list(Gs.edges())
        [(1, 2), (1, 4), (1, 5)]

        The run_sampling() method can also receive user-defined sampling parameters.

        Set alpha parameter:

        >>> Gs = MCGS_model.run_sampling(G, 0.5, alpha=2)
        >>> list(Gs.nodes())
        [1, 2, 3, 5]
        >>> list(Gs.edges())
        [(1, 2), (1, 3), (1, 5)]

        Set beta parameter:

        >>> Gs = MCGS_model.run_sampling(G, 0.5, beta=4)
        >>> list(Gs.nodes())
        [1, 3, 4, 5]
        >>> list(Gs.edges())
        [(1, 3), (1, 4), (1, 5)]

        Set loss_weight parameter:

        >>> Gs = MCGS_model.run_sampling(G, 0.5, loss_weight=[1, 1, 1])
        >>> list(Gs.nodes())
        [1, 4, 5, 6]
        >>> list(Gs.edges())
        [(1, 4), (1, 5), (5, 6)]

        """
        # initialize sampling parameter settings
        self.__init_sampling_settings(rate, alpha, beta, loss_weight)

        # initialize original graph information
        self.__init_original_graph_info(G)

        # identify minority structures
        minority_structures = {}
        pivot_star_structures = self.__identify_pivot_and_star_structures(G)
        minority_structures.update(pivot_star_structures)
        rim_tie_structures = self.__identify_rim_and_tie_structures(G)
        minority_structures.update(rim_tie_structures)

        # sampling graph
        Gs = nx.Graph()

        # sample minority structures
        self.__minority_structure_sampling(G, Gs, minority_structures)

        # sample majority structures
        size = round(len(G.nodes()) * self.__rate)
        Gs_size = len(Gs.nodes())
        if Gs_size < size:
            self.__majority_structure_sampling(G, Gs, size)
        elif Gs_size > size:
            self.__delete_node_randomly(Gs, size)

        # induce subgraph from original graph based on the sampled nodes
        Gs = G.subgraph(Gs.nodes())

        return Gs

    def __init_sampling_settings(self, rate, alpha=1, beta=2, loss_weight=None):
        """Initialize sampling parameter settings."""
        self.__rate = rate
        self.__alpha = alpha
        self.__beta = beta
        if loss_weight is None:
            loss_weight = [1, 0, 0]
        self.__loss_weight = np.array(loss_weight) / sum(loss_weight)

    def __init_original_graph_info(self, G):
        """Initialize original graph information."""
        # dictionary recording the degree of each node
        self.__G_node_degree_dict = {item[0]: item[1] for item in G.degree()}
        # dictionary recording one-step neighbors of each node
        self.__G_neighbors_dict = {_: set(G.neighbors(_)) for _ in G.nodes()}
        # minority structure record set
        self.__G_minority_structures_set = set()

    def __identify_pivot_and_star_structures(self, G):
        """Identify super pivots and huge stars.

        Parameters
        ----------
        G : Graph object
            The original graph data object with nodes, edges and other graph attributes.

        Returns
        -------
        minority_structures : dictionary
            Minority structure dictionary, recording the nodes of pivots and stars
            in the original graph.

        """

        def __DFS(current_node, father_node=None, grandfather_node=None):
            """Detect triangle structures by depth first traversal."""
            seen.add(current_node)
            for other_node in self.__G_neighbors_dict[current_node]:
                if other_node == grandfather_node:
                    triangle_node_set.add(current_node)
                    triangle_node_set.add(father_node)
                    triangle_node_set.add(grandfather_node)
                elif other_node == father_node:
                    continue
                elif father_node and (other_node in self.__G_neighbors_dict[father_node]):
                    triangle_node_set.add(current_node)
                    triangle_node_set.add(father_node)
                    triangle_node_set.add(other_node)
                elif father_node and (other_node not in self.__G_neighbors_dict[father_node]):
                    continue
                else:
                    if other_node not in seen:
                        __DFS(other_node, current_node, father_node)
                    else:
                        __DFS(other_node, current_node)

        original_nodes = set(G.nodes())
        seen = set()
        triangle_node_set = set()
        minority_structures = {
            'pivot': [],
            'star': []
        }

        # detect triangle structures among nodes with top 5% degrees,
        # and identify super pivots
        sorted_node_by_degree = sorted(original_nodes,
                                       key=lambda x: self.__G_node_degree_dict[x],
                                       reverse=True)
        nodes_5 = list(sorted_node_by_degree)[:math.ceil(len(original_nodes) * 0.05)]
        for node in nodes_5:
            if (node not in seen) or (node not in triangle_node_set): __DFS(node)
        minority_structures['pivot'] = [_ for _ in nodes_5 if _ in triangle_node_set]

        # detect triangle structures among nodes with degrees higher than average,
        # and identify huge stars
        degree_average = np.average(list(self.__G_node_degree_dict.values()))
        nodes_average = set(
            filter(lambda item: self.__G_node_degree_dict[item] >= degree_average,
                   original_nodes))
        for node in nodes_average:
            if (node not in seen) or (node not in triangle_node_set): __DFS(node)
        minority_structures['star'] = [_ for _ in nodes_average if _ not in triangle_node_set]

        return minority_structures

    def __identify_rim_and_tie_structures(self, G):
        """Identify rims and ties.

        Parameters
        ----------
        G : Graph object
            The original graph data object with nodes and edges.

        Returns
        -------
        minority_structures : dictionary
            Minority structure dictionary, recording the nodes of rims and ties in
            the original graph.

        """
        # one-degree node set in the original graph
        one_degree_node_set = set(
            filter(lambda item: self.__G_node_degree_dict[item] == 1,
                   self.__G_node_degree_dict.keys()))

        # cut point set in original graph
        cut_points = set(nx.articulation_points(G))
        # induced subgraph from the original graph based on cut points
        cut_points_graph = G.subgraph(cut_points)

        # dictionary recording the degree of each cut point in the induced subgraph
        cut_point_degree_dict = {_[0]: _[1] for _ in cut_points_graph.degree()}
        # dictionary recording the neighbor nodes of each cut points in the induced subgraph
        cut_point_neighbors_records = {_: list(cut_points_graph.neighbors(_)) for _ in cut_points}
        # dictionary recording the neighbors of each cut point in the original graph
        cut_point_neighbors_records_in_G = {_: set(G.neighbors(_)) for _ in cut_points}

        chains_list = []
        parachute_set = set()

        # start from cut points with degrees lower than 1 as end points of chain structures, 
        # and traverse all cut points
        iter_nodes = list(filter(lambda item: cut_point_degree_dict[item] <= 1, cut_points))
        seen = set()
        while len(cut_points - seen) > 0:
            for node in iter_nodes:
                if node not in seen:
                    seen.add(node)
                    temp_chain = [node]
                    current_node = node

                    while 1:
                        # detect parachute structure nodes
                        if cut_point_neighbors_records_in_G[current_node] - set(
                                cut_point_neighbors_records[current_node]):
                            parachute_set.add(current_node)
                        neighbors = cut_point_neighbors_records[current_node]

                        # the loop exits when the current cut point does not have cut point neighbors,
                        # or else, move to next node and continue traversal
                        if not neighbors or len(set(neighbors) - seen) >= 2:
                            break
                        else:
                            other = neighbors[0]
                            for _ in neighbors:
                                if _ not in seen:
                                    other = _
                                    break
                            if other not in seen:
                                seen.add(other)
                                temp_chain.append(other)
                                current_node = other
                            else:
                                break

                    if len(temp_chain) > 1: chains_list.append(temp_chain)

            remaining_cut_points = cut_points - seen
            cut_point_degree_dict = {_[0]: _[1] for _ in G.subgraph(remaining_cut_points).degree()}
            iter_nodes = set(
                filter(lambda item: cut_point_degree_dict[item] <= 1,
                       remaining_cut_points))
            if not iter_nodes:
                for cut_point in remaining_cut_points:
                    if cut_point_neighbors_records_in_G[cut_point] - set(
                            cut_point_neighbors_records[cut_point]):
                        parachute_set.add(cut_point)
                break

        minority_structures = {
            'rim': [],
            'tie': []
        }

        # sort the chains according to their length
        chains_list.sort(key=lambda chain: len(chain), reverse=True)

        for chain_item in chains_list:
            # get the set of one-degree neighbor nodes of the end nodes of each chain
            a_one_nodes = cut_point_neighbors_records_in_G[chain_item[0]] & one_degree_node_set
            b_one_nodes = cut_point_neighbors_records_in_G[chain_item[-1]] & one_degree_node_set

            # if any end node of the chain has one and only one neighbor with degree 1 in G,
            # then the chain is a rim, else a tie.
            if len(a_one_nodes) == 1 or len(b_one_nodes) == 1:
                minority_structures['rim'].append(chain_item)
            else:
                minority_structures['tie'].append(chain_item)

        # sort the parachute-shaped rims according to degrees and add them to rim records
        parachute_sorted = sorted(parachute_set,
                                  key=lambda item: self.__G_node_degree_dict[item],
                                  reverse=True)
        minority_structures['rim'].extend(list(parachute_sorted))

        return minority_structures

    def __minority_structure_sampling(self, G, Gs, minority_structures):
        """Sample the neighborhood structures of minority structures."""
        # truncate minority structures according to the importance of them
        for key, value in minority_structures.items():
            truncation = math.ceil(len(value) * self.__rate / self.__alpha)
            minority_structures[key] = value[:truncation]

        # sample the minority structures and their one-step neighbors
        for value in minority_structures.values():
            for item in value:
                if type(item) == int:
                    # the single node structures are sampled by points
                    self.__G_minority_structures_set.add(item)
                    Gs.add_node(item)

                    # conduct random sampling of neighbors
                    neighbor_list = list(self.__G_neighbors_dict[item])
                    sample_size = math.ceil(
                        len(neighbor_list) * self.__rate / self.__beta)
                    sample_node_list = random.sample(neighbor_list, sample_size)
                else:
                    # the chain structures are sampled as a whole
                    self.__G_minority_structures_set |= set(item)
                    Gs.add_nodes_from(item)

                    # conduct random sampling of neighbors
                    node_list = []
                    for node in item:
                        node_list.extend(list(self.__G_neighbors_dict[node]))
                    neighbor_list = list(set(node_list) - set(item))
                    sample_size = math.floor(
                        len(neighbor_list) * self.__rate / self.__beta)
                    sample_node_list = random.sample(neighbor_list, sample_size)

                Gs.add_nodes_from(sample_node_list)

    def __majority_structure_sampling(self, G, Gs, size):
        """Sample the majority structures according to their loss evaluation."""
        G_nodes = set(G.nodes())
        Gs_nodes = set(Gs.nodes())
        candidate_node_set = G_nodes - Gs_nodes
        current_Gs_degrees = {_[0]: _[1] for _ in list(G.subgraph(Gs_nodes).degree())}

        # the difference between the node degrees in the sample and the original graph
        sub_degrees = {key: self.__G_node_degree_dict[key] - value for
                       (key, value) in current_Gs_degrees.items()}

        # the temp variable of MSE, to avoid second computation of the mean square
        # error of degrees between Gs and G
        temp_MSE = ((np.array(list(sub_degrees.values()))) ** 2).sum()

        # get the number of connected components
        subgraph = G.subgraph(Gs_nodes)
        temp_graph = nx.Graph()
        temp_graph.add_nodes_from(subgraph.nodes())
        temp_graph.add_edges_from(subgraph.edges())
        # the number of connected components in original graph
        G_connected_components = len(sorted(nx.connected_components(G)))
        # the number of connected components in current sampling graph
        temp_connected_components = len(sorted(nx.connected_components(temp_graph)))
        
        # the similarity of connectivity between Gs and G, and the initial value is 1
        temp_NCC = 1

        # the structural similarity between Gs and G
        temp_JI = 0
        for node in Gs_nodes:
            if self.__G_node_degree_dict[node] != 0:
                temp_JI += current_Gs_degrees[node] / self.__G_node_degree_dict[node]

        while len(Gs) < size:
            min_loss_record_dict = {
                'MSE': {
                    'value': temp_MSE,
                    'node': None,
                    'new_edge_nodes': []
                },
                'NCC': {
                    'value': temp_NCC,
                    'node': None,
                    'new_edge_nodes': []
                },
                'JI': {
                    'value': temp_JI,
                    'node': None,
                    'new_edge_nodes': []
                }
            }
            max_loss_record_dict = {
                'MSE': temp_MSE,
                'NCC': temp_NCC,
                'JI': temp_JI
            }

            # loss metrics record of each candidate node
            all_node_record_dict = {}

            for iter_node in candidate_node_set:
                # the neighbors of iter_node in Gs
                new_edge_nodes = list(Gs_nodes & self.__G_neighbors_dict[iter_node])

                # calculate and update the MSE and JI
                new_MSE = temp_MSE + (self.__G_node_degree_dict[iter_node] - len(new_edge_nodes)) ** 2
                if self.__loss_weight[2] == 0 or not self.__G_neighbors_dict[iter_node]:
                    new_JI = 0
                else:
                    new_JI = temp_JI + len(new_edge_nodes) / len(self.__G_neighbors_dict[iter_node])

                for other_node in new_edge_nodes:
                    new_MSE += - 2 * sub_degrees[other_node] + 1
                    if self.__loss_weight[2] == 0 or self.__G_node_degree_dict[other_node] == 0:
                        new_JI += 0
                    else:
                        new_JI += 1 / self.__G_node_degree_dict[other_node]

                # calculate the NCC
                if self.__loss_weight[1] == 0:
                    new_NCC = 0
                elif new_edge_nodes:
                    new_NCC = temp_connected_components / G_connected_components
                else:
                    new_NCC = (temp_connected_components + 1) / G_connected_components

                # compare each loss metrics and update the minimal or maximum loss record dictionary
                if self.__loss_weight[0] != 0:
                    if new_MSE <= min_loss_record_dict['MSE']['value']:
                        min_loss_record_dict['MSE']['value'] = new_MSE
                        min_loss_record_dict['MSE']['node'] = iter_node
                        min_loss_record_dict['MSE']['new_edge_nodes'] = new_edge_nodes
                    if new_MSE >= max_loss_record_dict['MSE']:
                        max_loss_record_dict['MSE'] = new_MSE
                if self.__loss_weight[1] != 0:
                    if new_NCC <= min_loss_record_dict['NCC']['value']:
                        min_loss_record_dict['NCC']['value'] = new_NCC
                        min_loss_record_dict['NCC']['node'] = iter_node
                        min_loss_record_dict['NCC']['new_edge_nodes'] = new_edge_nodes
                    if new_NCC >= max_loss_record_dict['NCC']:
                        max_loss_record_dict['NCC'] = new_NCC
                if self.__loss_weight[2] != 0:
                    if new_JI <= min_loss_record_dict['JI']['value']:
                        min_loss_record_dict['JI']['value'] = new_JI
                        min_loss_record_dict['JI']['node'] = iter_node
                        min_loss_record_dict['JI']['new_edge_nodes'] = new_edge_nodes
                    if new_JI >= max_loss_record_dict['JI']:
                        max_loss_record_dict['JI'] = new_JI

                # save the calculation result to the record
                all_node_record_dict[iter_node] = {
                    'new_edge_nodes': new_edge_nodes,
                    'MSE': new_MSE,
                    'NCC': new_NCC,
                    'JI': new_JI
                }

            # calculate the overall loss according to the loss weight,
            # and get the optimal candidate node
            select_node = None
            if self.__loss_weight[0] == 0 and self.__loss_weight[1] == 0:
                select_node = min_loss_record_dict['JI']['node']
            elif self.__loss_weight[0] == 0 and self.__loss_weight[2] == 0:
                select_node = min_loss_record_dict['NCC']['node']
            elif self.__loss_weight[1] == 0 and self.__loss_weight[2] == 0:
                select_node = min_loss_record_dict['MSE']['node']
            else:
                min_loss_function_value = 1
                for key_node, value in all_node_record_dict.items():
                    MSE_detal = max_loss_record_dict['MSE'] - min_loss_record_dict['MSE']['value']
                    NCC_detal = max_loss_record_dict['NCC'] - min_loss_record_dict['NCC']['value']
                    JI_detal = max_loss_record_dict['JI'] - min_loss_record_dict['JI']['value']
                    try:
                        if self.__loss_weight[0] == 0:
                            loss_MSE_value = 0
                        else:
                            loss_MSE_value = (value['MSE'] - min_loss_record_dict['MSE']['value']) / MSE_detal

                        if self.__loss_weight[1] == 0:
                            loss_NCC_value = 0
                        else:
                            loss_NCC_value = (value['NCC'] - min_loss_record_dict['NCC'][0]) / NCC_detal

                        if self.__loss_weight[2] == 0:
                            loss_JI_value = 0
                        else:
                            loss_JI_value = (value['JI'] - min_loss_record_dict['JI'][0]) / JI_detal

                        # calculate the loss function of the current node and update the minimal loss record
                        current_loss_function_value = loss_MSE_value * self.__loss_weight[0] + \
                                                      loss_NCC_value * self.__loss_weight[1] + \
                                                      loss_JI_value * self.__loss_weight[2]
                        if current_loss_function_value <= min_loss_function_value:
                            min_loss_function_value = current_loss_function_value
                            select_node = key_node
                    except:
                        continue

            if select_node is None:
                select_node = candidate_node_set.pop()
            else:
                candidate_node_set.remove(select_node)

            # update the degree difference of the candidate node
            sub_degrees[select_node] = self.__G_node_degree_dict[select_node] - len(
                all_node_record_dict[select_node]['new_edge_nodes'])
            for other_node in all_node_record_dict[select_node]['new_edge_nodes']:
                sub_degrees[other_node] -= 1

            Gs.add_node(select_node)
            Gs_nodes.add(select_node)

            # update the MSE, JI, NCC and the number of connected components between G and Gs
            temp_MSE = all_node_record_dict[select_node]['MSE']
            temp_JI = all_node_record_dict[select_node]['JI']
            temp_NCC = all_node_record_dict[select_node]['NCC']
            if not all_node_record_dict[select_node]['new_edge_nodes']:
                temp_connected_components += 1

    def __delete_node_randomly(self, Gs, size):
        """Delete extra nodes randomly, prioritizing majority structures."""
        Gs_nodes = set(Gs.nodes())
        count = len(Gs_nodes) - size

        no_minority_structure_nodes = list(Gs_nodes - self.__G_minority_structures_set)
        if len(no_minority_structure_nodes) >= count:
            for _ in range(count):
                random_node = random.choice(no_minority_structure_nodes)
                Gs.remove_node(random_node)
                Gs_nodes.remove(random_node)
                no_minority_structure_nodes.remove(random_node)
        else:
            Gs.remove_nodes_from(no_minority_structure_nodes)
            Gs_nodes -= set(no_minority_structure_nodes)
            count -= len(no_minority_structure_nodes)
            Gs_nodes = list(Gs_nodes)
            for _ in range(count):
                random_node = random.choice(Gs_nodes)
                Gs.remove_node(random_node)
                Gs_nodes.remove(random_node)

    def __init__(self, *args, **kwargs):
        pass