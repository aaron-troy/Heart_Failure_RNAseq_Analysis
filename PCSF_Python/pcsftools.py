"""
Tools for deploying the OmicsIntegrator package to perform PCSF analysis.

OmicsIntegrator: https://github.com/fraenkel-lab/OmicsIntegrator2
"""
import pandas as pd
import itertools
import networkx as nx
import graph
import random
from _datetime import datetime

def run_pcsf(prizes, ppi, params, time = False, add_attribs = True):
    """
    Flagship function for running PCSF using the Omics integrator package
    Args:
        prizes: pandas DataFrame; seed nodes for PCSF solution. Should have columns ['name','prize','score','gene']
        ppi: pandas DataFrame; protein-protein interaction network in adj list format.
        params: dict; keys: strings for fast PCSF parameters, values: various.
        time: bool; optional, default False; set True to print current time. Useful in the case of successive calls
        add_attribs: bool; optional, default True. set False to skip adding attributes to nodes in the graph.
        Attributes include multitple centrality measures for each node and Louvain cluster labels.

    Returns:
        networkx graph object; PCSF result.
    """

    # If timing is desired, print the current time
    if time:
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print("Current Time =", current_time)

    # Construct the PCSF graph object from the PPI adj list
    pcsfGraph = graph.Graph(ppi, params)

    # Prepare the prizes
    pcsfGraph.prepare_prizes(prizes)

    # Get the solution
    (verts, edgs) = pcsfGraph.pcsf()

    # Build into networkx object
    resultForest, _ = pcsfGraph.output_forest_as_networkx(verts, edgs)

    if add_attribs:

        # Add some attributes - multiple centrality measures and Louvain clustering
        nx.set_node_attributes(resultForest, {node: {'b_centrality' :b_centrality} for node ,b_centrality in
                                              nx.betweenness_centrality(resultForest).items()})

        nx.set_node_attributes(resultForest, {node: {'ev_centrality' :ev_centrality} for node ,ev_centrality in
                                              nx.eigenvector_centrality(resultForest).items()})

        nx.set_node_attributes(resultForest, {node: {'deg_centrality' :deg_centrality} for node ,deg_centrality in
                                              nx.degree_centrality(resultForest).items()})

        nx.set_node_attributes(resultForest, {node: {'louvain_clusters' :str(cluster)} for node ,cluster in
                                              community.best_partition(resultForest).items()})

    return resultForest


def _get_similarity(a: list, b: list):
    """
    Compute the Jaccard similarity between two lists of items
    Args:
        a: list;
        b: list;

    Returns:
        float; Jaccard similarity between the two lists

    """
    return len(set(a).intersection(set(b))) / len(set(a).union(set(b)))


def get_path_costs(seeds, adj_list, print_prog=False, sampled_proportion = 1):
    """
    Compute the shortest path (minimum cost) within a graph for all combinations of a passed list of graph nodes.
    Args:
        seeds: list of graph nodes for use as start and end of the shortest paths
        adj_list: dataframe representation of adjacency list for the graph in question. Columns 1 and 2 should be
        the head and tail of each edge, column 3 the edge weight / cost
        print_prog: bool, optional, default False; Set True to print percentage of path costs computed
        sampled_proportion: float <= 1, > 0, optional, default 1; Proportion of all possible start/end combinations
        to compute path costs for. If <1, randomly sample the combinations until the desired proportion sampled.

    Returns:
        dic; keys: tuples (start, end), values: path costs (float)
    """
    # Check sampled proporation is reasonable
    try:
        assert sampled_proportion > 0
        assert sampled_proportion <= 1
    except AssertionError:
        print("Sampled proportion must be in (0, 1]")
        return None

    # Construct graph
    G = nx.from_pandas_edgelist(adj_list, source=adj_list.columns[0], target=adj_list.columns[1], edge_attr=adj_list.columns[2])

    # Get all possible start/end combinations
    combs = list(itertools.combinations(seeds, 2))

    # Dict for path costs
    path_costs = {}

    i = 0

    # Number of samples
    samples = int(len(combs) * sampled_proportion)

    for i in range(samples):

        if print_prog:
            print("Progress:", 100*i / samples, "%")

        # Get next home/dest pair to compute the path cost for
        (home, dest) = combs.pop(random.randrange(len(combs)))

        # Tr
        try:
            path_costs[(home, dest)] = nx.dijkstra_path_length(G, home, dest, weight="cost")
        except NetworkXNoPath:
            # If home/dest are not connected, continue
            continue
        i = i + 1

    return path_costs


def convert_to_gene_symbol(string_ids : list, map_path: str):
    """
    Convert a list of STRING PPI ids to a list of corresponding gene symbols
    Args:
        string_ids: list; STRING IDs to be converted to gene symbols
        map_path: str; path to file containing the map from STRING ID to gene symbols

    Returns:
        list; corresponding gene symbols
    """

    # Read in the map from the passed path
    tmp = pd.read_csv(map_path, sep='\t')
    string2gene = dict(zip(tmp['STRING'], tmp['display name']))

    # Output
    gene_syms = []

    for s in string_ids:
        try:
            gene_syms.append(string2gene[s])
        except KeyError:
            # If an ID isn't in the map, just keep the old name
            gene_syms.append(s)

    return gene_syms


def convert_to_string_id(genes, map_path):
    """
    Converts a list of gene symbols to a list of corresponding STRING PPI IDs
    Args:
        genes: list; gene symbols to convert to STRING IDs
        map_path: str; path with file containing the map between STRING IDs and gene symbols

    Returns:
        list of str; corresponding STRING IDs

    """

    # Read in the map from the passed path
    tmp = pd.read_csv(map_path, sep='\t')
    gene2string = dict(zip(tmp['display name'], tmp['STRING']))

    # Output
    string_ids = []

    for g in genes:
        try:
            string_ids.append(gene2string[g])
        except KeyError:
            # If a gene isn't in the map, just keep the old name
            string_ids.append(g)

    return string_ids


def get_network_details(G: networkx, atts: list=[], sort_by: bool=False):
    """
    Function to produce a summary DF containing centrality measures and Louvain cluster annotation
    Args:
        G: input graph, a networkx object
        atts: list of attributes to retrieve for each node in the graph
        sort_by: optional; attribute to sort the results by

    Returns:
        pandas dataframe with attribute details for each node in the input graph
    """
    att_df = pd.DataFrame(index=G.nodes)

    for a in atts:
        att_df[a] = nx.get_node_attributes(G, a).values()

    if sort_by:
        assert sort_by in att_df.columns
        return att_df.sort_values(sortBy, ascending=False)
    else:
        return att_df


def get_graph_edit_distance(graph_dict:dict):
    """
     As an alternative to Jaccard similarity, compute the graph edit distance between all combinations of graphs
     This is VERY slow; Not for use with large graphs
    Args:
        graph_dict: dict; keys: hashable identifiers, values: networkx objects; graphs to compute the distance between
    Returns:
        dict; keys: hashable identifier for pair of graphs, values: float, graph edit distance
    """

    # Get all combinations
    combs = list(itertools.product(graph_dict.keys(), graph_dict.keys()))

    edit_distances = {}

    for k1, k2 in combs:
        edit_distances[(k1, k2)] = nx.graph_edit_distance(graph_dict[k1], graph_dict[k2])

    return edit_distances


def make_solution_similarity_matrix(res_dict: dict):
    """
    Compute the similarity matrix for all PCSF solution in passed dictionary
    Args:
        res_dict: dict; keys: hashable PCSF identifier, values: PCSF graph object, solution

    Returns:
        pandas dataframe; columns 1,2 identify the PCSF solutions compared. column 3 is the computed Jaccard similarity
    """

    # All combinations
    combs = list(itertools.product(res_dict.keys(), res_dict.keys()))

    # Output
    sim_mat = []

    for k1, k2 in combs:
        (G1, G2) = (res_dict[k1], res_dict[k2])
        sim_mat.append([k1, k2, _get_similarity(G1, G2)])

    sim_mat = pd.DataFrame(sim_mat, columns=['PCSF1', 'PCSF2', 'M']).pivot(index='PCSF1', columns='PCSF2', values='M')

    return simMat