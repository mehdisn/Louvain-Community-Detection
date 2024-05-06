import argparse
from itertools import product
import networkx as nx
from collections import defaultdict
import random

def modularity(G, partition):
    m = G.size(weight="weight")
    degrees = dict(G.degree(weight="weight"))
    Q = 0
    for community in partition:
        for u, v in product(community, repeat=2):
            try:
                w = G[u][v].get("weight", 1)
            except KeyError:
                w = 0
            if u == v:
                w *= 2
            Q += w - degrees[u] * degrees[v] / (2 * m)
    return Q / (2 * m)

class CommunityTracker:
    def __init__(self):
        self.node_to_community_map = None
        self.m = 0.0
        self.degrees = None
        self.self_loops = None
        self.community_degrees = None
        self.community_self_loops = None

    def initialize_network_statistics(self, G):
        self.node_to_community_map = {}
        self.m = G.size(weight="weight")
        self.degrees = {}
        self.self_loops = {}
        self.community_degrees = {}
        self.community_self_loops = {}
        for community, node in enumerate(G):
            self.node_to_community_map[node] = community
            degree = G.degree(node, weight="weight")
            self.degrees[node] = self.community_degrees[community] = degree
            self_loop = 0
            if G.has_edge(node, node):
                self_loop = G[node][node].get("weight", 1)
            self.community_self_loops[community] = self.self_loops[node] = self_loop

    def remove(self, node, community, incident_weight):
        self.community_degrees[community] -= self.degrees[node]
        self.community_self_loops[community] -= incident_weight + self.self_loops[node]
        self.node_to_community_map[node] = None

    def insert(self, node, community, incident_weight):
        self.community_degrees[community] += self.degrees[node]
        self.community_self_loops[community] += incident_weight + self.self_loops[node]
        self.node_to_community_map[node] = community

class Louvain:
    def __init__(self, G, verbose=False, randomized=False):
        self.verbose = verbose
        self.randomized = randomized
        self.tracker = CommunityTracker()
        self.original_graph = G
        self.coarse_grain_graph = G
        self.community_history = []
        self.iteration_count = 0
        self.finished = False
        self.community_map = None
        self.communities = None

    def iterate(self):
        self.iteration_count += 1
        if self.verbose:
            print("Iteration: ", self.iteration_count)
        modified = False
        improved = True
        G = self.coarse_grain_graph
        self.tracker.initialize_network_statistics(G)
        community_map = self.tracker.node_to_community_map

        while improved:
            improved = False

            nodes = G.nodes()
            if self.randomized:
                nodes = list(G.nodes())
                random.seed()
                random.shuffle(nodes)

            for node in nodes:
                best_delta_Q = 0.0
                old_community = community_map[node]
                new_community = old_community
                neighbour_communities = self.get_neighbour_communities(
                    G, node, community_map)
                old_incident_weight = neighbour_communities.get(
                    old_community, 0)
                self.tracker.remove(node, old_community, old_incident_weight)
                for community, incident_wt in neighbour_communities.items():
                    delta_Q = self.calculate_delta_Q(
                        G, node, community, incident_wt)
                    if delta_Q > best_delta_Q:
                        best_delta_Q = delta_Q
                        new_community = community

                new_incident_weight = neighbour_communities[new_community]
                self.tracker.insert(node, new_community, new_incident_weight)
                if self.verbose:
                    message = "Moved node {} from community {} to community {}"
                    print(message.format(node, old_community, new_community))

                if new_community != old_community:
                    improved = True
                    modified = True

        if modified:
            self.relabel_community_map(community_map)
            self.community_history.append(community_map)
            self.coarse_grain_graph = self.generate_coarse_grain_graph(
                G, community_map)
        else:
            self.finished = True

    def get_neighbour_communities(self, G, node, community_map):
        neighbour_communities = defaultdict(int)
        for neighbour in G[node]:
            if neighbour != node:
                neighbour_community = community_map[neighbour]
                w = G[node][neighbour].get("weight", 1)
                neighbour_communities[neighbour_community] += w
        return neighbour_communities

    def calculate_delta_Q(self, G, node, community, incident_weight):
        sigma_tot = self.tracker.community_degrees[community]
        k_i = self.tracker.degrees[node]
        k_i_in = incident_weight
        m = self.tracker.m

        delta_Q = 2 * k_i_in - sigma_tot * k_i / m
        return delta_Q

    def generate_coarse_grain_graph(self, G, community_map):
        new_graph = nx.Graph()
        for community in set(community_map.values()):
            new_graph.add_node(community)
        for u, v, w in G.edges(data="weight", default=1):
            c1 = community_map[u]
            c2 = community_map[v]
            new_weight = w
            if new_graph.has_edge(c1, c2):
                new_weight += new_graph[c1][c2].get("weight", 1)
            new_graph.add_edge(c1, c2, weight=new_weight)
        return new_graph

    def relabel_community_map(self, community_map):
        community_labels = set(community_map.values())
        relabelled_communities = {j: i for i, j in enumerate(community_labels)}
        for node in community_map:
            community_map[node] = relabelled_communities[community_map[node]]

    def invert_community_map(self, community_map):
        inverted_community_map = defaultdict(list)
        for node in community_map:
            inverted_community_map[community_map[node]].append(node)
        return list(inverted_community_map.values())

    def generate_community_map(self, community_history):
        community_map = {node: node for node in self.original_graph}
        for node in community_map:
            for iteration in community_history:
                community_map[node] = iteration[community_map[node]]
        return community_map

    def run(self, output_file):
        i=0
        while not self.finished:
            i=i+1
            self.iterate()
            self.community_map = self.generate_community_map(self.community_history)
            self.communities = self.invert_community_map(self.community_map)
            mod = modularity(self.original_graph, self.communities)

            with open(output_file, "a") as f:
                f.write("Round " + str(i) + ":" + "\n")    
                for index, community in enumerate(self.communities):
                    f.write("Community " + str(index) + " " + "nodes: " + str(community) + "\n")
                f.write("Modularity: " + str(mod) + "\n\n")

        if self.verbose:
            print("Finished in {} iterations".format(self.iteration_count))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="Input file name")
    parser.add_argument("output_file", help="Output file name")
    args = parser.parse_args()
    output = args.output_file + "_with_Modularity"
    
    G = nx.Graph()
    with open(args.input_file) as f:
        for line in f:
            split_line = line.split()
            if len(split_line) == 2:
                u, v = map(int, split_line)
                G.add_edge(u, v)
            elif len(split_line) == 3:
                u, v, w = map(int, split_line)
                G.add_edge(u, v, weight=w)

    louvain = Louvain(G, verbose=False, randomized=False)
    louvain.run(output)
    partitions = louvain.communities

    with open(args.output_file, "w") as f:
        for index, community in enumerate(partitions):
            f.write("Community " + str(index) + ": " + str(community) + "\n")
    