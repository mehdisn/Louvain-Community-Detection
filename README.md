# Louvain Community Detection

This Python script implements the Louvain community detection algorithm for detecting communities in networks. The Louvain algorithm is a popular method for identifying communities in large networks based on modularity optimization.

## Requirements

- Python 3.x
- NetworkX

## Usage

Run the script using the following command:

```bash
python louvain.py input_file output_file
```
Replace input_file with the name of the file containing the network data and output_file with the desired name of the output file.

## Input File Format

The input file should contain the edges of the graph, with each line representing an edge. The format depends on whether the graph is weighted or unweighted:

-   For an unweighted graph: Each line should contain two integers representing the endpoints of an edge.
-   For a weighted graph: Each line should contain three integers representing the endpoints of an edge followed by the weight of the edge.

Example of an unweighted graph:

```bash
1 2
2 3
3 4
...
```

Example of a weighted graph:

```bash
1 2 5
2 3 3
3 4 7
...
```

## Output
The script will output the communities detected in the network to the specified output file. Each line corresponds to a community and lists the nodes belonging to that community. Additionally, the modularity of the partition is calculated and appended to the output file.

## Algorithm Details
The Louvain algorithm works by iteratively optimizing the modularity of the network partition. It starts with each node in its own community and iteratively merges communities to maximize modularity. The process continues until no further improvement can be made.

## Example
Suppose you have a network described in the input file network.txt. You can run the Louvain algorithm using the following command:

```bash
python louvain.py network.txt output.txt
```
This will generate an output file output.txt containing the detected communities and their modularity.
