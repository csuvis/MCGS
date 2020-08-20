# Mino-Centric Graph Sampling (MCGS)
The MCGS is a new graph sampling approach focusing on minority structures in graphs. 
It is designed to preserve minority structures in sampled graphs while balancing the
preservation of minority and majority structures and suppressing the generation of
new minority structures.

## Installation
The MCGS is implemented in python3 with dependency packages of networkx and numpy.  
You can install them by following commands:
```
$ pip install networkx
$ pip install numpy
```

## Simple example
Use MCGS module to directly sample an graph:
```
>>> import networkx as nx
>>> from MCGS import MCGS
>>> MCGS_model = MCGS()
>>> G = nx.Graph()
>>> G.add_nodes_from([1, 2, 3, 4, 5, 6, 7])
>>> G.add_edges_from([(1, 2), (1, 3), (1, 4), (1, 5), (5, 6), (5, 7), (6, 7)])
>>> Gs = MCGS_model.run_sampling(G, rate=0.5)
>>> list(Gs.nodes())
[1, 2, 4, 5]
>>> list(Gs.edges())
[(1, 2), (1, 4), (1, 5)]
```

## Tips
There are some several built-in graph data sets provided for users to use and 
visualize analysis, and you can execute the example program called usage_example.py 
by following command:
```
$ python usage_example.py
```

For more information about MCGS and the usage, please see MCGS.py.

## Reference Paper
+ Ying Zhao, Haojin Jiang, Qi'an Chen, Yaqi Qin, Huixuan Xie, Yitao Wu, Shixia Liu,
Zhiguang Zhou, Jiazhi Xia, and Fangfang Zhou. Preserving Minority Structures in
Graph Sampling[J]. *IEEE Transactions on Visualization and Computer Graphics*,
2021 (IEEE VIS 2020 VAST TVCG Track).