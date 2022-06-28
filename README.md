# pathfinder
Finds directed paths in directed graphs

To compile simply download the repo and run 
```
make
```

To run use:
```
./pathfinder --vertices vertices.npy --threads t --size n --pathdist pathdist.npy --indices indices.npy --indpty indices.npy --out outfile
```
Where the inputs are:
* vertices: address of numpy file of vertices to be considered, pathfinder will find all paths between pairs of elements from this list, dtype uint32
* threads: number of parallel threads
* size: number of vertices in graph
* pathdist: Address of numpy file containing path distance between the elements to consider, all paths of length pathdist[i*n+j] will be returned from vertices[i] to vertices[j], where n is number of vertices, must be 1 dimension and have length n^2, so can be obtained from an nxn numpy array using numpy flatten, dtype uint32
* indices: Address of numpy file of indices of csr format of graph, dtype uint32
* indptr: Address of numpy file of indptr of csr format of graph, dtype uint32
* out: Output address will save a file for each thread

IMPORTANT! The numpy files must have dtype uint32, if 64bit numpy files are inputted no error is displayed but things will go wrong.

Note graphs vertices are indexed starting at 0.
