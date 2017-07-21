## Exact Geodesics on Meshes

The project deals with finding exact geodesics on a polygonal meshes. Geodesics usually refer to shortest distances on a surface. We implemented a window propagation algorithm based on the <a href="https://github.com/srajangarg/geodesics/blob/master/surazhsky_geodesics.pdf" target="_blank">paper</a> by V. Surazhsky et al. We also juxtapose the exact geodesic against the shortest path along the mesh edges.<br><br>

The algorithm tries to 'unfold' the polygonal mesh from the source all the way to different destinations. It does so by maintaining windows on each edge, from which the source can be seen using some unfolding. These windows are further propagated (in a Djikstra's fashion) along the mesh by updating and merging with existing windows. Although the worst case time complexity of the algorithm is O(N<sup>2</sup>logN), for all practical purposes it runs in sub-quadratic time.

## Sample Results

Red lines represent the exact geodesic whereas the green path represents the shortest path via the mesh edges.

Simple

<img src="/images/simple.png" width="95%">

Bunny

<img src="/images/bunny.png" width="95%">

Beloved Homer

<img src="/images/homer1.png" width="95%"> 
<img src="/images/homer2.png" width="95%">
