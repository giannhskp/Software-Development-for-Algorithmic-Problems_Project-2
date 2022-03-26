# Project 2

Neighbor Search and Clustering for Time-Series and Vectors using Locality-sensitive hashing and Randomized Projection to Hypercube.

Vectors comparison is performed using L2 (Euclidean distance) metric.

Time series comparison is performed using Discrete Frechet or Continuous Frechet metric.

# Part 1 - Neighbor Search
In the first part of the project we implemented 3 different approaches:
  - Time series are represented as vectors in the Euclidean space. The distance between two vectors is computed using L2 (Euclidean distance) metric. The methods used for the Nearest Neighbor Search are: LSH for vectors and Projection to Hypercube (same as Project 1).
  - Time series are represented as a polygonal curves in R^2 space. The distance between two time-series is computed using Discrete Frechet metric (dynamic programming implemention). The method used for the Nearest Neighbor Search is LSH for discrete Time-Series, providing Nearest Neighbor Search, K Nearest Neighbors Search and Range Search. The hashing procedure of LSH algorithm is: snapping to 2D grid and elimination of dupicate snapped points followed by the Classic vector LSH procedures.
  - Time series are represented as a polygonal curves in R space, by abolishing the time dimension. The distance between two time-series is computed using Continuous Frechet metric. The method used for the Nearest Neighbor Search is LSH for continuous Time-Series, providing Nearest Neighbor Search, K Nearest Neighbors Search and Range Search. Filtering is applied to time-series in order to reduce their complexity. The hashing procedure of LSH algorithm is: filtering, snapping to 1D grid, minima-maxima elimination followed by the Classic vector LSH procedures. 
