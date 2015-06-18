This readme is for the N-Dimensional Interpolation program it is accompanying.  This routine will be able to use multiple interpolation methods in order to ascertain values for a set of points that can be related to a previously acquired set of training points. 

Currently, the interpolations schemes available are:

	Linear Nearest Neighbor (LN) - This method find a plane which exists from a set number of neighbors to the points to be solved.  The number of neighbors used is based directly on the dimensions of the problem. 

	Weight Nearest Neighbor (WN) - This method uses a selected number of found neighbors for each point to be solved and finds the value based on weightings related to nieghbor distance. The effect of distance on the weighting can also be adjusted.

	Cosine Nearest Neighbor (CN) - This method utilizes the cosine function in order to provide a smooth interpolation result for the point to be solved.  The number of neighbors used can be adjusted.

	Hermite Nearest Neighbor (HN) - This method, similar to the CN method, provides a smooth interpolation but uses the more complex hermite function.  The neighbors can also be adjusted in this method although it has a minimum of 5 neighbors. 
