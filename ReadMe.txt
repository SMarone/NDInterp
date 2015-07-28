==================== \\\ NDimInterp.py File Information /// ====================

#\ Author - Stephen Marone /\ Adviser - Justin Gray /\ Started June 8th, 2015 /#


================================== Objective ===================================

This readme is for the N-Dimensional Interpolation program it is accompanying called NDimInterp.py.  This routine will be able to use multiple interpolation methods in order to ascertain values for a user specified set of points (prediction points) that can be related to a previously acquired set of points with known values (training points). 

The program is written in python and takes advantage of the numpy, scipy, and matplotlib libraries throughout it.

================================= Nomenclature =================================

Training Points - Points which are accompanied by known dependent dimension values.

Prediction Points - Points in which the dependent dimension value is to be solved for.

Neighbor - A point which values are the closest or some of the closest that can be found to another specified point.  If set to 0, this number will be set to a default value.

KD-Tree - A sorting routine that enables the storage of points into tree 'leaves'.  Calling these leaves and then the respective point location in the program reduces cost greatly.

Leaves - The number of sections a KD-Tree will cut the training data into. If set to 0, this number will be set to a default value.

Tension - A fine adjustment in the hermite function that can change the 'tightness' of a function fit. Default is 0.

Bias - A fine adjustment in the hermite function that can change the attention to the neighbor greater than or less than the currently analyzed point.  Default is 0.

Distance Effect - The effect of the distance of the neighbors in the weighted nearest neighbors calculaton.  Default is 2.

Comp - This determines the type of function to use for the CR scheme.  -1 and -2 are general function types while 0 through 4 select a function based on the dimensionality of the problem.  Default is 2.

================================ I/O Structures ================================

The structures for the basic setup and query of the interpolation methods is desribed below.  If a section has <> around it, that can be changed to alter the results or rename a variable.

Setup Interpolation Method:

<Interpolator> = <InterpScheme>Interp(<trainpoints>, <trainvals>, <numberofleaves>,**)

	Interpolator - A function that is ready to be queried.
	InterpScheme - The scheme you wish to use (LN, WN, CR, etc.).
	trainpoints - The independent dimension values of the training points in the array shape of (<#oftrainingpoints>,<#ofindependentdimensions>).
	trainvals - The dependent dimension values of the training points in the array shape of (<#oftrainingpoints>,1).
	numberofleaves - The number of leaves to use in the KD-Tree.  Integer value.

** Note that the CR interpolation has optional arguments of the number of neighbors and the comp.

Query Interpolation Method for Dependent Values:

<predvals> = <Interpolator>(<predpoints>, **)

	predvals - The dependent dimension values of the prediction points that were solved for.  They are in the array shape of (<#ofpredictionpoints>,1).
	Interpolator - Described previously.
	predpoints - The independent dimension values of the prediction points that were queried for.  They are in the array shape of (<#ofpredictionpoints>,<#ofindependentdimensions>)

** Note that all interpolators except for LN and CR have optional arguments.  These are listed below and described in the nomenclature section:
WN - number of neighbors, distance effect
CN - number of neighbors
HN - number of neighbors, tension, bias, tight

Query Interpolation Method for Gradient Values:

<predgradients> = <Interpolator>.gradient(<predpoints>, **)

	predgradients - The gradient for each prediction point in the array shape of (<#ofpredictionpoints>,<#ofindependentdimensions>) since the gradient has a value per independent dimension.
	Interpolator - Described previously.
	predpoints - Described previously.

** Note that the gradient call works similarly to the normal call and can take the same optional arguments listed previously per interpolation scheme.

================= Available Test Problem and Distribution Types ================

Data Distributions:
	cart  - Cartesian or uniformly spaced distribution
	rand  - Random value distribution
	LH    - Linear Hypercube distribution

2 Dimensional Problems:
	2D3O  - 3rd order Legendre Polynomial. Recommended range of [-500,500].
	2D5O  - 5th order Legendre Polynomial. Recommended range of [-500,500].

3 Dimensional Problems:
	Plane - Simple plane. Dependency larger on 2nd dimension.
	PW    - Piecewise function of 3 planes. All intersect at [0,0].
	Crate - Egg crate function. Variable Peaks. Recommended range of [-512,512] 

5 Dimensional Problems:
	5D2O  - 2nd order test problem.
	5D4O  - 4th order test problem.

======================== Available Interpolation Schemes =======================

	Linear Nearest Neighbor (LN) - This method determines a plane (or function when only 2 dimensions exist) which exists from a set number of neighbors which contains to the dependent to be solved.  The number of neighbors used is based directly on the dimensions of the problem.  This interpolation is very quick, but does not mold well to high dimensional or complex problems.

	Weight Nearest Neighbor (WN) - This method uses a selected number of found neighbors for each point to be solved and finds the value based on weightings related to the inverse of the nieghbor's distance. The effect of distance on the weighting can also be adjusted.  This interpolation can be costly but may, on average, be my favorite interpolation. 

	Cosine Nearest Neighbor (CN) - This method utilizes the cosine function in order to provide a smooth interpolation result for the point to be solved.  The number of neighbors used can be adjusted.  Low dimensional problems prove difficult for this interpolation.

	Hermite Nearest Neighbor (HN) - This method, similar to the CN method, provides a smooth interpolation but uses the more complex hermite function.  The neighbors can also be adjusted in this method, although it has a minimum of 5 neighbors. This interpolation works well with complex problems and, if it's tight setting is set to True, can work decently with low dimensional problem sets.  It's ability to follow polynomials everywhere except for the range extremes is very notable.

	Compactly Supported Radial Basis Function (CR) - This method, although it technically uses a nearest neighbors method and KD-Tree like the other interpolations, has a different initialization when it receives training points.  Although this setup takes more time, the CR can adjust to complex functions well.  It uses a radial distance and a prescribed function (called a comp, there are multiple options for this) in order to designate a weighting matrix.  This matrix can provide interpolated values for predicted points with a query composed of the prescribed comp and the weighting matrix.

===================================== Tips =====================================

-This code only solves for 1 dependent dimension at a time.

-When the code refers to dimensions, it is equal to the number of independent dimensions + 1 dependent dimension.

-The comp in the CR interpolation specifies the order of the radial basis function to be used.  <-2> uses an 11th order, <-1> uses a 9th order, and any value from <0> to <4> uses an order equal to <floor((dimensions-1)/2) + (3*comp) +1>.
