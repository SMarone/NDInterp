import numpy as np
from  scipy import spatial
from scipy.sparse import csc_matrix
import scipy.sparse.linalg as spsl 
import scipy.interpolate as interp
from matplotlib import cm
import matplotlib.pyplot as plt
import time
import math
from mpl_toolkits.mplot3d import Axes3D
import copy as cp
from matplotlib.backends.backend_pdf import PdfPages
from contextlib import contextmanager
import sys, os

##   Written by Stephen Marone
##     Intern at the NASA GRC  
## Project Start - June 8th, 2015
## Program Version 12

## This is an n-dimensional interpolation code.  More information can be found
# at the end of the code along with the lines which run it.

# Local Methods ================================================================
'''
The Methods and functions in this section of the code are here in order to 
have an easy reference to change the test problems in this code.  This 
entire section can be removed if these interpolation routines are being 
integrated into another code.
'''


@contextmanager
def suppress_stdout():
	## This function just suppresses outputs
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout

def Solve(points, funct):		
	## Library of different problems to test interpolation
	if (funct == 'Crate'):
		## Egg Crate Function
		x = points[:,0]
		y = points[:,1]
		z = -x * np.sin(np.sqrt(abs(x - y - 47))) - \
		    (y+47) * np.sin(np.sqrt(abs(x/2 + y + 47)))
	elif (funct == 'PW'):
		## Piecewise Function with 3 Planes - made up
		x = points[:,0]
		y = points[:,1]
		dtype = type(points[0,0])
		z = np.empty(len(x), dtype=dtype)
		count = 0
		for i in range(len(x)):
			if (x[i] > 0):
				if (y[i] > 0):
					z[count] = x[i] + 8*y[i]
					count +=1
				else:
					z[count] = 5*x[i] + y[i]
					count +=1
			else:
				z[count] = 5*x[i] + 6*y[i]
				count +=1
	elif (funct == 'Plane'):
		## Simple Plane Function, made up
		x = points[:,0]
		y = points[:,1]
		z = 3*x+24*y
	elif (funct == '5D2O'):
		## Completely made up
		v = points[:,0]
		w = points[:,1]
		x = points[:,2]
		y = points[:,3]
		z = 15*(v**2) - w - 4*x + 3*(y*y)
		#z = 15*(v**2) - w - 4*x + np.sin(np.pi*(y*y))
	elif (funct == '5D4O'):
		## More bs stuff I made up
		v = points[:,0]
		w = points[:,1]
		x = points[:,2]
		y = points[:,3]
		z = 9*(v**4) + w*w - 7*x + 2*(y*y)
		#z = 9*(v**4) + w*w - np.sin(np.pi*x) + np.cos(np.pi*(y*y))
	elif (funct == '2D3O'):
		## 3rd Order Legendre Polynomial
		x = points[:,0]/500. ## The 500s ensures it is the normal range 
		z = 0.5*(5*(x**3)-3*x)*500
	elif (funct == '2D5O'):
		## 5rd Order Legendre Polynomial
		x = points[:,0]/500. ## The 500s ensures it is the normal range 
		z = (1/8.)*(63*(x**5) - 70*(x**3) + 15*x)*500
	else:
		print 'Function type not found.' #Srsly, read some instructions
		print
		return
	
	return z


# Local Data Class =============================================================
'''
The Methods and functions in this section of the code are here in order to 
easily setup and test data that might resemble good cases the interpolators
will run.  This entire section is unnecesary when wrapping the interpolation
classes into another code.
'''


class N_Data(object):
	## Main data type with N dimension.  Only needed if you want to use the 
	# example functions, check for error, or plot.
	def __init__(self, mini, maxi, numpts, funct='PW', dtype='rand', dims=3):
		## For funct possibilities, refer to the Solve function above
		## Dtype: cart for cartesian, LH for latin h-cube, rand for random
		self.funct = funct
		self.dtype = dtype
		self.dims = dims
		
		print "Creating Independent Data"
		print "-Function Type:", funct
		print "-Data Distribution Type:", dtype
		print "-Quantity of Points:", numpts
		print "-Dimensions:", dims
		for d in range(dims-1):
			print "-Range of dimension %s: [%s,%s]" % ((d+1), mini[d], maxi[d])
		print
		## Setup independent data for each function via distribution
		if (dtype == 'rand'):
			self.points = np.random.rand(numpts,dims-1) * \
										(maxi-mini) + mini
		
			if (dims == 2):
				self.points = np.sort(self.points, axis=0)
			'''
			if (funct == 'PW'):
				## This ensures points exist on discontinuities for piecewise
				self.points[0,0] = mini/2
				self.points[1,0] = mini
				self.points[2,0] = 0
				self.points[3,0] = maxi/2
				self.points[4,0] = maxi
				self.points[5:10,0] = 0
				self.points[0:5,1] = 0
				self.points[5,1] = mini/2
				self.points[6,1] = mini
				self.points[7,1] = maxi/2
				self.points[8,1] = maxi
			'''
		elif (dtype == 'LH'):
			points = np.zeros((numpts,dims-1), dtype="float")
			for d in range(dims-1):
				## Latin Hypercube it up
				points[:,d] = np.linspace(mini[d],maxi[d],numpts)
				np.random.shuffle(points[:,d])
			self.points = points
			
		elif (dtype == 'cart'):
			if (dims == 3):
				stp = (maxi-mini)/(np.sqrt(numpts)-1)
				x = np.arange(mini[0],maxi[0],stp[0])
				y = np.arange(mini[1],maxi[1],stp[1])
				x, y = np.meshgrid(x, y)
				self.points = np.transpose(np.array([x.flatten(), y.flatten()]))
			elif (dims == 2):
				stp = float(maxi[0]-mini[0])/float(numpts+1)
				x = np.arange((mini[0]+stp),maxi[0],stp)
				self.points = np.transpose(np.array([x.flatten()]))
			else:
				print 'Can not currently do cartesian in any dimension \
						except for 2 and 3. Blame meshgrid.'
				raise SystemExit
		else:
			print 'Distribution type not found.'
			raise SystemExit

	def AssignDep(self, values, gradient=[0]):
		## Takes numpy array of values and stores inside the class variable
		print '**Assigning Dependent Data**'
		print
		self.values = values[:,np.newaxis]
		self.gradient=gradient

	def CreateDep(self):
		## This creates the dependant from the specified function type
		print '**Creating Dependent Data from Problem Library**'
		print
		z = Solve(self.points,self.funct)	
		self.values = z[:,np.newaxis]

	def PlotResults(self, title, pltfile='None'):
		## Plots the value results of a 2 or 3 dimensional problem
		## Ensure points and values have been created
		try:
			self.values
		except NameError:
			print 'Values have not been set.'
			print
			return

		fig = plt.figure()
		if (self.dims == 3):
			x = self.points[:,0]
			y = self.points[:,1]
			ax = fig.gca(projection='3d')
			surf = ax.tricontour(x.flatten(), y.flatten(), 
								 self.values.real.flatten(), 
								 500, rstride=1, cstride=1, 
								 cmap=cm.coolwarm, linewidth=0, 
								 antialiased=False)
			fig.colorbar(surf, shrink=0.5, aspect=5)
			fig.suptitle('%s Results' % title, 
						 fontsize=14, fontweight='bold')
			ax.set_xlabel('X - Independent Variable')
			ax.set_ylabel('Y - Independent Variable')
			ax.set_zlabel('Z - Dependent Variable')
			ax.set_zticklabels([])

			## These adjustments are to improve view of the 3D problems when
			# the maximum and minimum are set to 500 and -500 for each 
			# independent dimension respectively.
			'''
			if (self.funct == 'PW'):
				ax.set_zlim([-5000,5000])
			elif (self.funct == 'Crate'):
				ax.set_zlim([-1000,1000])
			elif (self.funct == 'Plane'):
				ax.set_zlim([-300000,300000])
			'''
		elif (self.dims == 2):
			x = self.points[:,0].flatten()
			y = self.values[:,0].flatten()
			plt.plot(x, y, 'ko-')
			fig.suptitle('%s Results' % title, 
							fontsize=14, fontweight='bold')
			plt.xlabel('Independent Variable')
			plt.ylabel('Dependent Variable')
			plt.xlim(min(x), max(x))
		else:
			print 'Unable to plot this dimension of points.'
			print
			return
		if (pltfile != 'None'):
			pltfile.savefig()

	def PlotPoints(self, title, pltfile='None'):
		## Dim-1 plot of point locations for visualization
		try:
			## Ensure points have been created
			self.points
		except NameError:
			print 'Points have not been set.'
			print
			return

		if (self.dims == 3):
			fig = plt.figure()
			plt.scatter(self.points[:,0], self.points[:,1])
			fig.suptitle('%s Point Locations' % title, 
							fontsize=14, fontweight='bold')
			plt.xlabel('X Value')
			plt.ylabel('Y Value')
		elif (self.dims == 4):
			fig = plt.figure()
			ax = fig.gca(projection='3d')
			ax.scatter(self.points[:,0], self.points[:,1], 
						self.points[:,2], c=c, marker=m)
			fig.suptitle('%s Point Locations' % title, 
							fontsize=14, fontweight='bold')
			ax.set_ylabel('X Value')
			ax.set_ylabel('Y Value')
			ax.set_zlabel('Z Value')
		elif (self.dims == 2):
			zeros = np.zeros((len(self.points[:,0])), dtype='float')
			fig = plt.figure()
			plt.scatter(self.points[:,0], zeros)
			fig.suptitle('%s Point Locations' % title, 
							fontsize=14, fontweight='bold')
			plt.xlabel('X Value')
			plt.ylabel('')
			return
		else:
			print 'Unable to plot this dimension of points.'
			print
			return
		if (pltfile != 'None'):
			pltfile.savefig()

	def FindError(self, title, pltfile='None', plot=True, 
				  check=False, step=0.00000001):
		## Find and plot error of the found values
		try:
			## Ensure points and values have been created
			self.values
		except NameError:
			print 'Values have not been set.'
			return

		vals = self.values.flatten()
		loc = 1+np.arange(len(vals), dtype='float')
		zc = Solve(self.points,self.funct)	
		error = abs((zc-vals)/(zc+0.000000000001)) * 100
		avg = np.average(error)
		
		if check:
			for i in range(len(error)):
				if (error[i] > 100):
					print '**High Percent Error found of', error[i], '**'
					print '-Location:', xc[i], yc[i]
					print '-Calculated Value:', zc[i]
					print '-Found Value:', vals[i] 
		
		print '%s Error Information' % title
		print '-Average Result Percent Error:', avg
		print '-Max Result Percent Error:', np.max(error)

		if plot:
			fig = plt.figure()
			plt.plot(loc, error, 'k')
			fig.suptitle('%s Error Per Point' % title, 
						  fontsize=14, fontweight='bold')
			plt.xlabel('Point Identifier')
			plt.ylabel('Percent Error')
			plt.ylim((0,((avg+0.00000000001)*2)))
			if (pltfile != 'None'):
				pltfile.savefig()

		self.error = error
		self.actualvalues = zc
		if (np.all(self.gradient) != 0):
			## Use complex step to check nearby points with found gradient
			g = np.empty((len(vals),self.dims-1), dtype='float')
			points = np.zeros((len(vals),self.dims-1), dtype='complex')
			points += self.points

			for D in range(self.dims-1):
				points[:,D] += (step*1j)
				zs = Solve(points,self.funct)
				g[:,D] = zs.imag/ step
				points[:,D] = points[:,D].real

			gerror = abs((g-self.gradient)/(g+0.00000000000000001)) * 100.
			gerror = np.sum(gerror, axis=1)/2
			gavg = np.average(gerror)

			if check:
				for i in range(len(error)):
					if (gerror[i] > 100):
						print '**High Percent Error found of', gerror[i], '**'
						print '-Location:', points[i,:]
						print '-Calculated Value:', g[i,:]
						print '-Found Value:', self.gradient[i,:] 
			
			print '-Average Actual Gradient Percent Error:', gavg
			print '-Max Actual Gradient Percent Error:', np.max(gerror)
	
			if plot:
				if (self.dims == 2):
					col = 1
					rows = 1
				else:
					col = 2
					rows = int(self.dims/2)
				gbfig, gbax = plt.subplots(col,rows)
				gbfig.suptitle('%s Gradient wrt Independent' % title, 
							  fontsize=14, fontweight='bold')
				gbax = np.reshape(gbax,((self.dims-1)))
				for D in range(self.dims-1):
					ming = min(self.gradient[:,D])
					ming -= 0.05*abs(ming)
					maxg = max(self.gradient[:,D])
					maxg += 0.05*abs(ming)
					gbax[D].plot(points[:,D].real,g[:,D],'bo',
								 points[:,D].real,self.gradient[:,D],'rx')
					gbax[D].set_xlabel('Independent %s Value' % (D+1))
					gbax[D].set_ylabel('Gradient Value')
					gbax[D].set_ylim([ming,maxg])
				gbax[D].legend(['Complex Step in Function','Interpolation'],
						loc=4, prop={'size':8})
				if (pltfile != 'None'):
					pltfile.savefig()
				gfig = plt.figure()
				plt.plot(loc, gerror, 'k')
				gfig.suptitle('%s Actual Gradient Error Per Point' % title, 
							  fontsize=14, fontweight='bold')
				plt.xlabel('Point Identifier')
				plt.ylabel('Percent Error')
				plt.ylim((0,((gavg+0.0000000000001)*2)))
				if (pltfile != 'None'):
					pltfile.savefig()
	
			self.gaerror = gerror
		print

	def CheckGrad(self, title, Interp, itype, pltfile='None', plot=True,
				 check=False, step=0.00000001, N=5, DistEff=3, 
				 tension=0, bias=0, tight=False):
		loc = 1+np.arange(len(self.values), dtype='float')
		## Find and plot error of the found values
		
		if (np.all(self.gradient) == 0):
			## If the gradient is empty, just leave.  Ain't got time for that.
			print 'Gradient does not have a solution.'
			print
			return

		PrdiPts = np.empty((len(self.points[:,0]),self.dims-1), dtype='complex')
		PrdiPts.real = self.points
		gradi = np.empty((len(self.values),self.dims-1), dtype='float')
		vals = self.values.reshape((len(self.values)))
		for D in range(self.dims-1):
			PrdiPts[:,:].imag = 0.
			PrdiPts[:,D] = self.points[:,D] + (step*1.0j)
			with suppress_stdout():
				## Redo interpolation with complex addition but suppress output
				if (itype == 'LN'):
					zi = Interp(PrdiPts)
				elif (itype == 'WN'):
					zi = Interp(PrdiPts, N, DistEff)
				elif (itype == 'CN'):
					zi = Interp(PrdiPts, N)
				elif (itype == 'HN'):
					zi = Interp(PrdiPts, N, tension, bias, tight)
				elif (itype == 'CR'):
					zi = Interp(PrdiPts)
			gradi[:,D] = zi.imag/ step
		## Compare the difference in complex addition with the gradient answer
		gerror = abs((gradi-self.gradient)/(gradi+0.000000000000000001)) * 100.
		gerror = np.sum(gerror, axis=1)/(self.dims-1)
		gavg = np.average(gerror)

		if check:
			for i in range(len(error)):
				if (gerror[i] > 100):
					print '**High Percent Error found of', gerror[i], '**'
					print '-Location:', PrdiPts[i,:]
					print '-Calculated Value:', gradi[i,:]
					print '-Found Value:', self.gradient[i,:] 
			
		print '%s Gradient Check' % title
		print '-Average Gradient Percent Error:', gavg
		print '-Max Gradient Percent Error:', np.max(gerror)
		print

		if plot:
			if (self.dims == 2):
				col = 1
				rows = 1
			else:
				col = 2
				rows = int(self.dims/2)
			gbfig, gbax = plt.subplots(col,rows)
			gbfig.suptitle('%s Gradient wrt Independent' % title, 
						  fontsize=14, fontweight='bold')
			gbax = np.reshape(gbax,((self.dims-1)))
			for D in range(self.dims-1):
				ming = min(self.gradient[:,D])
				ming -= 0.05*abs(ming)
				maxg = max(self.gradient[:,D])
				maxg += 0.05*abs(ming)
				gbax[D].plot(PrdiPts[:,D].real,gradi[:,D],'bo',
							 PrdiPts[:,D].real,self.gradient[:,D],'rx')
				gbax[D].set_xlabel('Independent %s Value' % (D+1))
				gbax[D].set_ylabel('Gradient Value')
				gbax[D].set_ylim([ming,maxg])
			gbax[D].legend(['Complex Step in Interp.','Interpolation'], 
					loc=4, prop={'size':8})
			if (pltfile != 'None'):
				pltfile.savefig()
			gfig = plt.figure()
			plt.plot(loc, gerror, 'k')
			gfig.suptitle('%s Gradient Error Per Point' % title, 
						  fontsize=14, fontweight='bold')
			plt.xlabel('Point Identifier')
			plt.ylabel('Percent Error')
			plt.ylim((0,((gavg+0.0000000000001)*2)))
			if (pltfile != 'None'):
				pltfile.savefig()
	
		self.gcerror = gerror
		
	def PlotAll(self, title, Interp, itype, pltfile='None', erplot=True, 
			   check=False, ptplot=False, step=0.00000001, neighs=5, 
			   DistEff=3, tension=0, bias=0, tight=False):
		## This routine runs through all checks and plotting options/
		self.PlotResults(title,pltfile)
		if ptplot:
			self.PlotPoints(title,pltfile)
		self.FindError(title,pltfile,erplot,check,step)
		self.CheckGrad(title,Interp,itype,pltfile,erplot,check,step,
					   neighs,DistEff,tension,bias,tight)


# Interpolation Classes ========================================================
'''
The classes inside this section are the main parts of this code.  Below are
the interpolation methods all individually defined and ready to use.  Most
of the classes are children the the first class below, NNInterpBase. All of
the methods use a KD-Tree to sort the data; this improves the cost of each
query done with them.  To query for dependent answers, call methods are set.
To query for gradients, self.gradient methods are set.
'''


class NNInterpBase(object):
	## The shared aspects of all nearest neighbor interpolators
	def __init__(self, TrnPts, TrnVals, NumLeaves=2):
		## TrainPts and TrainVals are the known points and their
		# respective values which will be interpolated against.		
		## Grab the mins and normals of each dimension
		self.tpm = np.amin(TrnPts, axis=0)
		self.tpr = (np.amax(TrnPts, axis=0) - self.tpm)
		self.tvm = np.amin(TrnVals, axis=0)
		self.tvr = (np.amax(TrnVals, axis=0) - self.tvm)

		## This prevents against colinear data (range = 0)
		self.tpr[self.tpr == 0] = 1
		self.tvr[self.tvr == 0] = 1

		## Normalize all points
		self.tp = (TrnPts - self.tpm) / self.tpr
		self.tv = (TrnVals - self.tvm) / self.tvr

		## Set important variables
		self.dims = len(TrnPts[0,:]) + 1
		self.ntpts = len(TrnPts[:,0])

		'''
		print 'Interpolation Input Values'
		print '-KDTree Leaves:', NumLeaves
		print '-Number of Training Points:', self.ntpts
		print
		'''

		## Make training data into a Tree
		leavesz = math.ceil(self.ntpts/(NumLeaves+0.00000000000001))
		KData = spatial.cKDTree(self.tp,leafsize=leavesz)
		self.KData = KData

class LNInterp(NNInterpBase):
	def main(self, PrdPts, nppts, dims, nloc):
			## Extra Inputs for Finding the normal are found below
			## Number of row vectors needed always dimensions - 1
			r = np.arange(dims-1, dtype="int")
			
			## Diagonal counters are found here to speed up finding each normal
			da = np.zeros((dims, dims-1), dtype="int")
			db = np.zeros((dims, dims-1), dtype="int")
			for i in xrange(dims):
			    da[i] = (r+1+i) % dims
			    db[i] = (i-r-1) % dims
			## Planar vectors need both dep and ind dimensions
			trnd = np.concatenate((self.tp[nloc,:] ,
			                       self.tv[nloc].reshape(nppts,dims,1)),
			                       axis=2)
			## Set the prddata to a temporary array so it has the same
			# dimensions as the training data with the dependant as 0 
			prdtemp = np.concatenate((PrdPts, \
									np.zeros((nppts,1),dtype='float')), 1)
			
			nvect = np.empty((nppts,(dims-1),dims), dtype='float')
			for n in range(dims-1): ## Go through each neighbor
			    ## Creates array[neighbor, dimension] from NN results
			    nvect[:,n,:] = trnd[:,(n+1),:] - trnd[:,n,:]
			
			## Nvec used below to finish a hodge star product 
			# (in 3D essentially a cross product).
			normal = np.prod(nvect[:,r,da], axis=2) - \
					 np.prod(nvect[:,r,db], axis=2)
				
			## Use the point of the closest neighbor to 
			# solve for pc - the constant of the n-dimensional plane.
			pc = np.einsum('ij, ij->i',trnd[:,0,:],normal)
			return normal, prdtemp, pc

	def main2D(self, PrdPts, nppts, nloc):
			## Need to find a tangent instead of a normal, y=mx+b
			d = self.tp[nloc[:,1],0] - self.tp[nloc[:,0],0]
			d[abs(d) < 0.0000000000000001] = 100000000000000000000
			m = (self.tv[nloc[:,1],0] - self.tv[nloc[:,0],0]) / d
			b = self.tv[nloc[:,0],0] - (m * self.tp[nloc[:,0],0])	
			return m, b

	def __call__(self, PredPoints):
		## This method uses linear interpolation by defining a plane with
		# a set number of nearest neighbors to the predicted
		PrdPts = (PredPoints - self.tpm) / self.tpr
		nppts    = len(PrdPts[:,0])
		## Linear interp only uses as many neighbors as it has dimensions
		dims = self.dims
		## Find them neigbors.  Linear only uses #neighs=#dims
		## KData query takes (data, #ofneighbors) to determine closest 
		# training points to predicted data
		ndist, nloc = self.KData.query(PrdPts.real ,self.dims)
		'''
		print 'Linear Plane Nearest Neighbors (LN) KDTree Results'
		print '-Number of Predicted Points:', nppts
		print '-Number of Neighbors per Predicted Point:', dims
		print '-Nearest Neighbor Distance:', np.min(ndist)
		print '-Farthest Neighbor Distance:', np.max(ndist)
		print
		'''
	
		## Need to ensure there are enough dimensions to find the normal with
		if (self.dims > 2):
			normal, prdtemp, pc = self.main(PrdPts, nppts, self.dims, nloc)
			## Set all prdz from values on plane
			prdz = (np.einsum('ij, ij->i',prdtemp,normal)-pc)
			## Check to see if there are any colinear points and replace them
			n0 = np.where(normal[:,-1] == 0)
			prdz[n0] = self.tv[nloc[n0,0]]
			## Finish computation for the good normals
			n = np.where(normal[:,-1] != 0)
			prdz[n] = prdz[n]/-normal[n,-1]

		else:
			m, b = self.main2D(PrdPts, nppts, nloc)
			prdz = (m * PrdPts[:,0]) + b
		
		predz = (prdz * self.tvr) + self.tvm
		return predz

	def gradient(self, PredPoints):
		## Extra method to find the gradient at each location of a set of 
		# supplied predicted points.
		PrdPts = (PredPoints - self.tpm) / self.tpr
		nppts    = len(PrdPts[:,0])
		gradient = np.zeros((nppts, self.dims-1), dtype="float")
		## Linear interp only uses as many neighbors as it has dimensions
		dims = self.dims
		## Find them neigbors
		ndist, nloc = self.KData.query(PrdPts.real ,self.dims)
		## Need to ensure there are enough dimensions to find the normal with
		if (self.dims > 2):
			normal, prdtemp, pc = self.main(PrdPts, nppts, self.dims, nloc)
			if (np.any(normal[:,-1]) == 0):
				return gradient
			gradient = -normal[:,:-1]/normal[:,-1:]

		else:
			gradient[:,0], b = self.main2D(PrdPts, nppts, nloc)

		grad = gradient * (self.tvr/self.tpr)
		return grad

class WNInterp(NNInterpBase):
	## Weighted Neighbor Interpolation
	def __call__(self, PredPoints, N=5, DistEff=3):
		PrdPts = (PredPoints - self.tpm) / self.tpr
		nppts    = len(PrdPts[:,0])
		## Find them neigbors
		## KData query takes (data, #ofneighbors) to determine closest 
		# training points to predicted data
		ndist, nloc = self.KData.query(PrdPts.real ,N)
		'''
		print 'Weighted Nearest Neighbors (WN) KDTree Results'
		print '-Number of Predicted Points:', nppts
		print '-Number of Neighbors per Predicted Point:', N
		print '-Nearest Neighbor Distance:', np.min(ndist)
		print '-Farthest Neighbor Distance:', np.max(ndist)
		print
		'''

		## Setup problem
		vals = 0.0
		pdims = self.dims - 1
		dimdiff = np.subtract(PrdPts.reshape(nppts,1,pdims), self.tp[nloc,:])

		for D in range(pdims):
			if (PrdPts[0,D].imag > 0):
				## KD Tree ignores imaginary part, muse redo ndist if complex 
				ndist = np.sqrt(np.sum((dimdiff**2), axis=2))
				vals += 0.0j
		## Find the weighted neighbors per defined formula for distance effect
		part = ndist**DistEff
		sd = np.sum(1.0/part, axis=1)
		vals += self.tv[nloc] 
		wt = np.sum(vals[:,:,0]/part, axis=1)  
		predz = ((wt/sd) * self.tvr) + self.tvm

		return predz

	def gradient(self, PredPoints, N=5, DistEff=3):
		PrdPts = (PredPoints - self.tpm) / self.tpr
		nppts    = len(PrdPts[:,0])
		## Find them neigbors
		ndist, nloc = self.KData.query(PrdPts.real ,N)
		## Setup problem
		vals = 0.0
		pdims = self.dims - 1
		dimdiff = np.subtract(PrdPts.reshape(nppts,1,pdims), self.tp[nloc,:])
		## Find the weighted neighbors per defined formula for distance effect
		part = ndist**DistEff
		sd = np.sum(1.0/part, axis=1)
		vals += self.tv[nloc] 
		## Calculate necessaties for gradient
		dsd = np.sum((dimdiff/(ndist**(DistEff+2)).reshape(nppts,N,1)),
						axis=1).reshape(nppts,1,pdims)
		## Reshape a couple things for matrix operations	
		sumdist = sd.reshape(nppts,1,1)
		ndist = ndist.reshape(nppts,N,1)
		## Solve for gradient
		const = (DistEff/(sd*sd)).reshape(nppts,1)
		gradient = (const*np.sum((vals/(ndist**DistEff)) * \
				(dsd-(sumdist*dimdiff/(ndist*ndist))), axis=1))
		grad = gradient * (self.tvr/self.tpr)

		return grad

class CNInterp(NNInterpBase):
	## Cosine Neighbor Interpolation
	def main(self, PrdPts, nppts, D, N, neighvals):	
		orgneighs = np.empty((nppts,N,self.dims), dtype='float')
		anppts = np.arange(nppts)
		for t in range(nppts): 
			## I want this gone but can't figure out sorting
			orgneighs[t,:,:] = neighvals[t,neighvals[t,:,D].argsort(),:]
		## Make some temporary variables
		tprd = PrdPts.reshape(nppts,1,(self.dims-1))
		podiff = np.subtract(tprd[:,:,D], orgneighs[:,:,D])+0.00000000000001
		## Gives the indices of the closest neighbors to each prdpt per dim.
		cnd = np.argmin((np.abs(podiff)), axis=1)
		## If the closest neighbor per dim is a smaller value, add 1.
		cnd += np.ceil((podiff[anppts,cnd].real) / \
				(np.abs(podiff[anppts,cnd].real*2)))
		## Stay in the range!
		cnd[cnd == 0] = 1
		cnd[cnd == N] = N-1

		## Mark upper and lower closest neighbors
		zu = orgneighs[anppts,cnd,-1]
		zl = orgneighs[anppts,(cnd-1),-1]
		## Need to ensure no colinear data causes a divide by 0
		ddiff = (orgneighs[anppts,cnd,D]-orgneighs[anppts,(cnd-1),D])
		ddiff[ddiff < 0.00000000001] = podiff[:,(cnd-1)].real*2.
		ddiff = 1./ddiff
		diff = podiff[anppts,(cnd-1)] * ddiff

		return zu, zl, diff, ddiff

	def __call__(self, PredPoints, N=5):
		PrdPts = (PredPoints - self.tpm) / self.tpr
		nppts    = len(PrdPts[:,0])
		## Find them neigbors
		## KData query takes (data, #ofneighbors) to determine closest 
		# training points to predicted data
		ndist, nloc = self.KData.query(PrdPts.real ,N)
		
		'''
		print 'Cosine Nearest Neighbors (CN) KDTree Results'
		print '-Number of Predicted Points:', nppts
		print '-Number of Neighbors per Predicted Point:', N
		print '-Nearest Neighbor Distance:', np.min(ndist)
		print '-Farthest Neighbor Distance:', np.max(ndist)
		print
		'''
		## This sorts the rows by the values in the column specified. 
		## The neighvals dim is [neigh#, dimension] 
		neighvals = np.concatenate((self.tp[nloc,:] ,
									self.tv[nloc,:]),
		 							axis=2)
		tprdz = np.zeros((nppts), dtype='complex')
		for D in range(self.dims-1):
			zu, zl, diff, ddiff = self.main(PrdPts, nppts, D, N, neighvals)		
			## Take info from main function and insert into the cosine funct.
			mu2 = (1-np.cos(np.pi*diff))/2
			tprdz += zu * mu2 + zl * (1-mu2)

		## Average found solutions and remove normalization
		predz = ((tprdz/(D+1)) * self.tvr) + self.tvm
		return predz

	def gradient(self, PredPoints, N=5):
		PrdPts = (PredPoints - self.tpm) / self.tpr
		nppts    = len(PrdPts[:,0])
		gradient = np.zeros((nppts, self.dims-1), dtype="float")
		## Find them neigbors
		ndist, nloc = self.KData.query(PrdPts.real ,N)
		## This sorts the rows by the values in the column specified. 
		## The neighvals dim is [neigh#, dimension] 
		neighvals = np.concatenate((self.tp[nloc,:] ,
									self.tv[nloc,:]),
		 							axis=2)

		for D in range(self.dims-1):
			zu, zl, diff, ddiff = self.main(PrdPts, nppts, D, N, neighvals)		
			## Take info from main function and solve for grad of cosine funct.
			gradient[:,D] = (np.pi/2)*(zu-zl)*ddiff* \
					np.sin(np.pi*diff.real)/(self.dims-1) 

		## Remove normalization
		grad = gradient * (self.tvr/self.tpr)
		return grad


class HNInterp(NNInterpBase):
	## Hermite Neighbor Interpolation
	def HermFunctRes(self, y, mu, tension, bias):
		## Should  have 4 neighbor values (y) and percent distance along line 
		# from y1 and y2 to get the wanted value (mu).  Bias and tension 
		# are optional values.  This one returns the dependent values.
		mu2 = mu * mu
		mu3 = mu2 * mu
	
		m0  = (y[:,1]-y[:,0])*(1+bias)*(1-tension)/2
		m0 += (y[:,2]-y[:,1])*(1-bias)*(1-tension)/2
		m1  = (y[:,2]-y[:,1])*(1+bias)*(1-tension)/2
		m1 += (y[:,3]-y[:,2])*(1-bias)*(1-tension)/2

		a0 =  2*mu3 - 3*mu2 + 1
		a1 =    mu3 - 2*mu2 + mu
		a2 =    mu3 -   mu2
		a3 = -2*mu3 + 3*mu2

		return (a0*y[:,1] + a1*m0 + a2*m1 + a3*y[:,2])
		
	def HermFunctGrad(self, y, mu, dmu, tension, bias):
		## Should  have 4 neighbor values (y) and percent distance along line 
		# from y1 and y2 to get the wanted value (mu).  Bias and tension 
		# are optional values. This one returns the gradient values.
		mu2 = mu * mu
	
		m0  = (y[:,1]-y[:,0])*(1+bias)*(1-tension)/2
		m0 += (y[:,2]-y[:,1])*(1-bias)*(1-tension)/2
		m1  = (y[:,2]-y[:,1])*(1+bias)*(1-tension)/2
		m1 += (y[:,3]-y[:,2])*(1-bias)*(1-tension)/2
		
		## Derivatives of the 'a's
		b0 =  6*mu2 - 6*mu
		b1 =  3*mu2 - 4*mu + 1
		b2 =  3*mu2 - 2*mu
		b3 = -6*mu2 + 6*mu

		return ((dmu*(b0*y[:,1]+b1*m0+b2*m1+b3*y[:,2]))/(self.dims-1)).real

	def main(self, PrdPts, nppts, D, N, u, l, neighvals):
		orgneighs = np.empty((nppts,N,self.dims), dtype='float')
		anppts = np.arange(nppts)
		y = np.empty((nppts,4), dtype='float')
		for t in range(nppts): 
			## I want this gone but can't figure out sorting
			orgneighs[t,:,:] = neighvals[t,neighvals[t,:,D].argsort(),:]
		## Make some temporary variables
		tprd = PrdPts.reshape(nppts,1,(self.dims-1))
		podiff = np.subtract(tprd[:,:,D], orgneighs[:,:,D])+0.000000000001
		## Gives the indices of the closest neighbors to each prdpt per dim.
		cnd = np.argmin((np.abs(podiff)), axis=1)
		## If the closest neighbor per dim is a smaller value, add 1.
		cnd += np.ceil(podiff[anppts,cnd].real / \
				np.abs((podiff[anppts,cnd].real)*2))
		## Stay in the range!
		cnd[cnd <= 1] = 2
		cnd[cnd >= (N-1)] = N-2
		
		## Find location of value in min and max of neighvals to be used
		ddiff = (orgneighs[anppts,(cnd+u),D]-orgneighs[anppts,(cnd-l),D])
		## Avoid colinear causing /0
		ddiff[ddiff < 0.000000001] = podiff[:,(cnd-l)].real*2.
		ddiff = 1./ddiff
		diff = podiff[anppts,(cnd-l)] * ddiff 

		for n in range(4):
			## Only need 4 inputs.  Would like to remove this sometime
			y[:,n] = orgneighs[anppts,(cnd-2+n),-1]

		return y, diff, ddiff

	def __call__(self, PredPoints, N=5, tension=0, bias=0, tight=False):
		## Traindata has n dimensions, prddata has n-1 dimensions because last
		# dimension is the dependent one we are solving for.  
		## N neighbors will be found but only 4 will be used per dimension
		PrdPts = (PredPoints - self.tpm) / self.tpr
		nppts = len(PrdPts[:,0])
		tprdz = np.zeros((nppts), dtype='complex')
		u = 1
		l = 2
		if tight:
			## Added this because tight works much better with smaller problems
			# but very badly with larger problems
			u -= 1
			l -= 1
		if (N < 5):
			## Hermite requires at least 5 neighbors
			N = 5
		## Find them neigbors
		## KData query takes (data, #ofneighbors) to determine closest 
		# training points to predicted data
		ndist, nloc = self.KData.query(PrdPts.real ,N)

		'''
		print 'Hermite Nearest Neighbors (HN) KDTree Results'
		print '-Number of Predicted Points:', nppts
		print '-Number of Neighbors per Predicted Point:', N
		print '-Nearest Neighbor Distance:', np.min(ndist)
		print '-Farthest Neighbor Distance:', np.max(ndist)
		print
		'''
		## This sorts the rows by the values in the column specified. 
		## The neighvals dim is [neigh#, dimension] 
		neighvals = np.concatenate((self.tp[nloc,:] ,
									self.tv[nloc,:]),
		 							axis=2)

		for D in range(self.dims-1):
			y, diff, ddiff = self.main(PrdPts, nppts, D, N, u, l, neighvals)
			## Find necessary values and put into the Hermite function
			tprdz += self.HermFunctRes(y,diff,tension,bias)

		## Average found solutions and remove normalization
		predz = ((tprdz/(D+1)) * self.tvr) + self.tvm
		return predz

	def gradient(self, PredPoints, N=5, tension=0, bias=0, tight=False):
		## Traindata has n dimensions, prddata has n-1 dimensions because last
		# dimension is the dependent one we are solving for.  
		## N neighbors will be found but only 4 will be used per dimension
		PrdPts = (PredPoints - self.tpm) / self.tpr
		nppts = len(PrdPts[:,0])
		gradient = np.zeros((nppts, self.dims-1), dtype="float")
		u = 1
		l = 2
		if tight:
			## Added this because tight works much better with smaller problems
			# but very badly with larger problems
			u -= 1
			l -= 1
		if (N < 5):
			## Hermite requires at least 5 neighbors
			N = 5
		## Find them neigbors
		ndist, nloc = self.KData.query(PrdPts.real ,N)
		## This sorts the rows by the values in the column specified. 
		## The neighvals dim is [neigh#, dimension] 
		neighvals = np.concatenate((self.tp[nloc,:] ,
									self.tv[nloc,:]),
		 							axis=2)

		for D in range(self.dims-1):
			y, diff, ddiff = self.main(PrdPts, nppts, D, N, u, l, neighvals)
			## Find necessary variables and throw into the Hermite funct.
			gradient[:,D] =	self.HermFunctGrad(y,diff,ddiff,tension,bias)

		## Remove normalization
		grad = gradient * (self.tvr/self.tpr)
		return grad

class CRInterp(object):
	## Compactly Supported Radial Basis Function
	def FindR(self, npp, T, loc): 
		R = np.zeros((npp, self.ntpts), dtype="complex")
		## Choose type of CRBF R matrix
		if (self.comp == -1):
			## Comp #1 - a 
			Cf = (1.-T)**5
			Cb  = (8.+(40.*T)+(48.*T*T)+ \
					  (72.*T*T*T)+(5.*T*T*T*T))
		elif (self.comp == -2):
			## Comp #2
			Cf = (1.-T)**6.
			Cb = (6.+(36.*T)+(82.*T*T)+(72.*T*T*T)+ \
						(30.*T*T*T*T)+(5.*T*T*T*T*T))
		elif (self.comp == -3):
			## Comp #3
			Cf = T**0.
			Cb = np.sqrt((T*T)+1.)
		else:
			## The above options did not specify a dimensional requirement
			# in the paper found but the rest are said to only be guaranteed
			# as positive definite iff the dimensional requirements are.
			# Because of this, the user can select 0 through 4 to adjust to a 
			# level of order trying to be attained.
			if (self.dims <= 2):
				if (self.comp == 0):
					# This starts the dk comps, here d=1, k=0
					Cf = 1.-T
					Cb = np.ones((len(T[:,0]),len(T[0,:])), dtype="complex")
				elif (self.comp == 1):
					Cf = (1.-T)*(1.-T)*(1.-T)/12.
					Cb = 1.+(3.*T)
				elif (self.comp == 2):
					Cf = ((1.-T)**5)/840.
					Cb = 3.+(15.*T)+(24.*T*T)
				elif (self.comp == 3):
					Cf = ((1.-T)**7)/151200.
					Cb = 15.+(105.*T)+(285.*T*T)+(315.*T*T*T)
				elif (self.comp == 4):	
					Cf = ((1.-T)**9)/51891840.
					Cb = 105.+(945.*T)+(3555.*T*T)+(6795.*T*T*T)+\
							(5760.*T*T*T*T)
			elif (self.dims <= 4):
				if (self.comp == 0):
					Cf = (1.-T)
					Cb = (1.-T)
				elif (self.comp == 1):
					Cf = ((1.-T)**4)/20.
					Cb = 1.+(4.*T)
				elif (self.comp == 2):
					Cf = ((1.-T)**6)/1680.
					Cb = 3.+(18.*T)+(35.*T*T)
				elif (self.comp == 3):
					Cf = ((1.-T)**8)/332640.
					Cb = 15.+(120.*T)+(375.*T*T)+(480.*T*T*T)
				elif (self.comp == 4):		
					Cf = ((1.-T)**10)/121080960.
					Cb = (105.+(1050.*T)+(4410.*T*T)+(9450.*T*T*T)+\
							(9009.*T*T*T*T))
			elif (self.dims <= 6):
				if (self.comp == 0):
					Cf = (1.-T)*(1.-T)
					Cb = (1.-T)
				elif (self.comp == 1):
					Cf = ((1.-T)**5)/30.
					Cb = 1.+(5.*T)
				elif (self.comp == 2):
					Cf = ((1.-T)**7)/3024.
					Cb = 3.+(21.*T)+(48.*T*T)
				elif (self.comp == 3):
					Cf = ((1.-T)**9)/665280.
					Cb = 15.+(135.*T)+(477.*T*T)+(693.*T*T*T)
				elif (self.comp == 4):		
					Cf = ((1.-T)**11)/259459200.
					Cb = (105.+(1155.*T)+(5355.*T*T)+(12705.*T*T*T)+\
							(13440.*T*T*T*T))
			else:
				## Although not listed, this is ideally for 8 dim or less
				if (self.comp == 0):
					Cf = (1.-T)*(1.-T)
					Cb = (1.-T)*(1.-T)
				elif (self.comp == 1):
					Cf = ((1.-T)**6)/42.
					Cb = 1.+(6.*T)
				elif (self.comp == 2):
					Cf = ((1.-T)**8)/5040.
					Cb = 3.+(24.*T)+(63.*T*T)
				elif (self.comp == 3):
					Cf = ((1.-T)**10)/1235520.
					Cb = 15.+(150.*T)+(591.*T*T)+(960.*T*T*T)
				elif (self.comp == 4):
					Cf = ((1.-T)**12)/518918400.
					Cb = (105.+(1260.*T)+(6390.*T*T)+(16620.*T*T*T)+\
							(19305.*T*T*T*T))

		for i in xrange(npp):
			R[i,loc[i,:-1]] = Cf[i,:] * Cb[i,:]
		
		return R
	
	def FinddR(self, PrdPts, ploc, pdist):
		T  = (pdist[:,:-1]/pdist[:,-1:])
		## Solve for the gradient analytically
		## The first quantity needed is dRp/dt
		if (self.comp == -1):
			frnt = (1.-T)**4
			dRp  = frnt*((-5.*(8. + (40.*T) + (48.*T*T) + \
								(72.*T*T*T) + (5.*T*T*T*T))) + \
						((1.-T) * (40. + (96.*T) + \
								(216.*T*T) + (20.*T*T*T))))
		elif (self.comp == -2):
			frnt = (1.-T)**5
			dRp  = frnt*((-6. * (6. + (36.*T) + \
								(82.*T*T) + (72.*T*T*T) + \
								(30.*T*T*T*T) + (5.*T*T*T*T*T))) + \
						((1.-T)) * (36. + (164.*T) + \
								(216.*T*T) + (120.*T*T*T) + \
								(25.*T*T*T*T)))
		elif (self.comp == -3):
			dRp = T/np.sqrt((T*T)*1.)
		else:
			## Start dim dependent comps, review first occurance for more info
			if (self.dims <= 2):
				if (self.comp == 0):
					# This starts the dk comps(Wendland Functs), here d=1, k=0
					dRp  = -1.
				elif (self.comp == 1):
					dRp  = -T*(1.-T)*(1.-T)
				elif (self.comp == 2):
					frnt = ((1.-T)**4)/-20.
					dRp  = frnt*(T+(4.*T*T))
				elif (self.comp == 3):
					frnt = ((1.-T)**6)/-1680.
					dRp  = frnt*((3.*T)+(18.*T*T)+(35.*T*T*T))
				elif (self.comp == 4):		
					frnt = ((1.-T)**8)/-22176.
					dRp  = frnt*(T+(8.*T*T)+(25.*T*T*T)+(32.*T*T*T*T))
			elif (self.dims <= 4):
				if (self.comp == 0):
					dRp  = -2.*(1-T)
				elif (self.comp == 1):
					dRp  = -T*(1.-T)*(1.-T)*(1.-T)
				elif (self.comp == 2):
					frnt = ((1.-T)**5)/-30.
					dRp  = frnt*(T+(5.*T*T))
				elif (self.comp == 3):
					frnt = ((1.-T)**7)/-1008.
					dRp  = frnt*(T+(7.*T*T)+(16.*T*T*T))
				elif (self.comp == 4):		
					frnt = ((1.-T)**9)/-221760.
					dRp  = frnt*((5.*T)+(45.*T*T)+(159.*T*T*T)+(231.*T*T*T*T))
			elif (self.dims <= 6):
				if (self.comp == 0):
					dRp  = -3.*(1.-T)*(1.-T)	
				elif (self.comp == 1):
					dRp  = -T*((1.-T)**4)
				elif (self.comp == 2):
					frnt = ((1.-T)**6)/-42.
					dRp  = frnt*(T+(6.*T*T))
				elif (self.comp == 3):
					frnt = ((1.-T)**8)/-1680.
					dRp  = frnt*(T+(8.*T*T)+(21.*T*T*T))
				elif (self.comp == 4):		
					frnt = ((1.-T)**10)/-411840.
					dRp  = frnt*((5.*T)+(50.*T*T)+(197.*T*T*T)+(320.*T*T*T*T))
			else:
				## Although not listed, this is ideally for 8 dim or less
				if (self.comp == 0):
					dRp  = -4.*(1.-T)*(1.-T)*(1.-T)
				elif (self.comp == 1):
					dRp  = -T*((1.-T)**5)
				elif (self.comp == 2):
					frnt = ((1.-T)**7)/-56.
					dRp  = frnt*(T+(7.*T*T))
				elif (self.comp == 3):
					frnt = ((1.-T)**9)/-7920.
					dRp  = frnt*((3.*T)+(27.*T*T)+(80.*T*T*T))
				elif (self.comp == 4):
					frnt = ((1.-T)**11)/-720720.
					dRp  = frnt*((5.*T)+(55.*T*T)+(239.*T*T*T)+(429.*T*T*T*T))
		## Now need dt/dx
		xpi = np.subtract(PrdPts, self.tp[ploc[:,:-1],:])
		xpm = PrdPts - self.tp[ploc[:,-1:],:]
		dtx = (xpi-(T*T*xpm))/ \
				(pdist[:,-1:,:]*pdist[:,-1:,:]*T)
		## The gradient then is the summation across neighs of w*df/dt*dt/dx
		return np.sum((dRp * dtx * self.weights[ploc[:,:-1]]), axis=1)

	def __init__(self, TrnPts, TrnVals, NumLeaves=2, N=5, comp=2):
		## TrainPts and TrainVals are the known points and their
		# respective values which will be interpolated against.
		## Grab the min and range of each dimension
		self.tpm = np.floor(np.amin(TrnPts, axis=0))
		self.tpr = np.ceil(np.amax(TrnPts, axis=0) - self.tpm)
		## This prevents against colinear data (range = 0)
		self.tpr[self.tpr == 0] = 1
		## Normalize all points
		self.tp = (TrnPts - self.tpm) / self.tpr
		self.tv = TrnVals
		self.dims = len(TrnPts[0,:]) + 1
		self.ntpts = len(TrnPts[:,0])
		## Comp is an arbitrary value that pics a function to use
		self.comp = comp
		## The order of each comp can be found from the dims and its value
		order = int(np.floor((self.dims-1)/2) + (3*comp) + 1)
		if (comp < 0):
			## The two comps that do not follow the function above
			if (comp == -1):
				order = 9
			else:
				order = 11
		## Uses a KD Tree to set the distances, locations, and values 
		# of the closest neighbors to each point in the 
		'''
		print 'Interpolation Input Values'
		print '-Comp Equation Used:', comp
		print '-Resulting order:', order 
		print '-KDTree Leaves:', NumLeaves
		print '-Number of Training Points:', self.ntpts
		print
		'''
		## Make into training data into a Tree
		leavesz = math.ceil(self.ntpts/NumLeaves)
		## Start by creating a KDTree around the training points
		KData = spatial.cKDTree(self.tp,leafsize=leavesz)
		## For weights, first find the training points radial neighbors
		tdist, tloc = KData.query(self.tp, N, 0.0, 2)
		Tt  = tdist[:,:-1]/tdist[:,-1:]
		## Next determine weight matrix
		Rt = self.FindR(self.ntpts, Tt, tloc) 
		weights = (spsl.spsolve(csc_matrix(Rt), self.tv))[:,np.newaxis]
		## KData query takes (data, #ofneighbors) to determine closest 
		# training points to predicted data
		self.N = N
		self.KData = KData
		self.weights = weights

	def __call__(self, PredPoints):
		PrdPts = (PredPoints - self.tpm) / self.tpr
		nppts = len(PrdPts[:,0])
		## Setup prediction points and find their radial neighbors 
		pdist, ploc = self.KData.query(PrdPts.real, self.N)
		## Check if complex step is being run
		if (np.any(PrdPts[0,:].imag) > 0):
			dimdiff = np.subtract(PrdPts.reshape(nppts,1,(self.dims-1)), 
									self.tp[ploc,:])
			## KD Tree ignores imaginary part, muse redo ndist if complex 
			pdist = np.sqrt(np.sum((dimdiff*dimdiff), axis=2))
		
		## Take farthest distance of each point
		Tp  = pdist[:,:-1]/pdist[:,-1:]
		'''
		print 'Compactly Supported Radial Basis Function (CR) KDTree Results'
		print '-Number of Predicted Points:', nppts
		print '-Number of Radial Neighbers per Predicted Point:', self.N
		print '-Minimum Compact Support Domain Radius:', min(pdist[:,-1:])
		print '-Maximum Compact Support Domain Radius:', max(pdist[:,-1:])
		print
		'''
		Rp = self.FindR(nppts, Tp, ploc)
		predz = np.dot(Rp, self.weights).reshape(nppts)
		return predz

	def gradient(self, PredPoints):
		PrdPts = (PredPoints - self.tpm) / self.tpr
		#PrdPts = PredPoints
		nppts = len(PrdPts[:,0])
		## Setup prediction points and find their radial neighbors 
		pdist, ploc = self.KData.query(PrdPts, self.N)
		## Find Gradient
		grad = self.FinddR(PrdPts[:,np.newaxis,:], ploc, 
				pdist[:,:,np.newaxis])/self.tpr

		return grad


## Check and Run - Open Methods ================================================
'''
The functions in this section are unnecessary, but can be useful even if this
code is wrapped into another program.  Some may not be able to be run without 
a valid problem type.
'''


def CheckInputs(neighbors, trnpoints, NumLeaves, maximum, minimum, 
		problem, dimensions, DistanceEffect, tension, bias):
	## The first checks the inputs for default 0s or trips that may 
	# cause errors so that it can fix them.	
	## Organize inputs
	if (dimensions == 0):
		if((problem == '2D3O') or (problem == '2D5O')):
			dimensions = 2
		elif((problem == 'Plane') or (problem == 'PW') or (problem == 'Crate')):
			dimensions = 3
		elif((problem == '5D2O') or (problem == '5D4O')):
			dimensions = 5
		else:
			print 'Problem type %s does not exist.' % problem
			raise SystemExit	
	
	if (neighbors == 0):
		neighbors = max(int(trnpoints*dimensions/1500), min(10, trnpoints))
		print 'Number of neighbors has been changed to %s' % neighbors
		print
	
	if (NumLeaves == 0):
		NumLeaves = max(int(neighbors*5), 2)
		print 'Number of KD-Tree leaves has been changed to %s' % NumLeaves
		print
	
	if (NumLeaves > (trnpoints/neighbors)):
		## Can not have more leaves than training data points. Also 
		# need to ensure there will be enough neighbors in each leaf
		NumLeaves = max((trnpoints // (neighbors*1.5), 2))
		print 'Number of KD-Tree leaves has been changed to %s' % NumLeaves
		print
		if (NumLeaves > trnpoints):
			print 'The problem has too few training points for'
			print 'its number of dimensions.'
			raise SystemExit
	
	if (len(maximum) != (dimensions-1)):
		print "Due to being incorrectly inputted, the program is now"
		print "adjustng the maximum and minimum settings."
		print
		mini = minimum[0]
		maxi = maximum[0]
		minimum = np.zeros((dimensions-1), dtype="int")
		maximum = np.zeros((dimensions-1), dtype="int")
		for D in xrange(dimensions-1):
			minimum[D] = mini
			maximum[D] = maxi

	if ((trnpoints*dimensions) < 100):
		tight = True # Default algorithm had only true, change if bad 
		# results with neighs << trnpoints or a high dimension (>3)
	else:
		tight = False

	return problem, dimensions, neighbors, DistanceEffect, tension, \
			bias, tight, NumLeaves, maximum, minimum

def RunTestCode():
	print 
	print '^---- N Dimensional Interpolation ----^'
	print 
	
	## Check and shorten inputs for ease of typing:
	# pr=problem, dim=dimensions, N=neighbors, DE=DistanceEffect, 
	# t=tension, b=bias, tt=tight, NL=NumLeaves, maxi=maximum, mini=minimum
	pr, dim, N, DE, t, b, tt, NL, maxi, mini = \
			CheckInputs(neighbors, trnpoints, NumLeaves, maximum, minimum, 
					problem, dimensions, DistanceEffect, tension, bias)
	## Can override tt here
	
	## Create the Independent Data
	train = N_Data(mini, maxi, trnpoints, pr, trndist, dim)
	actul = N_Data(mini, maxi, (trnpoints*AbyT), pr, trndist, dim)
	predL = N_Data(mini, maxi, prdpoints, pr, prddist, dim)
	predW = cp.deepcopy(predL)
	predC = cp.deepcopy(predL)
	predH = cp.deepcopy(predL)
	predR = cp.deepcopy(predL)
	#predRbf = cp.deepcopy(predL)
	
	## Set Dependents for Training Data and Plot
	train.CreateDep()
	actul.CreateDep()
	
	# Setup Interpolation Methods around Training Points
	p0 = time.time()
	trainLNInt = LNInterp(train.points, train.values, NL)
	p1 = time.time()
	trainWNInt = WNInterp(train.points, train.values, NL)
	p2 = time.time()
	trainCNInt = CNInterp(train.points, train.values, NL)
	p3 = time.time()
	trainHNInt = HNInterp(train.points, train.values, NL)
	p4 = time.time()
	trainCRInt = CRInterp(train.points, train.values, NL, N, comp)
	p5 = time.time()
	#rbfi = interp.Rbf(train.points[:,0], train.points[:,1], train.values[:,0], function='multiquadric')
	#prbf = time.time()
	
	print
	print '^---- Running Interpolation Code ----^'
	print
	
	# Perform Interpolation on Predicted Points
	t0 = time.time()
	prdzL = trainLNInt(predL.points)
	tp5 = time.time()
	prdgL = trainLNInt.gradient(predL.points)
	t1 = time.time()
	prdzW = trainWNInt(predW.points, N, DE)
	t1p5 = time.time()
	prdgW = trainWNInt.gradient(predW.points, N, DE)
	t2 = time.time()
	prdzC = trainCNInt(predC.points, N)
	t2p5 = time.time()
	prdgC = trainCNInt.gradient(predC.points, N)
	t3 = time.time()
	prdzH = trainHNInt(predH.points, N, t, b, tt)
	t3p5 = time.time()
	prdgH = trainHNInt.gradient(predH.points, N, t, b, tt)
	t4 = time.time()
	prdzR = trainCRInt(predR.points)
	t4p5 = time.time()
	prdgR = trainCRInt.gradient(predR.points)
	t5 = time.time()
	#prdzrbf = rbfi(predRbf.points[:,0],predRbf.points[:,1])
	#trbf = time.time()
	
	## Assign all Dependents to Check Plots and Error
	predL.AssignDep(prdzL, prdgL)
	predW.AssignDep(prdzW, prdgW)
	predC.AssignDep(prdzC, prdgC)
	predH.AssignDep(prdzH, prdgH)
	predR.AssignDep(prdzR, prdgR)
	#predRbf.AssignDep(prdzrbf, [0])

	if CheckResults:

		## Plot All Results, Point Locations, and Error.  Since this will 
		# probably vary greatly per analysis, this function was not commented
		# a lot and was just written without much care for aesthetics or cost.
		print
		print '^---- Checking Results ----^'
		print
		
		if (dim == 3):
			actul.PlotResults('Actual Data')
			train.PlotResults('Training Data', pp)
		
			predL.PlotAll('LN Predicted Data', trainLNInt, 'LN', 
					pltfile='None', erplot=erplot, check=False, ptplot = False,
					step=step, neighs=N, DistEff=DE, tension=t, bias=b, tight=tt)
			predW.PlotAll('WN Predicted Data', trainWNInt, 'WN',
					pltfile='None', erplot=erplot, check=False, ptplot = False,
					step=step, neighs=N, DistEff=DE, tension=t, bias=b, tight=tt)
			predC.PlotAll('CN Predicted Data', trainCNInt, 'CN',
					pltfile='None', erplot=erplot, check=False, ptplot = False,
					step=step, neighs=N, DistEff=DE, tension=t, bias=b, tight=tt)
			predH.PlotAll('HN Predicted Data', trainHNInt, 'HN',
					pltfile='None', erplot=erplot, check=False, ptplot = False,
					step=step, neighs=N, DistEff=DE, tension=t, bias=b, tight=tt)
			predR.PlotAll('CR Predicted Data', trainCRInt, 'CR',
					pltfile='None', erplot=erplot, check=False, ptplot = False,
					step=step, neighs=N, DistEff=DE, tension=t, bias=b, tight=tt)
		#	predRbf.PlotResults('RBF Predicted Data', pltfile='None')
		#	predRbf.FindError('RBF Predicted Data', pltfile='None', plot=True, 
		#				  check=False, step=0.00000001)
		elif (dim == 2):
			fig = plt.figure()
			ax = plt.subplot(111)	
			Ax = actul.points[:,0].flatten()
			Ay = actul.values[:,0].flatten()
			Tx = train.points[:,0].flatten()
			Ty = train.values[:,0].flatten()
			Lx = predL.points[:,0].flatten()
			Ly = predL.values[:,0].flatten()
			Wx = predW.points[:,0].flatten()
			Wy = predW.values[:,0].flatten()
			Cx = predC.points[:,0].flatten().real
			Cy = predC.values[:,0].flatten().real
			Hx = predH.points[:,0].flatten().real
			Hy = predH.values[:,0].flatten().real
			Rx = predR.points[:,0].flatten().real
			Ry = predR.values[:,0].flatten().real
			xmi = np.min(Ax)
			xma = np.max(Ax)
			ymi = np.min(Ay)
			yma = np.max(Ay)
			ax.plot(Tx, Ty, 'bo', label='Training')
			ax.plot(Ax, Ay, 'k-', linewidth=2., label='Actual')
			ax.plot(Lx, Ly, 'm--', linewidth=2., label='Linear')
			ax.plot(Wx, Wy, '--', linewidth=2., color=[0.,.38,1.], label='Weighted')
			ax.plot(Cx, Cy, 'g--', linewidth=2., label='Cosine')
			ax.plot(Hx, Hy, 'r--', linewidth=2., label='Hermite')
			ax.plot(Rx, Ry, '--', linewidth=2., color=[1.,.48,0.], label='CS-RBF')
			box = ax.get_position()
			ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
			ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
			ax.set_title('%s Predicted Data Results' % problem,
							fontsize=16, fontweight='bold')
			ax.set_xlabel('Independent Variable', fontsize=14)
			ax.set_ylabel('Dependent Variable', fontsize=14)
			plt.tick_params(labelsize=14)
			plt.xlim(xmi, xma)
			plt.ylim(ymi, yma)
		
			if (allplot == True):
				fig = plt.figure()
				title = 'Linear'
				
				Px = predL.points[:,0].flatten()
				Py = predL.values[:,0].flatten()
				plt.plot(Ax, Ay, 'k-', Px, Py, 'r-', Tx, Ty, 'bo')
				fig.suptitle('%s Data Results' % title, 
								fontsize=14, fontweight='bold')
				plt.xlabel('Independent Variable')
				plt.ylabel('Dependent Variable')
				plt.xlim(xmi, xma)
				plt.ylim(ymi, yma)
				plt.legend(['Actual', 'Predicted', 'Training'], loc='upper left')
				
				fig = plt.figure()
				title = 'Weighted'
				
				Px = predW.points[:,0].flatten()
				Py = predW.values[:,0].flatten()
				plt.plot(Ax, Ay, 'k-', Px, Py, 'r-', Tx, Ty, 'bo')
				fig.suptitle('%s Data Results' % title, 
								fontsize=14, fontweight='bold')
				plt.xlabel('Independent Variable')
				plt.ylabel('Dependent Variable')
				plt.xlim(xmi, xma)
				plt.ylim(ymi, yma)
				plt.legend(['Actual', 'Predicted', 'Training'], loc='upper left')
			
				fig = plt.figure()
				title = 'Cosine'
				
				Px = predC.points[:,0].flatten().real
				Py = predC.values[:,0].flatten().real
				plt.plot(Ax, Ay, 'k-', Px, Py, 'r-', Tx, Ty, 'bo')
				fig.suptitle('%s Data Results' % title, 
								fontsize=14, fontweight='bold')
				plt.xlabel('Independent Variable')
				plt.ylabel('Dependent Variable')
				plt.xlim(xmi, xma)
				plt.ylim(ymi, yma)
				plt.legend(['Actual', 'Predicted', 'Training'], loc='upper left')
				
				fig = plt.figure()
				title = 'Hermite'
				
				Px = predH.points[:,0].flatten().real
				Py = predH.values[:,0].flatten().real
				plt.plot(Ax, Ay, 'k-', Px, Py, 'r-', Tx, Ty, 'bo')
				fig.suptitle('%s Data Results' % title, 
								fontsize=14, fontweight='bold')
				plt.xlabel('Independent Variable')
				plt.ylabel('Dependent Variable')
				plt.xlim(xmi, xma)
				plt.ylim(ymi, yma)
				plt.legend(['Actual', 'Predicted', 'Training'], loc='upper left')
			
				fig = plt.figure()
				title = 'CRBF'
				
				Px = predR.points[:,0].flatten().real
				Py = predR.values[:,0].flatten().real
				plt.plot(Ax, Ay, 'k-', Px, Py, 'r-', Tx, Ty, 'bo')
				fig.suptitle('%s Data Results' % title, 
								fontsize=14, fontweight='bold')
				plt.xlabel('Independent Variable')
				plt.ylabel('Dependent Variable')
				plt.xlim(xmi, xma)
				plt.ylim(ymi, yma)
				plt.legend(['Actual', 'Predicted', 'Training'], loc='upper left')
				
			predL.FindError('LN Predicted Data', 'None', erplot, False)
			predL.CheckGrad('LN Predicted Data',trainLNInt,'LN',
							pltfile='None',plot=erplot,check=False,
							step=step,N=N,DistEff=DE,tension=t,bias=b,tight=tt)
			
			predW.FindError('WN Predicted Data', 'None', erplot, False)
			predW.CheckGrad('WN Predicted Data',trainWNInt,'WN',
							pltfile='None',plot=erplot,check=False,
							step=step,N=N,DistEff=DE,tension=t,bias=b,tight=tt)
		
			predC.FindError('CN Predicted Data', 'None', erplot, False)
			predC.CheckGrad('CN Predicted Data',trainCNInt,'CN',
							pltfile='None',plot=erplot,check=False,
							step=step,N=N,DistEff=DE,tension=t,bias=b,tight=tt)
		
			predH.FindError('HN Predicted Data', 'None', erplot, False)
			predH.CheckGrad('HN Predicted Data',trainHNInt,'HN',
							pltfile='None',plot=erplot,check=False,
							step=step,N=N,DistEff=DE,tension=t,bias=b,tight=tt)
			
			predR.FindError('CR Predicted Data', 'None', erplot, False)
			predR.CheckGrad('CR Predicted Data',trainCRInt,'CR',
							pltfile='None',plot=erplot,check=False,
							step=step,N=N,DistEff=DE,tension=t,bias=b,tight=tt)
		else:
			
			predL.FindError('LN Predicted Data', 'None', erplot, False)
			predL.CheckGrad('LN Predicted Data',trainLNInt,'LN',
							pltfile='None',plot=erplot,check=False,
							step=step,N=N,DistEff=DE,tension=t,bias=b,tight=tt)
			predW.FindError('WN Predicted Data', 'None', erplot, False)
			predW.CheckGrad('WN Predicted Data',trainWNInt,'WN',
							pltfile='None',plot=erplot,check=False,
							step=step,N=N,DistEff=DE,tension=t,bias=b,tight=tt)
			predC.FindError('CN Predicted Data', 'None', erplot, False)
			predC.CheckGrad('CN Predicted Data',trainCNInt,'CN',
							pltfile='None',plot=erplot,check=False,
							step=step,N=N,DistEff=DE,tension=t,bias=b,tight=tt)
			predH.FindError('HN Predicted Data', 'None', erplot, False)
			predH.CheckGrad('HN Predicted Data',trainHNInt,'HN',
							pltfile='None',plot=erplot,check=False,
							step=step,N=N,DistEff=DE,tension=t,bias=b,tight=tt)
			predR.FindError('CR Predicted Data', 'None', erplot, False)
			predR.CheckGrad('CR Predicted Data',trainCRInt,'CR',
							pltfile='None',plot=erplot,check=False,
							step=step,N=N,DistEff=DE,tension=t,bias=b,tight=tt)
		
		print '<< Run Times >>'
		print
		print '-LN Interpolator Setup:', (p1-p0)
		print '-LN Interpolator Value Query:', (tp5-t0)
		print '-LN Interpolator Gradient Query:', (t1-tp5)
		print
		print '-WN Interpolator Setup:', (p2-p1)
		print '-WN Interpolator Value Query:', (t1p5-t1)
		print '-WN Interpolator Gradient Query:', (t2-t1p5)
		print
		print '-CN Interpolator Setup:', (p3-p2)
		print '-CN Interpolator Value Query:', (t2p5-t2)
		print '-CN Interpolator Gradient Query:', (t3-t2p5)
		print
		print '-HN Interpolator Setup:', (p4-p3)
		print '-HN Interpolator Value Query:', (t3p5-t3)
		print '-HN Interpolator Gradient Query:', (t4-t3p5)
		print
		print '-CR Interpolator Setup:', (p5-p4)
		print '-CR Interpolator Value Query:', (t4p5-t4)
		print '-CR Interpolator Gradient Query:', (t5-t4p5)
		print
		'''	
		print '-RBF Interpolator Setup:', (prbf-p5)
		print '-RBF Interpolator Value Query:', (trbf-t5)
		print
		
		with open("RBFTimes.txt", "a") as efile:
			efile.write("Number of Training Points:"+str(len(train.points[:,0]))+"\n")   
			efile.write("CSRBF Setup:"+str(p5-p4)+"\n")       
			efile.write("CSRBF Query:"+str(t4p5-t4)+"\n")       
			efile.write("RBF Setup:"+str(prbf-p5)+"\n")       
			efile.write("RBF Query:"+str(trbf-t5)+"\n")       
			efile.write("\n")
			
		with open("RunError.txt", "a") as efile:
			efile.write(predL.funct+ "\n")       
			efile.write("Run"+ "\n")
			efile.write(              "LN"+ "\n")
			efile.write(str(np.average(predL.error))  + "\n") 
			efile.write(str(    np.max(predL.error))  + "\n")
			efile.write(str(np.average(predL.gaerror))+ "\n")
			efile.write(str(    np.max(predL.gaerror))+ "\n")
			efile.write(str(np.average(predL.gcerror))+ "\n")
			efile.write(str(    np.max(predL.gcerror))+ "\n")
			efile.write(              "WN"+ "\n")
			efile.write(str(np.average(predW.error))  + "\n") 
			efile.write(str(    np.max(predW.error))  + "\n")
			efile.write(str(np.average(predW.gaerror))+ "\n") 
			efile.write(str(    np.max(predW.gaerror))+ "\n") 
			efile.write(str(np.average(predW.gcerror))+ "\n") 
			efile.write(str(    np.max(predW.gcerror))+ "\n") 
			efile.write(              "CN"+ "\n")
			efile.write(str(np.average(predC.error))  + "\n") 
			efile.write(str(    np.max(predC.error))  + "\n")
			efile.write(str(np.average(predC.gaerror))+ "\n") 
			efile.write(str(    np.max(predC.gaerror))+ "\n") 
			efile.write(str(np.average(predC.gcerror))+ "\n") 
			efile.write(str(    np.max(predC.gcerror))+ "\n") 
			efile.write(              "HN"+ "\n")
			efile.write(str(np.average(predH.error))  + "\n") 
			efile.write(str(    np.max(predH.error))  + "\n")
			efile.write(str(np.average(predH.gaerror))+ "\n") 
			efile.write(str(    np.max(predH.gaerror))+ "\n") 
			efile.write(str(np.average(predH.gcerror))+ "\n") 
			efile.write(str(    np.max(predH.gcerror))+ "\n") 
			efile.write(              "RN"+ "\n")
			efile.write(str(np.average(predR.error))  + "\n") 
			efile.write(str(    np.max(predR.error))  + "\n")
			efile.write(str(np.average(predR.gaerror))+ "\n")
			efile.write(str(    np.max(predR.gaerror))+ "\n")
			efile.write(str(np.average(predR.gcerror))+ "\n")
			efile.write(str(    np.max(predR.gcerror))+ "\n")
			efile.write("\n")
		'''
		#pp.close()


## More Information (As Promised) ==============================================
'''
		 \/\/ Nomenclature \/\/

LN-LinearNearest Neighbor,   WN-Weighted Nearest Neighbor, 
CN-Cosine Nearest Neighbor,  HN-Hermite Nearest Neighbor
prd - Predicted Points: Points which you want the values at
trn - Training Points:  Points with known values 
Ind - Independent: Dimensions which should be known
Dep - Dependent:   Dimensions which are a function of Ind

	    \/\/ Method Breakdown \/\/

# Local Methods ================================================================

df suppress_stdout():
df Solve(points, funct):		
	rtrn z

		\/\/ Class Breakdown \/\/

# Local Data Class =============================================================

clss N_Data(object):
	df __init__(self, mini, maxi, numpts, funct='PW', dtype='rand', dims=3):
	df AssignDep(self, values, gradient=[0]):
	df CreateDep(self):
	df PlotResults(self, title, pltfile='None'):
	df PlotPoints(self, title, pltfile='None'):
	df FindError(self, title, pltfile='None', plot=True, 
				  check=False, step=0.00000001):
	df CheckGrad(self, title, Interp, itype, pltfile='None', plot=True,
				 check=False, step=0.00000001, N=5, DistEff=3, 
				 tension=0, bias=0, tight=False):
	df PlotAll(self, title, Interp, itype, pltfile='None', erplot=True, 
			   check=False, ptplot = False, step=0.00000001, neighs=5, 
			   DistEff=3, tension=0, bias=0, tight=False):

# Interpolation Classes ========================================================

clss NNInterpBase(object):
	df __init__(self, TrnPts, TrnVals, NumLeaves=2):

clss LNInterp(NNInterpBase):
	df main(self, PrdPts, nppts, dims, nloc):
		rtrn normal, prdtemp, pc
	df main2D(self, PrdPts, nppts, nloc):
		rtrn m, b
	df __call__(self, PredPoints):
		rtrn predz
	df gradient(self, PredPoints):
		rtrn grad

clss WNInterp(NNInterpBase):
	df __call__(self, PredPoints, N=5, DistEff=3):
		rtrn predz
	df gradient(self, PredPoints, N=5, DistEff=3):
		rtrn grad

clss CNInterp(NNInterpBase):
	df main(self, PrdPts, nppts, D, N, neighvals):	
		rtrn zu, zl, diff, ddiff
	df __call__(self, PredPoints, N=5):
		rtrn predz
	df gradient(self, PredPoints, N=5):
		rtrn grad

clss HNInterp(NNInterpBase):
	df HermFunctRes(self, y, mu, tension, bias):
		rtrn (a0*y[:,1] + a1*m0 + a2*m1 + a3*y[:,2])
	df HermFunctGrad(self, y, mu, dmu, tension, bias):
		rtrn ((dmu*(b0*y[:,1]+b1*m0+b2*m1+b3*y[:,2]))/(self.dims-1)).real
	df main(self, PrdPts, nppts, D, N, u, l, neighvals):
		rtrn y, diff, ddiff
	df __call__(self, PredPoints, N=5, tension=0, bias=0, tight=False):
		rtrn predz
	df gradient(self, PredPoints, N=5, tension=0, bias=0, tight=False):
		rtrn grad

clss CRInterp(object):
	df FindR(self, npp, T, loc): 
		rtrn R
	df FinddR(self, PrdPts, ploc, pdist):
		rtrn np.sum((dRp * dtx * self.weights[ploc[:,:-1]]), axis=1)
	df __init__(self, TrnPts, TrnVals, N=5, NumLeaves=2, comp=2):
	df __call__(self, PredPoints):
		rtrn predz
	df gradient(self, PredPoints):
		rtrn grad

## Check Inputs and Results - Open Methods =====================================

df CheckInputs(neighbors, trnpoints, NumLeaves, maximum, minimum, 
		problem, DistanceEffect, tension, bias, tight):
	rtrn problem, dimensions, neighbors, DistanceEffect, tension, \
			bias, tight, NumLeaves, maximum, minimum

df RunTestCode():

Note - some vowels removed to ensure vim friendliness.
'''


## Running Code ================================================================

## Problem Inputs 
dimensions = 0 # If 0, the dimensions is chosen from the problem type
minimum = np.array([-500, -500]) # Minimum value for independent range
maximum = np.array([500, 500]) # Maximum value for independent range
trndist = 'cart' # Can be rand, LH, or cart (only for 2D and 3D)
prddist = 'rand' # Can be rand, LH, or cart (only for 2D and 3D)
problem = '2D3O' # Problem type, options seen in organize inputs loop below

trnpoints = 100 # Number of training pts, min of 5 because of Hermite lims
prdpoints = 3000 # Number of prediction points

neighbors = 0 # KD-Tree neighbors found, default 0, min 2
DistanceEffect = 2 # Effect of distance of neighbors in WN, default 2
tension = 0 # Hermite adjustable, loose is -ive, tight fit is +ive, default 0
bias = 0 # Attention to each closest neighbor in hermite, default 0
comp = 4 # Type of CRBF used in the CR interpolation, default 4
NumLeaves = 0 # Leaves of KD Tree, default 0

## Extra Inputs
CheckResults = True #All results will not be checked if this is false
step = 0.00000001 # Complex step step size
allplot = False #For 2D, plot separate graphs for each interp
erplot = False #Choose to plot error
AbyT = 100 # Multilpier for actual to training data points
ptplot = False #For plotting the predicted point locations
# Variable for saved file of multiple plots
#pp = PdfPages('ND_Interpolation_Plots.pdf') 
pp = 'None'

RunTestCode()

plt.show()                              

'''
Side note: PrdPts are predicted points technically, although it's not really 
that simple. They are better named pride points.  Everyone needs pride points,
but their value is unknown usually and can only be determined when related to
those around you. Usually taken from Corey Watt.
'''
