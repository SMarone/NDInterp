import numpy as np
from  scipy import spatial
from matplotlib import cm
import matplotlib.pyplot as plt
import time
import math
from mpl_toolkits.mplot3d import Axes3D

##   Written by Stephen Marone
##     Intern at the NASA GRC  
## Project Start - June 8th, 2015

## This is an n-dimensional interpolation code.  More information can be found
# at the end of the code along with the lines which run it.

# Setting up the original problem ==============================================

class N_Data(object):
	## Main data type with N dimension.  Only needed if you want to use the 
	# example functions, check for error, or plot.
	def __init__(self, mini, maxi, numpts, funct='PW', dtype='rand', dims=3):
		## The funct is Crate or PW for Piecewise.
		## The dtype is cart for cartesian or rand for random
		self.funct = funct
		self.dtype = dtype
		self.dims = dims
		
		print "Creating Independent Data"
		print "-Function Type:", funct
		print "-Data Distribution Type:", dtype
		print "-Quantity of Points:", numpts
		print "-Dimensions:", dims
		print "-Range of each dimension: [%s,%s]" % (mini, maxi)
		print
		## Setup independent data for each function via distribution
		if (dtype == 'rand'):
			self.points = np.random.rand(numpts,dims-1) * \
										(maxi-mini) + mini
			if (funct == 'PW'):
				## This ensures points exist on discontinuities
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

		elif (dtype == 'LH'):
			points = np.zeros((numpts,dims-1), dtype="float")
			for d in range(dims-1):
				## Latin Hypercube it up
				points[:,d] = np.linspace(mini,maxi,numpts)
				np.random.shuffle(points[:,d])
			self.points = points
			
		elif (dtype == 'cart'):
			if (dims == 3):
				stp = (maxi-mini)/(np.sqrt(numpts)-1)
				x = np.arange(mini,(maxi+stp),stp)
				y = np.arange(mini,(maxi+stp),stp)
				x, y = np.meshgrid(x, y)
				self.points = np.transpose(np.array([x.flatten(), y.flatten()]))
			else:
				print 'Can not currently do cartesian in any dimension \
						except for 3'
				raise SystemExit
			'''
			## Could not figure out why np would not do a meshgrid in N dims
			## Defaults for crate: mini=-512,maxi=512,numpts=81
			stp = (maxi-mini)/(np.sqrt(numpts)-1)
			xm = np.meshgrid(*[np.arange(mini,(maxi+stp),stp) for \
					d in range(dims)], sparse=False, indexing='ij')
			'''
		else:
			print 'Distribution type not found.'
			raise SystemExit

	def AssignDep(self, values):
		## Takes numpy array of values and stores inside the class variable
		self.values = values[:,np.newaxis]

	def CreateDep(self):
		## This not only creates z but the data subtype as well
		
		if (self.funct == 'Crate'):
  			x = self.points[:,0]
			y = self.points[:,1]
			print 'Setting up Crate Function Dependents'
			print
			## Egg Crate Function
			z = -x * np.sin(np.sqrt(abs(x - y - 47))) - \
			    (y+47) * np.sin(np.sqrt(abs(x/2 + y + 47)))
		elif (self.funct == 'PW'):
			x = self.points[:,0]
			y = self.points[:,1]
			print 'Setting up Piecewise Function Dependents'
			print
			z = np.empty(len(x))
			## Piecewise Function with 3 Planes
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
		else:
			print 'Function type not found.'
			raise SystemExit
		
		self.values = z[:,np.newaxis]

	def PlotResults(self, title):
		## Plots the value results of a 3 dimensional problem

		## Ensure points and values have been created
		try:
			self.values
		except NameError:
			print 'Values have not been set.'
			return

		x = self.points[:,0]
		y = self.points[:,1]
		fig = plt.figure()
		ax = fig.gca(projection='3d')
		surf = ax.tricontour(x.flatten(), y.flatten(), 
							 self.values.flatten(), 
							 100, rstride=1, cstride=1, 
							 cmap=cm.coolwarm, linewidth=0, 
							 antialiased=False)
		fig.colorbar(surf, shrink=0.5, aspect=5)
		fig.suptitle('%s Results' % title, 
					 fontsize=14, fontweight='bold')
		ax.set_xlabel('X - Independent Variable')
		ax.set_ylabel('Y - Independent Variable')
		ax.set_zlabel('Z - Dependent Variable')

	def PlotPoints(self, title):
		## Dim-1 plot of point locations for visualization

		## Ensure points have been created
		try:
			self.points
		except NameError:
			print 'Points have not been set.'
			return

		fig = plt.figure()
		if (self.dims == 3):
			plt.scatter(self.points[:,0], self.points[:,1])
			fig.suptitle('%s Point Locations' % title, 
							fontsize=14, fontweight='bold')
			plt.xlabel('X Value')
			plt.ylabel('Y Value')
		elif (self.dims == 4):
			ax = fig.gca(projection='3d')
			ax.scatter(self.points[:,0], self.points[:,1], 
						self.points[:,2], c=c, marker=m)
			fig.suptitle('%s Point Locations' % title, 
							fontsize=14, fontweight='bold')
			ax.set_xlabel('X Value')
			ax.set_ylabel('Y Value')
			ax.set_zlabel('Z Value')

	def FindError(self, title, plot=True, check=False):
		## Find and plot error of the found values
		
		## Ensure points and values have been created
		try:
			self.values
		except NameError:
			print 'Values have not been set.'
			return

		vals = self.values.flatten()
		xc = self.points[:,0]
		yc = self.points[:,1]
		loc = 1+np.arange(len(vals), dtype='float')
		
		if (self.funct == 'Crate'):
			zc = -xc * np.sin(np.sqrt(abs(xc - yc - 47))) - \
					(yc+47) * np.sin(np.sqrt(abs(xc/2 + yc + 47)))
		elif (self.funct == 'PW'):
			# Plot data against values from crate problem
			zc = np.empty(len(xc))
			# Make Piecewise Function with 3 Planes
			count = 0
			for i in range(len(xc)):
				if (xc[i] > 0):
					if (yc[i] > 0):
						zc[count] = xc[i] + 8*yc[i]
						count +=1
					else:
						zc[count] = 5*xc[i] + yc[i]
						count +=1
				else:
					zc[count] = 5*xc[i] + 6*yc[i]
					count +=1		
		else:
			print 'No function available for comparison.'
			print
			return
		error = abs((zc-vals)/(zc+0.000000000001)) * 100
		avg = np.average(error)
		
		if check:
			for i in range(len(error)):
				if (error[i] > 100):
					print '**High Percent Error found of', error[i], '**'
					print '-Location:', xc[i], yc[i]
					print '-Calculated Value:', zc[i]
					print '-Found Value:', vals[i] 
					print
		
		print '%s Error Information' % title
		print '-Average Percent Error:', avg
		print '-Max Percent Error:', np.max(error)
		print

		if plot:
			fig = plt.figure()
			plt.plot(loc, error, 'k')
			fig.suptitle('%s Error Per Point' % title, 
						  fontsize=14, fontweight='bold')
			plt.xlabel('Point Identifier')
			plt.ylabel('Percent Error')
			plt.ylim((0,((avg+0.00000000001)*2)))

		self.error = error
		self.actualvalues = zc

	def PlotAll(self,title, check=False):
		self.PlotResults(title)
		self.PlotPoints(title)
		self.FindError(title, True, check)

# Interpolation Classes ========================================================


class InterpBase(object):
	## The shared aspects of all interpolators
	def __init__(self, TrnPts, TrnVals, NumLeaves=8):
		## TrainPts and TrainVals are the known points and their
		# respective values which will be interpolated against.
		# The IType is the type of interpolation that the user wants.
		self.tp = TrnPts
		self.tv = TrnVals
		self.dims = len(TrnPts[0,:]) + 1
		self.ntpts = len(TrnPts[:,0])
		self.nl = NumLeaves

	def FindNeighs(self, PrdPts, N=5): 
		## Uses a KD Tree to set the distances, locations, and values 
		# of the closest neighbors to each point in the 
		nppts = len(PrdPts[:,0])
		prdz = np.zeros(nppts, dtype="float")
		numtrn = self.ntpts
		numleaves = self.nl

		prdz   =   np.zeros((nppts), dtype="float")
		gradient = np.zeros((nppts, self.dims-1), dtype="float")

		if (numleaves > (numtrn/N)):
			## Can not have more leaves than training data points. Also 
			# need to ensure there will be enough neighbors in each leaf
			numleaves = numtrn // (N*1.5)
			print 'Number of leaves has been changed to %s' % numleaves
			if (numleaves > numtrn):
				print 'The problem has too few training points for'
				print 'its number of dimensions.'
				raise SystemExit
		
		print 'Interpolation Input Values'
		print '-KDTree Leaves:', numleaves
		print '-Number of Neighbors per Predicted Point:', N
		print '-Number of Training Points:', numtrn
		print '-Number of Predicted Points:', nppts
		print
		## Make into training data into a Tree
		leavesz = math.ceil(numtrn/numleaves)
		KData = spatial.KDTree(self.tp,leafsize=leavesz)

		## KData query takes (data, #ofneighbors) to determine closest 
		# training points to predicted data
		ndist, nloc = KData.query(PrdPts ,N)
		return prdz, gradient, nppts, ndist, nloc

class LNInterp(InterpBase):
	def __call__(self, PrdPts):
		## This method uses linear interpolation by defining a plane with
		# a set number of nearest neighbors to the predicted
		dims = self.dims

		## Find them neigbors
		prdz, gradient, nppts, ndist, nloc = self.FindNeighs(PrdPts, N=dims)

		print 'Linear Plane Nearest Neighbors (LN) KDTree Results'
		print '-Nearest Neighbor Distance:', np.min(ndist)
		print '-Farthest Neighbor Distance:', np.max(ndist)
		print
		## Extra Inputs for Finding the normal are found below
		
		## Number of row vectors needed always dimensions - 1
		r = np.arange(dims-1, dtype="int")
		
		## Diagonal counters are found here to speed up finding each normal
		da = np.zeros((dims, dims-1), dtype="int")
		db = np.zeros((dims, dims-1), dtype="int")
		for i in xrange(dims):
		    da[i] = (r+1+i) % dims
		    db[i] = (i-r-1) % dims
		'''
		t2 = time.time()
		for t in range(nppts):  ## Go through each predict point 
			nvect = np.empty((0,dims), float)
			## Planar vectors need both dep and ind dimensions
			trnd = np.concatenate((self.tp[nloc[t,:],:] ,
									self.tv[nloc[t,:]]),
			 							axis=1)
			for n in range(dims-1): ## Go through each neighbor
				## Creates array[neighbor, dimension] from NN results
				nvect = np.append(nvect, [trnd[(n+1),:] - trnd[n,:]], 0)
	
			## It is known that found array will be size [dims-1, dims]
			# and can be used below to in an exterior product.
			normal = np.prod(nvect[r,da], axis=1) - \
					 np.prod(nvect[r,db], axis=1)

			## The pc is the constant of the n dimensional plane. 
			## It uses the point which is from the closest neighbor
			pc = np.dot(trnd[0,:],normal)

			## Set the prddata to a temporary array so it has the same
			# dimensions as the training data with z as 0 for now
			prdtemp = np.concatenate((PrdPts[t,:], np.array([0])), 1)
			
			if (normal[-1] == 0):
				print 'Error, dependent variable has infinite answers.'
				raise SystemExit


			gradient[t,:] = -normal[:-1]/normal[-1]
			prdz[t] = (np.dot(prdtemp,normal)-pc)/-normal[-1]
		t3 = time.time()
		'''

		## Planar vectors need both dep and ind dimensions
		trnd = np.concatenate((self.tp[nloc,:] ,
		                       self.tv[nloc].reshape(nppts,dims,1)),
		                       axis=2)
		## Set the prddata to a temporary array so it has the same
		# dimensions as the training data with z as 0 for now
		prdtemp = np.concatenate((PrdPts, np.zeros((nppts,1),float)), 1)
		
		nvect = np.empty((nppts,(dims-1),dims), float)
		for n in range(dims-1): ## Go through each neighbor
		    ## Creates array[neighbor, dimension] from NN results
		    nvect[:,n,:] = trnd[:,(n+1),:] - trnd[:,n,:]
		
		## It is known that found array will be size [dims-1, dims]
		# and can be used below to in an exterior product.
		normal = np.prod(nvect[:,r,da], axis=2) - \
				 np.prod(nvect[:,r,db], axis=2)
		
		## The pc is the constant of the n dimensional plane.
		## It uses the point which is from the closest neighbor
		pc = np.einsum('ij, ij->i',trnd[:,0,:],normal)
		
		if (np.any(normal[-1]) == 0):
		    print 'Error, dependent variable has infinite answers.'
		    raise SystemExit
		
		gradient = -normal[:,:-1]/normal[:,-1:]
		prdz = (np.einsum('ij, ij->i',prdtemp,normal)-pc)/-normal[:,-1]
		#t4 = time.time()
		#print 'Time old:', (t3-t2)
		#print 'Time new:', (t4-t3)
		#prddict['%s' % prddata[t,:]] = prdz[t] #from when I wanted a dict
		return prdz
	
class WNInterp(InterpBase):
	## Weighted Neighbor Interpolation
	def __call__(self, PrdPts, N=5, DistEff=5):

		prdz, gradient, nppts, ndist, nloc = self.FindNeighs(PrdPts, N=N)

		print 'Weighted Nearest Neighbors (WN) KDTree Results'
		print '-Nearest Neighbor Distance:', np.min(ndist)
		print '-Farthest Neighbor Distance:', np.max(ndist)
		print
		'''
		for t in range(nppts):  ## Go through each predict point 
			sumdist = np.sum(1/(ndist[t,:]**DistEff))
			wt = np.sum(self.tv[nloc[t,:]].T/(ndist[t,:]**DistEff))
			prdz[t] = wt/sumdist
		'''
		part = ndist**DistEff
		sumdist = np.sum(1/part, axis=1)
		vals = self.tv[nloc]
		wt = np.sum(vals[:,:,0]/part, axis=1)  
		prdz = wt/sumdist
		pdims = self.dims - 1

		dimdiff = np.subtract(PrdPts.reshape(nppts,1,pdims), self.tp[nloc,:])
		gradient = np.sum(((DistEff*vals * dimdiff) / \
						   (ndist.reshape(nppts,N,1)**2)) * 
						   (1-(1/(part.reshape(nppts,N,1) * \
								  sumdist.reshape(nppts,1,1))))
									, axis=1)
		return prdz

class CNInterp(InterpBase):
	## Cosine Neighbor Interpolation
	def __call__(self, PrdPts, N=5, tension=0, bias=0):

		prdz, gradient, nppts, ndist, nloc = self.FindNeighs(PrdPts, N=N)
		
		print 'Cosine Nearest Neighbors (CN) Interpolation KDTree Results'
		print '-Nearest Neighbor Distance:', np.min(ndist)
		print '-Farthest Neighbor Distance:', np.max(ndist)
		print
		'''
		t1 = time.time()	
		for t in range(nppts):  ## Go through each predict point 
			if (np.all(abs(PrdPts[t,:] - self.tp[nloc[t,0],:])) > 0.000001):
				## This sorts the rows by the values in the column specified. 
				## The neighvals dim is [neigh#, dimension] 

				neighvals = np.concatenate((self.tp[nloc[t,:],:] ,
											self.tv[nloc[t,:],:]),
				 							axis=1)
				orgneighs = neighvals[neighvals[:,0].argsort(),:]
				
				i = 0
				while True:
					i += 1
					if (i == N-1):
						i -= 1
						break
					elif (PrdPts[t,0] < orgneighs[i,0]):
						break
	
				diff = ((PrdPts[t,0] - orgneighs[i-1,0]) / 
				        (orgneighs[i,0]-orgneighs[i-1,0]))
				mu2 = (1-np.cos(diff*np.pi))/2
				prdz[t] = orgneighs[i-1,-1]*(1-mu2) + orgneighs[i,-1]*mu2
				#This is not as accurate but removes the while loop
				#diff = ((PrdPts[t,0] - orgneighs[0,0]) / 
				#        (orgneighs[1,0]-orgneighs[0,0]))
				#prdz[t] = y[0]*(1-diff) + y[1]*diff
				
				for D in range(self.dims-2):
					## Dim started above to initiate prdz, then continues to finish
					dim = D+1
					orgneighs = neighvals[neighvals[:,dim].argsort(),:]
					
					i = 0
					while True:
						i += 1
						if (i == N-1):
							i -= 1
							break
						elif (PrdPts[t,dim] < orgneighs[i,dim]):
							break
	
					## Place of prddata along mid neighbors
					diff = ((PrdPts[t,dim] - orgneighs[i-1,dim]) / 
					        (orgneighs[i,dim]-orgneighs[i-1,dim]))
					mu2 = (1-np.cos(diff*np.pi))/2
					prdz[t] += orgneighs[i-1,-1]*(1-mu2) + orgneighs[i,-1]*mu2
					#This is not as accurate but removes the while loop
					#diff = ((PrdPts[t,dim] - orgneighs[0,dim]) / 
				    #	    (orgneighs[1,dim]-orgneighs[0,dim]))
					#prdz[t] += y[0]*(1-diff) + y[1]*diff
					
					prdz[t] = prdz[t]/2 ## Average with previous
			else:
				# Closest neighbor overlaps predicted point
				#print 'That just happened'
				#print (PrdPts[t,:-1] - self.tp[nloc[t,0],:-1])
				prdz[t] = self.tv[nloc[t,0]]
		t2 = time.time()
		'''
		## This sorts the rows by the values in the column specified. 
		## The neighvals dim is [neigh#, dimension] 

		neighvals = np.concatenate((self.tp[nloc,:] ,
									self.tv[nloc,:]),
		 							axis=2)
		tprdz = np.zeros((1,nppts), dtype='float')
		orgneighs = np.empty((nppts,N,self.dims), float)
		anppts = np.arange(nppts)
		for D in range(self.dims-1):
			#orgneighs = neighvals[np.arange(nppts), 
			#					  neighvals[:,:,0].argsort(),
			#					  np.array([[np.arange(self.dims),]*N,]*nppts)]
			for t in range(nppts): 
				## I want this gone but can't figure out sorting
				orgneighs[t,:,:] = neighvals[t,neighvals[t,:,D].argsort(),:]
			## Make some temporary variables
			tprd = PrdPts.reshape(nppts,1,(self.dims-1))
			podiff = np.subtract(orgneighs[:,:,D], tprd[:,:,D])
			## Gives the indices of the closest neighbors to each prdpt per dim.
			cnd = np.argmin((np.abs(podiff)), axis=1)
			## If the closest neighbor per dim is a smaller value, add 1.
			cnd += np.ceil(-podiff[anppts,cnd] / \
					np.abs(podiff[anppts,cnd]*2))
			## Stay in the range!
			cnd[cnd == 0] = 1
			cnd[cnd == N] = N-1
			
			diff = -podiff[anppts,(cnd-1)] / \
			       (orgneighs[anppts,cnd,D] - orgneighs[anppts,(cnd-1),D])
			mu2 = (1-np.cos(diff*np.pi))/2
			tprdz += orgneighs[anppts,(cnd-1),-1] * (1-mu2) + \
					 orgneighs[anppts,cnd,-1] * mu2

		prdz = tprdz/(D+1)
		#t3 = time.time()
		#print 'Old Time:', (t2-t1)
		#print 'New Time:', (t3-t2)
		return prdz

class HNInterp(InterpBase):
	## Hermite Neighbor Interpolation
	def HermFunct(self, y, mu, tension, bias):
		## Should  have 4 neighbor values (y) and percent distance along line 
		# from y1 and y2 to get the wanted value (mu).  Bias and tension 
		# are optional values.
	
		mu2 = mu * mu
		mu3 = mu2 * mu
	
		m0  = (y[1]-y[0])*(1+bias)*(1-tension)/2
		m0 += (y[2]-y[1])*(1-bias)*(1-tension)/2
		m1  = (y[2]-y[1])*(1+bias)*(1-tension)/2
		m1 += (y[3]-y[2])*(1-bias)*(1-tension)/2
		a0 =  2*mu3 - 3*mu2 + 1
		a1 =    mu3 - 2*mu2 + mu
		a2 =    mu3 -   mu2
		a3 = -2*mu3 + 3*mu2

		return (a0*y[1] + a1*m0 + a2*m1 + a3*y[2])

	def HermFunctArr(self, y, mu, tension, bias):
		## Should  have 4 neighbor values (y) and percent distance along line 
		# from y1 and y2 to get the wanted value (mu).  Bias and tension 
		# are optional values.
	
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
	
	def __call__(self, PrdPts, N=8, tension=0, bias=0):
		## Traindata has n dimensions, prddata has n-1 dimensions because last
		# dimension is the dependent one we are solving for.  
		## N neighbors will be found but only 4 will be used per dimension

		prdz, gradient, nppts, ndist, nloc = self.FindNeighs(PrdPts, N=N)

		print 'Hermite Nearest Neighbors (HN) Interpolation KDTree Results'
		print '-Nearest Neighbor Distance:', np.min(ndist)
		print '-Farthest Neighbor Distance:', np.max(ndist)
		print
		'''
		t1 = time.time()
		for t in range(nppts):  ## Go through each predict point 
			if (np.all(abs(PrdPts[t,:] - self.tp[nloc[t,0],:])) > .001):
				## This sorts the rows by the values in the column specified. 
				## The neighvals dim is [neigh#, dimension] 

				neighvals = np.concatenate((self.tp[nloc[t,:],:] ,
											self.tv[nloc[t,:],:]),
				 							axis=1)
				orgneighs = neighvals[neighvals[:,0].argsort(),:]
				
				i = 1
				while True:
					i += 1
					if (i == N-2):
						i -= 1
						break
					elif (PrdPts[t,0] < orgneighs[i,0]):
						break
	
				diff = ((PrdPts[t,0] - orgneighs[i-1,0]) / 
				        (orgneighs[i+2,0]-orgneighs[i-1,0]))
				prdz[t] = self.HermFunct(orgneighs[(i-2):(i+2),-1],
										 diff, tension, bias)
				# This is not as accurate but removes the while loop
				#diff = ((PrdPts[t,0] - orgneighs[0,0]) / 
				#        (orgneighs[3,0]-orgneighs[0,0]))
				#prdz[t] = self.HermFunct(orgneighs[:4,-1],
				#						 diff, tension, bias)
				
				for D in range(self.dims-2):
					## Dim started above to initiate prdz, then continues to finish
					dim = D+1
					orgneighs = neighvals[neighvals[:,dim].argsort(),:]
					
					i = 1
					while True:
						i += 1
						if (i == N-2):
							i -= 1
							break
						elif (PrdPts[t,dim] < orgneighs[i,dim]):
							break
	
					## Place of prddata along mid neighbors
					#diff = ((PrdPts[t,dim] - orgneighs[i-1,dim]) / 
					#        (orgneighs[i+2,dim]-orgneighs[i-1,dim]))
					#prdz[t] += self.HermFunct(orgneighs[(i-2):(i+2),-1],
					#						  diff, tension, bias)
					# This is not as accurate but removes the while loop
					#diff = ((PrdPts[t,dim] - orgneighs[0,dim]) / 
					# 	    (orgneighs[3,dim]-orgneighs[0,dim]))
					#prdz[t] += self.HermFunct(orgneighs[:4,-1],
					#						  diff, tension, bias)
					#
					prdz[t] = prdz[t]/2 ## Average with previous
			else:
				# Closest neighbor overlaps predicted point
				#print 'That just happened'
				#print (PrdPts[t,:-1] - self.tp[nloc[t,0],:-1])
				prdz[t] = self.tv[nloc[t,0]]
		t2 = time.time()	
		'''
		## This sorts the rows by the values in the column specified. 
		## The neighvals dim is [neigh#, dimension] 

		neighvals = np.concatenate((self.tp[nloc,:] ,
									self.tv[nloc,:]),
		 							axis=2)
		tprdz = np.zeros((1,nppts), dtype='float')
		orgneighs = np.empty((nppts,N,self.dims), float)
		anppts = np.arange(nppts)
		y = np.empty((nppts,4), float)
		for D in range(self.dims-1):
			#orgneighs = neighvals[np.arange(nppts), 
			#					  neighvals[:,:,0].argsort(),
			#					  np.array([[np.arange(self.dims),]*N,]*nppts)]
			for t in range(nppts): 
				## I want this gone but can't figure out sorting
				orgneighs[t,:,:] = neighvals[t,neighvals[t,:,D].argsort(),:]
			## Make some temporary variables
			tprd = PrdPts.reshape(nppts,1,(self.dims-1))
			podiff = np.subtract(orgneighs[:,:,D], tprd[:,:,D])
			## Gives the indices of the closest neighbors to each prdpt per dim.
			cnd = np.argmin((np.abs(podiff)), axis=1)
			## If the closest neighbor per dim is a smaller value, add 1.
			cnd += np.ceil(-podiff[anppts,cnd] / \
					np.abs(podiff[anppts,cnd]*2))
			## Stay in the range!
			cnd[cnd <= 1] = 2
			cnd[cnd >= (N-1)] = N-2
		
			## Find location of value in min and max of neighvals to be used
			diff = -podiff[anppts,(cnd-2)] / \
			       (orgneighs[anppts,(cnd+1),D] - orgneighs[anppts,(cnd-2),D])
			for n in range(4):
				## Only need 4 inputs.  Would like to remove this sometime
				y[:,n] = orgneighs[anppts,(cnd-2+n),-1]
			tprdz += self.HermFunctArr(y, diff, tension, bias)
			#tprdz = self.HermFunctArr(orgneighs[anppts,np.array([np.arange((cnd-2),N),]*nppts),-1], diff,tension,bias)

		prdz = tprdz/(D+1)
		#t3 = time.time()
		#print 'Old Time:', (t2-t1)
		#print 'New Time:', (t3-t2)
		return prdz

## More Information (As Promised) ==============================================

print 
print '^---- N Dimensional Interpolation ----^'
print 


'''

		 \/\/ Nomenclature \/\/

LN-LinearNearest Neighbor,   WN-Weighted Nearest Neighbor, 
CN-Cosine Nearest Neighbor,  HN-Hermite Nearest Neighbor
prd - Predicted Points: Points which you want the values at
trn - Training Points:  Points with known values 
Ind - Independent: Dimensions which should be known
Dep - Dependent:   Dimensions which are a function of Ind


		\/\/ Class Breakdown \/\/

clss N_Data(object)
df __init__(self, mini, maxi, numpts, funct='PW', dtype='rand', dims=3):
df AssignDep(self, val)
df CreateDep(self)
df PlotResults(self, title)
df PlotPoints(self, title)
df FindError(self, title, plot=True, check=False)
df PlotAll(self, title, check=False)

clss InterpBase(object)
df __init__(self, TrnPts, TrnVals, NumLeaves=8)
df FindNeighs(self, PrdPts, N=5)
	rturn nppts, ndist, nloc

clss LNInterp(InterpBase)
df __call__(self, PrdPts)
	rturn prdz
	
clss WNInterp(InterpBase)
df __call__(self, PrdPts, N=8, DistEff=5)
	rturn prdz

clss CNInterp(InterpBase)
df __call__(self, PrdPts, N=5, tension=0, bias=0)
	rturn prdz

clss HNInterp(InterpBase)
df HermFunct(y, mu, tension, bias)
	rturn (a0*y[1] + a1*m0 + a2*m1 + a3*y[2])
df __call__(self, PrdPts, N=8, tension=0, bias=0)
	return prdz


Note - some vowels removed to ensure vim friendliness.
'''

## Running Code ================================================================


## Create the Independent Data
train = N_Data(-500, 500, 500000)
pred = N_Data(-500, 500, 1000, 'PW', 'LH')

## Set Dependents for Training Data and Plot
train.CreateDep()
train.PlotResults('Training Data')

## Setup Interpolation Methods around Training Points
trainLNInt = LNInterp(train.points, train.values, NumLeaves=100)
trainWNInt = WNInterp(train.points, train.values, NumLeaves=100)
trainCNInt = CNInterp(train.points, train.values, NumLeaves=100)
trainHNInt = HNInterp(train.points, train.values, NumLeaves=100)

print '^---- Running Interpolation Code ----^'
print

## Perform Interpolation on Predicted Points
t0 = time.time()
pred.AssignDep(trainLNInt(pred.points))
t1 = time.time()

## Plot All Results, Point Locations, and Error
pred.PlotAll('LN Predicted Data')

## Perform Interpolation on Predicted Points
t2 = time.time()
pred.AssignDep(trainWNInt(pred.points))
t3 = time.time()

## Plot All Results, Point Locations, and Error
pred.PlotAll('WN Predicted Data')

## Perform Interpolation on Predicted Points
t4 = time.time()
pred.AssignDep(trainCNInt(pred.points))
t5 = time.time()

## Plot All Results, Point Locations, and Error
pred.PlotAll('CN Predicted Data')

## Perform Interpolation on Predicted Points
t6 = time.time()
pred.AssignDep(trainHNInt(pred.points))
t7 = time.time()

## Plot All Results, Point Locations, and Error
pred.PlotAll('HN Predicted Data')

print 'Run Times'
print '-LN Interpolator:', (t1-t0)
print '-WN Interpolator:', (t3-t2)

print '-CN Interpolator:', (t5-t4)
print '-HN Interpolator:', (t7-t6)

plt.show()

## Side note: PrdPts are predicted points technically although it's not really 
# that simple. They are better named pride points.  Everyone needs pride points,
# but their value is unknown usually and can only be determined when related to
# those around you.
