from matplotlib import pylab 
import numpy as np

linecolors = [ 'red', 'green', 'blue', 'orange', 'magenta', 'cyan',
   'saddlebrown', 'skyblue', 'olivedrab', 'yellowgreen', 'black' ]

def listtounstruct(alpha, data):
	ti = 0
	ilength = len(data[alpha])
	jlengths = np.zeros((ilength), dtype='int')
	for i in range(ilength):
		jlengths[i] = len(data[alpha][i])
	numpts = np.sum(jlengths)
	uspts = np.empty((numpts, 3), dtype='float')
	usvls = np.empty((numpts, 3, 1), dtype='float')
	
	for i in range(ilength):
		for j in range(jlengths[i]):
			# Independent Variables in mapdata: 0-alpha, 1-speed, 2-Rline
			uspts[ti,:] = [mapdata[alpha][i][j][0],   \
					mapdata[alpha][i][j][1], mapdata[alpha][i][j][2]]
			# Dependent Variables in mapdata: 3-Wcorr, 4-PR, 5-eff 
			usvls[ti,:,0] = [mapdata[alpha][i][j][3], \
					mapdata[alpha][i][j][4], mapdata[alpha][i][j][5]]
			ti += 1

	return uspts, usvls, numpts, ilength, jlengths

def unstructtolist(uspts, usvls, numpts, ilength, jlengths):
	ldata = []
	i_list = []
	j_list = []
	ti = 0
	for i in range(ilength):
		for j in range(jlengths[i]):
			i_list.append([uspts[ti,0], uspts[ti,1], uspts[ti,2], \
							usvls[ti,0], usvls[ti,1], usvls[ti,2]])
		ldata.append(ilist)

	return ldata

def interp_compressor():
	# ==========================================================================
	#  STORE TRAINING DATA IN INTERPOLATION METHOD AND CREATE PREDICTED DATA 
	# ==========================================================================

	# Training points are of shape [alpha, point number, independent dimension].
	# Training values are if shape [alpha, point number, dependent dimension].
	# Prediction points are an array of shape like Training points
	# Prediction values are shaped like the training values and are the result
	# of the interpolation.

	numtrnpts = 0
	jstop = numpy.empty((len(mapdata[alpha])+1), dtype='int')
	jstop[0] = 0
	for i in range ( 0,len(mapdata[alpha])):
		numtrnpts += len(mapdata[alpha][i])
		jstop[i+1] = numtrnpts
	trainpts = numpy.empty((numtrnpts, 2), dtype='float')
	trainvls = numpy.empty((numtrnpts, 3, 1), dtype='float')
	speed = numpy.empty((len(mapdata[alpha])), dtype='float')
	mini = [0, 0]
	maxi = [0, 0]
	ti = 0

	# Alpha is being kept constant for each set of trainpts

	for i in range( 0,len( mapdata[alpha] ) ):
		speed[i] = mapdata[alpha][i][0][1]
		for j in range( 0,len( mapdata[alpha][i] ) ):
			# Independent Variables in mapdata: 0-alpha, 1-speed, 2-Rline
			trainpts[ti,:] = [mapdata[alpha][i][j][1], mapdata[alpha][i][j][2]]
			# Dependent Variables in mapdata: 3-Wcorr, 4-PR, 5-eff 
			trainvls[ti,:,0] = [mapdata[alpha][i][j][3], \
					mapdata[alpha][i][j][4], mapdata[alpha][i][j][5]]
			ti += 1

	# Min and maxes for each ind var needed to set up predicted pointsi
	mini[:] = [min(trainvls[:,0,0]), min(trainvls[:,1,0])]
	maxi[:] = [max(trainvls[:,0,0]), max(trainvls[:,1,0])]

	# ==========================================================================
	#  DETERMINE RESULTS FOR PREDICTED DATA
	# ==========================================================================

	# Create prediction points
	prdpts = N_Data(mini, maxi, prdpoints, 'None', prddist, 3)
	## Setup & Perform Interpolation on Predicted Pts around Training Pts
	prdvls = np.empty((prdpoints, 3), dtype='float') # 2nd dim - num of dep vars
	prdgrd = np.empty((prdpoints, 2, 3), dtype='float') # 3 ind & 3 dep var resp

	for dep in range(3): # Go through interp for each dependent variable
		if (problemtype == 'LN'):
			trainInt = LNInterp(trainpts, trainvls[:,dep,:], NumLeaves)
			prdvls[:,dep], prdgrd[:,:,dep] =  \
					trainInt(prdpts.points)
		elif (problemtype == 'WN'):
			trainInt = WNInterp(trainpts, trainvls[:,dep,:], NumLeaves)
			prdvls[:,dep], prdgrd[:,:,dep] =  \
					trainInt(prdpts.points, neighbors, DistanceEffect)
		elif (problemtype == 'CN'):
			trainInt = CNInterp(trainpts, trainvls[:,dep,:], NumLeaves)
			prdvls[:,dep], prdgrd[:,:,dep] =  \
					trainInt(prdpts.points, neighbors)
		elif (problemtype == 'HN'):
			trainInt = HNInterp(trainpts, trainvls[:,dep,:], NumLeaves)
			prdvls[:,dep], prdgrd[:,:,dep] =  \
					trainInt(prdpts.points, neighbors, tension, bias, tight)
	return jstop, speed, trainpts, trainvls, prdpts, prdvls
	# THIS ENDS THE DEFINITION OF interp_compressor() 


def plot_compressor():
# ==========================================================================
#  each compressor plot consists of the following elements based on map type:
#								 BETA			   RLINE
#	   a regular plot of	  beta lines		 speed lines
#	   a contour plot of	  efficiency		  efficiency
#	   a contour plot of	 speed lines			 R lines
#	   a scatter plot of	  op. points		  op. points
# ==========================================================================
	if scaled_map == 1:
		WC  = vals[:,0]*scalars[0]
		PR  = (vals[:,1]-1.)*scalars[1] + 1.
		EFF = vals[:,2]*scalars[2]
		if maptype == 'BETA':
			 NC = points[:,1]*scalars[3]
		if maptype == 'RLINE':
			 RL = points[:,1]*scalars[3]
	else:
		WC  = vals[:,0]
		PR  = vals[:,1]
		EFF = vals[:,2]
		if maptype == 'BETA':
			 NC = points[:,1]
		if maptype == 'RLINE':
			 RL = points[:,1]

	print WC.shape
	print PR.shape
	print EFF.shape

		# ======================================================================
		#  PLOT BETA LINES OR SPEED LINES
		# ======================================================================
	x = WC.flatten()
	y = PR.flatten()
	for i in range(len(speed)):
		if maptype == 'BETA':
			if show_lines == 1:
				if overlay == 0:
					pylab.plot( x, y, linewidth=0.5, linestyle = '--',
							color=linecolors[0] )
					pylab.text( x[len(x)-1], y[len(y)-1],
							'beta '+str( points[:,1] ),
							ha='left', va='bottom',
							fontsize='8', color='red', rotation=35 )
				else:
					pylab.plot( x, y, linewidth=0.5, linestyle = '--',
							color=linecolors[alpha] )
		if maptype == 'RLINE':
			if overlay == 0:
				pylab.plot( x[jstop[i]:jstop[i+1]], y[jstop[i]:jstop[i+1]],
						color='blue' )
				pylab.text( x[jstop[i]], y[jstop[i]],
						'speed '+str( speed[i] ), ha='right',
						va='bottom', fontsize='8', color='blue', rotation=-35 )
			else:
				pylab.plot( x, y, color=linecolors[alpha] )


	if overlay == 0:
		# ======================================================================
		#  PLOT EFFICIENCY CONTOURS ONLY ON NON-OVERLAID MAPS
		# ======================================================================
		if Veff != []:
			if filled == 0:
				pylab.contour(WC,PR,EFF,Veff)
				pylab.colorbar(ticks=Veff)
			else:
				pylab.contourf(WC,PR,EFF,Veff)
				pylab.colorbar(ticks=Veff)
		else:
			if filled == 0:
				pylab.contour(WC,PR,EFF)
			else:
				pylab.tricontour(WC.flatten(),PR.flatten(),EFF.flatten(),
						500, restride=1, cstride=1)

		# ======================================================================
		#  PLOT SPEED CONTOURS FOR BETA MAPS
		# ======================================================================
		if maptype == 'BETA':
			if Vspd != []:
				CS = pylab.contour(WC,PR,NC,Vspd,colors='blue')
				pylab.clabel(CS,Vspd,colors='blue',fmt='%3.0f')
			else:
				pylab.contour(WC,PR,NC,colors='blue')

		# ======================================================================
		#  PLOT RLINE CONTOURS FOR RLINE MAPS
		# ======================================================================
		if maptype == 'RLINE' and show_lines == 1:
			if Vspd != []:
				CS = pylab.contour(WC,PR,RL,Vspd, linewidths=0.5,
						colors='green')
				pylab.clabel(CS,Vspd,colors='green',fmt='%1.2f')
	#	   else:
		#	   pylab.tricontour(WC.flatten(),PR.flatten(),RL.flatten(), 
	#				   500, restride=1, cstride=1, 
	#				   linewidths=0.5, colors='green')

	else:
		# ======================================================================
		#  PLOT SPEED CONTOURS FOR BETA MAPS
		# ======================================================================
		if maptype == 'BETA':
			if Vspd != []:
				CS = pylab.contour(WC,PR,NC,Vspd,colors=linecolors[alpha])
				pylab.clabel(CS,Vspd,colors=linecolors[alpha],fmt='%3.0f')
			else:
				pylab.contour(WC,PR,NC,colors=linecolors[alpha])

		# ======================================================================
		#  PLOT RLINE CONTOURS FOR RLINE MAPS
		# ======================================================================
		if maptype == 'RLINE' and show_lines == 1:
			if Vspd != []:
				CS = pylab.contour(WC,PR,RL,Vspd,colors=linecolors[alpha])
				pylab.clabel(CS,Vspd,colors=linecolors[alpha],fmt='%1.2f')
			else:
				pylab.contour(WC,PR,RL,colors=linecolors[alpha])


	# ==========================================================================
	#  PLOT THE OPERATING POINTS
	# ==========================================================================
	pntx=[]  # array of corrected flow values for all of the saved points
	pnty=[]  # array of PR values for all of the saved points
	for p in range(0,len(points_data)):
		pntx.append( points_data[p][component][4] )
		pnty.append( points_data[p][component][5] )


	# Seidel attempt at Operating Point data labels  
	for p in range(0,len(points_data)):
		pointlabel= points_data[p][component][1]
		if (p/4.) == int(p/4.):
			yfact=1.08
			pylab.text( pntx[p], pnty[p]*yfact, pointlabel,
					fontsize='8', color='black', rotation=0 )
			pylab.plot( pntx[p], pnty[p], 'ks' , ms=4.0 )
		elif (p/3.) == int(p/3.):
			yfact=1.04
			pylab.text( pntx[p], pnty[p]*yfact, pointlabel,
					fontsize='8', color='blue', rotation=0 )
			pylab.plot( pntx[p], pnty[p], 'bs' , ms=4.0 )
		elif (p/2.) == int(p/2.):
			yfact=0.96
			pylab.text( pntx[p], pnty[p]*yfact, pointlabel,
					fontsize='8', color='green', rotation=0 )
			pylab.plot( pntx[p], pnty[p], 'gs' , ms=4.0 )
		else:
			yfact=0.92
			pylab.text( pntx[p], pnty[p]*yfact, pointlabel,
					fontsize='8', color='white', rotation=0 )
			pylab.plot( pntx[p], pnty[p], 'ws' , ms=4.0 )

	# ==========================================================================
	#  SET PLOT AXES IF SUPPLIED, LABEL AXES AND TITLE 
	# ==========================================================================
	if axes != []:
		pylab.axis( axes )

	pylab.xlabel('Wcorr')
	pylab.ylabel('PR')
	if multiple_alphas == 1 and overlay == 1:
		pylab.title(mapname+' MAP: alpha = all'+title)
	else:
		pylab.title(mapname+' MAP: alpha = '+str(mapdata[alpha][0][0][0])+title)

	# THIS ENDS THE DEFINITION OF plot_compressor()
	
