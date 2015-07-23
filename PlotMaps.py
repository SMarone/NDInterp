from matplotlib import pylab 
import numpy as np

problemtype = 'WN'

linecolors = [ 'red', 'green', 'blue', 'orange', 'magenta', 'cyan',
	'saddlebrown', 'skyblue', 'olivedrab', 'yellowgreen', 'black' ]

def listtounstruct(data, numpts, ilength, jlengths):
	ti = 0
	uspts = np.empty((numpts, 3), dtype='float')
	usvls = np.empty((numpts, 3, 1), dtype='float')
	
	for i in range(ilength):
		for j in range(jlengths[i]):
			# Independent Variables in mapdata: 0-alpha, 1-speed, 2-Rline
			uspts[ti,:]   = [data[i][j][0], data[i][j][1], data[i][j][2]]
			# Dependent Variables in mapdata: 3-Wcorr, 4-PR, 5-eff 
			usvls[ti,:,0] = [data[i][j][3], data[i][j][4], data[i][j][5]]
			ti += 1

	return uspts, usvls

def unstructtolist(uspts, usvls, numpts, ilength, jlengths):
	ldata = []
	i_list = []
	ti = 0
	for i in range(ilength):
		for j in range(jlengths[i]):
			i_list.append([uspts[ti,0], uspts[ti,1], uspts[ti,2], \
						usvls[ti,0,0], usvls[ti,1,0], usvls[ti,2,0]])
			ti += 1
		ldata.append(i_list)
		i_list = []

	return ldata

def plot_compressor(mapname, maptype, mapdata, scaled_map, scalars, overlay, 
					axes, filled, show_lines, Vspd, Veff): 
	# ==========================================================================
	# each compressor plot consists of the following elements based on map type:
	#								BETA					RLINE 
	#		 a regular plot of		beta lines			speed lines 
	#		 a contour plot of		efficiency			 efficiency 
	#		 a contour plot of	  speed lines				 R lines 
	#		 a scatter plot of		op. points			 op. points 
	# ========================================================================== 
	WC=[] 
	PR=[] 
	EFF=[] 
	NC=[] 
	RL=[] 
	 
	# determination of maximum number of j values 
	numJ=[] 
	for i in range( 0,len( mapdata[alpha] ) ): 
		numJ.append( len(mapdata[alpha][i]) ) 
	maxJ = max( numJ ) 
 
	#  for each alpha value, create PYTHON arrays used for plotting 
	for i in range( 0,len( mapdata[alpha] ) ):		 
		w=[]	 # array of speed values or Rline values 
		x=[]	 # array of corrected flow values 
		y=[]	 # array of PR values 
		z=[]	 # array of efficiency values 
 
		for j in range( 0,len( mapdata[alpha][i] ) ):			 
			if scaled_map == 1: 
				w.append( mapdata[alpha][i][j][2]*scalars[3] ) 
				x.append( mapdata[alpha][i][j][3]*scalars[0] ) 
				y.append(( mapdata[alpha][i][j][4]-1.)*scalars[1] + 1. ) 
				z.append( mapdata[alpha][i][j][5]*scalars[2] ) 
			else: 
				w.append( mapdata[alpha][i][j][2] )
				x.append( mapdata[alpha][i][j][3] )
				y.append( mapdata[alpha][i][j][4] )
				z.append( mapdata[alpha][i][j][5] )


		# if original map data is non-square, add duplicate data
		if len( mapdata[alpha][i] ) != maxJ:
			print "WARNING: for ", mapname, \
					": map data list is not square, adding data."
			for xtra in range( 0, maxJ-len( mapdata[alpha][i] ) ):
				if scaled_map == 1:
					w.append( mapdata[alpha][i][j][2]*scalars[3] )
					x.append( mapdata[alpha][i][j][3]*scalars[0] )
					y.append(( mapdata[alpha][i][j][4]-1.)*scalars[1] + 1. )
					z.append( mapdata[alpha][i][j][5]*scalars[2] )
				else:
					w.append( mapdata[alpha][i][j][2] )
					x.append( mapdata[alpha][i][j][3] )
					y.append( mapdata[alpha][i][j][4] )
					z.append( mapdata[alpha][i][j][5] )

		WC.append(x)
		PR.append(y)
		EFF.append(z)
		if maptype == 'BETA':
			NC.append(w)
		if maptype == 'RLINE':
			RL.append(w)

		#  PLOT BETA LINES OR SPEED LINES
		if maptype == 'BETA':
			if show_lines == 1:
				if overlay == 0:
					pylab.plot( x, y, linewidth=0.5, linestyle = '--', 
							color=linecolors[0] )
					pylab.text( x[len(x)-1], y[len(y)-1],
					'beta '+str( mapdata[alpha][i][j][1] ), ha='left', 
					va='bottom',
					fontsize='8', color='red', rotation=35 )
				else:
					pylab.plot( x, y, linewidth=0.5, linestyle = '--', 
							color=linecolors[alpha] )

		if maptype == 'RLINE':
			if overlay == 0:
				pylab.plot( x, y, color='blue' )
				pylab.text( x[0], y[0],
				'speed '+str(mapdata[alpha][i][j][1]), ha='right', va='bottom',
				fontsize='8', color='blue', rotation=-35 )
			else:
				pylab.plot( x, y, color=linecolors[alpha] )


	# turn these arrays into numpy arrays 
	WC=np.array(WC)
	PR=np.array(PR)
	EFF=np.array(EFF)
	NC=np.array(NC)
	RL=np.array(RL)

	if overlay == 0:
		#  PLOT EFFICIENCY CONTOURS ONLY ON NON-OVERLAID MAPS
		if Veff != []:
			if filled == 0:
				pylab.contour(WC,PR,EFF,Veff)
				cb = pylab.colorbar(ticks=Veff)
				cb.ax.set_ylabel('Veff')
			else:
				pylab.contourf(WC,PR,EFF,Veff)
				cb = pylab.colorbar(ticks=Veff)
				cb.ax.set_ylabel('Veff')
		else:
			if filled == 0:
				pylab.contour(WC,PR,EFF)
				cb = pylab.colorbar()
				cb.ax.set_ylabel('Eff')
			else:
				pylab.contourf(WC,PR,EFF)
				cb = pylab.colorbar()
				cb.ax.set_ylabel('Eff')

		#  PLOT SPEED CONTOURS FOR BETA MAPS
		if maptype == 'BETA':
			if Vspd != []:
				CS = pylab.contour(WC,PR,NC,Vspd,colors='blue')
				pylab.clabel(CS,Vspd,colors='blue',fmt='%3.0f')
			else:
				pylab.contour(WC,PR,NC,colors='blue')

		#  PLOT RLINE CONTOURS FOR RLINE MAPS
		if maptype == 'RLINE' and show_lines == 1:
			if Vspd != []:
				CS = pylab.contour(WC,PR,RL,Vspd, linewidths=0.5, 
						colors='green')
				pylab.clabel(CS,Vspd,colors='green',fmt='%1.2f')
			else:
				pylab.contour(WC,PR,RL, linewidths=0.5, colors='green')

	else:
		#  PLOT SPEED CONTOURS FOR BETA MAPS
		if maptype == 'BETA':
			if Vspd != []:
				CS = pylab.contour(WC,PR,NC,Vspd,colors=linecolors[alpha])
				pylab.clabel(CS,Vspd,colors=linecolors[alpha],fmt='%3.0f')
			else:
				pylab.contour(WC,PR,NC,colors=linecolors[alpha])

		#  PLOT RLINE CONTOURS FOR RLINE MAPS
		if maptype == 'RLINE' and show_lines == 1:
			if Vspd != []:
				CS = pylab.contour(WC,PR,RL,Vspd,colors=linecolors[alpha])
				pylab.clabel(CS,Vspd,colors=linecolors[alpha],fmt='%1.2f')
			else:
				pylab.contour(WC,PR,RL,colors=linecolors[alpha])


	#  PLOT THE OPERATING POINTS
	pntx=[]	 # array of corrected flow values for all of the saved points
	pnty=[]	 # array of PR values for all of the saved points
	for p in range(0,len(points_data)):
		pntx.append( points_data[p][component][4] )
		pnty.append( points_data[p][component][5] )


	# Seidel attempt at Operating Point data labels  
	for p in range(0,len(points_data)):
		pointlabel= points_data[p][component][1]
		if (p/4.) == int(p/4.):
			yfact=1.08
			pylab.text( pntx[p], pnty[p]*yfact, pointlabel, fontsize='8', 
					color='black', rotation=0 )
			pylab.plot( pntx[p], pnty[p], 'ks' , ms=4.0 )
		elif (p/3.) == int(p/3.):
			yfact=1.04
			pylab.text( pntx[p], pnty[p]*yfact, pointlabel, fontsize='8', color='blue', rotation=0 )
			pylab.plot( pntx[p], pnty[p], 'bs' , ms=4.0 )
		elif (p/2.) == int(p/2.):
			yfact=0.96
			pylab.text( pntx[p], pnty[p]*yfact, pointlabel, fontsize='8', color='green', rotation=0 )
			pylab.plot( pntx[p], pnty[p], 'gs' , ms=4.0 )
		else:
			yfact=0.92
			pylab.text( pntx[p], pnty[p]*yfact, pointlabel, fontsize='8', color='white', rotation=0 )
			pylab.plot( pntx[p], pnty[p], 'ws' , ms=4.0 )

	#  SET PLOT AXES IF SUPPLIED, LABEL AXES AND TITLE 
	if axes != []:
		pylab.axis( axes )

	pylab.xlabel('Wcorr')
	pylab.ylabel('PR')
	if multiple_alphas == 1 and overlay == 1:
		pylab.title( mapname + ' MAP: alpha = all' )
	else:
		pylab.title( mapname + ' MAP: alpha = ' + \
				str( mapdata[alpha][0][0][0] ) )

	# THIS ENDS THE DEFINITION OF plot_compressor()

def interpcheck():
	# Training points are shape [point number, independent dimension].
	# Training values are shape [point number, dependent dimension].
	trnpts, trnvls = \
			listtounstruct(mapdata[alpha], numpts, ilength, jlengths)
	## Make prediction points slightly off of training points
	prdpts = trnpts + [0.0, 0.025, 0.1]
	## Setup & Perform Interpolation on Pred Pts around Train Pts
	## Prediction points  and vals shaped like training 
	## Prediction gradient is shape [numpts, indeps, deps]
	prdvls = np.empty((numpts, 3, 1), dtype='float') 
	prdgrd = np.empty((numpts, 3, 3), dtype='float')
	## Go through the chosen interp for each dependent variable
	for dep in range(3): 
		if (problemtype == 'LN'):
			trainInt = LNInterp(trnpts, trnvls[:,dep,:], NumLeaves)
			prdvls[:,dep, 0], prdgrd[:,:,dep] =  \
					trainInt(prdpts)
		elif (problemtype == 'WN'):
			trainInt = WNInterp(trnpts, trnvls[:,dep,:], NumLeaves)
			prdvls[:,dep, 0], prdgrd[:,:,dep] =  \
					trainInt(prdpts, neighbors, DistanceEffect)
		elif (problemtype == 'CN'):
			trainInt = CNInterp(trnpts, trnvls[:,dep,:], NumLeaves)
			prdvls[:,dep, 0], prdgrd[:,:,dep] =  \
					trainInt(prdpts, neighbors)
		elif (problemtype == 'HN'):
			trainInt = HNInterp(trnpts, trnvls[:,dep,:], NumLeaves)
			prdvls[:,dep, 0], prdgrd[:,:,dep] =  \
					trainInt(prdpts, neighbors, tension, bias, tight)
		elif (problemtype == 'CR'):
			trainInt = CRInterp(trnpts, trnvls[:,dep,:], neighbors, NumLeaves, comp)
			prdvls[:,dep, 0], prdgrd[:,:,dep] =  \
					trainInt(prdpts)
	prddata.append( \
			unstructtolist( prdpts, prdvls, numpts, ilength, jlengths))
	plot_compressor(mapname, maptype, prddata, scaled_map, scalars, overlay, 
		axes, filled, show_lines, Vspd, Veff)
	## END interpcheck DEFINITION

#======= RUN CODE ========
## Read list of component map files and operating points
execfile("../MapData/mapCompList.pyth.txt")
execfile("../MapData/mapOpPoints.pyth.txt")
execfile("NDimInterp.py")
pylab.spectral()

for component in range(len(component_list)-1):

	alpha = 0
	axes = []
	Veff = []
	Vspd = []
	prddata = []

	## Read this component data
	execfile(component_list[component])
	
	ilength = len(mapdata[alpha])
	jlengths = np.zeros((ilength), dtype='int')
	for i in range(ilength):
		jlengths[i] = len(mapdata[alpha][i])
	numpts = np.sum(jlengths)


	## Start setting up plotting
	if (maptype == 'RLINE') or (maptype == 'BETA'):
		for alpha in range(len(mapdata)):
			if (alpha!= 0):
				pylab.figure()
			mapname = 'Training Data'
			plot_compressor(mapname, maptype, mapdata, scaled_map, scalars,
					overlay, axes, filled, show_lines, Vspd, Veff)
			pylab.figure()
			mapname = 'Predicted Data'
			interpcheck()


pylab.show()
