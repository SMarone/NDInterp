#=============================================================================
#  PYTHON SCRIPT FOR PLOTTING TURBOMACHINERY MAPS
# =============================================================================
from matplotlib import pylab
# import griddata no longer necessary
import numpy

# set some colors for 3-D maps (standard colors are b:blue, g:green, r:red, c:cyan, m:magenta, y:yellow, k:black, w:white
linecolors = [ 'red', 'green', 'blue', 'orange', 'magenta', 'cyan', 
   'saddlebrown', 'skyblue', 'olivedrab', 'yellowgreen', 'black' ]

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
#                                 BETA               RLINE
#       a regular plot of      beta lines         speed lines
#       a contour plot of      efficiency          efficiency
#       a contour plot of     speed lines             R lines
#       a scatter plot of      op. points          op. points
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
	#		else:
		#		pylab.tricontour(WC.flatten(),PR.flatten(),RL.flatten(), 
	#					500, restride=1, cstride=1, 
	#					linewidths=0.5, colors='green')

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


def plot_turbine():
   WC=[]
   PR=[]
   EFF=[]
   NC=[]

   # determination of maximum number of j values
   numJ=[]
   for i in range( 0,len( mapdata[alpha] ) ):
      numJ.append( len(mapdata[alpha][i]) )
   maxJ = max( numJ )

   # ==========================================================================
   #  for each alpha value, create PYTHON arrays used for plotting
   # ==========================================================================
   for i in range( 0,len( mapdata[alpha] ) ):
      w=[]    # array of PR values
      x=[]    # array of corrected flow values
      y=[]    # array of efficiency values
      z=[]    # array of speed values
   
      for j in range( 0,len( mapdata[alpha][i] ) ):
         if scaled_map == 1:
            w.append( (mapdata[alpha][i][j][2]-1.)/scalars[1] + 1. )
            x.append( mapdata[alpha][i][j][3]*scalars[0] )
            y.append( mapdata[alpha][i][j][4]*scalars[2] )
            z.append( mapdata[alpha][i][j][1]*scalars[3] )
         else:
            w.append( mapdata[alpha][i][j][2] )
            x.append( mapdata[alpha][i][j][3] )
            y.append( mapdata[alpha][i][j][4] )
            z.append( mapdata[alpha][i][j][1] )

      # if original map data is non-square, add duplicate data
      if len( mapdata[alpha][i] ) != maxJ:
         print "WARNING: for ", mapname, ": map data list is not square, adding data."
         for xtra in range( 0, maxJ-len( mapdata[alpha][i] ) ):
            if scaled_map == 1:
               w.append( (mapdata[alpha][i][j][2]-1.)/scalars[1] + 1. )
               x.append( mapdata[alpha][i][j][3]*scalars[0] )
               y.append( mapdata[alpha][i][j][4]*scalars[2] )
               z.append( mapdata[alpha][i][j][1]*scalars[3] )
            else:
               w.append( mapdata[alpha][i][j][2] )
               x.append( mapdata[alpha][i][j][3] )
               y.append( mapdata[alpha][i][j][4] )
               z.append( mapdata[alpha][i][j][1] )


      WC.append(x)
      EFF.append(y)
      PR.append(w)
      NC.append(z)


      # =======================================================================
      #  PLOT SPEED LINES (DEFAULT VALUES OR USER-INPUT)
      #  note: different indent levels is correct
      # =======================================================================
      if overlay == 0:
         if Vspd == []:
            pylab.plot( w, x, color='blue' )
            pylab.text( w[len(w)-1], x[len(x)-1], 
            'speed '+str(mapdata[alpha][i][j][1]), ha='left', va='center', 
            fontsize='8', color='blue', rotation=0 )
      else:
         if Vspd == []:
            pylab.plot( w, x, color=linecolors[alpha] )
            pylab.text( w[len(w)-1], x[len(x)-1], 
            'speed '+str(mapdata[alpha][i][j][1]), ha='left', va='center', 
            fontsize='8', color=linecolors[alpha], rotation=0 )

   # user input speed values
   if overlay == 0:
      if Vspd != []:
         CS = pylab.contour(PR,WC,NC,Vspd,colors='blue')
         pylab.clabel(CS,Vspd,colors='blue',fmt='%3.0f')
   else:
      if Vspd != []:
         CS = pylab.contour(PR,WC,NC,Vspd,colors=linecolors[alpha])
         pylab.clabel(CS,Vspd,colors='blue',fmt='%3.0f')


   # turn these arrays into numpy arrays
   WC=numpy.array(WC)
   PR=numpy.array(PR)
   EFF=numpy.array(EFF)
   NC=numpy.array(NC)
   
   # pylab.contourf(PR,WC,EFF)

   if overlay == 0:
      # =======================================================================
      #  PLOT EFFICIENCY CONTOURS ONLY ON NON-OVERLAID MAPS
      # =======================================================================
      if Veff != []:
         if filled == 0:
            pylab.contour(PR,WC,EFF,Veff)
            pylab.colorbar(ticks=Veff)
         else:
            pylab.contourf(PR,WC,EFF,Veff)
            pylab.colorbar(ticks=Veff)
      else:
         if filled == 0:
            pylab.contour(PR,WC,EFF)
         else: 
            pylab.contourf(PR,WC,EFF)


   # ==========================================================================
   #  PLOT THE OPERATING POINTS
   # ==========================================================================
   pnty=[]    # array of corrected flow values for all of the saved points
   pntx=[]    # array of PR values for all of the saved points
   for p in range(0,len(points_data)):
      pntx.append( points_data[p][component][5] )
      pnty.append( points_data[p][component][4] )

   # Seidel attempt at Operating Point data labels  
   for p in range(0,len(points_data)):
      pointlabel= points_data[p][component][1]
      if (p/4.) == int(p/4.): 
         yfact=1.06
         pylab.text( pntx[p], pnty[p]*yfact, pointlabel, fontsize='8', color='black', rotation=0 )
         pylab.plot( pntx[p], pnty[p], 'ks' , ms=4.0 )   # 'ko' = black oval, 'bs' = blue square, ms=markerSize
      elif (p/3.) == int(p/3.):
         yfact=1.02
         pylab.text( pntx[p], pnty[p]*yfact, pointlabel, fontsize='8', color='blue', rotation=0 )
         pylab.plot( pntx[p], pnty[p], 'bs' , ms=4.0 )   # 'ko' = black oval, 'bs' = blue square, ms=markerSize
      elif (p/2.) == int(p/2.):
         yfact=0.98
         pylab.text( pntx[p], pnty[p]*yfact, pointlabel, fontsize='8', color='green', rotation=0 )
         pylab.plot( pntx[p], pnty[p], 'gs' , ms=4.0 )   # 'ko' = black oval, 'bs' = blue square, ms=markerSize
      else:
         yfact=0.94
         pylab.text( pntx[p], pnty[p]*yfact, pointlabel, fontsize='8', color='white', rotation=0 )
         pylab.plot( pntx[p], pnty[p], 'ws' , ms=4.0 )   # 'ko' = black oval, 'bs' = blue square, ms=markerSize
         
   #pylab.plot( pntx, pnty, 'bo' )   # 'ko' = black oval, 'bs' = blue square



   # ==========================================================================
   #  SET PLOT AXES IF SUPPLIED, LABEL AXES AND TITLE 
   # ==========================================================================
   if axes != []:
      pylab.axis( axes )

   pylab.ylabel('Wcorr')
   pylab.xlabel('PR')
   if multiple_alphas == 1 and overlay == 1:
      pylab.title(mapname+' MAP: alpha = all')
   else:
      pylab.title(mapname+' MAP: alpha = '+str(mapdata[alpha][0][0][0]))

   # THIS ENDS THE DEFINITION OF plot_turbine()



# =============================================================================
#  READ THE LIST OF COMPONENT MAP FILES AND THE OPERATING POINT DATA
# =============================================================================
execfile("../MapData/mapCompList.pyth.txt")
execfile("../MapData/mapOpPoints.pyth.txt")
execfile("NDimInterp.py")
pylab.spectral()

# Interpolation Inputs
prdpoints = 1000  # Number of prediction points
prddist = 'LH' # Can be rand, LH, or cart (only for 2D and 3D)
dim = 4 # Num of ind vars plus 1 dep var for each interp
problemtype = 'WN' # Can be LN, WN, CN, or HN.  View ReadMe for more info
neighbors = 20 # KD-Tree neighbors found, default ~ 1/1000 trnpoints, min 2
NumLeaves = 100 # Leaves of KD Tree, default of about 1 per 500 training points
DistanceEffect = 2 # FOR WN ONLY: Effect of distance of neighbors, default 2
tension = 0 # FOR HN ONLY: loose is -ive, tight fit is +ive, default 0
bias = 0 # FOR HN ONLY: Attention to each closest neighbor in hermite, default 0
tight = True # FOR HN ONLY: Default algorithm had only true, change if bad 
             # results with neighs << trnpoints or a high dimension (>3)

# =============================================================================
#  CREATE THE PLOTS FOR EACH COMPONENT MAP IN TURN
# =============================================================================
for component in range(0,len(component_list)-1):
   alpha = 0
   axes=[]
   Veff=[]
   Vspd=[]
   vals = []
   points = []
   # ==========================================================================
   #  READ THE DATA AND OPTIONS FOR THIS MAP
   # ==========================================================================
   execfile( component_list[component] )

   # ==========================================================================
   #  for 2-D maps: one figure (the first alpha)
   #  for 3-D maps: one figure with all alphas and no efficiency contours OR
   #     a separate figure for each alpha (based on OVERLAY option)
   # ==========================================================================
   if maptype == 'RLINE' or maptype == 'BETA':
      if multiple_alphas == 0:
         pylab.figure()
         jstop, speed, trainpts, trainvls, prdpts, prdvls = interp_compressor()
         del vals 
         del points 
         print trainvls
         print trainpts
         vals = trainvls
         points = trainpts
         title = 'Training Points'
         plot_compressor()

         pylab.figure()
         del vals 
         del points 
         vals = prdvls
         points = prdpts.points
         print prdvls
         print prdpts.points
         title = 'Predicted Points'
         plot_compressor()
   
      else:
         if overlay == 1:
            pylab.figure()
            for alpha in range( 0,len( mapdata ) ):
               jstop, speed, trainpts, trainvls, prdpts, prdvls = \
					   interp_compressor()
               del vals 
               del points 
               vals = trainvls
               points = trainpts
               title = 'Training Points'
               plot_compressor()
               # add a text label for each alpha in the upper left corner
               # axes may be different for each alpha if not specified - SMJ
               ax=pylab.axis()
               xloc = ax[0] + ( ax[1] - ax[0] )*0.05
               yloc = ax[3] - ( ax[3] - ax[2] )*0.05*(alpha+1)
               pylab.text( xloc, yloc, ('alpha='+str( mapdata[alpha][0][0][0] )), color=linecolors[alpha] )

            pylab.figure()
            for alpha in range( 0,len( mapdata ) ):
               jstop, speed, trainpts, trainvls, prdpts, prdvls = \
					   interp_compressor()
               del vals 
               del points 
               vals = prdvls
               points = prdpts.points
               pylab.figure
               title = 'Predicted Points'
               plot_compressor
   
               # add a text label for each alpha in the upper left corner
               # axes may be different for each alpha if not specified - SMJ
               ax=pylab.axis()
               xloc = ax[0] + ( ax[1] - ax[0] )*0.05
               yloc = ax[3] - ( ax[3] - ax[2] )*0.05*(alpha+1)
               pylab.text( xloc, yloc, ('alpha='+str( mapdata[alpha][0][0][0] )), color=linecolors[alpha] )
   
         else:
            for alpha in range( 0,len( mapdata ) ):
               pylab.figure()
               jstop, speed, trainpts, trainvls, prdpts, prdvls = \
					   interp_compressor()
               del vals 
               del points
               vals = trainvls
               points = trainpts
               title = 'Training Points'
               plot_compressor()

               pylab.figure()
               del vals 
               del points 
               vals = prdvls
               points = prdpts.points
               title = 'Predicted Points'
               plot_compressor()


   if maptype == 'TURBone' or maptype == 'TURBtwo':
      if multiple_alphas == 0:
         pylab.figure()
         plot_turbine()

      else:
         if overlay == 1:
            pylab.figure()
            for alpha in range( 0,len( mapdata ) ):
               plot_turbine()

               # add a text label for each alpha in the upper left corner
               # NOTE: axes may be different for each alpha if not specified - SMJ
               ax=pylab.axis()
               xloc = ax[0] + ( ax[1] - ax[0] )*0.05
               yloc = ax[3] - ( ax[3] - ax[2] )*0.05*(alpha+1)
               pylab.text( xloc, yloc, ('alpha='+str( mapdata[alpha][0][0][0] )), color=linecolors[alpha] )

         else:
            for alpha in range( 0,len( mapdata ) ):
               pylab.figure()
               plot_turbine()

# =============================================================================
#  CREATE MAPS FROM THE INTERPOLATED DATA
# =============================================================================

# =============================================================================
#  DISPLAY ALL OF THE CREATED PLOTS
# =============================================================================
pylab.show()

