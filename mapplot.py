#=============================================================================
#  PYTHON SCRIPT FOR PLOTTING TURBOMACHINERY MAPS
# =============================================================================
from matplotlib import pylab
# import griddata no longer necessary
import numpy

# set some colors for 3-D maps (standard colors are b:blue, g:green, r:red, c:cyan, m:magenta, y:yellow, k:black, w:white
linecolors = [ 'red', 'green', 'blue', 'orange', 'magenta', 'cyan', 
   'saddlebrown', 'skyblue', 'olivedrab', 'yellowgreen', 'black' ]

def plot_compressor():
   # ==========================================================================
   #  each compressor plot consists of the following elements based on map type:
   #                                 BETA               RLINE
   #       a regular plot of      beta lines         speed lines
   #       a contour plot of      efficiency          efficiency
   #       a contour plot of     speed lines             R lines
   #       a scatter plot of      op. points          op. points
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

   # ==========================================================================
   #  for each alpha value, create PYTHON arrays used for plotting
   # ==========================================================================
   for i in range( 0,len( mapdata[alpha] ) ):
      w=[]    # array of speed values or Rline values
      x=[]    # array of corrected flow values
      y=[]    # array of PR values
      z=[]    # array of efficiency values

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
         print "WARNING: for ", mapname, ": map data list is not square, adding data."
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

      # =======================================================================
      #  PLOT BETA LINES OR SPEED LINES
      # =======================================================================
      if maptype == 'BETA':
         if show_lines == 1:
            if overlay == 0:
               pylab.plot( x, y, linewidth=0.5, linestyle = '--', color=linecolors[0] )
               pylab.text( x[len(x)-1], y[len(y)-1], 
               'beta '+str( mapdata[alpha][i][j][1] ), ha='left', va='bottom', 
               fontsize='8', color='red', rotation=35 )
            else:
               pylab.plot( x, y, linewidth=0.5, linestyle = '--', color=linecolors[alpha] )

      if maptype == 'RLINE':
         if overlay == 0:
            pylab.plot( x, y, color='blue' )
            pylab.text( x[0], y[0], 
            'speed '+str(mapdata[alpha][i][j][1]), ha='right', va='bottom', 
            fontsize='8', color='blue', rotation=-35 )
         else:
            pylab.plot( x, y, color=linecolors[alpha] )


   # turn these arrays into numpy arrays 
   WC=numpy.array(WC)
   PR=numpy.array(PR)
   EFF=numpy.array(EFF)
   NC=numpy.array(NC)
   RL=numpy.array(RL)

   if overlay == 0:
      # =======================================================================
      #  PLOT EFFICIENCY CONTOURS ONLY ON NON-OVERLAID MAPS
      # =======================================================================
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
            pylab.contourf(WC,PR,EFF)

      # =======================================================================
      #  PLOT SPEED CONTOURS FOR BETA MAPS
      # =======================================================================
      if maptype == 'BETA':
         if Vspd != []:
            CS = pylab.contour(WC,PR,NC,Vspd,colors='blue')
            pylab.clabel(CS,Vspd,colors='blue',fmt='%3.0f')
         else:
            pylab.contour(WC,PR,NC,colors='blue')

      # =======================================================================
      #  PLOT RLINE CONTOURS FOR RLINE MAPS
      # =======================================================================
      if maptype == 'RLINE' and show_lines == 1:
         if Vspd != []:
            CS = pylab.contour(WC,PR,RL,Vspd, linewidths=0.5, colors='green')
            pylab.clabel(CS,Vspd,colors='green',fmt='%1.2f')
         else:
            pylab.contour(WC,PR,RL, linewidths=0.5, colors='green')

   else:
      # =======================================================================
      #  PLOT SPEED CONTOURS FOR BETA MAPS
      # =======================================================================
      if maptype == 'BETA':
         if Vspd != []:
            CS = pylab.contour(WC,PR,NC,Vspd,colors=linecolors[alpha])
            pylab.clabel(CS,Vspd,colors=linecolors[alpha],fmt='%3.0f')
         else:
            pylab.contour(WC,PR,NC,colors=linecolors[alpha])

      # =======================================================================
      #  PLOT RLINE CONTOURS FOR RLINE MAPS
      # =======================================================================
      if maptype == 'RLINE' and show_lines == 1:
         if Vspd != []:
            CS = pylab.contour(WC,PR,RL,Vspd,colors=linecolors[alpha])
            pylab.clabel(CS,Vspd,colors=linecolors[alpha],fmt='%1.2f')
         else:
            pylab.contour(WC,PR,RL,colors=linecolors[alpha])

   
   # ==========================================================================
   #  PLOT THE OPERATING POINTS
   # ==========================================================================
   pntx=[]    # array of corrected flow values for all of the saved points
   pnty=[]    # array of PR values for all of the saved points
   for p in range(0,len(points_data)):
      pntx.append( points_data[p][component][4] )
      pnty.append( points_data[p][component][5] )
      
      
   # Seidel attempt at Operating Point data labels  
   for p in range(0,len(points_data)):
      pointlabel= points_data[p][component][1]
      if (p/4.) == int(p/4.): 
         yfact=1.08
         pylab.text( pntx[p], pnty[p]*yfact, pointlabel, fontsize='8', color='black', rotation=0 )
         pylab.plot( pntx[p], pnty[p], 'ks' , ms=4.0 )   # 'ko' = black oval, 'bs' = blue square, ms=markerSize
      elif (p/3.) == int(p/3.):
         yfact=1.04
         pylab.text( pntx[p], pnty[p]*yfact, pointlabel, fontsize='8', color='blue', rotation=0 )
         pylab.plot( pntx[p], pnty[p], 'bs' , ms=4.0 )   # 'ko' = black oval, 'bs' = blue square, ms=markerSize
      elif (p/2.) == int(p/2.):
         yfact=0.96
         pylab.text( pntx[p], pnty[p]*yfact, pointlabel, fontsize='8', color='green', rotation=0 )
         pylab.plot( pntx[p], pnty[p], 'gs' , ms=4.0 )   # 'ko' = black oval, 'bs' = blue square, ms=markerSize
      else:
         yfact=0.92
         pylab.text( pntx[p], pnty[p]*yfact, pointlabel, fontsize='8', color='white', rotation=0 )
         pylab.plot( pntx[p], pnty[p], 'ws' , ms=4.0 )   # 'ko' = black oval, 'bs' = blue square, ms=markerSize
      
   #pylab.plot( pntx, pnty, 'ks' , ms=4.0 )   # 'ko' = black oval, 'bs' = blue square, ms=markerSize
   
           
            
   # plot single point
   # pylab.plot( [points_data[p][component][4]] , [points_data[p][component][5]], 'ko', ms=10.0 )


   # ==========================================================================
   #  SET PLOT AXES IF SUPPLIED, LABEL AXES AND TITLE 
   # ==========================================================================
   if axes != []:
      pylab.axis( axes )

   pylab.xlabel('Wcorr')
   pylab.ylabel('PR')
   if multiple_alphas == 1 and overlay == 1:
      pylab.title( mapname + ' MAP: alpha = all' )
   else:
      pylab.title( mapname + ' MAP: alpha = ' + str( mapdata[alpha][0][0][0] ) )

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
      pylab.title( mapname + ' MAP: alpha = all' )
   else:
      pylab.title( mapname + ' MAP: alpha = ' + str( mapdata[alpha][0][0][0] ) )

   # THIS ENDS THE DEFINITION OF plot_turbine()


# =============================================================================
#  READ THE LIST OF COMPONENT MAP FILES AND THE OPERATING POINT DATA
# =============================================================================
execfile("../MapData/mapCompList.pyth.txt")
execfile("../MapData/mapOpPoints.pyth.txt")
pylab.spectral()

# =============================================================================
#  CREATE THE PLOTS FOR EACH COMPONENT MAP IN TURN
# =============================================================================
for component in range(0,len(component_list)-1):

   alpha = 0
   axes=[]
   Veff=[]
   Vspd=[]
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
         plot_compressor()
   
      else:
         if overlay == 1:
            pylab.figure()
            for alpha in range( 0,len( mapdata ) ):
               plot_compressor()
   
               # add a text label for each alpha in the upper left corner
               # NOTE: axes may be different for each alpha if not specified - SMJ
               ax=pylab.axis()
               xloc = ax[0] + ( ax[1] - ax[0] )*0.05
               yloc = ax[3] - ( ax[3] - ax[2] )*0.05*(alpha+1)
               pylab.text( xloc, yloc, ('alpha='+str( mapdata[alpha][0][0][0] )), color=linecolors[alpha] )
   
         else:
            for alpha in range( 0,len( mapdata ) ):
               pylab.figure()
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
#  DISPLAY ALL OF THE CREATED PLOTS
# =============================================================================
pylab.show()


