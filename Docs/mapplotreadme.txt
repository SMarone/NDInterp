//
//------------------------------------------------------------------------
//                                                                       |
//   File Name:     maptest_PYTHON.mdl                                   |
//   Date(s):       May 1, 2008                                          |
//   Author:        Scott Jones                                          |
//                                                                       |
//   Description:   Model used to test/debug the mapplot() function      |
//                  Tested with NPSS_1.6.4 - Rev: AE                     |
//                                                                       |
//------------------------------------------------------------------------
//
// GENERAL NOTES - PLOTTING COMPRESSOR AND TURBINE MAPS IN NPSS
// 
// This file demonstrates how to add the capability of plotting compressor and 
// turbine maps during an NPSS run.  The plots are generated using the PYTHON
// code, which must be installed.  To incorporate map-plotting functionality 
// into an existing NPSS model, do the following:
//
// 1. ADD THE FILES
//    Add the NPSS interpreted language file, "mapplot.fnc", via a #include to 
//    your NPSS model.  Add the PYTHON code script, "mapplot.py", to your 
//    working directory or make sure it's in your path.
//
//
// 2. CHANGE ANY PLOT DEFAULTS (OPTIONAL)
//    There are several plotting options directly available to the user by 
//    declaring specific variable names and values at the map level for each
//    particular component.  Available options are:
//
//    integer SCALED_MAP: this determines whether the plotted map will be of 
//        the actual unscaled tabular data or of the scaled tabular data.  The 
//        map scalars are saved in the individual map data files.  The default
//        is TRUE.
//
//    integer MULTIPLE_ALPHAS: if set to FALSE, the data corresponding to the 
//        first alpha value only will be plotted.  If set to TRUE, the data for
//        all alpha values will be plotted.  The default is FALSE.
//
//    integer OVERLAY (used only when MULTIPLE_ALPHAS is TRUE): if set to FALSE
//        each alpha value will be shown on its own separate plot.  If TRUE,
//        all alphas will be shown on the same plot with a "grid" of speed and
//        R-lines/beta-lines of a particular color.  Efficiency contours will 
//        not be shown for clarity.  The default is FALSE.
//
//    integer FILLED: if TRUE, filled color spaces are used to mark the edges
//        of efficiency contours; otherwise colored lines are used.  The 
//        default is TRUE.
//
//    integer SHOW_LINES: if TRUE, R-lines/beta-lines are shown on the map;
//        otherwise they are not.  Applicable for compressor maps only.  The 
//        default is TRUE.
//
//    real array AXES[]: if this variable is not declared, the map axes will be
//        determined by PYTHON.  If this variable is declared (with four
//        values), the map will have its axes equal to XMIN, XMAX, YMIN, YMAX. 
//        The default is this variable is undeclared.
//
//    real array CONTOUR_SPD[]: if this variable is not declared, lines of
//        constant speed (beta-line maps) or lines of constant R-value (R-line
//        maps) will be shown at values determined by PYTHON.  If this variable
//        is declared, these lines will be at the values set in this array.  
//        The default is this variable is undeclared.
//
//    real array CONTOUR_EFF[]: if this variable is not declared, lines of
//        constant adiabatic efficiency will be shown at values determined by 
//        PYTHON.  If this variable is declared, these lines will be at the 
//        values set in this array.  The default is this variable is undeclared.
//
//
// 3. CALL THE FUNCTIONS
//    The file "mapplot.fnc" declares two functions.  The first, 
//    'createMapDataFiles()', creates several files containing the names of the 
//    model's turbomachinery components and their data.  This function is best
//    called after the model DESIGN POINT has been run, and is only called once.
//    The second function, 'saveOpPoint("")', saves the CURRENT compressor and 
//    turbine operating points to a file; this does not have to be done for 
//    every run point; do it only for those operating points of interest.  The
//    'saveOpPoint("")' function takes a string argument: this descriptor is 
//    set to "NOSAVE" on the first call if the user does not want any operating
//    points to be plotted; otherwise the descriptor must be set to "DONE" on 
//    the last call to this function.  All other times the user can leave the 
//    descriptor blank or set it to something meaningful such as "ADP", "top of 
//    climb", "RTO", etc.
//
//
// 4. PLOT THE MAPS
//    To plot the maps at a particular point in an NPSS run, do a system call
//    to the PYTHON script file: system( "mapplot.py" ); 
//    To plot the maps outside of NPSS using PYTHON, simply input "mapplot.py" 
//    in the command prompt window.  This assumes that the map data files have 
//    been previously generated.
//
//
//  GENERAL NOTES
//
//  - for the code to work best, the map data must be "square": there must be
//    the same number of "R-line" values for each speed value (R-line maps) or
//    the same number of speed values for each beta value (beta-line maps).  If 
//    this requirement is not met, extra data values will be added during 
//    plotting and the message
//       "WARNING: for mapname: map data list is not square, adding data."
//    will be issued.
//
//  - when overlaying maps of multiple alpha, it is best to specify the axes to 
//    be used; otherwise each alpha may have its own axes prior to being shown
//    on a single plot and labels can be affected.
//
//  - when specifying contours (particularly of R-line), the "edge", or min and
//    max values, should be incremented slightly (e.g., a R-line of 1.0001 
//    instead of 1.0).  This is due to the fact that machine binary 
//    representation of numbers may not be exact.
//
//  - the input data required by the PYTHON script is saved in several files:
//    "mapCompList.txt" contains a list of all turbomachinery components in 
//        the model.
//    "mapDataNAME.txt", where NAME is the name of a particular turbomachinery 
//        component, contains that component's tabular map data, its scale 
//        factors, and the option variables for that plot.
//    "mapOpPoints.txt" contains the saved operating points for all of the 
//        components.
//

setThermoPackage("GasTbl");

//------------------------------------------------------------------------
//                          user-defined elements
//------------------------------------------------------------------------
#include <mapplot.fnc>

//------------------------------------------------------------------------
//                    user-defined tables and functions
//------------------------------------------------------------------------


//------------------------------------------------------------------------
//                           output data viewers 
//------------------------------------------------------------------------
#include <maptest.view> 


//------------------------------------------------------------------------
//                            model definition 
//------------------------------------------------------------------------

Element FlowStart Inlet { 
   Pt = 15.000; 
   Tt = 550.00; 
   W = 100.00; 
} 


Element Compressor Comp1 { 
// #include <../ncp01.map>
   #include <Maps\FAN01.map>

 //S_map.alpha = 0.00;
 //S_map.RlineMap = 2.0000;
 //S_map.NcDes = 1.000;
   S_map.PRdes = 2.0000;
   S_map.effDes = 0.8500;
} 

Element Compressor Comp2 { 
// #include <../ncp23.map>
   #include <Maps\HPC02.map>

 //S_map.alpha = 0.00;
 //S_map.RlineMap = 1.8000;
 //S_map.NcDes = 1.000;
   S_map.PRdes = 2.500;
   S_map.effDes = 0.8800;
} 

Element Compressor Comp3 { 
// #include <../ncp03.map>
   #include <Maps\HPC02.map>

 //S_map.alpha = 0.00;
 //S_map.RlineMap = 2.0000;
   S_map.NcDes = 0.900;
   S_map.PRdes = 7.000;
   S_map.effDes = 0.8350;
} 


Element FuelStart FuelIn { 
   LHV = 18500.;
} 


Element Burner Burner { 
   effBase = 1.000;
   switchHotLoss = "input";
   switchBurn = "TEMPERATURE";
   dPqPBase = 0.0500;
   TtCombOut = 2800.0;
} 


Element Turbine Turb1 { 
// #include <ncp04.map>
   #include <Maps\HPT02.map>

   S_map.effDes = 0.8800;
   S_map.parmMap = 3.500;  // this is actual PR, and will change to about 3.69
 //S_map.parmGeomDes  = 1.0;
   S_map.parmMapDes = 5.500; //3.688;   // this sets map PR scalar
 //S_map.parmNcDes = 100.0;
} 

Element Turbine Turb2 { 
// #include <ncp17.map>
   #include <Maps\HPT02.map>

   S_map.effDes = 0.8700;
   S_map.parmMap = 2.000;  // this is actual PR, and will change to about 2.12
   S_map.parmGeom = 1.00;
 //S_map.parmGeomDes  = 1.0;
   S_map.parmMapDes = 2.000;     // this sets map PR scalar
 //S_map.parmNcDes = 100.0;
} 


Element Shaft Shaft1 { 
   ShaftInputPort COMP3, TURB1;
   Nmech = 20000; 
} 

Element Shaft Shaft2 { 
   ShaftInputPort COMP1, COMP2, TURB2;
   Nmech = 7000; 
} 

Element FlowEnd Exit {  } 


//------------------------------------------------------------------------
//                           component linkages
//------------------------------------------------------------------------
linkPorts( "Inlet.Fl_O"      , "Comp1.Fl_I"        , "FL1"  ) ;
linkPorts( "Comp1.Fl_O"      , "Comp2.Fl_I"        , "FL2"  ) ;
linkPorts( "Comp2.Fl_O"      , "Comp3.Fl_I"        , "FL3"  ) ;
linkPorts( "Comp3.Fl_O"      , "Burner.Fl_I"       , "FL4"  ) ;
linkPorts( "Burner.Fl_O"     , "Turb1.Fl_I"        , "FL5"  ) ;
linkPorts( "Turb1.Fl_O"      , "Turb2.Fl_I"        , "FL6"  ) ;
linkPorts( "Turb2.Fl_O"      , "Exit.Fl_I"         , "FL7"  ) ;
linkPorts( "FuelIn.Fu_O"     , "Burner.Fu_I"       , "Fuel" ) ;

//------------------------------------------------------------------------
//                           shaft connections
//------------------------------------------------------------------------
linkPorts( "Comp1.Sh_O", "Shaft2.COMP1",  "comp1wk"  ) ;
linkPorts( "Comp2.Sh_O", "Shaft2.COMP2",  "comp2wk"  ) ;
linkPorts( "Comp3.Sh_O", "Shaft1.COMP3",  "comp3wk"  ) ;
linkPorts( "Turb1.Sh_O", "Shaft1.TURB1",  "turb1wk"  ) ;
linkPorts( "Turb2.Sh_O", "Shaft2.TURB2",  "turb2wk"  ) ;


//------------------------------------------------------------------------
//                           run design point 
//------------------------------------------------------------------------
solver.maxJacobians = 25;
solver.maxIterations = 50;
setOption( "switchDes", "DESIGN" );
autoSolverSetup(); 

run(); page.display(); 



//------------------------------------------------------------------------
//       DECLARE MAP PLOTTING VARIBLES (BEFORE CREATING DATA FILES)
//  note: this can also be done within the component instantiations, above
//------------------------------------------------------------------------
Comp1.S_map.S_eff { 
 //int SCALED_MAP = TRUE;
   int MULTIPLE_ALPHAS = FALSE;
 //int OVERLAY = FALSE;
 //int FILLED = TRUE;
 //int SHOW_LINES = TRUE;
 //real AXES[] = { 30., 110., 0.6, 2.6 };
 //real CONTOUR_SPD[] = { 1.001, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.199 };
 //real CONTOUR_EFF[] = { 0.00, 0.40, 0.50, 0.60, 0.70, 0.80, 0.85, 0.90, 0.92, 0.940, 0.960 };
} 

Comp2.S_map.S_eff { 
   int SCALED_MAP = TRUE;
   int MULTIPLE_ALPHAS = TRUE;
   int OVERLAY = TRUE;
 //int FILLED = TRUE;
   int SHOW_LINES = TRUE;
   real AXES[] =  { 15., 75., 0.5, 3.5 };
 //real CONTOUR_SPD[] = { 1.001, 1.2, 1.4, 1.6, 1.8, 1.999 };
 //real CONTOUR_EFF[] = { 0.00, 0.20, 0.40, 0.60, 0.80, 1.00 };
} 

Comp3.S_map.S_eff { 
   int SCALED_MAP = TRUE;
   int MULTIPLE_ALPHAS = TRUE;
 //int OVERLAY = TRUE;
 //int FILLED = TRUE;
 //int SHOW_LINES = TRUE;
 //real AXES[] = {10., 35., 0.0,14.0 };
 //real CONTOUR_SPD[] = { 1.001, 1.10, 1.2, 1.30, 1.4, 1.499 };
 //real CONTOUR_EFF[] = { 0.00, 0.40, 0.50, 0.60, 0.70, 0.80, 0.85, 0.90, 0.92, 0.940, 0.960 };
} 

Turb1.S_map.S_eff { 
   int SCALED_MAP = TRUE;
   int MULTIPLE_ALPHAS = FALSE;
 //int OVERLAY = TRUE;
 //int FILLED = TRUE;
   real AXES[] = { 2.0, 8.0, 6.90, 7.10 };
 //real CONTOUR_SPD[] = { 1.001, 1.10, 1.2, 1.30, 1.4, 1.499 };
   real CONTOUR_EFF[] = { 0.75, 0.78, 0.80, 0.82, 0.84, 0.86, 0.87, 0.88, 0.89 };
} 

Turb2.S_map.S_eff { 
   int SCALED_MAP = TRUE;
   int MULTIPLE_ALPHAS = TRUE;
   int OVERLAY = TRUE;
 //int FILLED = TRUE;
   real AXES[] = { 1.0, 6.0, 10.0, 40.0 };
   real CONTOUR_SPD[] = { 4525., 5650., 6680. };
 //real CONTOUR_EFF[] = { 0.75, 0.78, 0.80, 0.82, 0.84, 0.86, 0.87, 0.88, 0.89 };
} 


//------------------------------------------------------------------------
//               CREATE THE MAP DATA FILES IN PYTHON FORMAT 
//------------------------------------------------------------------------
createMapDataFiles();
saveOpPoint("");


//------------------------------------------------------------------------
//                         run off-design point 
//------------------------------------------------------------------------
solver.maxJacobians = 50;
solver.maxIterations = 100;
setOption( "switchDes", "OFFDESIGN" );
autoSolverSetup(); 

run(); page.display(); 
saveOpPoint("");


//------------------------------------------------------------------------
//                           run more points 
//------------------------------------------------------------------------
Inlet.W = 80.0;
Inlet.Pt = 16.0;
Inlet.Tt = 575.0; 
Burner.TtCombOut = 3000.0; 

run(); page.display(); 
saveOpPoint("DONE");


//------------------------------------------------------------------------
//        PERFORM A SYSTEM CALL TO RUN PYTHON WITH THIS INPUT FILE
//  note: this file will access the map data files you created above
//------------------------------------------------------------------------
system("python mapplot.py");

