ONE-DIMENSIONAL METHANE ENGINE COOLING ANALYSIS (OMECA)
-------------------------------------------------------

This is the repository of the OMECA program. It contains the following files:
- CEA.py 
    Python file that reads CEA output files
- CIRAcalibrated.py, CIRAflux.py, CIRAvalidation.py
    Python files that contain OMECA calculations of the Hyprob engine (see Ch. 3 of thesis)
- copper-design.py
    Python file that contains OMECA calculation of generic 10 kN engine with copper wall 
    (see Ch. 3 of thesis)
- nozzleContour.csv
    Nozzle contour definition of generic 10 kN engine
- thermoClass.py
    Python file that contains additional definitions. Most notably, thermodynamic model of methane
- Thesis Luka Denies
    The MSc thesis that contains further information and comparison of OMECA to CFD results
    
There are two additional directories:
- dataCea  
    Contains output files of CEA for the Hyprob and generic 10 kN engines
- dataCira
    Contains data for Hyprob engine, digitised from papers by Pizzarelli et al.


