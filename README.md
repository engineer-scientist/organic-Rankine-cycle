# organic-Rankine-cycle
Contains programs for computational analysis of organic Rankine cycle (ORC) with different working fluids (pure and mixed), expander inlet temperatures, and pump outlet pressures. Some programs for ORC without economizer, and others for ORC with economizer. (Programs are in MATLAB, and invoke REFPROP for physical properties of substances used as working fluids.)

Work includes systematic analysis of 75 pure substances, and sample analysis of 16 sample mixtures, as ORC working fluids.

Output graphs are also included here.

-----------------------------------------------------------------------------------------------------------------

*Units*

In all MATLAB programs, the physical quantities have the following units (as per the refpropm.m function of NIST).
Pressure: kPa. 
Temperature: K. 
Specific enthalpy or specific work: J / kg. 
Specific entropy: J / (K kg). 
Mass flow rate: kg / s. 
Power: W. 
Economizer effectiveness is dimensionless.

Pure substances (for use as working fluids) are mentioned as per their REFPROP file names (for example, ‘CYCLOHEX’ for cyclohexane).

*Nomenclature of variables in MATLAB programs*

In these MATLAB programs, variables that store values of physical properties usually have the following name format: *name_unit*.
Examples:
• ambient_T_K is the variable storing the value of ambient temperature, in kelvins. 
• h3_Jperkg is the variable storing the value of specific enthalpy of working fluid at state 3, in joules per kilogram. 
  
Short forms:
RC: Rankine cycle.
WF: working fluid.
T: temperature.
P: pressure.
h: specific enthalpy.
s: specific entropy.
W: work.
