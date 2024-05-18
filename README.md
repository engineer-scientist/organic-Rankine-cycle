# organic-Rankine-cycle

Organic Rankine cycle (ORC) is a technology for producing work by utilizing heat sources of temperature around 100 ⁰C to 200 ⁰C (e.g. solar thermal energy, geothermal energy, biomass energy, industrial waste heat, etc.).

This repository contains programs for computational analysis of ORC with different working fluids (pure and mixed), expander inlet temperatures, and pump outlet pressures. Some programs are of ORC without economizer, and others are of ORC with economizer. (Programs are in MATLAB, and invoke REFPROP for physical properties of substances used as working fluids.)

*Inputs* include variable parameters like working fluid, expander inlet temperature, pump outlet pressure; and constant parameters like ambient temperature and component efficiencies.
*Outputs* are net work produced and overall efficiency.

This project includes systematic analysis of 75 pure substances, and sample analysis of 16 sample mixtures, as ORC working fluids.

Output graphs (3-dimensional and 2-dimensional) and ORC schematics are in the "image-folder".

This work is a part of my master's degree research project in Mechanical Engineering, at Indian Institute of Science (IISc).
I also did case study on 100 kW ORC power plant in IISc Challakere campus, and construction of lab-scale experimental setup in IISc Bangalore campus.

-----------------------------------------------------------------------------------------------------------------

*Units*

In all MATLAB programs, the physical quantities have the following units (as per the refpropm.m function of NIST).
- Pressure: kPa, 
- Temperature: K, 
- Specific enthalpy or specific work: J / kg, 
- Specific entropy: J / (K kg), 
- Mass flow rate: kg / s, 
- Power: W, 
- Economizer effectiveness is dimensionless.

Pure substances (for use as working fluids) are mentioned as per their REFPROP file names (for example, ‘CYCLOHEX’ for cyclohexane).

*Nomenclature of variables in MATLAB programs*

In these MATLAB programs, variables that store values of physical properties usually have the following name format: *name_unit*.
Examples: 
ambient_T_K is the variable storing the value of ambient temperature, in kelvins. 
h3_Jperkg is the variable storing the value of specific enthalpy of working fluid at state 3, in joules per kilogram. 
  
Short forms:
- RC: Rankine cycle, 
- WF: working fluid, 
- T: temperature, 
- P: pressure, 
- h: specific enthalpy, 
- s: specific entropy, 
- W: work.
