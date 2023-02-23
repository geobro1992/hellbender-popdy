# hellbender-popdy
Hellbender Population Dynamics Model

Contributors: G. C. Brooks and H. K. Kindsvater

This repository contains code to conduct a population viability analysis for hellbenders. A brief description of each scripts contents is provided below.

## hellbenders_vital_rates_and_deterministic_projections.R
This file contains the numbers and their associated sources used to parameterize hellbender vital rate functions, as well as a deterministic demographic model. The code can be used to simulate population projections under various levels of nest failure. 

## hellbenders_nest_availability_and_fate.R
This file contains the foundational demographic model from the first file and extends it to include stochasticity. Both juvenile and adult mortality can vary randomly across years. The simulations quantify extinction risk across a range of nest availability and nest failure rates.  

## hellbenders_nest_fate_and_recruitment.R
This file contains the foundational demographic model from the first file and extends it to include stochasticity. Both juvenile and adult mortality can vary randomly across years. The simulations quantify extinction risk across a range of nest failure rates and the degree of density dependence in larval recruitment.  
