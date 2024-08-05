
/*************************************************************************
Creates TWOWAY FE weights following Chaisemartin and D'Haultfoeuille 2020 AER
************************************************************************/

clear all
set more off

global input "../input"
global output "../output"

cap ssc install twowayfeweights

/*******************************************************************************

Step 1: load data

******************************************************************i*************/

use $input/psmjcleaninterview.dta, clear


/******************************************************************************

Run Parker style regression and first difference regression

******************************************************************************/

*Parker replication
ivreghdfe  d_l_totexp2 rbtindicator age d_num_adults d_perslt18   if insample [weight=finlwt21], cluster(cuid) absorb(intdate)  

*Change on change (Note: THIS May be WHAT THEY ARE RUNNING!)
cap gen d_rbtindicator = rbtindicator-l3.rbtindicator
ivreghdfe  d_l_totexp2 d_rbtindicator age d_num_adults d_perslt18   if insample [weight=finlwt21], cluster(cuid) absorb(intdate)  

*Add group fe (Note this assigns last rebate date for households that receive multiple rebates)
gen _rbtgroup = intdate if rbtindicator == 1
egen rbtgroup = max(_rbtgroup), by(cuid) 
replace rbtgroup = 600 if rbtgroup == . //asigns max value to never treated
drop _rbtgroup

ivreghdfe  d_l_totexp2 d_rbtindicator age d_num_adults d_perslt18   if insample [weight=finlwt21], cluster(cuid) absorb(intdate rbtgroup)  

*Remove controls (to correspond with Chaisemartin weights)
ivreghdfe  d_l_totexp2 d_rbtindicator    if insample [weight=finlwt21], cluster(cuid) absorb(intdate rbtgroup)  

/******************************************************************************
Estimate weights using twowayfeweights command
*****************************************************************************/


twowayfeweights d_l_totexp2 rbtgroup intdate d_rbtindicator if insample, weight(finlwt21) type(fdS) path(../output/cdhweights.dta)



preserve
use ../output/cdhweights.dta, clear

format Group_TWFE %tm
format Time_TWFE %tm

sort weight
restore

twowayfeweights totexp2 rbtgroup intdate rbtindicator if insample, weight(finlwt21) type(feTR) path(../output/cdhweights.dta)
