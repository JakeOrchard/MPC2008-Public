**** REBATE2008_COUNTER_MON.DO

***
*** STUDIES THE COUNTERFACTUAL OF NO STIMULUS, USING BRODA-PARKERPSJM 2013 AER ESTIMATES OF TIMING OF CONSUMPTION RESPONSE TO THE 2008 TAX REBATE

** ORIGINALLY FOR GOV JEP, UPDATED FOR TAYLOR DISCUSSION, THEN MPC AUGUST 2021

** Valerie Ramey, revised September 27, 2021

***************************************************************************************************


drop _all
clear all

set more 1

capture log close

set scheme s1color

/*
********************************************************************************
* I. BEA DATA IMPORT - Already done and saved, so can leave commented out until need revised data
********************************************************************************

* First import file containing both the rebate and the deflator for nd + sv
import excel rebates.xlsx, first sheet("rebate")
gen mdate = m(1959m1) + _n-1
tsset mdate, m

sort mdate
tempfile rebate
save `rebate'

clear

*Download data from FRED
set fredkey b2d9da12107485ec086c5691d41d626c

import fred DSPI PCE PCEND PCES PCEDG DNRGRC1M027SBEA PMSAVE PSAVERT PCEPI  ///
  DNDGRG3M086SBEA DSERRG3M086SBEA DDURRG3M086SBEA PCEPILFE DNRGRG3M086SBEA DFXARG3M086SBEA ///
  UMCSENT
gen mdate = mofd(daten)
tsset mdate, m
order mdate
drop daten datestr

drop if mdate<m(1959m1)

merge 1:1 mdate using `rebate'

rename DSPI ndisp_income
rename PCE ncons
rename PCEND ncnd
rename PCES ncsv
rename PCEDG ncdur
rename DNRGRC1M027SBEA ncnrg
rename PMSAVE nsaving
rename PSAVERT nsavingrt
rename PCEPI pcons
rename DNDGRG3M086SBEA pcnd
rename DSERRG3M086SBEA pcsv
rename DDURRG3M086SBEA pcdur
rename PCEPILFE pcxnrgfd
rename DFXARG3M086SBEA pcfood
rename DNRGRG3M086SBEA pcnrg

rename UMCSENT umcsent

gen ncndsv = ncnd + ncsv
gen ncxnrg = ncons - ncnrg

label var ndisp_income "nominal disposable income"
label var ncons "nominal consumption expenditures"
label var ncnd "nominal nondurable consumption expenditures"
label var ncsv "nominal services consumption expenditures"
label var ncdur "nominal durables consumption expenditures"
label var ncnrg "nominal energy goods and services consumption expenditures"
label var nsaving "nominal personal saving"
label var nsavingrt "nominal personal saving rate"
label var pcons "price deflator for consumption expenditures"
label var pcnd "price deflator for nondurable consumption expenditures"
label var pcsv "price deflator for services consumption expenditures"
label var pcndsv "price deflator for nondurable+services consumption expenditures"
label var pcdur "price deflator for durables consumption expenditures"
label var pcxnrgfd "price deflator consumption excluding food and energy"
label var pcfood "price deflator for consumption, food"
label var pcnrg "price deflator for energy goods and services consumption"
label var ncxnrg "nominal consumption less energy goods and services"
label var umcsent "U of Michigan Consumer Sentiment"
label var nrebate "nominal rebate"

save bea_consumption.dta, replace
*/

********************************************************************************
* II.  IMPORT SIMULATION IRFS FROM NK MODEL, MANIPULATE AND SAVE
********************************************************************************

*  Variables in excel file in order: rebate, g, y, c, co, cr, t, to, tr, h, ho, hr, w, iv, u, rk, r, pi, k, b
*  see .mod programs for definitions

*Rule-of-thumb (ROT) fraction = gamma = (Excel sheet # - 1)/10

* Choose NK model specification
* NK model possibilities:  baseline, taylor_rule; also wlags variations

local model baseline

forvalues i = 1/10 {

  import excel ///
  rebatesim_`=`i'-1'=A gsim_`=`i'-1'=B ysim_`=`i'-1'=C csim_`=`i'-1'=D cosim__`=`i'-1'=E ///
  crsim_`=`i'-1'=F tsim__`=`i'-1'=G tosim_`=`i'-1'=H trsim_`=`i'-1'=I hsim_`=`i'-1'=J ///
  hosim_`=`i'-1'=K hrsim_`=`i'-1'=L wsim_`=`i'-1'=M ivsim_`=`i'-1'=N usim_`=`i'-1'=O ///
  rksim_`=`i'-1'=P rsim_`=`i'-1'=Q pisim_`=`i'-1'=R ksim_`=`i'-1'=S bsim_`=`i'-1'=T ///
      using `model'_irfs.xlsx, sheet("Sheet`i'")
	  
  gen h = _n - 2  /* h = -1 is steady state, h = 0 is month of shock */

  
  * Normalize by steady-state value (exclude rates r rk pi)
 * foreach var in g y c cr t to tr h ho hr w iv u k b {
   foreach var in y c {
  	  scalar `var'sim_ss = `var'sim_`=`i'-1'
      replace `var'sim_`=`i'-1' = `var'sim_`=`i'-1'/`var'sim_ss
  }
 
  
  label var rebatesim_`=`i'-1' "simulated tax rebate"
  label var ysim_`=`i'-1' "simulated real output (=gap), gamma = .`=`i'-1'"
  label var csim_`=`i'-1' "simulated real consumption, gamma = .`=`i'-1'"
  label var pisim_`=`i'-1' "simulated inflation, gamma = .`=`i'-1'"
  
  tempfile junk`i'
  save `junk`i''
  clear

}

use `junk1'

forvalues i = 2/10 {
	merge 1:1 h using `junk`i'', nogen
}

gen mdate = m(2008m1) + _n-1
tsset mdate, m

/* Create simulated price level p from inflation pi, then annualize gross inflation 
   and gross nominal interest rate */
   
forvalues i = 0/9 {
	gen psim_`i' = 1 if mdate==m(2008m1)
	replace psim_`i' = pisim_`i'*L.psim_`i' if mdate>=m(2008m2)
	    label var psim_`i' "simulated price, gamma = .`i'"
	  gen pisim_ann_`i' = pisim_`i'^12
	   label var pisim_ann_`i' "annualized gross inflation rate, gamma = .`i'"	    
	  gen rsim_ann_`i' = rsim_`i'^12
	    label var rsim_ann_`i' "annualized gross nominal interest rate, gamma = .`i'" 
}


tempfile junk
save `junk'

clear

********************************************************************************
* III.  USE BEA DATA AND PRODUCE GRAPHS OF CONSUMPTION, ETC. DURING 2008
********************************************************************************

use ../input/bea_consumption.dta

foreach var in cons rebate {
  replace n`var' = n`var'/12 /*put on monthly actual basis*/
  gen r`var' = n`var'*100/pcons  /*real*/
}

foreach var in cnd csv cdur cndsv {
  replace n`var' = n`var'/12 /*put on monthly actual basis*/
  gen r`var' = n`var'*100/p`var'  /*real*/
}



label var ncons "nominal consumption"
label var rcons "real consumption"
/*
tw (scatter ncons mdate, c(l) clp(l) ms(i) clw(medthick)) ///
   (scatter rcons mdate, c(l) clp(-) ms(i) clw(medthick) yaxis(2)) ///
   if mdate>=m(2007m1) & mdate<=m(2009m12), ///
   xtitle("month") ysize(4) xsize(7) scale(1.1) name(cons) ///
   title("Nominal and Real NIPA Personal Consumption Expenditures", size(medsmall))

tw scatter nrebate mdate ///
   if mdate>=m(2008m1) & mdate<=m(2008m12),  c(l l l l l l) clp(l - - - - -) ms(d i i i i i) ///
   clw(medthick medthick medthick medthick medthick medthick)  ytitle("billions, monthly rates") ///
   xtitle("month") ysize(4) xsize(7) scale(1.1) name(rebate) ///
   title("2008 Rebates (nominal)", size(medsmall)) 
*/

/*
foreach var in cons cnd csv cdur cnrg cxnrgfd cfood {
	gen lp`var' = ln(p`var')
 tw scatter lp`var' mdate ///
   if mdate>=m(2008m1) & mdate<=m(2008m12),  c(l l l l l l) clp(l - - - - -) ms(d i i i i i) ///
   clw(medthick medthick medthick medthick medthick medthick)  ytitle("log price index") ///
   xtitle("month") ysize(4) xsize(7) scale(1.1) name(lp`var') ///
   title("2008 log price `var'", size(medium)) 
}

*graph combine lpcons lpcxnrgfd lpcnrg, col(1) ysize(8) xsize(6) iscale(0.6)

graph combine lpcnd lpcsv lpcdur , col(1) ysize(8) xsize(6) iscale(0.6)

label var lpcons "total consumption"
label var lpcxnrgfd "consumption ex. food and nrg"
label var lpcnrg "energy goods & services consumption"

/*
tw (scatter lpcons lpcxnrgfd mdate, c(l l) clp(l -) ms(d d) clw(medthick medthick)) ///
   (scatter lpcnrg mdate, yaxis(2) c(l l) clp(-) ms(d i) clw(medthick medthick) clc(purple) mc(purple)) ///
   if mdate>=m(2008m1) & mdate<=m(2008m12), ///
   xtitle("month") ysize(4) xsize(7) scale(1.1) name(combo1) ///
   title("2008 Log PCE Prices", size(medsmall)) 
*/

 tw scatter lpcons lpcxnrgfd mdate ///
   if mdate>=m(2008m1) & mdate<=m(2008m12),  c(l l l l l l) clp(l - - - - -) ms(d i i i i i) ///
   clw(medthick medthick medthick medthick medthick medthick)  ytitle("log price index") ///
   xtitle("month") ysize(4) xsize(7) scale(1.1) name(combo1) ///
   title("2008 Log PCE Prices", size(medsmall)) 
*/


********************************************************************************
* III.  MERGE BEA DATA AND MODEL IRFS AND CREATE COUNTERFACTUAL
********************************************************************************
merge 1:1 mdate using `junk', nogen

forvalues i = 0/9 {
    gen ncounter`i' = ncons/(psim_`i'*csim_`i')
    gen rcounter`i' = rcons/csim_`i'
	label var ncounter`i' "mpc=.`i'"
	label var rcounter`i' "mpc=.`i'"
}

   label var rcons "actual"
   label var ncons "actual"


tw scatter ncons ncounter1 ncounter2 ncounter3 ncounter4 ncounter5 ncounter6 ncounter7 mdate ///
   if mdate>=m(2008m1) & mdate<=m(2008m12),  c(l l l l l l l l) clp(l - - - - - - -) ///
   ms(d i i i i i i i) ///
   clw(medthick medthick medthick medthick medthick medthick medthick medthick)  ytitle("billions, monthly rates") ///
   xtitle("month") ysize(4) xsize(7) scale(1.1) name(ncounter) ylabel(770 790 810 830 850, grid) ///
   title("Rule-of-Thumb Households Consume 100% in Current Month", size(medsmall)) ///
   subtitle("`model'")
graph export ../output/mpc_counterfactual.eps, replace   
 
*  title("Actual and Counterfactual Nominal Consumption", size(medsmall)) ///
*   subtitle("(NK Model, assumes rule-of-thumb HH MPC = 1 in current month)", size(small))
 
 /*
local spec wlags

tw scatter ncons ncounter1 ncounter2 ncounter3 ncounter4 ncounter5 mdate ///
   if mdate>=m(2008m1) & mdate<=m(2008m12),  c(l l l l l l) clp(l - - - - -) ms(d i i i i i) ///
   clw(medthick medthick medthick medthick medthick medthick)  ytitle("billions, monthly rates") ///
   xtitle("month") ysize(4) xsize(7) scale(1.1) name(ncounter) ylabel(800 810 820 830 840 850, grid) ///
   title("Rule-of-Thumb Households Consume 2/3 in Month 0, 1/6 in 1, 1/6 in 2", size(medsmall)) ///

graph combine ncountercontemp ncounterwlags, col(1) xsize(5) ysize(7) name(nominal)



local spec contemp

tw scatter rcons rcounter1 rcounter2 rcounter3 rcounter4 rcounter5 mdate ///
   if mdate>=m(2008m1) & mdate<=m(2008m12),  c(l l l l l l) clp(l - - - - -) ms(d i i i i i) ///
   clw(medthick medthick medthick medthick medthick medthick)  ytitle("billions, monthly rates") ///
   xtitle("month") ysize(4) xsize(7) scale(1.1) name(rcounter) ylabel(850 860 870 880 890 900, grid) ///
   title("Rule-of-Thumb Households Consume 100% in Current Month", size(medsmall)) ///
   
local spec wlags

tw scatter rcons rcounter1 rcounter2 rcounter3 rcounter4 rcounter5 mdate ///
   if mdate>=m(2008m1) & mdate<=m(2008m12),  c(l l l l l l) clp(l - - - - -) ms(d i i i i i) ///
   clw(medthick medthick medthick medthick medthick medthick)  ytitle("billions, monthly rates") ///
   xtitle("month") ysize(4) xsize(7) scale(1.1) name(rcounter) ylabel(850 860 870 880 890 900, grid) ///
   title("Rule-of-Thumb Households Consume 2/3 in Month 0, 1/6 in 1, 1/6 in 2", size(medsmall)) ///

graph combine rcountercontemp rcounterwlags, col(1) xsize(5) ysize(7) name(real)
