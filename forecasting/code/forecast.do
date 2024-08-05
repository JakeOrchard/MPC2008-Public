**** FORECAST.DO

***
*** Creates forecasts for 2008 using a time series model

*** Required data: 

***      completedata_for_forecasting.dta  (created by build_forecast_data.do)
                   

** Valerie Ramey, revised April 14, 2022

********************************************************************************
version 16
drop _all
clear all

set more 1

capture log close

set scheme s1color

use ../output/completedata_for_forecasting.dta

**** SAMPLE PERIOD FOR ESTIMATION 
****   starts 1984 because of documented change in oil price effects on consumption

gen sample = mdate>=m(1984m1) & mdate<=m(2019m12)

********************************************************************************
* I. CREATE FORECASTS - A, B, C, and D baseline
********************************************************************************

********************************************************************************
*A. With Lehman Brothers dummy, exogenous oil - most pessimistic
********************************************************************************

local lagcontrols L(1/6).lrcons L(1/6).lrdisp_income L(1/6).lpcons L(1/6).ebp ///
   L(1/6).recession L(1/6).lrpoil
   
local lagcontrols_norecession L(1/6).lrcons L(1/6).lrdisp_income L(1/6).lpcons L(1/6).ebp ///
    L(1/6).lrpoil

reg lrcons `lagcontrols' ebp recession lrpoil L(0/6).lehman if sample

estimates store consreg

reg lrdisp_income `lagcontrols' ebp recession lrpoil L(0/6).lehman if sample

estimates store yreg

reg lpcons `lagcontrols' ebp recession lrpoil L(0/6).lehman if sample

estimates store preg

reg ebp `lagcontrols' recession lrpoil L(0/6).lehman if sample

estimates store ebpreg

forecast create mymodel
forecast estimates consreg
forecast estimates yreg
forecast estimates preg
forecast estimates ebpreg

forecast exogenous lrpoil recession lehman

forecast solve, begin(m(2008m2)) end(m(2009m12)) suffix(forA) simulate(betas, statistic(stddev,  suffix(forA_sd)) reps(100))

********************************************************************************
* B. No Lehman Brothers dummy, exogenous oil
********************************************************************************

reg lrcons `lagcontrols' ebp recession lrpoil if sample

estimates store consreg

reg lrdisp_income `lagcontrols' ebp recession lrpoil if sample

estimates store yreg

reg lpcons `lagcontrols' ebp recession lrpoil if sample
estimates store preg

reg ebp `lagcontrols' recession lrpoil if sample

estimates store ebpreg

forecast create mymodelb, replace
forecast estimates consreg
forecast estimates yreg
forecast estimates preg
forecast estimates ebpreg
*forecast estimates oilreg

forecast exogenous lrpoil recession

forecast solve, begin(m(2008m2)) end(m(2009m12)) suffix(forB) simulate(betas, statistic(stddev,  suffix(forB_sd)) reps(100))

********************************************************************************
* C. With Lehman Brothers dummy, endogenous oil
********************************************************************************

reg lrcons `lagcontrols' ebp recession lrpoil L(0/6).lehman if sample

estimates store consreg

reg lrdisp_income `lagcontrols' ebp recession lrpoil L(0/6).lehman if sample

estimates store yreg

reg lpcons `lagcontrols' ebp recession lrpoil L(0/6).lehman if sample

estimates store preg

reg ebp `lagcontrols' recession lrpoil L(0/6).lehman if sample

estimates store ebpreg

reg lrpoil `lagcontrols' ebp recession L(0/6).lehman if sample

estimates store oilreg

forecast create mymodelc, replace
forecast estimates consreg
forecast estimates yreg
forecast estimates preg
forecast estimates ebpreg
forecast estimates oilreg

forecast exogenous recession

forecast solve, begin(m(2008m2)) end(m(2009m12)) suffix(forC) simulate(betas, statistic(stddev,  suffix(forC_sd)) reps(100))

********************************************************************************
* D. No Lehman Brothers dummy, endogenous oil
********************************************************************************

reg lrcons `lagcontrols' ebp recession lrpoil if sample

estimates store consreg

reg lrdisp_income `lagcontrols' ebp recession lrpoil if sample

estimates store yreg

reg lpcons `lagcontrols' ebp recession lrpoil if sample

estimates store preg

reg ebp `lagcontrols' recession lrpoil if sample

estimates store ebpreg

reg lrpoil `lagcontrols' ebp recession if sample

estimates store oilreg

forecast create mymodeld, replace
forecast estimates consreg
forecast estimates yreg
forecast estimates preg
forecast estimates ebpreg
forecast estimates oilreg

forecast exogenous recession

forecast solve, begin(m(2008m2)) end(m(2009m12)) suffix(forD) simulate(betas, statistic(stddev,  suffix(forD_sd)) reps(100))


********************************************************************************
* E. No Recession, endogenous oil
********************************************************************************

reg lrcons `lagcontrols_norecession' ebp lrpoil if sample

estimates store consreg

reg lrdisp_income `lagcontrols_norecession' ebp lrpoil if sample

estimates store yreg

reg lpcons `lagcontrols_norecession' ebp lrpoil if sample

estimates store preg

reg ebp `lagcontrols_norecession' lrpoil if sample

estimates store ebpreg

reg lrpoil `lagcontrols_norecession' ebp if sample

estimates store oilreg

forecast create mymodele, replace
forecast estimates consreg
forecast estimates yreg
forecast estimates preg
forecast estimates ebpreg
forecast estimates oilreg

forecast exogenous recession

forecast solve, begin(m(2008m2)) end(m(2009m12)) suffix(forE) simulate(betas, statistic(stddev,  suffix(forE_sd)) reps(100))

********************************************************************************
* F. No Recession, exogenous oil
********************************************************************************

reg lrcons `lagcontrols_norecession' ebp lrpoil if sample

estimates store consreg

reg lrdisp_income `lagcontrols_norecession' ebp lrpoil if sample

estimates store yreg

reg lpcons `lagcontrols_norecession' ebp lrpoil if sample

estimates store preg

reg ebp `lagcontrols_norecession' lrpoil if sample

estimates store ebpreg

reg lrpoil `lagcontrols_norecession' ebp if sample

estimates store oilreg

forecast create mymodelf, replace
forecast estimates consreg
forecast estimates yreg
forecast estimates preg
forecast estimates ebpreg

forecast exogenous lrpoil recession

forecast solve, begin(m(2008m2)) end(m(2009m12)) suffix(forF) simulate(betas, statistic(stddev,  suffix(forF_sd)) reps(100))


********************************************************************************
* II. GRAPH FORECASTS
********************************************************************************

*Summarize real consumption in 2007Q4 for use later
preserve
gen qdate = qofd(dofm(mdate))
sum rcons if qdate == yq(2007,4)
local rcons_0704 = r(mean)
restore


drop if mdate<m(2008m1)

*Generate 1- and 2-sd error bands for real consumption
foreach var in lrcons lpcons {
  foreach model in A B C D E F {
	gen `var'lc1for`model' = `var'for`model' - `var'for`model'_sd 
	gen `var'lc2for`model' = `var'for`model' - 2*`var'for`model'_sd 
	gen `var'hc1for`model' = `var'for`model' + `var'for`model'_sd 
	gen `var'hc2for`model' = `var'for`model' + 2*`var'for`model'_sd 
  }
}
* Normalize starting value to 0

  scalar ebp_start = ebp
  replace ebp = ebp - ebp_start
  label var ebp "Actual"
  foreach model in A B C D E F {
  	scalar ebpfor`model'_start = ebpfor`model'
	replace ebpfor`model' = ebpfor`model' - ebpfor`model'_start
	label var ebpfor`model' "Forecast `model'"
  }
 
 *	scalar `var'_start = `var'
*	replace `var' = 100*(`var' - `var'_start)
*	  foreach model in A B C D {
*	  	scalar `var'for`model'_start = `var'for`model'
*		replace `var'for`model' = 100*(`var'for`model' - `var'for`model'_start)
*	  }

foreach var in lrcons lrdisp_income lpcons {
	
	scalar `var'_start = `var'
	replace `var' = exp(`var')
	 foreach model in A B C D E F {
	  	scalar `var'for`model'_start = `var'for`model'
		replace `var'for`model' = exp(`var'for`model')
     }
	 
	label var `var' "Actual"
	label var `var'forA "Pessmistic Forecast"
	label var `var'forB "Forecast B"
	label var `var'forC "Model C"
	label var `var'forD "Forecast D"
	label var `var'forE "Regular Forecast"
	label var `var'forF "Model D"
}	




foreach var in lrconslc1 lrconslc2 lrconshc1 lrconshc2 lpconslc2 {
	 foreach model in A B C D E F {
	  	scalar `var'for`model'_start = `var'for`model'
		replace `var'for`model' = exp(`var'for`model')
     }
}

preserve
keep if mdate>=m(2008m1) & mdate<=m(2008m12)

local var lrcons
	tw scatter `var' `var'forA   `var'forE mdate ///
	if mdate>=m(2008m1) & mdate<=m(2008m12), clc(black red blue green gray) ///
  c(l l l l l l) clp(l l - _ -) ms(i i i i i i) clw(medthick medthick medthick medthick medthick) ///
  xlabel(, valuelabel) xtitle(Forecast) ytitle("billions of $, monthly rate")  ylabel(810(10)850) ///
  scale(1.2)  yscale(titlegap(*10)) xlabel(`=ym(2008,1)'(3)`=ym(2008,10)')
  graph export ../output/fig_forecasts_cons.eps, replace
  
  local var lrcons
	tw (scatter `var' `var'forA  `var'lc2forA `var'hc2forA  mdate )///
	if mdate>=m(2008m1) & mdate<=m(2008m12), clc(black red  red red) ///
  c(l l l l l l) clp(l l - - -) ms(i i i i i i) clw(medthick medthick medthick medthick medthick) ///
  xlabel(, valuelabel) xtitle(Forecast) ytitle("billions of $, monthly rate")  ylabel(810(10)850) ///
  scale(1.2)  yscale(titlegap(*10)) xlabel(`=ym(2008,1)'(3)`=ym(2008,10)') legend(order(1 2 ))
  graph export ../output/fig_forecasts_cons_wCI.eps, replace
  
  
  
local var lrcons
	tw scatter `var' `var'forA   `var'forE `var'forC `var'forF mdate ///
	if mdate>=m(2008m1) & mdate<=m(2008m12), clc(black red blue green gray) ///
  c(l l l l l l) clp(l l - _ -) ms(i i i i i i) clw(medthick medthick medthick medthick medthick) ///
  xlabel(, valuelabel) xtitle(Forecast) ytitle("billions of $, monthly rate")  ylabel(810(10)850) ///
  scale(1.2)  yscale(titlegap(*10)) xlabel(`=ym(2008,1)'(3)`=ym(2008,10)')
  graph export ../output/fig_forecasts_cons_all.eps, replace  

local var lrcons
	tw scatter `var' `var'forA   `var'forF mdate ///
	if mdate>=m(2008m1) & mdate<=m(2008m12), clc(black red blue green gray) ///
  c(l l l l l l) clp(l l - _ -) ms(i i i i i i) clw(medthick) title("Real Consumption") name(`var')

local var lrdisp_income
	tw scatter `var' `var'forA   `var'forF mdate  ///
	if mdate>=m(2008m1) & mdate<=m(2008m12), clc(black red blue green gray) ///
  c(l l l l l l) clp(l l - _ -) ms(i i i i i i) clw(medthick) title("Real Disposable Income") name(`var')
  
local var lpcons
	tw scatter `var' `var'forA   `var'forF mdate  ///
	if mdate>=m(2008m1) & mdate<=m(2008m12), clc(black red blue green gray) ///
  c(l l l l l l) clp(l l - _ -) ms(i i i i i i) clw(medthick) title("Consumption Deflator") name(`var')
  
local var ebp
	tw scatter `var' `var'forA   `var'forF mdate  ///
	if mdate>=m(2008m1) & mdate<=m(2008m12), clc(black red blue green gray) ///
  c(l l l l l l) clp(l l - _ -) ms(i i i i i i) clw(medthick)  title("Gilchrist-Zakrajsek Risk Spread") name(`var')

grc1leg lrcons lrdisp_income lpcons ebp, cols(2) ysize(4) xsize(4) iscale(0.7) ///
  legendfrom(lrcons) name(forecasts)

  graph export ../output/fig_forecasts.eps, replace
restore 

preserve
 keep mdate lrconsforA lpconsforA  
gen rcons_0704 = `rcons_0704'
 save ../output/forecasts.dta, replace  
 restore 
 replace lrconslc2forA = lrconsforA if lrconslc2forA  ==.
 replace lpconslc2forA = lpconsforA if lpconslc2forA  ==.
 keep mdate lrconslc2forA lpconslc2forA
gen rcons_0704 = `rcons_0704'
 save ../output/forecasts_lowerCI.dta, replace 
 /*
********************************************************************************
* III. SUPPLEMENT 1: SPECIFICATION CHECK: Do any additional variables lead to more 
*         pessimistic paths than Forecast A?
******************************************************************************** 
 

*Adds variable `varx' to Forecast A base, assumes additional variable is exogenous

foreach varx in umcsent nstockprice ffr {

reg lrcons `lagcontrols' L(0/6).`varx' ebp recession lrpoil L(0/2).lehman if sample

estimates store consreg

reg lrdisp_income `lagcontrols' L(0/6).`varx' ebp recession lrpoil L(0/2).lehman if sample

estimates store yreg

reg lpcons `lagcontrols' L(0/6).`varx' ebp recession lrpoil L(0/2).lehman if sample

estimates store preg

reg ebp `lagcontrols' L(0/6).`varx' recession lrpoil L(0/2).lehman if sample

estimates store ebpreg

reg `varx' `lagcontrols' L(1/6).`varx' recession lrpoil L(0/2).lehman if sample

estimates store `varx'reg

forecast create mymodel`varx', replace
forecast estimates consreg
forecast estimates yreg
forecast estimates preg
forecast estimates ebpreg
forecast estimates `varx'reg

forecast exogenous lrpoil recession lehman

forecast solve, begin(m(2008m2)) end(m(2009m12)) suffix(for`varx') 

tw scatter lrcons lrconsforA lrconsfor`varx' mdate ///
	if mdate>=m(2008m1) & mdate<=m(2008m12), ///
	 c(l l l) clp(l l - ) ms(i i i ) clw(medthick) title("Log Real Consumption, `varx' added", size(medsmall)) name(`varx')

}

graph combine umcsent nstockprice ffr
*/

/*
********************************************************************************
* IV. SUPPLEMENT 2: FORECAST E - EXCLUDE RECESSION DUMMY FROM FORECAST D SPECIFICATION
*       TO SHOW HOW OPTIMISTIC THE FORECAST IS WITHOUT THE RECESSION DUMMY 	
********************************************************************************
* E. No Lehman Brothers dummy, endogenous oil, no recession dummy
********************************************************************************

reg lrcons `lagcontrols_norecession' ebp lrpoil if sample

estimates store consreg

reg lrdisp_income `lagcontrols_norecession' ebp lrpoil if sample

estimates store yreg

reg lpcons `lagcontrols_norecession' ebp lrpoil if sample

estimates store preg

reg ebp `lagcontrols_norecession' lrpoil if sample

estimates store ebpreg

reg lrpoil `lagcontrols_norecession' if sample

estimates store oilreg

forecast create mymodele, replace
forecast estimates consreg
forecast estimates yreg
forecast estimates preg
forecast estimates ebpreg
forecast estimates oilreg

forecast exogenous recession

forecast solve, begin(m(2008m2)) end(m(2009m12)) suffix(forE) 

**************************************************************
*/
/*
********************************************************************************
* V. SUPPLEMENT 3: OIL PRICE GRAPHS
********************************************************************************

*graph combine lrcons lrdisp_income lpcons ebp, col(2) 

keep mdate rcons ncons rcndsv ncndsv rcdur ncdur lrcons lrdisp_income lpcons ebp lrpoil ///

   
label var lrconsforA "Lehman dummy, exogenous oil"
label var lrconsforB "No Lehman dummy, exgenous oil"
label var lrconsforC "Lehman dummy, endogenous real oil prices"
label var lrconsforD "No Lehman dummy, endogenous real oil prices"


label var lrpoilforC "Forecast C"
label var lrpoilforD "Forecast D"
label var lrpoil "Actual"
local var lrpoil
	tw scatter `var' `var'forC `var'forD mdate if mdate>=m(2008m1) & mdate<=m(2008m12), ///
	clc(black red blue green gray) ///
  c(l l l ) clp(l - _) ms(i i i i i) clw(medthick) title("Log Real Oil Prices") name(`var')

save forecasts.dta, replace
*/
