**** JPS_COUNTER.DO

***
*** COUNTERFACTUALS USING JOHNSON, PARKER, SOULELES ESTIMATES OF MPC ON 2001 REBATE

** REQUIRES:  baseline_2001_irfs.xlsx
**            JPS_consumption_rebate.xlsx

**            Creates freddata.dta and bea_jps_rebate_combined.dta

***************************************************************************************************

drop _all
clear all

set more 1

capture log close

set scheme s1color

capture log close
log using jps_counterfactual_results.log, replace

*global path "C:\Users\vrame\Dropbox\MPC - Valerie"
*global path "D:\Dropbox\MPC - Valerie"

********************************************************************************
* I.  IMPORT SIMULATION IRFS FROM NK MODEL, MANIPULATE AND SAVE AS SIM2001.DTA
********************************************************************************

*  Variables in excel file in order: rebate, g, y, c, co, cr, t, to, tr, h, ho, hr,
*   w, iv, u, rk, r, pi, k, b
*  see .mod programs for definitions

* Rule-of-thumb (ROT) fraction = 0.375


import delimited  using ../output/baseline_2001_irfs.csv, clear

 rename (v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 ) (rebatesim ysim csim cosim crsim tsim tosim trsim hsim hosim hrsim wsim ivsim usim rksim rsim pisim ksim bsim)
	  
  gen h = _n - 2  /* h = -1 is steady state, h = 0 is month of shock */
  
  drop if h>122
  
  * Ratio to steady-state value (exclude rates r rk pi)
  foreach var in y c co cr t to tr h ho hr w iv u k b {
  	  scalar `var'sim_ss = `var'sim
	  gen `var'sim_ss = `var'sim_ss
      gen `var'sim_ratio = `var'sim/`var'sim_ss
	  gen `var'sim_dev = `var'sim- `var'sim_ss
  } 
  
  label var rebatesim "simulated tax rebate"
  label var ysim "simulated real output (=gap), gamma = 0.375"
  label var csim "simulated real consumption, gamma = 0.375"
  label var pisim "simulated inflation, gamma = 0.375"


gen mdate = m(2001m5) + _n-1
tsset mdate, m

keep mdate csim csim_ss csim_ratio
tempfile sim2001
save `sim2001'
clear

********************************************************************************
* II. FRED DATA IMPORT AND CONVERT TO MONTHLY - Uncomment if you want to update the FRED data
********************************************************************************

*A. Download data from FRED

    freduse DSPI UNRATE PCE PCEND PCES PCEDG DNRGRC1M027SBEA PMSAVE PSAVERT PCEPI ///
	  CPIAUCSL DNDGRG3M086SBEA DSERRG3M086SBEA DDURRG3M086SBEA PCEPILFE DNRGRG3M086SBEA ///
	  DFXARG3M086SBEA UMCSENT FEDFUNDS WTISPLC GS3M GS10 CUSR0000SETA01 TERMCBAUTO48NS
	  
  gen mdate = mofd(daten)
  tsset mdate, m
  order mdate

drop daten 

rename DSPI ndisp_income
rename UNRATE ur
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
rename FEDFUNDS ffr
rename WTISPLC npoil
rename GS3M tbill3m
rename GS10 tbond10y
rename CPIAUCSL pcpi

gen ncndsv = ncnd + ncsv
gen ncxnrg = ncons - ncnrg

* B. Convert NIPA variables to monthly

  foreach var in ndisp_income ncons ncnd ncsv ncdur ncndsv ncxnrg ncnrg nsaving {
  	 replace `var' = `var'/12
  }
 

********************************************************************************
* III.  DEFLATORS AND REAL VALUES - CREATE NONDUR + SERVICES, RENORMALIZE DEFLATORS
********************************************************************************

********************************************************************************
* A. Create nondur + services deflator using Whelan's method
********************************************************************************

gen rcnd = 100*ncnd/pcnd
gen rcsv = 100*ncsv/pcsv

gen chain_curr_p = (rcnd*pcnd + rcsv*pcsv)/(L.rcnd*pcnd + L.rcsv*pcsv)
gen chain_lag_p = (rcnd*L.pcnd + rcsv*L.pcsv)/(L.rcnd*L.pcnd + L.rcsv*L.pcsv)

gen rcndsv = ncndsv if mdate==m(2012m7) /* Current BEA is in 2012 $ */

*dynamically generate values after 2012m7
replace rcndsv = L.rcndsv*sqrt(chain_curr_p*chain_lag_p) if mdate>m(2012m7)

*dynamically generate values before 2012m7, must be done manually b/c Stata gen only goes forward
gen t = _n


summ t if mdate==m(2012m7)

local t_base = r(mean)

gen backwardt = `t_base' - t

forvalues i = 1/`t_base' {
	
	replace rcndsv = F.rcndsv/sqrt(F.chain_curr_p*F.chain_lag_p) if backwardt == `i'

}

gen pcndsv = 100*ncndsv/rcndsv /* Create chained deflator */

drop rcnd rcsv rcndsv /* will create them later with other similar variables */

********************************************************************************
* B. Renormalize PCE price indexes so that they = 1 in 2001m5 (Tax cuts passed June 7, 2001)
********************************************************************************
* Renormalize PCE price indexes so that they = 1 in 2001m5, so that real = nominal

  foreach var in cons cnd csv cndsv cdur {
	
	egen p`var'_2001m5 = max(cond(mdate == m(2001m5), p`var', .))
	replace p`var' = p`var'/p`var'_2001m5
	label var p`var' "Price deflator for `var', = 1 in 2001m5"
}

********************************************************************************
* C. CREATE REAL CONSUMPTION
********************************************************************************
  gen rcons = ncons/pcons
    label var rcons "real total PCE (based on PCE deflator)"
  gen rcndsv = ncndsv/pcndsv
    label var rcndsv "consumption of nondur + services, divided by joint deflator"

********************************************************************************
* IV. Label FRED data
********************************************************************************

label var ndisp_income "nominal disposable income, monthly"
label var ur "unemployment rate"
label var ncons "nominal total consumption expenditures, monthly"
label var ncnd "nominal nondurable consumption expenditures, monthly"
label var ncsv "nominal services consumption expenditures, monthly"
label var ncdur "nominal durables consumption expenditures, monthly"
label var ncnrg "nominal energy goods and services consumption expenditures, monthly"
label var nsaving "nominal personal saving"
label var nsavingrt "nominal personal saving rate"
label var pcons "price deflator for consumption expenditures"
label var pcnd "price deflator for nondurable consumption expenditures"
label var pcsv "price deflator for services consumption expenditures"
label var pcdur "price deflator for durables consumption expenditures"
label var pcxnrgfd "price deflator consumption excluding food and energy"
label var pcfood "price deflator for consumption, food"
label var pcnrg "price deflator for energy goods and services consumption"
label var ncxnrg "nominal consumption less energy goods and services"
label var umcsent "U of Michigan Consumer Sentiment"
label var ffr "effective federal funds rate"
label var npoil "spot crude oil price, W. Texas Intermediate"
label var tbill3m "3-month treasury bill yield"
label var tbond10y "10-year treasury bond yield"
label var pcpi "consumer price index all urban all items"

label var pcndsv "price deflator for nondurables + services consumption expenditures"

label data "FRED data for 2001 rebate paper, monthly rates"
save freddata.dta, replace
clear


********************************************************************************
* V. IMPORT BEA JPS CATEGORY DATA AND REBATES
********************************************************************************
import excel "../input/JPS_consumption_rebate.xlsx", sheet("monthly") first

gen mdate = m(1959m1) + _n-1

tsset mdate, m

label var ncndur_jpscat "nominal consumption, JPS nondurable categories, monthly"
label var nrebate "nominal rebate, monthly"

merge 1:1 mdate using freddata.dta, nogen

sort mdate

*gen rcndsv = ncndsv/pcndsv

gen rcndur_jpscat = ncndur_jpscat/pcndsv
gen rrebate = nrebate/pcndsv

label data "Merge of freddata.dta with JPS_consumption_rebate.xlsx"
save bea_jps_rebate_combined.dta, replace

********************************************************************************
* VI. CREATE COUNTERFACTULAS
********************************************************************************
  merge 1:1 mdate using `sim2001', nogen


* A. MICRO COUNTERFACTUAL

     scalar mpc = 0.375 /* JPS Table 2 2SLS estimate of mpc */
     scalar mpc_tab4_contemp = 0.386 /* JPS Table 4 2SLS contemporaneous mpc */
     scalar mpc_tab4_lag = 0.273 /* JPS Table 4 2SLS implied 1-quarter lag mpc */

     * Assume 1/3rd each month of quarter

    gen induced_cf3mo = mpc*(1/3)*(rrebate + L.rrebate + L2.rrebate)
    gen rjpscat_cf3mo = rcndur_jpscat - induced_cf3mo
    gen rjpscat_cftab4 = rcndur_jpscat - mpc_tab4_contemp*(1/3)*(rrebate + L.rrebate + L2.rrebate) ///
                     - mpc_tab4_lag*(1/3)*(L3.rrebate + L4.rrebate + L5.rrebate)

    label var rcndur_jpscat "Actual"
    label var rjpscat_cf3mo "Counterfactual"
    label var rjpscat_cftab4 "Counterfactual (table 4)"


* B. GE COUNTERFACTUAL - USE MODEL SIMULATION RESULTS

    * Need scaling factor because jps categories use only part of consumption
	 *  Actual total C/Y = 0.67, model calibrated c/y = 0.62, cjps/y = 0.357
	 * cjps is 53% of total PCE in the data
	 
    local cfactor = 1  // 1 now that scaling is taken care of in dynare
    gen induced_ge = `cfactor'*(csim - csim_ss)*rcndur_jpscat/csim_ss

    gen rge_cf = rcndur_jpscat - induced_ge

    label var rge_cf "Counterfactual: Macro"

    gen dc = rcndur_jpscat*(csim-csim_ss)/csim_ss

    summ rrebate dc if mdate>=m(2001m6) & mdate<=m(2002m3)


* Main graphs

/*
* working graph
tw scatter rcndur_jpscat rjpscat_cf3mo rjpscat_cftab4 mdate if mdate>=m(2001m3) & mdate<m(2002m6), ///
  c(l l l) clp(l - _) clw(medthick medthick medthick) clc(black purple orange) mc(black purple orange) ///
  ms(i i i)  ytitle("billions of $, monthly rate") xtitle("month") ysize(4) xsize(6) ///
  title("JPS Real Nondurable Consumption - Actual vs. Counterfactual", size(medium)) name(jpscat_cf) ///
  subtitle("(Assumes equal spending across months of quarter)", size(medsmall)) 
*/

* graphs to upload to overleaf

* Micro counterfactual graph
tw scatter rcndur_jpscat rjpscat_cf3mo mdate if mdate>=m(2001m3) & mdate<m(2002m6), ///
  c(l l l) clp(l - _) clw(medthick medthick medthick) clc(black purple) mc(black purple) ///
  ms(i i i)  ytitle("billions of $, monthly rate") xtitle("month") ysize(4) xsize(6) name(jpscat_cf) 
*  title("JPS Real Nondurable Consumption - Actual vs. Counterfactual", size(medium)) 
*  subtitle("(Assumes equal spending across months of quarter)", size(medsmall)) 
   graph export ../output/fig_jpscat_cf.eps, replace

label var rjpscat_cf3mo "Counterfactual: Micro"

* Micro and GE counterfactual graph
tw scatter rcndur_jpscat rjpscat_cf3mo rge_cf mdate if mdate>=m(2001m3) & mdate<m(2002m6), ///
  c(l l l) clp(l - l) clw(medthick medthick medthick) clc(black purple green) mc(black purple green) ///
  ms(i i i)  ytitle("billions of $, monthly rate") xtitle("month") ysize(4) xsize(6) name(jpscat_gecf) 
 * title("JPS Real Nondurable Consumption - Actual vs. Counterfactuals", size(medium)) ///
 * subtitle("(Assumes equal spending across months of quarter)", size(medsmall)) 

 graph export ../output/fig_jpscat_gecf.eps, replace
/*
********************************************************************************
* VII. COMPARE TO HISTORICAL JPS CONSUMPTION DROPS
********************************************************************************

gen d3lrjpscat_cf3mo = 100*ln(rjpscat_cf3mo/L3.rjpscat_cf3mo)

gen d3lrcndur_jpscat = 100*ln(rcndur_jpscat/L3.rcndur_jpscat)

list mdate d3lrjpscat_cf3mo if mdate>=m(2001m5) & mdate<m(2002m3)

gsort d3lrcndur_jpscat

list mdate d3lrcndur_jpscat in 1/50
*/


********************************************************************************
* VIII. COMPARE TO Forecasted values
********************************************************************************
merge 1:1 mdate using ../input/forecasts2001 , nogen


label var rcndur_jpscatforA "Pessimistic Forecast"
label var rcndur_jpscatforD "Baseline Forecast"

* Micro and GE counterfactual graph with forecast
tw (scatter rcndur_jpscat rjpscat_cf3mo rge_cf rcndur_jpscatforA rcndur_jpscatforD mdate if mdate>=m(2001m3) & mdate<m(2002m6), ///
  c(l l l l l) clp(l - l l l) clw(medthick medthick medthick medthick medthick) clc(black purple green red blue) mc(black purple green red blue) ///
  ms(i i i i i) ) , ytitle("billions of $, monthly rate") xtitle("month") ysize(4) xsize(6) name(jpscat_gecf_for) 
 * title("JPS Real Nondurable Consumption - Actual vs. Counterfactuals", size(medium)) ///
 * subtitle("(Assumes equal spending across months of quarter)", size(medsmall)) 

 graph export ../output/fig_jpscat_gecf_for.eps, replace
 

********************************************************************************
* IX. Cumlative drop in consumption compared to actual
********************************************************************************
 
 
*Cumulative Log-Difference: counterfatual and forecasts to actual 

foreach var in rge_cf rjpscat_cf3mo rcndur_jpscatforA{
	
	gen d_`var' = ((`var') - (rcndur_jpscat)  )
	gen cum_d_`var' = sum(d_`var')
	*Percent change based on amount of total consumption over period
	gen _total_consumption = rcndur_jpscat 
	replace _total_consumption = 0 if mdate <= m(2001m5)
	gen total_consumption = sum(_total_consumption)
	gen pd_`var' = 100*cum_d_`var'/_total_consumption
// 	stop
	drop *total_consumption 
// 	stop
}

merge 1:1 mdate using ../input/forecasts2001_errorbands, nogen 

*Center error band around forecast 
replace EH_cumrcndur = EH_cumrcndur + pd_rcndur_jpscatforA
replace EL_cumrcndur = EL_cumrcndur + pd_rcndur_jpscatforA


tw (rarea EH_cumrcndur EL_cumrcndur mdate if mdate>=m(2001m6) & mdate<m(2002m3), fcolor(grey) fintensity(30) lcolor(grey)) (scatter pd_rge_cf pd_rjpscat_cf3mo pd_rcndur_jpscatforA mdate if mdate>=m(2001m6) & mdate<m(2002m3), ///
  c(l l l ) clp(l - l ) clw(medthick medthick medthick ) clc(green purple red) mc(green purple  red) ///
  ms(i i i i) ) , ytitle("Percent Change Cumulative Consumption") xtitle("month") ysize(4) xsize(6) legend(order(1 "Largest Historical Forecast Errors" 2 "Counterfactual: Macro" 3 "Counterfactual: Micro" 4 "Pessimistic Forecast")) name(cumulative_difference) 
 * title("JPS Real Nondurable Consumption - Actual vs. Counterfactuals", size(medium)) ///
 * subtitle("(Assumes equal spending across months of quarter)", size(medsmall)) 

graph export ../output/fig_jpscat_cumulative_forerror.eps, replace


*Pointwise errors
 tw (rarea EHrcndur ELrcndur mdate if mdate>=m(2001m6) & mdate<m(2002m3), fcolor(grey) fintensity(30) lcolor(grey))  (scatter rcndur_jpscat rjpscat_cf3mo rge_cf rcndur_jpscatforC rcndur_jpscatforB mdate if mdate>=m(2001m3) & mdate<m(2002m6), ///
  c(l l l l l) clp(l - l l l) clw(medthick medthick medthick medthick medthick) clc(black purple green red blue) mc(black purple green red blue) ///
  ms(i i i i i) ) , ytitle("billions of $, monthly rate") xtitle("month") ysize(4) xsize(6) name(jpscat_gecf_for_witherrors) 
 * title("JPS Real Nondurable Consumption - Actual vs. Counterfactuals", size(medium)) ///
 * subtitle("(Assumes equal spending across months of quarter)", size(medsmall)) 

capture log close
