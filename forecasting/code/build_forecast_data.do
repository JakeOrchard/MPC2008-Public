**** BUILD_FORECAST_DATA.DO

***
*** Builds data set for forecasting model for 2008

*** Required data:
*       Download from FRED or use previous download saved in freddata_for_forecasting.dta 
*       Rebate data from rebates.xlsx
*       Auxiliary variables from auxiliary_forecasting_data.xlsx
*       Michigan survey data from clean_michigan_small.dta

***  Output:  completedata_for_forecasting.dta

** Valerie Ramey, revised April 14, 2022
*Update May 26,2023 J. Orchard -- Adds in ND series based on JPS 2006 definition

***************************************************************************************************
version 16
drop _all
clear all

set more 1

capture log close

set scheme s1color


********************************************************************************
* I. FRED DATA IMPORT - Uncomment if you want to update the FRED data
********************************************************************************

*Download data from FRED
set fredkey INSERT FREDKEY HERE

import fred DSPI UNRATE PCE PCEND PCES PCEDG DNRGRC1M027SBEA PMSAVE PSAVERT PCEPI  ///
  DNDGRG3M086SBEA DSERRG3M086SBEA DDURRG3M086SBEA PCEPILFE DNRGRG3M086SBEA DFXARG3M086SBEA ///
  UMCSENT FEDFUNDS WTISPLC GS3M GS10 
gen mdate = mofd(daten)
tsset mdate, m
order mdate
drop daten datestr

drop if mdate<m(1959m1)

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

*Merge in VIX index, which starts as a daily index 
preserve 
import fred VIXCLS, clear
gen mdate = mofd(daten)
collapse (mean) vix = VIXCLS, by(mdate)
tempfile vix
save `vix'
restore

merge 1:1 mdate using `vix', keep(1 3) nogen


gen ncndsv = ncnd + ncsv
gen ncxnrg = ncons - ncnrg

********************************************************************************
* II. Label FRED data
********************************************************************************

label var ndisp_income "nominal disposable income"
label var ur "unemployment rate"
label var ncons "nominal total consumption expenditures"
label var ncnd "nominal nondurable consumption expenditures"
label var ncsv "nominal services consumption expenditures"
label var ncdur "nominal durables consumption expenditures"
label var ncnrg "nominal energy goods and services consumption expenditures"
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
label var vix "CBOE VIX Index"

save ../output/freddata_for_forecasting.dta, replace

clear


********************************************************************************
* III. AUXILIARY DATA IMPORT - 
********************************************************************************

* Rebates

import delimited "../input/rebates.csv", encoding(ISO-8859-2)  
gen mdate = m(1959m1) + _n-1
tsset mdate, m

label var nrebate "nominal rebate"

sort mdate
tempfile rebate
save `rebate'

clear

* S&P 500, NBER recession dates, Gilchrist-Zakrajsek ebp, Ramey-Vine variables
import delimited "../input/auxiliary_forecasting_data.csv" 
gen mdate = m(1959m1) + _n-1
tsset mdate, m

gen npgasrv = npgas*rvfactor

label var nstockprice "S&P 500, nominal, from Shiller"
label var npgas "PCE deflator, gasoline and other motor fuel"
label var rvfactor "Ramey-Vine multiplicative factor for rationing costs"
label var npgasrv "Gas prices augmented with rationing costs"
label var umcsent_cargas "UM Consumer sentiment, bad time to buy car b/c gas price or rationing"
label var gz_spr "Gilchrist-Zakrajsek spread"
label var ebp "Gilchrist-Zakrajsek version of ebp"

sort mdate
tempfile auxiliary
save `auxiliary'

clear

* Michigan consumer survey data/rebates

use ../input/clean_michigan_small.dta

keep date PFE E12M DUR_R_ALL VEH_R_ALL

rename date mdate
tsset mdate, m

rename PFE pfe
rename E12M e12m
rename DUR_R_ALL dur_r_all
rename VEH_R_ALL veh_r_all

sort mdate
tempfile michigan
save `michigan'

clear

*JPS ND Monthly series
import excel "../input/JPS_consumption_rebate.xlsx", sheet("monthly") firstrow
gen mdate = m(1959m1) + _n-1

tsset mdate, m

label var ncndur_jpscat "nominal consumption, JPS nondurable categories, monthly"
label var nrebate "nominal rebate 2001, monthly"

rename nrebate nrebate2001
keep mdate nrebate2001 ncndur_jpscat
tempfile jpsnd
save `jpsnd'

use ../output/freddata_for_forecasting.dta

********************************************************************************
* IV. Create real nondur + services consumption and deflator using Whelan's method
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

label var pcndsv "price deflator for nondurables + services consumption expenditures"

drop rcnd rcsv rcndsv /* will create them later with other similar variables */

********************************************************************************
* V. Merge data and create variables for forecasting
********************************************************************************

merge 1:1 mdate using `rebate', nogen
merge 1:1 mdate using `auxiliary', nogen
merge 1:1 mdate using `michigan', nogen
merge 1:1 mdate using `jpsnd', nogen


gen term_spr = tbond10y - tbill3m
label var term_spr "10-year - 3-month spread"

* Dummy for Lehman Brothers

gen lehman = 0
replace lehman = 1 if mdate==m(2008m9)

label var lehman "dummy for Lehman Brothers bankruptcy"

*Dummy for 9/11
gen m911 = mdate == m(2001m9)

* Convert to monthly

foreach var in ndisp_income ncons ncnd ncsv ncndsv ncdur ncnrg {
  	replace `var' = `var'/12
  }

* take logs
foreach var in ndisp_income ncons ncnd ncsv ncndsv ncdur ncnrg pcons pcnd pcsv pcndsv pcdur ///
  pcnrg npoil nstockprice npgas npgasrv  ncndur_jpscat{
  	gen l`var' = ln(`var')
	label var l`var' "log of `var'"
  }

* create real consumption
foreach var in cons cnd csv cndsv cdur {
	su p`var' if mdate == ym(2008,1)
	replace p`var' = p`var' / r(mean) * 100
	gen r`var' = 100*n`var'/p`var'
	gen lr`var' = ln(r`var')
	
	label var r`var' "real `var'"
	label var lr`var' "log real `var'"
}



* create other real variables using pce deflator

foreach var in disp_income poil stockprice pgas pgasrv{
	gen lr`var' = ln`var' - lpcons
	label var lr`var' "log real `var'"
}


save ../output/completedata_for_forecasting.dta, replace










