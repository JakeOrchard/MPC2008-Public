**** BUILD_NARRATIVE_DATA.DO

***
*** COMBINES DATA FROM FRED WITH REBATE DATA TO CREATE THE 2008 NARRATIVE DATASET

** Valerie Ramey, revised April 24, 2022

**  REQUIRED FILES:

*     rebates.xlsx
*     other data is downloaded directly from FRED

***************************************************************************************************

drop _all
clear all

set more 1

capture log close

set scheme s1color

********************************************************************************
* I. DATA IMPORT AND LABEL 
********************************************************************************

* First import file containing both the rebate and the deflator for nd + sv
import excel ../input/rebates.xlsx, first sheet("rebate")
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
  UMCSENT MICH FEDFUNDS
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
rename MICH expected_infl_umich
rename FEDFUNDS ffr

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
label var pcdur "price deflator for durables consumption expenditures"
label var pcxnrgfd "price deflator consumption excluding food and energy"
label var pcfood "price deflator for consumption, food"
label var pcnrg "price deflator for energy goods and services consumption"
label var ncxnrg "nominal consumption less energy goods and services"
label var umcsent "U of Michigan Consumer Sentiment"
label var ffr "effective federal funds rate"
label var expected_infl_umich "U of Michigan Inflation Expectations"
label var nrebate "nominal rebate"

********************************************************************************
* IV. SAVE DATA
********************************************************************************

save ../output/narrative_data.dta, replace

********************************************************************************
* V. EXPORT REBATE AS FRACTION OF EXPENDITURE
********************************************************************************

* note: both variables annualized so correct ratio
use ../output/narrative_data.dta, clear
su ncons if mdate>=ym(2008,1) & mdate<=ym(2008,3)
gen rel_rebate = nrebate / r(mean)
keep mdate rel_rebate
keep if mdate>=ym(2008,1) & mdate<=ym(2008,12)
save ../output/rel_rebate.dta, replace

use ../output/narrative_data.dta, clear
su ncons if mdate>=ym(2001,5) & mdate<=ym(2002,6)
gen rel_rebate = nrebate / r(mean)
keep mdate rel_rebate
keep if mdate>=ym(2001,5) & mdate<=ym(2002,6)
save ../output/rel_rebate2001.dta, replace


