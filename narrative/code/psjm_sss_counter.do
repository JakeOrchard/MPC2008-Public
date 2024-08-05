**** PSJM_SSS_COUNTER.DO

***
*** UPDATES SAHM-SHAPIRO-SLEMROD MOTOR VEHICLE COUNTERFACTUAL USING PARKER ET AL. MPC

** REQUIRES: Update_of_Sahm-Shapiro-Slemrod-Motorvehicles.xlsx

***************************************************************************************************

drop _all
clear all

set more 1

capture log close

set scheme s1color

* IMPORT DATA - DATA ARE ANNUALIZED, BILLIONS OF $
import excel "../input/Update_of_Sahm-Shapiro-Slemrod-Motorvehicles.xlsx", sheet("main") first

gen mdate = m(2007m1) + _n-1

tsset mdate, m

replace rebate = 0 if rebate==.

* Convert to a monthly basis

gen nexpmv = npce_new_mv_saar/12
gen rebate = rebate_ar/12
gen induced = induced_saar/12
gen induced_nsa = induced_nsaar/12

replace induced = 0 if induced==. & month>=2007.99 & month<2008.25

gen nexpmvcf = nexpmv - induced

label var nexpmv "actual"
label var nexpmvcf "counterfactual"





* Alternative graphs with spending more distributed

gen induced_alt = 0.357*(rebate + L.rebate + L2.rebate)/3
gen induced_alt2 = 0.357*(rebate + L.rebate + L2.rebate + L3.rebate)/4
gen nexpmvcf_alt = nexpmv - induced_alt
gen nexpmvcf_alt2 = nexpmv - induced_alt2

replace induced_alt = 0 if induced_alt==. & month>=2007.99 & month<2008.25
replace induced_alt2 = 0 if induced_alt2==. & month>=2007.99 & month<2008.25

label var nexpmvcf_alt "3-month smooth"
label var nexpmvcf_alt2 "4-month smooth"


* Compare induced_nsa to induced - the differences are very small

// tw scatter induced induced_nsa mdate if mdate>=m(2008m1) & mdate<=m(2008m12), c(l l) clp(l -) clw(medthick medthick) ///
//   ms(i i)  ytitle("billions of $, monthly rate") xtitle("month") ///
//   ysize(4) xsize(7) scale(1.5) name(compare_sa_nsa) ///


* Main graph

tw scatter nexpmv nexpmvcf mdate if mdate>=m(2007m9) & mdate<=m(2009m6), ///
  c(l l) clp(l -) clw(medthick medthick) clc(black purple) mc(black purple) ///
  ms(i i)  ytitle("billions of $, monthly rate") xtitle("month") ///
  ysize(4) xsize(7) scale(1.5) name(mv_sss_counter) ///

graph export ../output/fig_sss_mv_counter.eps, replace

tw scatter nexpmv nexpmvcf nexpmvcf_alt nexpmvcf_alt2 mdate if mdate>=m(2007m9) & mdate<=m(2009m6), ///
   c(l l l l) clp(l - - - -) clw(medthick medthick medthick medthick) clc(black) ///
  ms(i i i i)  ytitle("billions of $, monthly rate") xtitle("month") ///
  ysize(4) xsize(7) scale(1.5)  ///

 
graph export ../output/fig_sss_mv_counter_alt.eps, replace
