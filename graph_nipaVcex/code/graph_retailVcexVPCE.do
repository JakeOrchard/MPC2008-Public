clear all
set more off




/******************************************************************************
*1. Pull in clean retail data and merge with PCE and CEX data
******************************************************************************/


use ../input/retail_sales_clean

*Fix timing
gen date = mofd(dofc(__date______))
format date %tm


*Rename variables
rename __sales____Retail_sales__total nrts_tot 


merge 1:1 date using ../output/pce_cex_series.dta, nogenerate


*Goods
// gen rPCEG = 100*PCEG_month/pcons_good
// gen rpceg = 100*nipa_total_good_pop_sa/pcons_good
gen rnrts_tot = nrts_tot / 1000 / pcons_good * 100

label var rnrts_tot "Real Retail Sales"
// label var PCEG "PCE for goods"

set scheme s1color

tw (scatter rnrts_tot date if mdate>m(2007m1) & mdate<=m(2009m12), c(l l) clp(l) clw(medthick) xline(`=ym(2008,5)', lp(-) lc(red)) ///
    ytitle("billions of $, monthly rate", size(medsmall)) clc(purple) mc(purple))  ///
    (scatter rPCEG date if mdate>m(2007m1) & mdate<=m(2009m12), c(l l) clp(l) clw(medthick) yaxis(2) clc(blue) mc(blue))


graph export ../output/fig_real_rts_pcegoods.pdf, as(pdf) replace
