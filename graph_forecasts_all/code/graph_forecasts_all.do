
set scheme s1color
set fredkey b2d9da12107485ec086c5691d41d626c

use ../input/forecasts.dta
sum rcons_0704
local rcons_0704 = r(mean) //Real Consumption 2007Q4 monthly rate


import fred PCEC96,  clear 
gen idate = date(datestr,"YMD")
gen qdate = qofd(idate)

collapse (mean) pcec96 = PCEC96, by(qdate) //Monthly rate
format qdate %tq
tsset qdate
keep qdate pcec96
tempfile rpce
save `rpce'

use qdate RCONS* using ../input/levels_SPF0704.dta, clear
merge 1:1 qdate using ../input/greenbook08.dta, nogen 
merge 1:1 qdate using `rpce', nogen keep(1 3)


*Goldman Sachs forecasts
tsset qdate

gen rcons_gs = `rcons_0704' if qdate==q(2007q4) /* convert to monthly rates to match our other graph magnitudes */
replace rcons_gs = L.rcons_gs if qdate==q(2008q1)
replace rcons_gs = 0.995^.25*L.rcons_gs if qdate==q(2008q2)
replace rcons_gs = 0.995^.25*L.rcons_gs if qdate==q(2008q3)
replace rcons_gs = 1.01^.25*L.rcons_gs if qdate==q(2008q4)

*Normalize Real and SPF spending to 2007Q4 

foreach i in "pcec96" "RCONS_min" "RCONS_med"{

	gen _`i'_074 = `i' if qdate == yq(2007,4)
	egen `i'_074 = max(_`i'_074)
	replace `i' = `rcons_0704'*`i'/`i'_074

}

*Convert Greenbook to Dollars

foreach i in "20071205" "20080123"{

	gen rcons_gb_`i' = `rcons_0704' if qdate==q(2007q4) /* convert to monthly rates to match our other graph magnitudes */
	replace rcons_gb_`i' = ((1+gRPCE_`i'/100)^.25)*L.rcons_gb_`i' if qdate==q(2008q1)
	replace rcons_gb_`i' = ((1+gRPCE_`i'/100)^.25)*L.rcons_gb_`i' if qdate==q(2008q2)
	replace rcons_gb_`i' = ((1+gRPCE_`i'/100)^.25)*L.rcons_gb_`i' if qdate==q(2008q3)
	replace rcons_gb_`i' = ((1+gRPCE_`i'/100)^.25)*L.rcons_gb_`i' if qdate==q(2008q4)
}



/***************************************************************************
Graphs Forecasts
****************************************************************************/

keep if qdate < yq(2009,1) & qdate > yq(2007,3)

*Labels variables
lab var rcons_gs "Goldman-Sachs"
lab var rcons_gb_20071205 "Greenbook, 12/2007"
lab var rcons_gb_20080123 "Greenbook, 01/2008"
lab var pcec96 "Actual RPCE"
lab var RCONS_min "SPF Min."
lab var RCONS_med "SPF Med."


format qdate %tq





tw (scatter rcons_gb_2007 qdate, c(l l) clp(l) clw(medthick)  ///
    ytitle("", size(medsmall)) clc(forest_green) mc(forest_green))  ///
    (scatter rcons_gb_2008 qdate, c(l l) clp(l) clw(medthick)  clc(green) mc(green)) ///
    (scatter RCONS_min qdate, c(l l) clp(l) clw(medthick)  clc(ltblue) mc(ltblue)) ///
    (scatter RCONS_med qdate, c(l l) clp(l) clw(medthick)  clc(blue) mc(blue)) ///
    (scatter rcons_gs qdate, c(l l) clp(l) clw(medthick)  clc(purple) mc(purple)) ///
    (scatter pcec96 qdate, c(l l) clp(l) clw(medthick)  clc(black) mc(black)) ///
    , name(rcons_forecast, replace) scale(1.2) ///
    xlabel(, valuelabel) xtitle(Forecast) ytitle("billions of $, monthly rate")  yscale(titlegap(*10)) graphregion(margin(2 6 2 2)) ylabel(810(10)850)


graph export ../output/all_forecasts.eps, replace
cap graph export ../output/all_forecasts.png,  replace width(2400)




