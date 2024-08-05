**** MOTOR_VEHICLE_ANALYSIS.DO

***
*** STUDIES MOTOR VEHICLES SALES IN 2008

** Valerie Ramey, revised November 19, 2021

***************************************************************************************************


drop _all
clear all

set more 1

capture log close

set scheme s1color

********************************************************************************
* I. DATA IMPORT 
********************************************************************************

import excel "../input/retail_pce_dat.xlsx", first sheet("motor_vehicles")

*CAUTION - check start date to make sure the follow line is correct

gen mdate = m(2001m5) + _n-1
tsset mdate, m

sort mdate

drop if mdate<m(2007m1) | mdate>m(2009m12)

* These are from the BEA Table 7.2.5S - SSAR, thousands of units
* Motor Vehicle (MV) = auto + light truck (lgtrk)

*convert to monthly rate

replace newunitsales_auto_cons = newunitsales_auto_cons/12
replace newunitsales_lgtrk_cons = newunitsales_lgtrk_cons/12
gen newunitsales_mv_cons = newunitsales_auto_cons + newunitsales_lgtrk_cons
gen newsales_mv_cons = (pavg_auto_cons*newunitsales_auto_cons + pavg_lgtrk_cons*newunitsales_lgtrk_cons)/1000
gen rnewsales_mv_cons = 100*newsales_mv_cons/pcons_newmv

tw scatter newunitsales_mv_cons month, c(l l) clp(l l) clw(medthick medthick) ///
  xline(2008.3333) ytitle("1,000s of units, monthly rate") title("Unit Sales of New MV to Consumers") ///
  subtitle("May 2008 indicated by vertical line")   name(newunitmv)
  
tw scatter newsales_mv_cons month, c(l l) clp(l l) clw(medthick medthick) xline(2008.3333) ///
   title("$ Sales of New MV to Consumers (based on avg. price)") ///
   subtitle("May 2008 indicated by vertical line") ytitle("millions of $, monthly rate") name(newsalesmv)	
   
tw scatter rnewsales_mv_cons month, c(l l) clp(l l) clw(medthick medthick) xline(2008.3333) ///
   ytitle("millions of 2012 $, monthly rate")  name(rnewsalesmv) ///
   title("Real $ Sales of New MV to Consumers") subtitle("May 2008 indicated by vertical line")  	

tw (scatter newunitsales_auto_cons month, c(l l) clp(l) clw(medthick) xline(2008.3333) ) ///
    (scatter newunitsales_lgtrk_cons month, c(l l) clp(l) clw(medthick) yaxis(2)), ///
	 title("Unit Sales of New MV to Consumers by Segment")  name(segment)

tw (scatter pavg_auto_cons month, c(l l) clp(l) clw(medthick) xline(2008.3333) ) ///
    (scatter pavg_lgtrk_cons month, c(l l) clp(l) clw(medthick) yaxis(2)), ///
	title("Average Price of a MV by Segment") subtitle("May 2008 indicated by vertical line")  name(psegment)
	
tw (scatter newunitsales_auto_cons month, c(l l) clp(l) clw(medthick) xline(2008.3333) ) ///
    (scatter newunitsales_auto_bus month, c(l l) clp(l) clw(medthick) yaxis(2)), ///
	 title("Unit Sales of New MV to Consumers vs. Business")  name(sector)
	 
tw (scatter pavg_auto_cons month, c(l l) clp(l) clw(medthick) xline(2008.3333) ) ///
    (scatter pavg_auto_bus month, c(l l) clp(l) clw(medthick) yaxis(2)), ///
	title("Average Price of a MV Consumer vs. Business") subtitle("May 2008 indicated by vertical line")  name(psector)