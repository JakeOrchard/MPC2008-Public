**** NARRATIVE_2001_GRAPHS.DO

***
*** STUDIES THE BEHAVIOR OF KEY MACRO SERIES IN 2001

** Valerie Ramey, revised March 30, 2023

** REQUIRED FILES:

*      bea_jps_rebate_combined.dta

***************************************************************************************************

drop _all
clear all

set more 1

capture log close

set scheme s1color

use bea_jps_rebate_combined.dta

********************************************************************************
*   NARRATIVE GRAPHS
********************************************************************************

* A. REBATE GRAPH

     tw scatter nrebate mdate if mdate>=m(2001m1) & mdate<=m(2001m12), ///
       c(l ) clp(l ) ms(d ) clw(medthick ) ytitle("billions of $, monthly rates") ///
       clc(dknavy) mc(dknavy) xtitle("month") name(rebate) scale(1.2) ysize(4) xsize(7)
	 
	 label var nrebate "Rebate, Billions of $"
	 label var ndisp_income "Disposable Income"
	  tw (scatter nrebate mdate, c(l) clp(-) clw(medthick) ms(i) clc(blue) mc(blue) ///
	 yaxis(1) ylabel(0 5 10 15 20 25, axis(1))) ///
	 (scatter ndisp_income mdate, yaxis(2) c(l) clp(l) clw(medthick) ms(i) clc(green) mc(green) ///
	 ylabel(640 645 650 655 660 665, axis(2)) ytitle("Disposable Income, Billions of $, Monthly", axis(2) )) ///
	 if mdate>=m(2001m1) & mdate<=m(2001m12), legend(cols(1) position(11) ring(0)) ///
      xtitle("month") ysize(4) xsize(6) scale(1) name(nrebatey_combo)
      graph export ../output/fig_ndisprebate01.eps, replace

      
