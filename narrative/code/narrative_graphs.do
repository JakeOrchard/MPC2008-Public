**** NARRATIVE_GRAPHS.DO

***
*** STUDIES THE BEHAVIOR OF KEY MACRO SERIES IN 2008

** Valerie Ramey, revised April 24, 2022
** Edited December, 13, 2022 by Jacob Orchard

** REQUIRED FILES:

*      narrative_data.dta, created by build_narrative_dat.do

**  STATA ADD-IN COMMANDS - must be installed before running this program

      * To install the command, first type the phrase below in the Stata command window
	  *    then run the program.
	  
      *  net install grc1leg,from( http://www.stata.com/users/vwiggins/)
***************************************************************************************************

drop _all
clear all

set more 1

capture log close

set scheme s1color

local yearlist `"2008" "2001"'

********************************************************************************
* I.  READ IN DATA, CREATE REAL AND MONTHLY BASIS VARIABLES
********************************************************************************
foreach year in "`yearlist'"{
	
	if "`year'" == "2008"{
		local basemonth = ym(2008,1)
		local rebatebegin = ym(2008,5)
		local graphrange "mdate>=m(2007m7) & mdate<=m(2009m3)"
		local graphrangeincome  "mdate>=m(2007m7) & mdate<=m(2009m6)"
		local graphrangerebate = "mdate>=m(2008m1) & mdate<=m(2008m12)"
		local yearlabel ""
		local axissize = 82
	}
	if "`year'" == "2001"{
		local basemonth = ym(2001,1)
		local rebatebegin = ym(2001,7)
		local graphrange "mdate>=m(2001m1) & mdate<=m(2001m12)"
		local graphrangeincome  "mdate>=m(2001m1) & mdate<=m(2001m12)"
		local graphrangerebate = "mdate>=m(2001m1) & mdate<=m(2001m12)"
		local yearlabel "01"
		local axissize = 60
	}
	
	use ../output/narrative_data.dta

	foreach var in cons rebate disp_income {
	  replace n`var' = n`var'/12 /*put on monthly actual basis*/
	    sum pcons if mdate == `basemonth'
	  gen r`var' = n`var'*r(mean)/pcons  /*real, using Jan 2008 base of 93.307*/
	}
	
	foreach var in cnd csv cdur cndsv cnrg {
	  replace n`var' = n`var'/12 /*put on monthly actual basis*/
	  gen r`var' = n`var'*100/p`var'  /*real*/
	}
	replace ncxnrg = ncxnrg/12

	gen infl = 1200*ln(pcons/L.pcons)

	label var infl "inflation, pce, annualized"

	label var ncons "nominal"
	label var rcons "real"
	label var ndisp_income "disposable income"
	label var rdisp_income "disposable income"
	

	********************************************************************************
	* II.  NARRATIVE GRAPHS
	********************************************************************************

	* A. REBATE GRAPH

	      tw scatter nrebate mdate if `graphrangerebate', ///
		   xline(`rebatebegin', lp(-) lc(red)) ///
		c(l ) clp(l ) ms(d ) clw(medthick ) ytitle("billions of $, monthly rates") ///
		clc(dknavy) mc(dknavy) xtitle("month") name(rebate) scale(1.2) ysize(4) xsize(7)
		   
		   graph export ../output/fig_rebates`yearlabel'.eps, replace 
	
	if "`year'" == "2001"{
		   tw (scatter nrebate mdate, c(l ) clp(l ) ms(d ) clw(medthick ) clc(dknavy) mc(dknavy)  yaxis(1) ylabel(0 5 10 15 20 25, axis(1))) (scatter ndisp_income mdate, c(l ) clp(l ) ms(d ) clw(medthick ) clc(green) mc(green) yaxis(2) ylabel(640 645 650 655 660 665, axis(2))) if `graphrangerebate', ///
		ytitle("rebate, billions of $", axis(1)) ytitle("disposable income, billions of $", axis(1)) ///
		 xtitle("month") name(rebate_ndisp) scale(1.2) ysize(4) xsize(7)
		   
		   graph export ../output/fig_rebate_ndisp`yearlabel'.eps, replace 
	}

	* B. DISPOSABLE INCOME AND CONSUMPTION GRAPHS - COMBINED REAL AND NOMINAL

	   * REAL AND NOMINAL CONSUMPTION GRAPHS

	     label var ncons "nominal"
	     label var rcons "real"
	     label var ndisp_income "nominal"
	     label var rdisp_income "real"

	  * Graph 1: Disposable income (real & nominal); Graph 2: Consumption (real & nominal)
		sum ndisp_income  if `graphrangeincome'
		local min_ndisp = 868
		local max_ndisp = `min_ndisp' + `axissize'
		local min_ndispscale = 870
		
	    tw scatter rdisp_income ndisp_income mdate if `graphrangeincome', title("Disposable Income") ///
	      xline(`rebatebegin', lp(-) lc(red) lw(thin))  ytitle("billions of $, monthly rates") ///
	      c(l l ) clp(l -) clw(medthick medium) clc(blue green) mc(blue green) ms(d i) ///
	      xtitle("month")  yscale(range(`min_ndisp' `max_ndisp')) ylabel(`min_ndispscale'(20)`max_ndisp') legend(size(small)) name(disp_y)
	      
	      if "`year'" == "2008"{
		local min_ncons = 808
		local max_ncons = `min_ncons' + `axissize'
		local min_nconsscale = 810
	      }
	      if "`year'" == "2001"{
		local min_ncons = 560
		local max_ncons = `min_ncons' + `axissize'
		local min_nconsscale = 560
		
	      }
	      
	   
	   tw scatter rcons ncons mdate if `graphrangeincome', title("Consumption Expenditure") ///
	      xline(`rebatebegin', lp(-) lc(red) lw(thin))  ytitle("billions of $, monthly rates") ///
	      c(l l ) clp(l -) clw(medthick medium) clc(blue green) mc(blue green) ms(d i) ///
	      xtitle("month")  yscale(range(`min_ncons' `max_ncons')) ylabel(`min_nconsscale'(20)`max_ncons') legend(size(small)) name(cons)  

		 * Create combined graph with a common legend - must pre-install grc1leg
		 
	//     grc1leg disp_y cons, cols(2) ysize(3) xsize(8) iscale(1.1) legendfrom(cons)
		graph combine disp_y cons, cols(2) ysize(3.5) xsize(8) iscale(1.1) 

		graph export ../output/fig_cy_combo`yearlabel'.eps, replace

		
		*Real Disposable Income and Consumption Expenditures
		label var rdisp_income "Real Disposable Income"
		label var rcons "Real Consumption Expenditure"
		*2008
		if "`year'" == "2008"{
		tw (scatter rcons mdate if `graphrangeincome', c(l) clp(-) clw(thick) ms(i) clc(blue) mc(blue) ///
			yaxis(1) ylabel(815 825 835 845 855 865 875, axis(1)) ) ///
			(scatter rdisp_income mdate if `graphrangeincome', yaxis(2) c(l) clp(l) clw(medthick) ms(i) clc(green) mc(green) ///
			ylabel(875 885 895 905 915 925 935, axis(2))) ///
			if mdate>=m(2007m7) & mdate<=m(2009m6), ///
			xline(`rebatebegin', lp(-) lc(red)) ///
			legend(cols(1) ring(0) position(2)) ///
			xtitle("Month") ysize(4) xsize(6) scale(1) name(rycons_combo)  
		}
		 if "`year'" == "2001"{
		 	tw (scatter rcons mdate if `graphrangeincome', c(l) clp(-) clw(thick) ms(i) clc(blue) mc(blue) ///
			yaxis(1) ylabel(570 580 590 600 610, axis(1)) ) ///
			(scatter rdisp_income mdate if `graphrangeincome', yaxis(2) c(l) clp(l) clw(medthick) ms(i) clc(green) mc(green) ///
			ylabel(630 640 650 660 670, axis(2))) ///
			if `graphrangeincome', ///
			xline(`rebatebegin', lp(-) lc(red)) ///
			legend(cols(1) ring(0) position(2)) ///
			xtitle("Month") ysize(4) xsize(5) scale(1) name(rycons_combo)  
			
			
		 }
			graph export ../output/fig_cy_pres`yearlabel'.eps, replace
			
		
		* CONSUMPTION BY CATEGORY GRAPHS

	   label var rcnd "nondurables"
	   label var rcdur "durables"
	   label var rcsv "services"

	   tw scatter rcnd mdate if `graphrange', ///
	     xline(`rebatebegin', lp(-) lc(red)) c(l l l) clp(l - _) clw(medthick ) ///
	     clc(navy) mc(navy) ytitle("billions of $, monthly rate") xtitle("month") ///
	     legend(size(small)) title("Nondurables") name(rcndgraph)
	   
	   tw scatter rcdur mdate if `graphrange', ///
	     xline(`rebatebegin', lp(-) lc(red)) c(l l l) clp(l - _) clw(medthick ) ///
		 clc(dkgreen) mc(dkgreen) ytitle("billions of $, monthly rate") xtitle("month") ///
		 title("Durables") legend(size(small)) name(rcdurgraph)
	   
	   tw scatter rcsv mdate if `graphrange', ///
	      xline(`rebatebegin', lp(-) lc(red)) c(l l l) clp(l - _) clw(medthick ) ///
		  clc(maroon) mc(maroon) ytitle("billions of $, monthly rate") xtitle("month")  ///
	      title("Services") legend(size(small)) name(rcsvgraph)

	   graph combine rcndgraph rcdurgraph rcsvgraph, col(2) ysize(4) xsize(7) iscale(0.75) ///
	      name(PCE_by_type)

	  graph export ../output/fig_pce_type`yearlabel'.eps, replace
	  
	  if "`year'" == "2001"{
		 	tw (scatter rcndur_jpscat mdate if `graphrangeincome', c(l) clp(-) clw(thick) ms(i) clc(blue) mc(blue) ///
			yaxis(1) ylabel(300 310 320 330 340, axis(1)) ) ///
			(scatter rdisp_income mdate if `graphrangeincome', yaxis(2) c(l) clp(l) clw(medthick) ms(i) clc(green) mc(green) ///
			ylabel(630 640 650 660 670, axis(2))) ///
			if `graphrangeincome', ///
			xline(`rebatebegin', lp(-) lc(red)) ///
			legend(cols(1) ring(0) position(2)) ///
			xtitle("Month") ysize(4) xsize(5) scale(1) name(rycjps_combo)  
			
			graph export ../output/fig_jpscy_pres`yearlabel'.eps, replace
			
			
			
		 }
			

	* D. CONSUMER PRICE INDEX GRAPHS

	  foreach var in cons cnd csv cdur cnrg cxnrgfd cfood {
		su p`var' if mdate==`basemonth'
		gen p`var'norm = p`var' / r(mean) * 100
	  
		 gen lp`var' = ln(p`var')
		
	      tw scatter lp`var' mdate if `graphrangerebate', ///
		c(l ) clp(l) clw(medthick ) ytitle("log price index") xtitle("month") /// 
			ysize(4) xsize(7) scale(1.1) name(lp`var') ///
		title("`year' log price `var'", size(medium)) 
	  }

	   label var lpcons "total"
	   label var lpcxnrgfd "excl. food, energy"
	   label var lpcnrg "energy goods & services"
	   
	   label var pconsnorm "total"
	   label var pcxnrgfdnorm "excl. food, energy"
	   label var pcnrgnorm "energy goods & services"


	 tw scatter pconsnorm pcxnrgfdnorm mdate if `graphrange', ///
	     xline(`rebatebegin', lp(-) lc(red) lw(thin)) ytitle("price index, jan `year' = 100") ///
	     c(l l ) clp(l l) clw(medthick medthick) clc(navy dkgreen) mc(navy dkgreen) ms(d i) ///
	     xtitle("month") ysize(4) xsize(6) scale(1.2) legend(size(small)) name(figpnorm)


	   tw scatter pcnrgnorm mdate if `graphrange', ///
	     xline(`rebatebegin', lp(-) lc(red) lw(thin)) ytitle("price index, jan `year' = 100") ///
	     c(l l ) clp(l -) clw(medthick medium) clc(sienna) mc(sienna) ms(d i) ///
	     xtitle("month") ysize(4) xsize(6) scale(1.2) legend(on) legend(size(small)) name(figpnrgnorm)
	   
	   graph combine figpnorm figpnrgnorm, col(2) ysize(3) xsize(8) iscale(1.1) name(price_deflator_combo)

	  graph export ../output/fig_pconsnorm_combo`yearlabel'.eps, replace

	  *E. FEDERAL FUNDS RATES
	  
	   	* create ex-post real rate
	   	gen rffr = ffr - 1200*D.lpcons      /* uses PCE deflator */
	   	gen rffr2 = ffr - 1200*D.lpcxnrgfd  /* uses PCE deflator excluding energy */

		* create ex-ante real rate
		gen rffr_umich = ffr - expected_infl_umich

	    tw scatter ffr mdate if `graphrange', ///
	      xline(`rebatebegin', lp(-) lc(red) lw(thin)) ytitle("percent") ///
	      c(l l ) clp(l -) clw(medthick medium) clc(dkgreen) mc(dkgreen) ms(d i) ///
	      xtitle("month") ysize(4) xsize(6) scale(1.2) name(ffr) title("Nominal")
	 
	    tw scatter rffr_umich mdate if `graphrange', ///
	      xline(`rebatebegin', lp(-) lc(red) lw(thin)) ytitle("percent") ///
	      c(l l ) clp(l -) clw(medthick medium) clc(sienna) mc(sienna) ms(d i) ///
	      xtitle("month") ysize(4) xsize(6) scale(1.2) name(rffr) title("Real (ex ante)")
	 
	  graph combine ffr rffr, col(2) ysize(3) xsize(8) iscale(1.1) name(fedfunds)

	graph export ../output/fig_ffrante_combo`yearlabel'.eps, replace
	clear all
}

