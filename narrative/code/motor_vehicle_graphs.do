**** MOTOR_VEHICLE_GRAPHS.DO

***
*** CREATES GRAPHS OF MOTOR VEHICLE SALES IN 2008

** Valerie Ramey, revised April 24, 2022

** REQUIRED DATA:

    * Motor_vehicle_data.xlsx

***************************************************************************************************


drop _all
clear all

set more 1

capture log close

set scheme s1color


********************************************************************************
* I. DATA IMPORT AND TRANSFORMATION
********************************************************************************

  import excel ../input/Motor_vehicle_data.xlsx, first sheet("motor_vehicles")

  *CAUTION - check start date in Excel file to make sure the follow line is correct

  gen mdate = m(2001m5) + _n-1
  tsset mdate, m

  sort mdate

  drop if mdate<m(2007m1) | mdate>m(2009m12)

  * These are from the BEA Table 7.2.5S - SAAR, thousands of units
  * Motor Vehicle (MV) = auto + light truck (lgtrk)

  * Unit sales in 1,000s, SAAR
  * Average price is in $

  *convert to monthly rate

  replace newunitsales_auto_cons = newunitsales_auto_cons/12
  replace newunitsales_lgtrk_cons = newunitsales_lgtrk_cons/12
  gen newunitsales_mv_cons = newunitsales_auto_cons + newunitsales_lgtrk_cons
  gen newsales_mv_cons = (pavg_auto_cons*newunitsales_auto_cons + pavg_lgtrk_cons*    newunitsales_lgtrk_cons)/1000000
  gen rnewsales_mv_cons = 100*newsales_mv_cons/pcons_newmv

  gen newsales_auto_cons = pavg_auto_cons*newunitsales_auto_cons/12

  * put prices in thousands
  replace pavg_auto_cons = pavg_auto_cons/1000
  replace pavg_lgtrk_cons = pavg_lgtrk_cons/1000

  label var newunitsales_auto_cons "autos"
  label var newunitsales_lgtrk_cons "light trucks"

  label var pavg_auto_cons "autos"
  label var pavg_lgtrk_cons "light trucks"

********************************************************************************
* II. GRAPH DATA
********************************************************************************

   tw scatter newunitsales_mv_cons mdate, c(l l) clp(l l) clw(medthick medthick) ///
     clc(dkgreen) mc(dkgreen)  xline(`=ym(2008,5)', lp(-) lc(red)) ///
     ytitle("1,000s of units, monthly rate") title("Unit Sales (1,000s)") xtitle("month") name(newunitmv)
  
   tw scatter newsales_mv_cons mdate, c(l l) clp(l l) clw(medthick medthick) ///
     xline(`=ym(2008,5)', lp(-) lc(red)) title("Sales ($ billions)") clc(maroon) ///
	 mc(maroon) ytitle("billions of $, monthly rate") xtitle("month") name(newsalesmv)	

   * A. UNIT SALES VS. DOLLAR SALES
   
   * manuscript
   *  graph combine newunitmv newsalesmv , col(2) ysize(3) xsize(8) iscale(1.3) name(combomv)	

   * slides
   graph combine newunitmv newsalesmv , col(2) ysize(3) xsize(8) iscale(1.1) name(combomv)	
	 

graph export ../output/fig_mv_units_combo.eps, replace	

	 
   * B. 4 GRAPH COMBO - UNIT VS. DOLLAR, SEGMENT

   tw (scatter newunitsales_auto_cons mdate, c(l l) clp(l) clw(medthick) ///
       xline(`=ym(2008,5)', lp(-) lc(red)) clc(navy) mc(navy) ) ///
      (scatter newunitsales_lgtrk_cons mdate, c(l l) clp(l) clw(medthick) yaxis(2) ///
	  clc(orange) mc(orange)), title("Unit Sales by Segment (1,000s)") name(segment)

   tw (scatter pavg_auto_cons mdate, c(l l) clp(l) clw(medthick) ///
      xline(`=ym(2008,5)', lp(-) lc(red)) clc(navy) mc(navy) ) ///
      (scatter pavg_lgtrk_cons mdate, c(l l) clp(l) clw(medthick) yaxis(2) ///
	  clc(orange) mc(orange)), title("Average Price by Segment ($1,000s)") name(psegment)

   graph combine newunitmv newsalesmv segment psegment, col(2) ysize(4) xsize(6) ///
      iscale(.6) name(combo_by_segment)

   * C. MORE SEGMENT DETAIL - NOT CURRENTLY USED
   
   tw scatter subcompact compact intermediate fullsize luxury mdate, c(l l l l l) ///
     clp(l l l l l l) clw(medthick medthick medthick medthick medthick) ///
	 xline(`=ym(2008,5)', lp(-) lc(red)) ytitle("units, monthly rate")  name(autosegments) ///
     title("Sales of Autos by Segment") subtitle("May 2008 indicated by vertical line") 
   
   tw scatter small_pickup large_pickup small_van large_van cross_utility utility mdate, ///
      c(l l l l l l) clp(l l l l l l) clw(medthick medthick medthick medthick medthick medthick) ///
	  xline(`=ym(2008,5)', lp(-) lc(red)) ytitle("units, monthly rate")  name(lgtrksegments) ///
      title("Sales of Light Trucks by Segment") subtitle("May 2008 indicated by vertical line") 

