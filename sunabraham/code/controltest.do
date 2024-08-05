
// egen totalrebatedates = sum(rbtindicator) if insample, by(cuid)

// ivreghdfe d_totexp2 rbtindicator if insample & totalrebatedates<=1, absorb(intdate)  


// ivreghdfe d_totexp2 rbtindicator d_num_adults if insample & totalrebatedates<=1, absorb(intdate)  

// program drop eventstudyweights

cd "/Users/johanneswieland/Dropbox/MPC/Tasks/sunabraham/code"

use tempfile.dta, clear


gen mdate = mofd(dofc(INTDATE))
format mdate %tm

tsset CUID mdate, monthly

// egen totalrebatedates = sum(RBT_INDICATOR), by(CUID)

// gen pastrebate = (L3.RBT_INDICATOR==1) | (L6.RBT_INDICATOR==1) | (L9.RBT_INDICATOR==1)

local yvar "d_CARTKN"

local weight //"[weight=FINLWT21]"

local controls "AGE" //"AGE d_NUM_ADULTS d_PERSLT18"

local sample "if INSAMPLE " //  & totalrebatedates>0


ivreghdfe d_TOTEXP2 RBT_INDICATOR `controls'  `weight' `sample', cluster(CUID) absorb(mdate) 

areg d_TOTEXP2 RBT_INDICATOR `controls'   `weight' `sample', cluster(CUID) absorb(mdate) 



foreach var of varlist  CARTKN { //TOTEXP2

	reg d_`var' RBT_INDICATOR  `controls' i.mdate  `weight'  `sample', cluster(CUID)
	
// 	reg d_`var' RBT_INDICATOR L3.RBT_INDICATOR `controls' i.mdate  `weight'  `sample' , cluster(CUID)

//  	reg d_`var' RBT_INDICATOR L3.RBT_INDICATOR  L3.`var' `controls' i.mdate  `weight'  `sample', cluster(CUID)
	
	qui reg d_`var' i.mdate `controls' `weight' `sample'
	predict yvarres, r
	

	qui reg `var' i.mdate `controls' `weight' `sample'
	predict yvarreslvl, r
	gen yvarreslvllag = L3.yvarreslvl
	
// 	by RBT_INDICATOR, sort: su yvarres `sample' & mdate==ym(2008,7) `weight', detail
// 	by RBT_INDICATOR, sort: su yvarres `sample' & mdate==ym(2008,8) `weight', detail
// 	by RBT_INDICATOR, sort: su yvarres `sample' & mdate==ym(2008,9) `weight', detail
	by RBT_INDICATOR, sort: su yvarres `sample' & mdate==ym(2008,10) `weight', detail
	by RBT_INDICATOR, sort: su yvarreslvl `sample' & mdate==ym(2008,10) `weight', detail
	by RBT_INDICATOR, sort: su yvarreslvllag `sample' & mdate==ym(2008,10) `weight', detail
	
	
	
// 	su yvarres `sample' & mdate==ym(2008,10) `weight', detail
// 	sort CUID mdate
	
// 	qui reg `controls' i.mdate `weight' `sample'
// 	predict controlsrest, r
// // 	by RBT_INDICATOR, sort: su controlsrest `sample' & mdate==ym(2008,10) `weight', detail
// // 	by RBT_INDICATOR, sort: su controlsrest `sample' `weight', detail
// // 	i.mdate#i.RBT_INDICATOR

	
// 	reg yvarres i.mdate i.mdate#i.RBT_INDICATOR  controlsrest `sample' `weight'
	
// 	reg yvarres ind1 i.group   controlsrest `sample' `weight', nocons
	
// // 	reg yvarres i.mdate   `sample' `weight'
// // 	reg controlsrest i.mdate   `sample' `weight'
	
	qui drop yvarres* // controlsrest

}
// areg CARTKN RBT_INDICATOR `controls' i.mdate  `weight', cluster(CUID) absorb(CUID)
/*
*/


/*
areg RBT_INDICATOR `controls' if totalrebatedates<=1, absorb(mdate) 
predict testvar2, r
su testvar2, detail

hdfe d_TOTEXP2 RBT_INDICATOR, absorb(i.mdate `controls') gen(hdfe) tol(0.0000000001)
su hdfeRBT_INDICATOR, detail

reg yvar testvar
reg hdfed_TOTEXP2 hdfeRBT_INDICATOR

// eventstudyweights RBT_INDICATOR, cohort(COHORT) rel_time(RELATIVE_TIME) controls(i.mdate `controls')
