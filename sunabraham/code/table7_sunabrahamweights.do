clear all
set more off

global input "../input"
global output "../output"

/*******************************************************************************

Step 1: load data

*******************************************************************************/

* unzip file
*unzipfile "../input/psmjsample", replace
*For windows users
unzipfile "$input/psmjsample", replace

* import file
import delimited using "psmjsample.csv"

* delete file
erase "psmjsample.csv"

/*******************************************************************************

Step 2: Set Time and XS dimension

*******************************************************************************/

* convert interview data
qui gen intdate2 = mofd(date(intdate, "YMD"))
qui drop intdate
qui gen intdate = intdate2
qui drop intdate2
format intdate %tm

* tsset
xtset cuid intdate

*Ever received rebate
egen rbtreceived = max(rbtindicator),by(cuid)

*Labels variables
label variable rbtamt "ESP"
label variable rbtindicator "I(ESP)"
label variable eft "EFT"

gen rbtindicatoriv = rbtindicator
gen rbtindicatorcheck = rbtindicator*check
gen rbtindicatoreft = rbtindicator*eft

gen rbtamtcheck = rbtamt*check
gen rbtamteft = rbtamt*eft

label var rbtreceived "Household Received Rebate"

gen carn = cartkn > 0
gen caru = cartku >0

*Cohort date and relative time
gen _rbtdate = intdate if rbtindicator == 1
egen rbtdate = max(_rbtdate), by(cuid)

gen rm = intdate - rbtdate

gen g0 = rm == 0
gen g3 = rm == 3
gen gl3 = rm == -3

*Control Cohort

gen controlcohort = rbtreceived == 0


/******************************************************************************
Run baseline regression and compute weights
******************************************************************************/

**

reghdfe cartkn rbtamt age d_num_adults d_perslt18   [w=finlwt21] if insamplelvl == 1 & rbtreceived == 1, absorb(rbtdate intdate)

eventstudyweights g0 [w=finlwt21] if insamplelvl == 1 & rbtreceived == 1, controls(age d_num_adults d_perslt18 i.intdate i.rbtdate) cohort(rbtdate) rel_time(rm) saveweights("../output/weights")

mat list e(weights)

preserve
import excel "../output/weights.xlsx", clear firstrow

keep g0 rbtdate rm
reshape wide g0, i(rm) j(rbtdate)

graph twoway line g0581 g0582 g0583 g0584 rm, xtitle("Months since rebate") ytitle("Weight in TWFE rbtamt coefficient") graphregion(fcolor(white)) scheme(sj) legend(order(1 "June" 2 "July" 3 "August" 4 "September") subtitle("Rebate Interview"))

graph export weights_table7.png, replace

restore


/******************************************************************************
Run interaction weighted estimator from Sun and Abraham 2020
******************************************************************************/

*Cars
eventstudyinteract cartkn  g0 [w=finlwt21] if insamplelvl == 1 , cohort(rbtdate) control_cohort(controlcohort) /*
          */  covariates(age d_num_adults d_perslt18 ) absorb(i.intdate) vce(cluster intdate)



eventstudyinteract cartkn  g0 g3 [w=finlwt21] if insamplelvl == 1 , cohort(rbtdate) control_cohort(controlcohort) /*
          */  covariates(age d_num_adults d_perslt18 ) absorb(i.intdate) vce(cluster intdate)
	
*All spending
eventstudyinteract totexp2  g0 [w=finlwt21] if insamplelvl == 1 , cohort(rbtdate) control_cohort(controlcohort) /*
          */  covariates(age d_num_adults d_perslt18 ) absorb(i.intdate) vce(cluster intdate)



eventstudyinteract ndexp  g0 g3 gl3 [w=finlwt21] if insample == 1 , cohort(rbtdate) control_cohort(controlcohort) /*
          */  covariates(age d_num_adults d_perslt18 ) absorb(i.intdate i.cuid) vce(cluster intdate)
	
  
	




