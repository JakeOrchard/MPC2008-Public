cd "/Users/johanneswieland/Dropbox/MPC/Tasks/sunabraham/code"


use "../output/satest.dta", clear

tsset HHind Date 

xtreg Y D i.Date, fe

eventstudyweights D, controls(i.HHind i.Date) cohort(Cohort) rel_time(RelTime) 

mat list e(weights)
