clear all
set more off

/*******************************************************************************

Step 0: Load Data and define interview schedule variable

*******************************************************************************/


use ../output/psmjsampleinterview_wlabels.dta, clear


/****************************************************************************
Month between rebate receipt and interview
*****************************************************************************/
gen ms_rebate = intdate - firstrbtdate
*replace ms_rebate = 0 if rbtamt == 0


gen onemonth = ms_rebate == 1 & rbtindicator
gen twomonth = ms_rebate == 2 & rbtindicator
gen threemonth = ms_rebate ==3 & rbtindicator
gen lag_fourmonth = ms_rebate == 4 & lag1rbtindicator
gen lag_fivemonth = ms_rebate == 5 & lag1rbtindicator
gen lag_sixmonth = ms_rebate == 6 & lag1rbtindicator

gen ms_rebateshort = ms_rebate
replace ms_rebateshort = 0 if ms_rebate > 3 | ms_rebate < 0

label var onemonth "    \hspace{\parindent} Interview 1 month post rebate"
label var twomonth "    \hspace{\parindent} Interview 2 months post rebate"
label var threemonth "    \hspace{\parindent} Interview 3 months post rebate"

label var lag_fourmonth "   \hspace{\parindent} Interview 4 months post rebate"
label var lag_fivemonth "  \hspace{\parindent}  Interview 5 months post rebate"
label var lag_sixmonth "  \hspace{\parindent}  Interview 6 months post rebate"



xi gen i.incdinterview 

/*******************************************************************************

Step 1: Locals for regressions

*******************************************************************************/

local dependentvars "d_other_nipa d_mv d_pce" 

local varorder " onemonth twomonth threemonth  lag_fourmonth lag_fivemonth lag_sixmonth"

local controlvars "age d_num_adults d_perslt18"
local extravar "_Iincdinte* c.lag_pce#ms_rebateshort c.lag_mv#ms_rebateshort"
local absorbvars "intdate"
local hhfe "No"
local leveltitle "changes"
local frequency "interview"
local versionlist `""  "lag_interaction" "lc"' 


*Label Dependent Variables
label var d_pce "PCE"
label var lag_mv "Lag Motor Vehicle"
label var d_mv "MV and Parts"
label var d_other_nipa "Other Spending"
label var lag_pce "Lag Total Expenditure"
label var lag1rbtindicator "Lag Rebate Indicator"


*Income indicator 
xi gen i.incdinterview 


*Drop multiple rebate households
egen numrbts = sum(rbtindicator), by(cuid)
drop if numrbts > 1
		
/*******************************************************************************

Step 2: Table 2 Replication

*******************************************************************************/

	foreach depvar in `dependentvars' {
local variablename "`: variable label `depvar''"

foreach rebatesample in "" "& everrbtindicator" {
	
	local rebatetitle = "Full Sample"
	local rebatesamplesave "a"
	if "`rebatesample'"!="" {
		local rebatesamplesave = "b"
		local rebatetitle = "Rebate Only Sample"
	}
	
	foreach version in "`versionlist'"{		
	local groups1
	
	
	

if "`version'" == ""{
		display "ivreghdfe  `depvar'  onemonth twomonth threemonth lag1rbtindicator `controlvars'  if insample `rebatesample' [weight=finlwt21], cluster(cuid) absorb(`absorbvars')  "
		 ivreghdfe  `depvar'   onemonth twomonth threemonth lag1rbtindicator `controlvars'  if insample `rebatesample' [weight=finlwt21], cluster(cuid) absorb(`absorbvars')  
		 
qui estadd local decilefe `"No"' 
qui estadd local lagexp `"No"' 

test onemonth = twomonth 
local tmptest_onetwo: di %3.2f r(p)
qui estadd local test_onetwo "`tmptest_onetwo'" 
test onemonth = threemonth 
local tmptest_onethree: di %3.2f r(p)
qui estadd local test_onethree "`tmptest_onethree'" 
}

if "`version'" == "lag_interaction"{
		display "ivreghdfe  `depvar'  onemonth twomonth threemonth  lag_fourmonth lag_fivemonth lag_sixmonth`controlvars'  if insample `rebatesample' [weight=finlwt21], cluster(cuid) absorb(`absorbvars')  "
		 ivreghdfe  `depvar' onemonth twomonth threemonth  lag_fourmonth lag_fivemonth lag_sixmonth `controlvars'  if insample `rebatesample' [weight=finlwt21], cluster(cuid) absorb(`absorbvars')  

qui estadd local decilefe `"No"' 
qui estadd local lagexp `"No"' 

test onemonth = twomonth 
local tmptest_onetwo: di %3.2f r(p)
qui estadd local test_onetwo "`tmptest_onetwo'" 
test onemonth = threemonth 
local tmptest_onethree: di %3.2f r(p)
qui estadd local test_onethree "`tmptest_onethree'" 


}

if "`version'" == "lc"{
		display "ivreghdfe  `depvar' onemonth twomonth threemonth  lag_fourmonth lag_fivemonth lag_sixmonth `controlvars' `extravar'  if insample `rebatesample' [weight=finlwt21], cluster(cuid) absorb(`absorbvars')  "
		 ivreghdfe  `depvar' onemonth twomonth threemonth  lag_fourmonth lag_fivemonth lag_sixmonth `controlvars' `extravar'  if insample `rebatesample' [weight=finlwt21], cluster(cuid) absorb(`absorbvars')  

qui estadd local decilefe `"Yes"' 
qui estadd local lagexp `"Yes"' 

test onemonth = twomonth 
local tmptest_onetwo: di %3.2f r(p)
qui estadd local test_onetwo "`tmptest_onetwo'" 
test onemonth = threemonth 
local tmptest_onethree: di %3.2f r(p)
qui estadd local test_onethree "`tmptest_onethree'" 

}

		 
		 /*if "`depvar'" == "d_totexp2"{
		 	stop
		 }*/
	
		qui eststo
		qui estadd local timefe `"Yes"'
		qui estadd local hhfe `"`hhfe'"'
		qui estadd local hhcontrols `"Yes"'
		
		 	 //`: variable label `depvar''
		
	}	
	
}
	* this code produces the regression table
	local nregs = 6
	local order `"`varorder'"' 
	local indicate
	local midrules1 `" \cmidrule(l{.75em}){2-4} \cmidrule(l{.75em}){5-7}"'
	local groups2 `" &   \multicolumn{3}{c}{ Full Sample} &  \multicolumn{3}{c}{Rebate Only} \\ "'  	 
	local groups `" "`groups2'" "`midrules1'"  "'
	local stats " test_onetwo test_onethree  timefe hhcontrols decilefe lagexp  N"
	local stats_fmt "%12.2f %12.2f %3s %3s %12.2f %12.0fc"
	local stats_label `" `"One month = Two month p-value"' `"One month = Three month p-value"'   `"Time Effects"' `"Household Controls"' "Income Decile FE" "Lag Expenditure Interacted Rebate Month" `"Observations"' "'
	local num_stats: word count `stats' 
	local layout
	forvalues l = 1/`num_stats' {
		local layout `"`layout' "\multicolumn{1}{c}{@}" "'
	}
	local keepvars `"`varorder'"' 
	local dropvars 
	local title `"`variablename' Response to Rebate:"' 
	local table_preamble `" "\begin{table}[!t] \centering \sisetup{table-format=3.2} \def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}" "\caption{`title'}" "\begin{tabularx}{\hsize}{@{\hskip\tabcolsep\extracolsep\fill}l*{`nregs'}{S}}" "\\" "\hline" "'
	local prehead `"prehead(`table_preamble' `groups')"'			
	local posthead `"posthead(`"\hline"' `"\multicolumn{`=`nregs'+1'}{l}{Right hand side variables:}\\"' `"\\"')"'
	local notes `"\$ \:^{*}\:p<0.1,\:\:^{**}\:p<0.05,\:\:^{***}\:p<0.01 \$ "'
	local filename `"table_randomtest_`depvar'"'
	local prefoot(" ")
	local postfoot `"postfoot(`"\hline\hline \end{tabularx} \begin{minipage}{\hsize} \rule{0pt}{9pt} \footnotesize `notes'  \end{minipage} \label{tab:`filename'} \end{table}"')"'
	
	
esttab * using "../output/`filename'.tex", rename(1.ms_rebatec#c.rbtindicator ms_rebate 2.ms_rebatec#c.rbtindicator ms_rebate) replace cells(b(star fmt(a2)) se(par fmt(a2) abs)) starlevels(\$^{*}$ 0.1 \$^{**}$ 0.05 \$^{***}$ 0.01) drop(`dropvars', relax) 	keep(`keepvars') indicate(`indicate') `prehead' `posthead' `postfoot' order(`order') label varlabel(`varlabels') stats(`stats', layout(`layout') fmt(`stats_fmt') 		labels(`stats_label')) collabels(,none) numbers nomtitles substitute(# `" X "' tabular* tabularx `"{1}{c}{("' `"{1}{L}{("') width(\hsize)
		
	eststo clear
		}	
	
	
	
