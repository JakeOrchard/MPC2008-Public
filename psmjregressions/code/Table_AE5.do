clear all
set more off

cap cd ../psmjregressions/code




/*******************************************************************************

Step 0: Load Data

*******************************************************************************/


use ../output/psmjsampleinterview_wlabels.dta, clear

drop nipa* psmj*





/*******************************************************************************

Step 1: Locals for regressions

*******************************************************************************/


local dependentvars "rbtamt d_pce" //rbtamt must go first
local varorder "rbtindicator lag_pce lag_mv"
local samplelist `"" "& everrbtindicator"'


local absorbvars "intdate"

local controlvars "age d_num_adults d_perslt18"

local versionlist `""' 
local extracontrollist `"" "_Iincdinte* lag_pce lag_mv"'

*Never treated and currently treated sample
gen rel_interview = intdate-firstrbtintdate
gen nevertreatsamp = everrbtind == 0 | rel_interview == 0


keep if insample 
egen numrbts = sum(rbtindicator), by(cuid)
keep if numrbts <2

xtset cuid intdate

*Label Dependent Variables
label var d_pce "PCE"
label var lag_mv "Lag Motor Vehicle"
label var d_mv "Motor Vehicle and Parts Spending"
label var d_other_nipa "Other Spending"
label var lag_pce "Lag Total Expenditure"
label var lag1rbtindicator "Lag Rebate Indicator"




/*******************************************************************************

Step 2: Table 2 Replication with BJS

*******************************************************************************/

xi gen i.incdinterview



foreach depvar in `dependentvars' {



local variablename "`: variable label `depvar''"
local regnumber = 1


	foreach rebatesample in "`samplelist'" {
		foreach extra in "`extracontrollist'"{
		
						
*DID Imputation Estimator 				
			
			display "did_imputation `depvar' cuid intdate firstrbtintdate [w=finlwt21] if insample `rebatesample' , cluster(cuid) fe(`absorbvars' `hhfe') controls(`controlvars' `extra') horizons(0) maxit(10000)"
			did_imputation `depvar' cuid intdate firstrbtintdate [w=finlwt21] if insample `rebatesample' , cluster(cuid) fe(`absorbvars' `hhfe') controls(`controlvars' `extra') horizons(0) autosample maxit(10000)
			local coef = _b[tau0]
			

		
		
		
		
		
		qui eststo
		
		*Extra controls indicator
		if "`extra'" == ""{
		qui estadd local decilefe `"No"'
	
		}
		if "`extra'" != ""{
			qui estadd local decilefe `"Yes"'
		}
		 
		local regnumber = `regnumber' + 1
		
		
		*Create implied MPC
		if "`depvar'" == "rbtamt"{
			local rebate`regnumber' = `coef'
		
		}
		else{
		
		local temp_mpc: di %3.2f `coef'/`rebate`regnumber''
		qui estadd local implied_mpc "`temp_mpc'"
		
	}
	
		}		
}

	* this code produces the regression table
	local nregs = 4
	
	
	local order `"`varorder'"' 
	local indicate
	local midrules1 `"\cmidrule(l{.75em}){2-`=`nregs'+1'}"'
	
	local groups1 `" &  \multicolumn{2}{c}{Full Sample} &  \multicolumn{2}{c}{Rebate Only Sample} "' 
	
	

	local groups2 `" &  \multicolumn{4}{c}{OLS}  &  \multicolumn{4}{c}{DID Imputation} \\"'  	 
	local groups `" "`groups1'" \\ "`midrules1'"  "'
	local stats " implied_mpc decilefe N"
	local stats_label `" `"Implied MPC"' `"Income Decile FE"' `"Observations"' "'
	local stats_fmt " %3.2f %3s  %12.0fc"
	
	
	if "`depvar'" == "rbtamt"{
	
		local stats "extravars N"
		local stats_fmt " %3s  %12.0fc"
		local stats_label `"`"Extra controls"' `"Observations"' "'
	}
	local num_stats: word count `stats' 
	local layout
	forvalues l = 1/`num_stats' {
		local layout `"`layout' "\multicolumn{1}{c}{@}" "'
	}
	local keepvars `"`varorder'"' 
	local dropvars 
	local filename `"Table_AE5"'
	local tablename `"Table_AE5"'
	
	local title `"Contemporaneous Household `variablename'  Response to Rebate: BJS Method"'
	local notes `"Notes: The dependent variable is the change in `variablename'.   Regressions include interview (time) fixed effects, and household level controls for age, change in number of adults, and change in number of children.   Standard errors, in parentheses, are clustered at the household level. Significance is indicated by: \$ \:^{*}\:p<0.1,\:\:^{**}\:p<0.05,\:\:^{***}\:p<0.01 \$. "'
	

	
	if "`depvar'" == "rbtamt"{
		local title `" First Stage: Rebate Amount Conditional on Rebate Receipt"'
	
		local notes `"Notes: The dependent variable is the dollar value of Econoimic Stimulus Payments (ESP) received by the household. Standard errors, in parentheses, are clustered at the household level. Significance is indicated by: \$ \:^{*}\:p<0.1,\:\:^{**}\:p<0.05,\:\:^{***}\:p<0.01 \$. All regressions include interview (time) fixed effects, as well as household level controls for age, change in number of adults, and change in number of children. Extra controls refer to additional controls for household income decile and lagged total spending. Rebate sample includes only households that receieve a rebate at some point during our sample period.  "'
	
	
	}
	
	
	
	local table_preamble `" "\begin{table}[!t] \centering \sisetup{table-format=3.2} \def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}" "\caption{`title'}" "\begin{tabularx}{\hsize}{@{\hskip\tabcolsep\extracolsep\fill}l*{`nregs'}{S}}"   "\hline" "'
	
		
		local postfoot `"postfoot(`"\hline\hline \end{tabularx} \begin{minipage}{\hsize} \rule{0pt}{9pt} \footnotesize `notes'  \end{minipage} \label{tab:`tablename'} \end{table}"')"'
	
	local prehead `"prehead(`table_preamble' `groups')"'			
	local posthead `"posthead(`"\hline"' `"\\"')"'
	
	local prefoot(" ")
	

	
	esttab * using "../output/`filename'.tex", rename(tau0 rbtindicator tau3 lag1rbtindicator __00000V rbtindicator) replace cells(b(star fmt(a2)) se(par fmt(a2) abs)) starlevels(\$^{*}$ 0.1 \$^{**}$ 0.05 \$^{***}$ 0.01) drop(`dropvars', relax) ///
		keep(`keepvars') indicate(`indicate') `prehead' `posthead' `postfoot' order(`order') label varlabel(`varlabels') stats(`stats', layout(`layout') fmt(`stats_fmt') ////
		labels(`stats_label')) collabels(,none) numbers nomtitles substitute(# `" X "' tabular* tabularx `"{1}{c}{("' `"{1}{L}{("') width(\hsize)
		
		
	
	eststo clear
	
	}
	
	

