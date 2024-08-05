clear all
set more off

/*******************************************************************************

Step 0: Load Data

*******************************************************************************/



foreach sample in "" "if firstrbt!=." {
	
	local groups2 "True MPC:"
	local midrules2
	local groups1 "Specification: & \multicolumn{3}{c}{Table \ref{tab:Table_E1_d_pce}, Column 1} & \multicolumn{3}{c}{Table \ref{tab:Table_E1_d_pce}, Column 4} "
	local regcount = 0
		
	if regexm("`sample'","firstrbt")==1 {
		local samplename "rebateonly"
		local sampletitle "Rebate Only Sample"
	}
	else {
		local samplename ""
		local sampletitle "Full Sample"
	}
// 		local midrules1 `"\cmidrule(l{.75em}){2-`=`nregs'/2+1'}\cmidrule(l{.75em}){`=`nregs'/2+2'-`=`nregs'+1'}"'
// 		local groups1 `" &  \multicolumn{`=`nregs'/2'}{c}{Full Sample} &  \multicolumn{`=`nregs'/2'}{c}{Rebate Only Sample} "' 

	foreach gamma in 3 5 9 { /**/
	
		local groups2 `"`groups2' &  \multicolumn{1}{c}{0.`gamma'}  "'
	
		use "../output/simulbaselinegammampc0`gamma'.dta", clear
		
		tsset Household Date
		
		/*******************************************************************************
		
		Step 1: Create Treatment Variable
		
		*******************************************************************************/
		
		qui gen Totexp = Durable + Nondurable
		qui gen d_Totexp = d_Durable + d_Nondurable
		
		qui egen firstrbt = max(Date * (Rebate==1)), by(Household)
		
		qui replace firstrbt = . if firstrbt == 0
		
		
		local regcount = `regcount' + 2
		
		label variable Rebate "Rebate Coeff."
		
		
		/*******************************************************************************
		
		Step 2a: OLS Regression
		
		*******************************************************************************/
		

		ivreghdfe  RebateAmount Rebate   `sample' , absorb(Date)  
		
		local rebateamtOLS = _b[Rebate]
		
		ivreghdfe  d_Totexp Rebate  `sample' , absorb(Date)  
		
		qui eststo, prefix(ols)
		
		local temp_mpcOLS: di %3.2f _b[Rebate]/`rebateamtOLS'
		qui estadd local implied_mpc "`temp_mpcOLS'"

// 		local groups3 `"`groups3' &  \multicolumn{1}{c}{TWFE}  "'
		
		/*******************************************************************************
		
		Step 2b: Our specification
		
		*******************************************************************************/
		

		ivreghdfe  RebateAmount Rebate L3.Rebate   `sample' , absorb(Date)  
		
		local rebateamtOLS = _b[Rebate]
		
		ivreghdfe  d_Totexp Rebate  L3.Rebate L3.Totexp L3.Durable `sample' , absorb(Date)  
		
		qui eststo, prefix(new)
		
		local temp_mpcOLS: di %3.2f _b[Rebate]/`rebateamtOLS'
		qui estadd local implied_mpc "`temp_mpcOLS'"

// 		local groups3 `"`groups3' &  \multicolumn{1}{c}{TWFE}  "'		
		
		
// 		/*******************************************************************************
		
// 		Step 2b: BJS Regression
		
// 		*******************************************************************************/
		
// 		cd ../../did_imputation
		
		
		
// 			qui did_imputation RebateAmount Household Date firstrbt `sample',  fe(Date) horizons(0) autosample maxit(10000)
			
// 			local rebateamt = _b[tau0]
			
// 			qui did_imputation Totexp Household Date firstrbt `sample',  fe(Date)  horizons(0) autosample maxit(10000)
			
// 			qui eststo
			
// 			local temp_mpc: di %3.2f _b[tau0]/`rebateamt'
// 			qui estadd local implied_mpc "`temp_mpc'"
// 			local groups2 `"`groups2' &  \multicolumn{2}{c}{0.`gamma'}  "' 
// 			local midrules2 `"`midrules2'\cmidrule(l{.75em}){`regcount'-`=`regcount'+1'}"'	 
			
// 			local groups3 `"`groups3' &  \multicolumn{1}{c}{BJS}  "' 
		
		
// 		cd ../model/code
		
		
	}
	
	foreach gamma in 3 5 9 { 
	
		local groups2 `"`groups2' &  \multicolumn{1}{c}{0.`gamma'}  "'
	}


	/*******************************************************************************
		
	Step 3: Table
		
	*******************************************************************************/
	
	* this code produces the regression table
	local nregs = `regcount'
		
	
	
	local varorder `"Rebate"'	
	local order `"`varorder'"' 
	local indicate
	
	
	local groups `" "`groups1'" \\ "`midrules2'" "`groups2'" \\ "'
	local stats " implied_mpc "
	local stats_fmt " %3.2f "
	local stats_label `" `"Implied MPC"'  "'
	
	
	local num_stats: word count `stats' 
	local layout
	forvalues l = 1/`num_stats' {
		local layout `"`layout' "\multicolumn{1}{c}{@}" "'
	}
	local keepvars `"`varorder'"' 
	local dropvars 
	local filename `"modelregression`samplename'"'
	local tablename `"`filename'"'
		
	local title `"\caption{Rebate Coefficient MPC Estimates in Model: `sampletitle'}"'
	local notes `"Notes: The dependent variable is the change in total expenditures from the previous interview. The rebate size in the model is \$950 for all households conditional on receiving a rebate. The rebate only sample includes only households that receieve a rebate at some point during our sample period.  "'
	
	
	local table_preamble `" "\begin{table}[!t] \centering \sisetup{table-format=1.2} \def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}" "`title'" "\begin{tabularx}{\hsize}{@{\hskip\tabcolsep\extracolsep\fill}l*{`nregs'}{S}}" "\multicolumn{`=`nregs'+1'}{c}{Dependent Variable: Change in Total Expenditures} \\" "\hline" "'
	local postfoot `"postfoot(`"\hline\hline \end{tabularx} \begin{minipage}{\hsize} \rule{0pt}{9pt} \footnotesize `notes'  \end{minipage} \label{tab:`tablename'} \end{table}"')"'
	local prehead `"prehead(`table_preamble' `groups')"'			
// 	local posthead `"posthead(`"\hline"' )"'
	
	local prefoot(" ")
	
		
		
		
			esttab ols* new* using "../output/`filename'.tex", rename(tau0 Rebate __00000V rbtindicator) replace cells(b(fmt(a2)) ) nostar drop(`dropvars', relax) ///
			keep(`keepvars') indicate(`indicate') `prehead' `posthead' `postfoot' order(`order') label varlabel(`varlabels') stats(`stats', layout(`layout') fmt(`stats_fmt') ////
			labels(`stats_label')) collabels(,none) nonumbers nomtitles substitute(# `" X "' tabular* tabularx `"{1}{c}{("' `"{1}{L}{("') width(\hsize)
			
	eststo clear
	
	

}
