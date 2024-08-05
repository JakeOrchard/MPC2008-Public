

clear all

use ../input/psmjsampleinterview.dta





/****************************************************************************
Define interview schedule
*****************************************************************************/

gen _intmonth = month(dofm(intdate))

gen schedule = mod(_intmonth+2,3)

gen rbtmo = month(dofm(firstrbtdate))

*Percent of stimulus check receiving households that receive the check in May, June or July
gen mjj = rbtmo == 5 | rbtmo == 6 | rbtmo == 7
tab mjj if everrbtindicator == 1 





/*****************************************************************************
Output Table
*************************************************************************/



local checklist `"" "& check"'


foreach check in "`checklist'"{

	estimates drop _all
	eststo clear
	local ncols = 4

	
	local title `"Distribution of CEX Interview Schedule"'

	local notes `"Notes: Data in column 1 come from the entire CEX Sample 2007-2009. Data in columns 2-4 come from our subsample.  "'

	*Randomness of schedule
	if "`check'" == ""{
	eststo col1: estpost tab  schedule //Overall CEX 2007-2009
	
	local panelname `" &  \multicolumn{`ncols'}{c}{Panel A: EFT and Check Recipients } "'
	local colnames `" & \multicolumn{1}{c}{Overall CEX}   & \multicolumn{1}{c}{May Cohort}  & \multicolumn{1}{c}{June Cohort}  & \multicolumn{1}{c}{July Cohort} "'
	local filename `"interview_randomness"'
	local table_preamble `" "\begin{table}[!t] \centering \sisetup{table-format=3.2} \def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}" "\caption{`title'}" "\begin{tabularx}{\hsize}{@{\hskip\tabcolsep\extracolsep\fill}l*{`ncols'}{S}}" "\\" "\hline\hline" "'
		local postfoot `"postfoot(`"\hline"')"'

	}
	
	else if  "`check'" == "& check"{
	
	local extra "extrac(1)" //Mising column
	
	local panelname `" &  \multicolumn{`ncols'}{c}{Panel B: Check Recipients Only } "'
	local colnames `" &   & \multicolumn{1}{c}{May Cohort}  & \multicolumn{1}{c}{June Cohort}  & \multicolumn{1}{c}{July Cohort} "'
	local filename `"interview_randomnessPB"'
	local table_preamble `"  "\\" "\hline\hline" "'
	local postfoot `"postfoot(`"\hline \end{tabularx} \begin{minipage}{\hsize} \rule{0pt}{9pt} \footnotesize `notes'  \end{minipage} \label{tab:`filename'} \end{table}"')"'

	}
	
	eststo col2: estpost tab schedule if rbtmo == 5 & insample `check' //May Cohort
	eststo col3: estpost tab schedule if  rbtmo == 6 & insample `check' //June Cohort
	eststo col4: estpost tab schedule if  rbtmo == 7 & insample `check' //July Cohort
	
	*Test if largest number is = .33
	cap drop testmay testjune testjuly
	gen testmay = schedule == 2 if rbtmo == 5 & insample `check' //May Cohort
	gen testjune = schedule == 0 if  rbtmo == 6 & insample `check' //June Cohort
	gen testjuly = schedule == 1 if  rbtmo == 7 & insample `check' //July Cohort
	
	reg testmay, robust
	test _b[_cons] == .33
	
	reg testjune, robust
	test _b[_cons] == .33
	
	
	reg testjuly, robust
	test _b[_cons] == .33
	

	
	
	local indicate
	local midrules1 `" \cmidrule(l{.75em}){2-3} \cmidrule(l{.75em}){4-5} \cmidrule(l{.75em}){6-7} \cmidrule(l{.75em}){8-9} "'
	local groups `" " `panelname' \\ `colnames'" \\ "'
	local dropvars "Total"


	local prehead `"prehead(`table_preamble' `groups')"'			
	local posthead `"posthead(`"\hline"' `"\multicolumn{1}{l}{Interview Schedule} & \\"' `"\\"')"'
	
	
	local prefoot(" ")

	esttab * using "../output/`filename'.tex", `extra' replace cells(pct(fmt(0) par("" "\%"))) drop(`dropvars', relax)  indicate(`indicate') `prehead' `posthead' `postfoot' label coeflabels(0 "Jan-Apr-Jul-Oct" 1 "Feb-May-Aug-Nov" 2 "Mar-Jun-Sep-Dec") nonumbers stats(`stats', layout(`layout') fmt(`stats_fmt') labels(`stats_label')) collabels(,none) nomtitles substitute(# `" X "' tabular* tabularx `"{1}{c}{("' `"{1}{L}{("') width(\hsize)
	estimates drop _all
		

		eststo clear
}




*Shell script to create one table from the two panels:

shell cat ../output/interview_randomnessPB.tex >> ../output/interview_randomness.tex 




/*****************************************************************************
Test if share is = .33
**************************************************************************/




