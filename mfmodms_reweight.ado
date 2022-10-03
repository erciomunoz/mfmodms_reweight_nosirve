*! version 1.0   Ercio Munoz 10/3/2022 
	
/* 
- Program to reweight the population given demographic projections (age & education) and sector shares using wentropy. 
- For education, it assumes that young cohorts (i.e., 20-30) in the future will have similar levels
of education as the current young individuals.
- Note that observations with missing values for education are dropped
- Note that calif must have 3 categories (1=unskilled, 2=semi-skilled, and 3=skilled)
*/

cap program drop mfmodms_reweight
program define mfmodms_reweight
version 16.0

syntax [anything], AGE(string) EDUcation(string) GENDER(string) HHSize(string) /// 
 ID(string) IWeights(string) COUNTRY(string) IYear(string) TYear(string) ///
 Generate(string) MATCH(string) POPDATA(string) INDUSTRY(string) GROWTH(string) ///
 EMPLOYMENT(string) VARIANT(string) [ OMIT ] 

	** Installing needed software **
*	cap ssc install csipolate
*	cap ssc install gtools
	
********************************************************************************
	** Bringing UN population projections and aggregating cohorts ** 
********************************************************************************
quietly	{
preserve
			
	use "`popdata'" if country=="`country'" & Variant=="`variant'", clear

	sort country cohort
	replace cohort="P0009" if inlist(cohort,"P0004","P0509") /* 1st cohort */
	replace cohort="P1019" if inlist(cohort,"P1014","P1519") /* 2nd cohort */
	replace cohort="P2029" if inlist(cohort,"P2024","P2529") /* 3rd cohort */
	replace cohort="P3039" if inlist(cohort,"P3034","P3539") /* 4th cohort */
	replace cohort="P4049" if inlist(cohort,"P4044","P4549") /* 5th cohort */
	replace cohort="P5059" if inlist(cohort,"P5054","P5559") /* 6th cohort */
	replace cohort="P6069" if inlist(cohort,"P6064","P6569") /* 7th cohort */				
	replace cohort="P70up" if inlist(cohort,"P7074","P7579","P8084","P8589","P9094","P9599","P100up") /* 8th cohort, 70 and up */
	gcollapse (sum) yf* ym*, by(country cohort)
			
	gen cohort_merge = _n
	drop country cohort
	sort cohort_merge
	tempfile popdata
	save `popdata', replace
/* `popdata' has 8 cohorts and one column for each year/gender */		
				
	gen country_merge = "`country'"
	collapse (sum) y*, by(country_merge)
					
	forval p = 1991/2100 {
		gen double y`p'1 = yf`p' + ym`p'
	}
				
	drop ym* yf*
		
	sort country_merge
	tempfile popdatatot
	save `popdatatot', replace
	/* `popdatatot' has total population by year */
restore
	}

********************************************************************************		
	** Alpha parameter for excluding cells with few observations and verifying data
********************************************************************************
quietly	{
	* Pre-defined option		
	local adjust_missing = 1 

	local stp = (`tyear' - `iyear')
							
	* Verifying that EDUcation has not missing values
	keep if `education'!=.
	levelsof `education', local(alledus)
	
	* Verifying that the target year selected is within range
	if `tyear' < 1991 | `tyear' > 2100 {
		noi di " "
		noi di in red "Target year is out of range"
		exit
	}	

	* Saving original data 
		tempfile base
		save "`base'"
	}

********************************************************************************
	** Generating target matrix based on education and age
********************************************************************************
quietly {
	* Generate age - cohorts
	tempvar cohort
	tempvar pop
	gen `cohort'=8
	label var `cohort' "Age Group"
	
	* Generate Age Groups
	forval x = 1/8 {
		replace `cohort' = `x' if `age' >=(`x'-1)*10 & `age' <= ((`x'-1)*10)+9
	}
			
	gen double `pop'=1
	
	collapse (sum) `pop'  [fw=round(`iweights')], by(`cohort' `education' `gender')

	bys `education': table `cohort' `gender', c(sum `pop')

	reshape wide `pop' , i(`cohort' `gender') j(`education') 
	reshape wide `pop'*, i(`cohort') j(`gender') 
		
	sort `cohort'
	mvencode _all, mv(0) // Change missing values to numeric values

	tempvar pop1	// Population for Males (gender=1)
	tempvar pop2    // Population for Women (gender=2)

	egen double `pop1' = rsum(`pop'*1) 
	egen double `pop2' = rsum(`pop'*2)
		
	replace `pop1'=. if `pop1'==0
	replace `pop2'=. if `pop2'==0
		
	* Calculate Shares
		
	foreach j of local alledus {
		gen double sh`j'1 = `pop'`j'1/`pop1'
		gen double sh`j'2 = `pop'`j'2/`pop2'
	}
			
	list `cohort' `pop1' `pop2' sh11 sh21 sh31 sh12 sh22 sh32
	sort `cohort'

	local steps = round((`tyear' - `iyear')/10)+2

	forval step = 1/`steps' {
		
		local yyyy = `iyear' + 10*(`step'-1)
		local yyy0 = `yyyy' - 10
		
		gen sh11_`yyyy'=. 
		gen sh21_`yyyy'=.
		gen sh31_`yyyy'=.

		gen sh12_`yyyy'=. 
		gen sh22_`yyyy'=.
		gen sh32_`yyyy'=.
			
		if `iyear'==`yyyy' {

			replace sh11_`yyyy' = sh11
			replace sh21_`yyyy' = sh21
			replace sh31_`yyyy' = sh31
				
			replace sh12_`yyyy' = sh12
			replace sh22_`yyyy' = sh22
			replace sh32_`yyyy' = sh32
		
		}
			
			/* Here we use the same shares from baseline to 3 youngest cohorts */
		if `iyear'!=`yyyy' {
			replace sh11_`yyyy'=sh11 if `cohort'<=3
			replace sh21_`yyyy'=sh21 if `cohort'<=3
			replace sh31_`yyyy'=sh31 if `cohort'<=3

			replace sh12_`yyyy'=sh12 if `cohort'<=3
			replace sh22_`yyyy'=sh22 if `cohort'<=3
			replace sh32_`yyyy'=sh32 if `cohort'<=3
	
			local count1 = ((`yyyy' - `iyear')/10)
			di `count1'
			
			/* Replacing share for cohort i with share of previous cohort getting older */
			forval i = 1/`count1' {
				replace sh11_`yyyy'=sh11_`yyy0'[_n-1] if `cohort'==3+`i'
				replace sh21_`yyyy'=sh21_`yyy0'[_n-1] if `cohort'==3+`i'
				replace sh31_`yyyy'=sh31_`yyy0'[_n-1] if `cohort'==3+`i'

				replace sh12_`yyyy'=sh12_`yyy0'[_n-1] if `cohort'==3+`i'
				replace sh22_`yyyy'=sh22_`yyy0'[_n-1] if `cohort'==3+`i'
				replace sh32_`yyyy'=sh32_`yyy0'[_n-1] if `cohort'==3+`i'
			}	
		
				
			local count2 = 8 - ((`yyyy' - `iyear')/10) + 3	
			local count3 =      ((`yyyy' - `iyear')/10)		
		
			forval i = 1/`count2' {
		
				replace sh11_`yyyy'=sh11_`iyear'[_n-`count3'] if `cohort'==3+`count3'+`i'
				replace sh21_`yyyy'=sh21_`iyear'[_n-`count3'] if `cohort'==3+`count3'+`i'
				replace sh31_`yyyy'=sh31_`iyear'[_n-`count3'] if `cohort'==3+`count3'+`i'

				replace sh12_`yyyy'=sh12_`iyear'[_n-`count3'] if `cohort'==3+`count3'+`i'
				replace sh22_`yyyy'=sh22_`iyear'[_n-`count3'] if `cohort'==3+`count3'+`i'
				replace sh32_`yyyy'=sh32_`iyear'[_n-`count3'] if `cohort'==3+`count3'+`i'

			}
				
		}
	}
		
		if `adjust_missing' == 1 {
			gen R1m = sh11==0|sh21==0|sh31==0
			gen R1f = sh12==0|sh22==0|sh32==0			
			
			forval step = 1/`steps' {
			
				local yyyy = `iyear' + 10*(`step'-1)

		 		replace sh11_`yyyy' = 0 if R1m==1 & sh11==0
				replace sh21_`yyyy' = 0 if R1m==1 & sh21==0
		 		replace sh31_`yyyy' = 0 if R1m==1 & sh31==0

			 	replace sh12_`yyyy' = 0 if R1f==1 & sh12==0
				replace sh22_`yyyy' = 0 if R1f==1 & sh22==0
			 	replace sh32_`yyyy' = 0 if R1f==1 & sh32==0

				egen _Tsh1_`yyyy' = rsum(sh*1_`yyyy')
				egen _Tsh2_`yyyy' = rsum(sh*2_`yyyy')

		 		replace sh11_`yyyy' = sh11_`yyyy'/_Tsh1_`yyyy'
	 			replace sh21_`yyyy' = sh21_`yyyy'/_Tsh1_`yyyy'
	 			replace sh31_`yyyy' = sh31_`yyyy'/_Tsh1_`yyyy'

		 		replace sh12_`yyyy' = sh12_`yyyy'/_Tsh2_`yyyy'
		 		replace sh22_`yyyy' = sh22_`yyyy'/_Tsh2_`yyyy'
	 			replace sh32_`yyyy' = sh32_`yyyy'/_Tsh2_`yyyy'

		 		drop _T*
		 	
		 		}	
		 	
			 drop R1m R1f
		 
		 }
		
		
		* Interpolation 
		
		tempvar last
		gen `last' = 1
		
		keep `pop1' `pop2' `cohort' sh11_`iyear' - `last'
					
		reshape long sh11_ sh21_ sh31_ sh12_ sh22_ sh32_, i(`cohort')
		rename _j year
		levelsof year, local(yrs)
		
		foreach var of local yrs {
			expand 10 if year==`var'
		}
		
		foreach var of varlist sh11_ - sh32_ {
		bys `cohort' year: replace `var'=. if _n>1
		}
		
		bys `cohort' year: replace year = year[_n-1]+1 if _n>=2
			drop if year>2100
		

	
		foreach var of varlist sh11_ - sh32_ {
			bys `cohort': ipolate `var' year, gen (e_`var')
			replace `var' = e_`var' if `var'==.
			drop e_`var'
		}

		
		levelsof year, local(allyears)
		
		reshape wide sh*, i(`cohort') j(year)
		
		gen cohort_merge = `cohort'
		sort cohort_merge
		
		merge cohort_merge using "`popdata'"
			
		tab _merge
		keep if _merge==3
		drop _merge 
		
		noi di "`allyears'"
		
	foreach y of local allyears {
	 
		tempvar sumpop`y'
		tempvar UNsumpop`y'
			
		if "`match'" == "UN" {
			egen double `sumpop`y'' = sum(yf`y' + ym`y')
			sum `sumpop`y''
			local lsumpop`y' = r(mean)
		}
			
		if "`match'" == "HH" {
		
			if `y'==`iyear' {
				egen double `sumpop`y'' = sum(`pop1' + `pop2')
				sum `sumpop`y''
				local lsumpop`y' = r(mean)
				noi di "Population in year: `y': `lsumpop`y''" 
			}
			
			if `y'> `iyear' {
				local z = `y' - 1	// Option for all years
				
				
				replace `pop1' = (`pop1')*(ym`y'/ym`z')		// pop1 is for males
				replace `pop2' = (`pop2')*(yf`y'/yf`z')		// pop2 is for females
			
				egen double `sumpop`y'' = sum(`pop1' + `pop2')
				sum `sumpop`y''
				local lsumpop`y' = r(mean)
				noi di "Population in year: `y': `lsumpop`y''" 
			}
		}			
		
		if "`match'" == "UN" {	
			replace sh11_`y' = (sh11_`y' * ym`y')/`lsumpop`y''
			replace sh21_`y' = (sh21_`y' * ym`y')/`lsumpop`y''
			replace sh31_`y' = (sh31_`y' * ym`y')/`lsumpop`y''
				
			replace sh12_`y' = (sh12_`y' * yf`y')/`lsumpop`y''
			replace sh22_`y' = (sh22_`y' * yf`y')/`lsumpop`y''
			replace sh32_`y' = (sh32_`y' * yf`y')/`lsumpop`y''
		}
			
		if "`match'" == "HH" {
			replace sh11_`y' = (sh11_`y' * `pop1')/`lsumpop`y''
			replace sh21_`y' = (sh21_`y' * `pop1')/`lsumpop`y''
			replace sh31_`y' = (sh31_`y' * `pop1')/`lsumpop`y''
				
			replace sh12_`y' = (sh12_`y'* `pop2')/`lsumpop`y''
			replace sh22_`y' = (sh22_`y'* `pop2')/`lsumpop`y''
			replace sh32_`y' = (sh32_`y'* `pop2')/`lsumpop`y''
		}
		
	}
		
	drop yf* ym* 
	
		forval a = `iyear'/`tyear' {	
		mkmat sh12_`a'                 	, matrix(A`a') 
		mkmat sh22_`a'   				, matrix(B`a')
		mkmat sh32_`a' if `cohort' <=8	, matrix(C`a')
		mkmat sh11_`a'                 	, matrix(D`a')
		mkmat sh21_`a' 					, matrix(E`a')
		mkmat sh31_`a' if `cohort' <=8 , matrix(F`a')
		}
		
		
		forval a = `iyear'/`tyear' {
			matrix define const`a' = [ A`a' \ B`a' \ C`a' \ D`a' \ E`a' \ F`a' ]
		}

		
*		matrix list const`tyear'
}				
	
********************************************************************************
	**  Applying the calibration command
********************************************************************************		
	use `base', clear

if "`omit'" == "" {
		
		foreach g in 2 1 {
			foreach e of local alledus {
				forval i = 1/8 {
					
					if `g' == 2 local gm = "f"	// Females
					if `g' == 1 local gm = "m"	// Males

				
					local j = (`i'-1)*10
				
					if `i'==1 {	
						gen double c`e'`gm'0`j' = (age >= (`i'-1)*10 & age<= ((`i'-1)*10)+9 & `gender'==`g' & calif==`e')
					}

					else {

						gen double c`e'`gm'`j' = (age >= (`i'-1)*10 & age<= ((`i'-1)*10)+9 & `gender'==`g' & calif==`e')
					}
				}
			}
		}
		
	qui replace `industry'=99 if missing(`industry')
	qui ta `industry', g(industry)	
	local r = `r(r)'-1
	qui svyset [pw=`iweights']
	qui svy: total industry1-industry`r'
	mat shares = e(b)'	
	mata : st_matrix("employment2", colsum(st_matrix("shares"))) /* Number of employed individuals at baseline */
	mat shares = hadamard(shares,`growth')
	mata : st_matrix("employment3", colsum(st_matrix("shares"))) /* Number of employed individuals after applying VA-elasticities */

	mat shares = shares/`lsumpop`tyear'' /* Target is employment as share of population */	
	
	mat shares = shares*( ( (employment2[1,1]/`lsumpop`iyear'')*`employment') / (employment3[1,1]/`lsumpop`tyear'')) /* Re-scaling the number to match employment levels */

	matrix const`tyear' = const`tyear'\shares
	
	* This is to generate the constrains for the maxentropy, first, at the household level. 
	gcollapse (mean) c1f00-c3m70 industry1-industry`r', by(`id' `hhsize' `iweights')
			
	noi di "Wentropy for country `country' in year `tyear'"
	noi di "The constraint matrix for year `tyear'"
	matrix list const`tyear'
	noi di "Wentropy c1f00-c3m70 industry1-industry`r', old(`iweights') new(newwgt) constraints(const`tyear') pop(`lsumpop`tyear'')"	
	
	wentropy c1f00-c3m70 industry1-industry`r', old(`iweights') new(newwgt) constraints(const`tyear') pop(`lsumpop`tyear'')
	
	qui gen newiwgt = newwgt/hhsize
	qui keep `id' newiwgt
	clonevar idh_merge = `id'
	sort idh_merge
	tempfile max`tyear'
	qui save `max`tyear'', replace
	
use `base', clear
	
	clonevar idh_merge = `id'
	sort idh_merge
	qui merge m:1 idh_merge using `max`tyear''
	qui tab _merge
	qui drop _merge
	sort idh_merge
	qui drop idh_merge
	qui gen double `generate'`tyear' = round(newiwgt)
	drop newiwgt
	
}
end
