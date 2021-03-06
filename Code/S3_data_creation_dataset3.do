
/*************************************************/
/* Data creation for Dataset 3                   */
/* Outcome based on A,L,W                        */
/* Missingness based on W or A,W                 */
/*************************************************/

set more off
capture program drop datagen
global cloudstor "C:\Users\z3312911\Cloudstor"
sysdir set PLUS "$cloudstor\ado\"

program datagen3

scalar drop _all

// Define scalars which are passed to the program by the simulation code
local numdata="`1'" // Number of simulation iterations

// Set intercepts for missing data mechanisms
// This varies the amount of missing data without changing the associations with other variables
	matrix missconsmod=(0.25,-3.85,-4.35)
	matrix missconssev=(0.50,-1.50,-1.95)

// Loop over number of simulation iterations (currently set at 500)	
forvalues dat=1/`numdata' {
	clear
	set obs 500
	gen id=_n

	// Define coefficients for each of the data generating mechanisms
	matrix lco=(2.0,1.5,0.2,0.5,-0.8,0.5)
	matrix aco=(1.5,1.0,0.5,1.0,0.5,0.5,-0.5,0.7,0.3)
	matrix marco=(0.5,0.5,1.0,1.5,0.5,0.1,-0.2,0.1)
	matrix outco=(0.5,1.0,0.5,1.0,0.5,0.5,0.8,-0.5)
	
	// Define correlations for baseline and latent variables
	matrix C1=(1,0.3,0.3 \ 0.3,1,0.3 \ 0.3,0.3,1)
	matrix C2=(1,0.9,0.8,0.7,0.6 \ 0.9,1,0.9,0.8,0.7 \ 0.8,0.9,1,0.9,0.8 \ 0.7,0.8,0.9,1,0.9 \ 0.6,0.7,0.8,0.9,1)

	// Create 3xbaseline, correlated, time-constant variables
	drawnorm oa ob oc, corr(C1) means(0 0 0) sds(1 1 1)
	// Create 2xtime-varying, correlated, unobserved (latent) variables
	drawnorm u0 u1 u2 u3 u4, corr(C2) means(0 0 0 0 0) sds(1 1 1 1 1)
	drawnorm v0 v1 v2 v3 v4, corr(C2) means(0 0 0 0 0) sds(1 1 1 1 1)
	
	// Define subject specific error terms for each variable
	gen el_ijk=rnormal(0,1)
	gen ew_ijk=rnormal(0,1)
	gen ea_ijk=rnormal(0,1)
	gen em_ijk=rnormal(0,1)
	gen ey_ijk=rnormal(0,1)
		
	reshape long u v, i(id) j(obs)
	xtset id obs
	
	// Initialise variables
	gen ll=rbinomial(1,0.25)
	gen lw=rbinomial(1,0.25)
	gen la=rbinomial(1,0.20)
	gen l=.
	gen w=.
	gen a=.
	gen y=.
	gen missmcarmod=0
	gen missmarmod1=0
	gen missmarmod2=0
	gen missmcarsev=0
	gen missmarsev1=0
	gen missmarsev2=0

	// Loop over waves to create variable at each wave (allowing variables to be associated with the previous wave)
	forvalues i=0/4 {
		
		replace ll=l1.l if obs==`i' & l1.l!=.
		replace lw=l1.w if obs==`i' & l1.w!=.
		replace la=l1.a if obs==`i' & l1.a!=.
		
		replace l = runiform()<invlogit(-3.7 + ///
		lco[1,1]*ll + lco[1,2]*la + lco[1,3]*obs + ///
		lco[1,4]*oa + lco[1,5]*ob + lco[1,6]*oc + 2.0*u + el_ijk) if obs==`i'
		
		replace w = runiform()<invlogit(-3.7 + ///
		lco[1,1]*lw + lco[1,3]*obs + ///
		lco[1,4]*oa + lco[1,5]*ob + lco[1,6]*oc + 2.0*v + ew_ijk) if obs==`i'
		
		replace a = runiform()<invlogit(-2.3 + ///
		aco[1,1]*la + aco[1,2]*l + aco[1,3]*ll + aco[1,6]*obs + ///
		aco[1,7]*oa + aco[1,8]*ob + aco[1,9]*oc + ea_ijk) if obs==`i'
		
		replace y = rnormal(5 + ///
		outco[1,1]*la + outco[1,2]*a + outco[1,3]*ll + outco[1,4]*l - outco[1,3]*lw - outco[1,4]*w + outco[1,5]*obs + ///
		outco[1,6]*oa + outco[1,7]*ob + outco[1,8]*oc + ///
		3.0*u + ey_ijk,1) if obs==`i'

		// Moderate missing data ~25%
		// MCAR is totally random based on scalar %
		replace missmcarmod = runiform()<missconsmod[1,1] if obs==`i' & obs>0
		// 2 versions of MAR, which differ based on whether missingness is associated with exposure
		replace missmarmod1 = runiform()<invlogit(missconsmod[1,2] + ///
		marco[1,3]*lw + marco[1,4]*w + ///
		marco[1,5]*obs - marco[1,6]*oa + marco[1,7]*ob + marco[1,8]*oc + 2.0*v + em_ijk) if obs==`i' & obs>0
		replace missmarmod2 = runiform()<invlogit(missconsmod[1,3] + ///
		marco[1,1]*la + marco[1,2]*a + marco[1,3]*lw + marco[1,4]*w + ///
		marco[1,5]*obs - marco[1,6]*oa + marco[1,7]*ob + marco[1,8]*oc + 2.0*v + em_ijk) if obs==`i' & obs>0
				
		// Severe missing data ~50%
		// MCAR is totally random based on scalar %
		replace missmcarsev = runiform()<missconssev[1,1] if obs==`i' & obs>0
		// 2 versions of MAR, which differ based on whether missingness is associated with exposure
		replace missmarsev1 = runiform()<invlogit(missconssev[1,2] + ///
		marco[1,3]*lw + marco[1,4]*w + ///
		marco[1,5]*obs - marco[1,6]*oa + marco[1,7]*ob + marco[1,8]*oc + 2.0*v + em_ijk) if obs==`i' & obs>0
		replace missmarsev2 = runiform()<invlogit(missconssev[1,3] + ///
		marco[1,1]*la + marco[1,2]*a + marco[1,3]*lw + marco[1,4]*w + ///
		marco[1,5]*obs - marco[1,6]*oa + marco[1,7]*ob + marco[1,8]*oc + 2.0*v + em_ijk) if obs==`i' & obs>0
	}

	drop if obs==0
	keep id y a l w la ll lw oa ob oc obs missmcarmod missmarmod1 missmarmod2 missmcarsev missmarsev1 missmarsev2
	order id obs oa ob oc ll lw la l w a y missmcarmod missmarmod1 missmarmod2 missmcarsev missmarsev1 missmarsev2
		
	save "$cloudstor\PhD\Paper 4 - Missing data simulation\Data\S3\Simulation `dat'.dta", replace
	
}
	
end

local numdata=1000

set seed 24543
datagen3 `numdata'

// Quick program to double-check the amount of missing data is correct
capture program drop getscalars
program getscalars

global i=$i+1
use "$cloudstor\PhD\Paper 4 - Missing data simulation\Data\S3\Simulation $i.dta", clear
	
	mean missmcarmod missmarmod1 missmarmod2 missmcarsev missmarsev1 missmarsev2 if obs>0
	scalar mmcarmod=_b[missmcarmod]
	scalar mmar1mod=_b[missmarmod1]
	scalar mmar2mod=_b[missmarmod2]
	scalar smcarsev=_b[missmcarsev]
	scalar smar1sev=_b[missmarsev1]
	scalar smar2sev=_b[missmarsev2]
	
end

global i=0
simulate mmcarmod=mmcarmod mmar1mod=mmar1mod mmar2mod=mmar2mod ///
smcarsev=smcarsev smar1sev=smar1sev smar2sev=smar2sev, reps(1000): getscalars
su, sep(3)


// Quick program to double-check the amount of missing data is correct
capture program drop testcollider
program testcollider

global i=$i+1
use "$cloudstor\PhD\Paper 4 - Missing data simulation\Data\S3\Simulation $i.dta", clear
	
	logit a w l la lw ll if missmarmod1!=1
	scalar wmarmod1=_b[w]
	logit a w l la lw ll if missmarmod2!=1
	scalar wmarmod2=_b[w]
	logit a w l la lw ll if missmarsev1!=1
	scalar wmarsev1=_b[w]
	logit a w l la lw ll if missmarsev2!=1
	scalar wmarsev2=_b[w]

	
end

global i=0
simulate wmarmod1=wmarmod1 wmarmod2=wmarmod2 wmarsev1=wmarsev1 wmarsev2=wmarsev2, reps(1000): testcollider
su, sep(3)

