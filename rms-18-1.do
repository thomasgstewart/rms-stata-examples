*-------------------------------------------------------------------------------
* rms 18.1
*
*-------------------------------------------------------------------------------
getvdata prostate

* Code missing as missing
replace age = . if age == .z
replace wt = . if wt == .z
replace ekg = . if ekg == .z
replace sz = . if sz == .z
replace sg = . if sg == .z
* Code time of death = .5 for folks that died in first month 
replace dtime = .5 if dtime == 0

* multiple imputation (will only do single imputation for demonstration purposes)
misstable pattern
mi set flong
mi register imputed age wt sz ekg sg
mi impute chained (pmm, knn(10)) age wt sz sg (ologit) ekg = i.rx dtime status i.pf i.hx sbp dbp hg i.bm ap, add(1) burnin(20)
drop if _mi_m == 0
mi unset, asis


* Variable transformations

* all cause mortality variable
gen acm = 1*(status >= 2)

* log of ap
gen logap = log(ap)

* Combine into mean arterial bp
gen map = (2*dbp + sbp)/3

* Combine ekg and hx
gen heart = hx + 1*(ekg >=3)

* keep variables which will be in the model
keep dtime acm rx age wt pf heart map hg sg sz logap bm


* Setup up Cox model data
stset dtime, failure(acm==1) 

* Kaplan-Meier curve (no covariate adjustment)
sts graph

* plotdata for partial effect plots after the model is fit
centile age wt pf heart map hg sg sz logap bm, centile(50)
* reference: plotdata, at(rx=1 age=73 wt=98 pf=1 heart=1 map=10 hg=14 sg=10 sz=11 logap=-.35 bm=0)
plotdata, at(rx=1; age=55/82; wt=98; pf=1; heart=1; map=10; hg=14; sg=10; sz=11; logap=-.35; bm=0)
plotdata, at(rx=1; age=73; wt=75/130; pf=1; heart=1; map=10; hg=14; sg=10; sz=11; logap=-.35; bm=0)
plotdata, at(rx=1; age=73; wt=98; pf=1; heart=1; map=10; hg=9(.5)17; sg=10; sz=11; logap=-.35; bm=0)

* Construct restricted cubic splines
mkspline_plotindicator age, nknots(4)
mkspline_plotindicator wt, nknots(3)
mkspline_plotindicator map, nknots(3)
mkspline_plotindicator hg, nknots(4)
mkspline_plotindicator sg, nknots(3)
mkspline_plotindicator sz, nknots(3)
mkspline_plotindicator logap, nknots(5)

* Fit Cox model
stcox rcs* i.rx pf heart i.bm if _plotindicator == 0, efron

* check proportional hazard assumption
cox_lhhat_ci, xb
egen float xbgroup = cut(time_xbhat) if _plotindicator == 0, group(4)
stphplot, by(xbgroup)

* Outliers
predict dfb_*, dfbeta
gen index = _n
twoway (scatter dfb_1 index)
twoway (scatter dfb_2 index)
* ... twoway (scatter dfb_26 index)

* Hypothesis Tests
* rx (association)
test (2.rx 3.rx 4.rx)
* age (association)
test (rcs_age_1 rcs_age_2 rcs_age_3)
* age (linearity)
test (rcs_age_2 rcs_age_3)
* pf (association)
test (pf)

* Measures of model performance
estat concordance

* Partial effect plots on log hazards scale
twoway ///
  (rarea time_xbhat_lb time_xbhat_ub age, sort fcolor(orange_red%50) lcolor(orange_red)) ///
  (line time_xbhat age, sort fcolor(orange_red%50) lcolor(orange_red) lwidth(thick)) ///
  if _plotindicator == 1 ///
  , ytitle(log relative hazard)
  
twoway ///
  (rarea time_xbhat_lb time_xbhat_ub wt, sort fcolor(orange_red%50) lcolor(orange_red)) ///
  (line time_xbhat wt, sort fcolor(orange_red%50) lcolor(orange_red) lwidth(thick)) ///
  if _plotindicator == 2 ///
  , ytitle(log relative hazard)
  
twoway ///
  (rarea time_xbhat_lb time_xbhat_ub hg, sort fcolor(orange_red%50) lcolor(orange_red)) ///
  (line time_xbhat hg, sort fcolor(orange_red%50) lcolor(orange_red) lwidth(thick)) ///
  if _plotindicator == 3 ///
  , ytitle(log relative hazard)
  
* Partial effect plots on 5 year probability scale
cox_lhhat_ci, spt(60) drop
twoway ///
  (line sp_60 age, sort fcolor(orange_red%50) lcolor(orange_red) lwidth(thick)) ///
  if _plotindicator == 1 ///
  , ytitle(P(5 year survival))

twoway ///
  (line sp_60 wt, sort fcolor(orange_red%50) lcolor(orange_red) lwidth(thick)) ///
  if _plotindicator == 2 ///
  , ytitle(P(5 year survival))
  
twoway ///
  (line sp_60 hg, sort fcolor(orange_red%50) lcolor(orange_red) lwidth(thick)) ///
  if _plotindicator == 3 ///
  , ytitle(P(5 year survival))

* Partial effect plots on median survival time
cox_lhhat_ci, qt(50 75) drop
twoway ///
  (line q50_time age, sort fcolor(orange_red%50) lcolor(orange_red) lwidth(thick)) ///
  if _plotindicator == 1 ///
  , ytitle(Median survival time (months))

twoway ///
  (line q50_time wt, sort fcolor(orange_red%50) lcolor(orange_red) lwidth(thick)) ///
  if _plotindicator == 2 ///
  , ytitle(Median survival time (months))
  
twoway ///
  (line q50_time hg, sort fcolor(orange_red%50) lcolor(orange_red) lwidth(thick)) ///
  if _plotindicator == 3 ///
  , ytitle(Median survival time (months))



* Generate a difference in prediction plot
* Fit Cox model
* Age effect comparing bm = 1 to bm = 0 (ref)
stcox rcs* i.rx pf heart i.bm c.(rcs_age*)#i.bm if _plotindicator == 0, efron

gen lb = .
gen ub = .
gen p = .
local N = _N

	forvalues i = 1/`N' {
		if _plotindicator[`i'] == 1 {
			local a1 = rcs_age_1[`i']
			local a2 = rcs_age_2[`i']
			local a3 = rcs_age_3[`i']
			lincom _b[1.bm#c.rcs_age_1]*`a1' + _b[1.bm#c.rcs_age_2]*`a2' + _b[1.bm#c.rcs_age_3]*`a3' +  _b[1.bm]
			replace ub = r(ub) if _n == `i'
			replace lb = r(lb) if _n == `i'
			replace p = r(estimate) if _n == `i'
		}
	}

twoway ///
  (line lb age, sort lcolor(blue) lwidth(thin) lpattern(dash)) ///
  (line ub age, sort lcolor(blue) lwidth(thin) lpattern(dash)) ///
  (line p age, sort lcolor(blue) lwidth(thick) lpattern(solid)) ///
  if _plotindicator == 1
  
  
// Hazard ratio plot for age (age = 60 is reference) 
gen lb2 = .
gen ub2 = .
gen p2 = .
local N = _N


        forvalues i = 1/`N' { 
                if _plotindicator[`i'] == 1 {
                        local a1 = rcs_age_1[`i']
                        local a2 = rcs_age_2[`i']
                        local a3 = rcs_age_3[`i']
			local ref1 = rcs_age_1[508]
			local ref2 = rcs_age_2[508]
			local ref3 = rcs_age_3[508]
                        lincom _b[c.rcs_age_1]*(`a1'-`ref1') + _b[c.rcs_age_2]*(`a2'-`ref2') + _b[c.rcs_age_3]*(`a3'-`ref3'), hr
                        replace ub2 = r(ub) if _n == `i'
                        replace lb2 = r(lb) if _n == `i'
                        replace p2 = r(estimate) if _n == `i'
                }
        }
replace lb2 = 1 if _n == 508
replace ub2 = 1 if _n == 508

twoway ///
  (line lb2 age, sort lcolor(blue) lwidth(thin) lpattern(dash)) ///
  (line ub2 age, sort lcolor(blue) lwidth(thin) lpattern(dash)) ///
  (line p2 age, sort lcolor(blue) lwidth(thick) lpattern(solid)) ///
  if _plotindicator == 1 ///
  , legend(order(3 "Hazard ratio (age 60 is reference"))
