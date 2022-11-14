//////////////
/*
Toy model as in De Groote (2022): A Dynamic Model of Effort Choice in High 
School.

I simulate data of a 2 period binary school choice model in which students also 
choose their level of effective effort. In the first period, students decide to 
go to high school and, if they do, set their optimal level of effort. Graduation 
is uncertain, but a function of effective effort. If they graduate, they can 
decide to enter college in period 2 and receive its total expected lifetime 
utility. Drop out from school is a terminal action and its utility is normalized 
to 0. 

The notation is explained below and follows the paper unless mentioned 
otherwise. 

I then estimate two models. First, I estimate a misspecified model: the pure 
discrete choice model. This is a model that does not allow for changes in effort. 
I then repeat the exercise for the correctly specified model.

Predicted values from the estimated models are in CAPITAL LETTERS. 

To illustrate the different impact of dynamic incentives, I run a counterfactual
in which only a 50% random sample of high school graduates is allowed to enter 
college. In the status quo everyone can enter.
*/
//////////////

//set version of STATA
version 15

//set path
/*Not needed because of sysuse*/

//open log
cap log close
*log using ...,replace

//seed
set seed 10

/////////////////////////////////////////////////////
//////////// CREATE DATASET /////////////////////////
/////////////////////////////////////////////////////

//load some data Stata has on individual characteristics
sysuse nlsw88,clear
gen white = (race==1)
keep idcode white south

//make dataset bigger
expand 10

drop idcode
gen idcode=_n

//make panel
expand 2
sort idcode
by idcode: gen t=_n
xtset idcode t

//draw four extreme value type 1 taste shocks
gen ev11=-ln(-ln(runiform())) 
gen ev01=-ln(-ln(runiform())) 
gen ev12=-ln(-ln(runiform())) 
gen ev02=-ln(-ln(runiform())) 

//draw 1 logistic performance shock
gen eta=-ln(-ln(runiform())) - (-ln(-ln(runiform())) )
by idcode: replace eta=eta[1]

//fixed cost of going to (high) school in period 1 (note: "fc" is "C^0" in paper)
gen fc = 1 - 1*white + 0.5*south 

//expected lifetime utility of going to school (college) in period 2
gen psi = 3 + 0.5*white - 0.5*south 
	
//expected lifetime utility of no school in period 1 -> 0
//expected lifetime utility of no school in period 2 -> 0

//marginal costs of effort (note: "mc" is "c" in paper) 
gen mc = exp(-1.8-0.2*white-1*south)


////solve in status quo scenario where everyone with HS degree can enter college
//effort
gen y=sqrt(0.95*ln(1+exp(psi))/mc)-1
sum y

//variable cost
gen vc = mc*y

//flow utility
gen u = -fc -vc

//success probability after period 1
gen phi = invlogit(ln(y))

//emax
gen emax = 0.5772 + phi*ln(1+exp(psi))
sum emax
drop phi

//schooling in period 1
gen school = (u+0.95*emax+ev11 > ev01) if t==1
sum school if t==1

//success
gen success = (ln(y)+eta)>0 if t==1 & school==1
replace success = 0 if t==1 & school==0

//schooling in period 2
replace school = (psi+ev12 > 0+ev02) if t==2 & L.school==1 & L.success==1
replace school=0 	if t==2 & (L.school==0 | (L.school==1 & L.success==0))

//summarize
sum mc y vc fc u psi
sum school success if t==1
sum school if t==2


////solve in counterfactual: 50% allowed to enter college
//effort
gen y_counter = sqrt(0.95*0.5*ln(1+exp(psi))/mc)-1
sum y_counter

//allow for corner solutions in counterfactuals
replace y_counter = max(0,y_counter) 						
sum y_counter

//variable cost
gen vc_counter = mc*y_counter

//flow utility
gen u_counter = -fc-vc_counter

//success prob
gen phi_counter=invlogit(ln(y_counter))

//emax
gen emax_counter = 0.5772 + phi_counter*0.5*ln(1+exp(psi))
drop phi_counter

//schooling in period 1
gen school_counter = (u_counter+0.95*emax_counter+ev11 > ev01) if t==1

//success
gen success_counter = (ln(y_counter)+eta)>0 if t==1 & school==1
replace success_counter = 0 if t==1 & school_counter==0

//schooling in period 2
replace school_counter = (psi + ev12 > 0+ev02) ///
					if t==2 & L.school_counter==1 & L.success_counter==1
replace school_counter=0 ///
				if t==2 & (L.school_counter==0 ///
				| (L.school_counter==1 & L.success_counter==0))

//50% chance to be allowed to enter college in period 2
gen unif=runiform() if t==2
replace school_counter=0 if t==2 & unif>0.5

//summarize
global tabstatopts varwidth(20) stats(mean) format(%9.3g) c(s)
tabstat mc y_counter vc_counter fc u_counter psi,$tabstatopts
tabstat school_counter success_counter if t==1,$tabstatopts
tabstat school_counter  if t==2,$tabstatopts



///////////////////////////////////////////////////////////////////////////
/////////////////// START ESTIMATION PURE DISCRETE CHOICE MODEL ///////////
//////////////////////////////////////////////////////////////////////////
		
//college stage (period 2)	
logit school white south if t==2 & L.school==1 & L.success==1	
est store college						
predict COLLEGE						

//success
global vars white south
logit success $vars if t==1 & school==1 	
*logit success $vars c.( $vars )#c.( $vars ) if t==1 & school==1 //more flexible spec. 	
predict PHI
predict Y,xb
replace Y=exp(Y)
		
//calculate emax
est restore college
predict PSI,xb
gen EMAX = 0.5772 + PHI*ln(1+exp(PSI))
	
//estimate high school stage (period 1)
constraint 1 EMAX=0.95
logit school white south EMAX if t==1,constraint(1)
est store hs_dynamic
predict HS_DYNAMIC_PR
	
//predictions
gen SCHOOL = HS_DYNAMIC if t==1
replace SCHOOL = HS_DYNAMIC*PHI*COLLEGE if t==2
gen SUCCESS = HS_DYNAMIC*PHI if t==1

//results
tabstat Y y SCHOOL school SUCCESS success if t==1,$tabstatopts
tabstat SCHOOL school if t==2,$tabstatopts

/*
CONCLUSION: Fitting the data in the status quo with the misspecified model is 
not the issue.
*/


/////// run counterfactual: 50% allowed in college

//effort exogenous
gen Y_COUNTER = Y
gen PHI_COUNTER = PHI
				
//prediction of the new emax
gen EMAX_COUNTER = 0.5772 + 0.5*PHI_COUNTER*ln(1+exp(PSI))
					
//counterfactual prediction of conditional value function in high school 
est restore hs_dynamic
gen V_COUNTER=_b[_cons]+_b[white]*white ///
					+_b[south]*south ///
					+_b[EMAX]*EMAX_COUNTER 
		
//predictions
gen SCHOOL_COUNTER = invlogit(V_COUNTER) if t==1
replace SCHOOL_COUNTER = invlogit(V_COUNTER)*0.5*PHI_COUNTER*COLLEGE if t==2
gen SUCCESS_COUNTER = SCHOOL_COUNTER*PHI_COUNTER  if t==1
										
//results
tabstat Y_COUNTER y_counter SCHOOL_COUNTER school_counter ///
		SUCCESS_COUNTER success_counter if t==1,$tabstatopts
tabstat SCHOOL_COUNTER school_counter if t==2,$tabstatopts

/*
CONCLUSION: because effort cannot be adjusted, we see overprediction of success 
and underprediction of schooling in period 1. Period 2 shows the resulting bias 
in college enrollment numbers: the overprediction of success is driving the bias 
upwards, while the underprediction of schooling period 1 is driving it downwards. 
The overprediction of success dominates, leading to a positive bias.
*/

//drop predictions
drop COLLEGE-SUCCESS_COUNTER



///////////////////////////////////////////////////////////////////////
/////////////////// START ESTIMATION WITH EFFORT CHOICE ///////////////
//////////////////////////////////////////////////////////////////////

//college stage (period 2)	
logit school white south ///
			if t==2 & L.school==1 & L.success==1			
est store college					
predict COLLEGE						

//success 
logit success $vars if t==1 & school==1 	
*logit success $vars c.( $vars )#c.( $vars ) if t==1 & school==1 //more flexible spec. 	
predict PHI
predict Y,xb
replace Y=exp(Y)
		
//calculate emax
est restore college
predict PSI,xb
gen EMAX = 0.5772 + PHI*ln(1+exp(PSI))
	
//mc
gen MC = 0.95*ln(1+exp(PSI))/((1+Y)^2)
gen VC = MC*Y
	
//estimate high school stage (period 1)
constraint 1 EMAX=0.95
constraint 2 VC=-1
logit school white south EMAX VC if t==1,constraint(1 2) col
est store hs_dynamic
predict HS_DYNAMIC
predict FC,xb
replace FC = -(FC-0.95*EMAX+VC)
	
//predictions
gen SCHOOL = HS_DYNAMIC if t==1
replace SCHOOL = HS_DYNAMIC*PHI*COLLEGE if t==2
gen SUCCESS = HS_DYNAMIC*PHI if t==1

//results
tabstat Y y SCHOOL school SUCCESS success if t==1,$tabstatopts
tabstat SCHOOL school if t==2,$tabstatopts

/*
CONCLUSION: Also a good fit of the data in the status quo.
*/


/////// run counterfactual: 50% allowed in college
		
//prediction of new effort level, vc and prob
gen Y_COUNTER = sqrt(0.95*0.5*ln(1+exp(PSI))/MC)-1
replace Y_COUNTER = max(0,Y_COUNTER) 

gen VC_COUNTER = MC*Y_COUNTER
gen PHI_COUNTER = invlogit(ln(Y_COUNTER))
		
//prediction of the new emax
gen EMAX_COUNTER = 0.5772 + 0.5*PHI_COUNTER*ln(1+exp(PSI))
					
//counterfactual prediction of conditional value function in high school 
est restore hs_dynamic
gen V_COUNTER=_b[_cons]+_b[white]*white ///
					+_b[south]*south ///
					+_b[EMAX]*EMAX_COUNTER ///
					+_b[VC]*VC_COUNTER
		
//predictions
gen SCHOOL_COUNTER = invlogit(V_COUNTER) if t==1
replace SCHOOL_COUNTER = invlogit(V_COUNTER)*0.5*PHI_COUNTER*COLLEGE if t==2
gen SUCCESS_COUNTER=SCHOOL_COUNTER*PHI_COUNTER if t==1
			
//results
tabstat Y_COUNTER y_counter SCHOOL_COUNTER school_counter ///
			SUCCESS_COUNTER success_counter if t==1,$tabstatopts
tabstat SCHOOL_COUNTER school_counter if t==2,$tabstatopts

drop COLLEGE-SUCCESS_COUNTER
	

/*
CONCLUSION: this estimation approach does allow effort to be adjusted in
a counterfactual, removing all biases.
*/

			
cap log close
	