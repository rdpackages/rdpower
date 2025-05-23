********************************************************************************
* Auxiliary functions for rdpow package
* Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
********************************************************************************
*!version 2.2 22-May-2025

version 13

********************************************************************************
*** Power function

capture mata: mata drop rdpower_powerfun()
mata:
function rdpower_powerfun(real scalar n,real scalar tau,real scalar stilde,real scalar z)
{
	x = 1-normal(sqrt(n)*tau/stilde+z)+normal(sqrt(n)*tau/stilde-z)
	return(x)
}
mata mosave rdpower_powerfun(), replace
end

********************************************************************************
*** Power function derivative w.r.t n

capture mata: mata drop rdpower_powerfun_dot()
mata:
function rdpower_powerfun_dot(real scalar n,real scalar tau,real scalar stilde,real scalar z)
{
	x = (normalden(sqrt(n)*tau/stilde-z)-normalden(sqrt(n)*tau/stilde+z))*tau/(2*stilde*sqrt(n))
	return(x)
}
mata mosave rdpower_powerfun_dot(), replace
end

********************************************************************************
*** Power function derivative w.r.t tau

capture mata: mata drop rdpower_powerfun_dot_tau()
mata:
function rdpower_powerfun_dot_tau(real scalar n,real scalar tau,real scalar stilde,real scalar z)
{
	x = (normalden(sqrt(n)*tau/stilde-z)-normalden(sqrt(n)*tau/stilde+z))*sqrt(n)/stilde
	return(x)
}
mata mosave rdpower_powerfun_dot_tau(), replace
end


********************************************************************************
*** Newton-Raphson for sample size calculation

capture mata: mata drop rdpower_powerNR()
mata:
void rdpower_powerNR(real scalar x0,real scalar tau,real scalar stilde,real scalar z,real scalar beta)
{
	tol = 1
	iter = 0
	while (tol>epsilon(1)){
		++iter
		k = 1	
		// check if derivative at x0 is too small:
		while (rdpower_powerfun_dot(x0,tau,stilde,z)<.00001){
			x0 = 1.2*x0*(rdpower_powerfun(x0,tau,stilde,z)<=beta) + 0.8*x0*(rdpower_powerfun(x0,tau,stilde,z)>beta)
			++iter
		}
		
		x1 = x0 - (rdpower_powerfun(x0,tau,stilde,z)-beta)/rdpower_powerfun_dot(x0,tau,stilde,z)
		
		// check if x1 is negative or too small:
		while (x1<2){
			x1 = x0 - k*(rdpower_powerfun(x0,tau,stilde,z)-beta)/rdpower_powerfun_dot(x0,tau,stilde,z)
			k = k/2
			++iter
		}
		if (x1==.){
			break
		}
		tol = abs(rdpower_powerfun(x1,tau,stilde,z)-beta)
		x0 = x1
	}
	b = rdpower_powerfun(x1,tau,stilde,z)
	x1 = ceil(x1)
	st_numscalar("m_rdpower",x1)
	st_numscalar("iter_rdpower",iter)
	st_numscalar("powercheck_rdpower",b)
}
mata mosave rdpower_powerNR(), replace
end

********************************************************************************
*** Newton-Raphson for MDE calculation

capture mata: mata drop rdpower_powerNR_mde()
mata:
void rdpower_powerNR_mde(real scalar n,real scalar tau0,real scalar stilde,real scalar z,real scalar beta)
{
	tol = 1
	iter = 0
	while (tol>epsilon(1)){
		++iter
		// check if derivative at tau0 is too small:
		while (rdpower_powerfun_dot_tau(n,tau0,stilde,z)<.00001){	
			tau0 = 1.5*tau0*(rdpower_powerfun(n,tau0,stilde,z)<=beta) + 0.5*tau0*(rdpower_powerfun(n,tau0,stilde,z)>beta)
			++iter
		}

		tau1 = tau0 - (rdpower_powerfun(n,tau0,stilde,z)-beta)/rdpower_powerfun_dot_tau(n,tau0,stilde,z)
		
		if (tau1==.){
			break
		}
		tol = abs(rdpower_powerfun(n,tau1,stilde,z)-beta)
		tau0 = tau1
	}
	b = rdpower_powerfun(n,tau1,stilde,z)
	st_numscalar("mde_rdpower",tau1)
	st_numscalar("iter_rdpower",iter)
	st_numscalar("powercheck_rdpower",b)
}
mata mosave rdpower_powerNR_mde(), replace
end

