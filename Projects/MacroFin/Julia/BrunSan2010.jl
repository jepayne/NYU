#=
	filename: BrunSan2010.jl
	author: Jonathan Payne

	cd("/Users/Juddy/GitHub/NYU/Projects/MacroFin/Julia")

	This file replicates Brunnermeir Sannikov in Julia

	This file computes the equilibrium in the setting of Brunnermeier and
	Sannikov under the assumptions that (1) experts can issue only debt and
	not equity and (2) the production set of experts is (a - iota, Phi(iota)
	- delta) and households, (a_ - iota, Phi(iota) - delta_), where the
	function [Phi iota] = investment(q) solves for Phi and iota

=#

#include("odeev.jl")
global qmax

using Roots





function rk4step(fdash, t, x, h)
    # Calculates the the nth step in the classical fourth order Runga-Kutta
    # procedure for solving a differential equation of the form:
    #   x'(t) = fdash(x, t),
    # Inputs:
    #   * fdash: user defined function
    #   * t: previous value of independent variable
    #   * x: previous value of dependent variable
    #   * h: increment
    # Outputs:
    #   * xnew: next step
    #   * tnew: new time
    k1 = h*fdash(x, t)
    k2 = h*fdash(x + k1/2, t + h/2)
    k3 = h*fdash(x + k2/2, t + h/2)
    k4 = h*fdash(x + k3, t + h)
    xnew = x + (k1 + 2*k2 + 2k3 + k4) / 6
    tnew = t + h
    return xnew, tnew
end

function rk4uniform(fdash, x0, tspan, N=100)
    # Approximate the solution of the initial value problem using the
    # classical fourth order Runge-Kutta method with uniform step sizes:
    #   * x'(t) = fdash(x, t),
    #   * x(t0) = x0
    # Inputs:
    #   * fdash: user defined function
    #   * t0: initial value of the independent variable
    #   * x0: initial value of the dependent variable(s) (column vector)
    #   * tf: final value of the independent variable
    #   * N: number of universally sized time steps to be taken to advance the
    #        the solution from t = t0 to t = tf
    # Outputs:
    #   * wi: vector / matrix containing the values of the approximate solution
    #         to the ode
    #   * ti: number of uniformly sized time steps to be taken to advance t
    t0, tf = tspan
    neqn = length(x0)
    ti = linspace(t0, tf, N+1)
    wi = zeros(neqn, N+1)
    wi[1:neqn, 1] = x0
    h = (tf - t0) / N
    for i = 1:N
        #t = t0; x = x0
        x0, t0 = rk4step(fdash, t0, x0, h)
        wi[1:neqn, i+1] = x0
    end
    return wi, ti
end

function odelay(func, y0, xspan,
                 events, TOLERANCE = 1e-6)#, fsolve_args=[], args...)
    #if fsolve_args == None
    #    fsolve_args = []
    #end
    N = length(xspan)
    x0 = xspan[1] # initial point

    X = [x0]
    Nf = length(y0)
    sol = zeros(N, Nf)
    sol[1,:] = float64(y0)
    TE, YE, IE = [], [], [] # to store where events occur

    # initial value of events
    e = zeros(length(events), N)
    for (i, event) in enumerate(events)
        e[i, 1], isterminal, direction = event(y0, x0)
    end

    # now we step through the integration
    for (i, x1) in enumerate(xspan[1:end-1])
        x2 = xspan[i + 1]
        f1 = sol[i, :]

        #fdash = func; x0 = f1'; tspan = [x1, x2]; N = 1
        f2, wz = rk4uniform(func, f1', [x1, x2], 1)#, args...)
        
        #X[i+1] = x2
        #println(X[i+1])
        X = push!(X, x2)
        sol[i+1,:] = f2[:, end]
        
        # check event functions. At each step we compute the event
        # functions, and check if they have changed sign since the
        # last step. If they changed sign, it implies a zero was
        # crossed.        
        for (j, event) in enumerate(events)
            e[j, i + 1], isterminal, direction = event(sol[i + 1], X[i + 1])
                
            if ((e[j, i + 1] * e[j, i] < 0) || abs(e[j, i + 1]) < TOLERANCE || abs(e[j, i]) < TOLERANCE)
                xLt = X[end]       # Last point
                fLt = sol[i+1,:]
                eLt = e[j, i+1]
        
                function objective(x)
                    # evaluate ode from xLT to x
                    txspan = [xLt, x]
                    tempsol, wyz = rk4uniform(func, fLt', txspan, length(txspan))#, args...)
                    tsol = tempsol[:, end]
                    val, isterminal, direction = event(tsol, x)
                    return val[1]
                end
                xZ = fzero(objective, xLt)#, **fsolve_args)  # this should be the
                                              # value of x that makes
                                              # the event zero
                # now evaluate solution at this point, so we can
                # record the function values here.
                txspan = [xLt, xZ]
                tempsol, wyz = rk4uniform(func, fLt', txspan, length(txspan))#, args...)
                fZ = tempsol[:,end]

                vZ, isterminal, direction = event(fZ, xZ)

                COLLECTEVENT = false
                if direction == 0
                    COLLECTEVENT = true
                elseif (e[j, i + 1] > e[j, i] ) && direction == 1
                    COLLECTEVENT = true
                elseif (e[j, i + 1] < e[j, i] ) && direction == -1
                    COLLECTEVENT = true
                end

                if COLLECTEVENT
                    push!(TE, xZ)
                    push!(YE, fZ)
                    push!(IE, j)

                    if isterminal
                        X[end] = xZ[1]
                        sol[i+1,:] = fZ[1]
                        return X, sol[1:i+1,:], TE, YE, IE
                    end
                end

            end
        end
    end
    # at the end, return what we have
    return X, sol, TE, YE, IE

end









# Function finds the optimal investment rate iota and rate of capital creation
# Phi given the price q
function investment(q)
	# Input:
	# 	q = price of capital
	# Output:
	# 	Phi, iota
	Phi = (q - 1)/kappa
	iota = Phi + kappa*Phi^2/2
	return Phi, iota
end


a      = 0.11  # productivity of expert
a_     = 0.05  # productivity of household
rho    = 0.06  # expert discount rate
r      = 0.05  # riskfree rate of return
sigma  = 0.025 # volatility of capital ito process
delta  = 0.03  # expert discount rate
delta_ = 0.08  # household discount rate
kappa  = 10    # Adjustment cost parameter, such that investment
                # ...iota = Phi + kappa Phi^2/2 is required
                # ...to build capital Phi, is defined by investment(q)


# Some preparation: determine first best Q, and the worst Q
# This algorithm solves the following fixed point problem iteratively:
# q = max_{iota(q)} (a - iota)/(r - (Phi(iota(q)) - delta))
QL = 0 # initial lower bound on price
QR = 10 # initial upper bound on prices
for iter = 1:50
    qmax = (QL + QR)/2
    Phi, iota = investment(qmax)
    value = (a - iota)/(r + delta - Phi)
    if iota > a # iota is outside the bounds 0 < iota < a
        QR = qmax # iota is too large because the price is too large
    elseif value > qmax # qmax is to the right of fixed point q*
        QL = qmax # price too large
    elseif r + delta < Phi # Adjustment cost greater than output
        println("first-best price is infinite") # warn that first-best price is
        # infinite, but allow computation to continue
        QR = qmax  # without allowing q to grow above the level where
        # Phi > r + delta 
    else
        QR = qmax
    end
end

# determine q_, the value of capital at eta = 0
QL = 0 # initial lower bound on price
QR = 10 # initial upper bound on prices
q_ = -1
for iter = 1:50
    q_ = (QL + QR)/2
    Phi, iota = investment(q_)
    value = (a_ - iota)/(r + delta_ - Phi)
    if iota > a_
        QR = q_ # price is too large
    elseif value < q_  # price still too large
        QR = q_
    else
        QL = q_
    end
end

# evntfcn
# ------------------------------
# [value,isterminal,direction] = evntfcn(eta, F) returns three
# variables used by ode45 to determine when to terminate integration
function event1(F, eta)
    value = qmax - F[3] # difference between qmax and q, 
    isterminal = true # terminate computation in all three cases
    direction = 0 # event occurs whether we get there from above or below
    return value, isterminal, direction
end

function event2(F, eta)
    value = F[2] # first derivative of theta
    isterminal = true
    direction = 0
    return value, isterminal, direction
end

function event3(F, eta)
    value = F[4] # first derivative of q
    isterminal = true
    direction = 0
    return value, isterminal, direction
end
events = [event1, event2, event3]

# fnct
# ------------------------------
# This funciton implements proposition II.4. That is, it runs an algorithm
# for finding psi(eta) and coefficients of the ODEs for q(eta) and theta(eta)
# given (eta, q(eta), q'(eta), theta(eta), theta'(eta))
#
# Inputs:
# eta = scalar
# f = [theta, theta', q, q'], a 4x1 vector
# r, rho, a, a_, delta, delta_, sigma = parameters
# 
# Outputs:
# fp = the derivative of f with respect to eta
# dyn = [psi, sigma_eta*eta, sigma_q, mu_eta*eta, mu_q, iota, leverage, rk, r_k]
#
# assumes that the production set of experts is (a - iota, Phi(iota) - delta)
# and households, (a_ - iota, Phi(iota) - delta_), where the function
# [Phi iota] = investment(q) solves for Phi and iota

function fnct(eta, f, r, rho, a, a_, delta, delta_, sigma)
    # search for psi between eta (lower bound) and
    # min(f(3)/f(4) + eta, 1) (upper bound)
    Phi, iota = investment(f[3])

    psi = 0
    sigma_eta_eta,sigma_q,sigma_theta,risk_premium,household_premium = 0,0,0,0,0
    psi_L = eta; psi_R = min(f[3]/f[4] + eta, 1)

    for n = 1:50
        psi = (psi_L + psi_R)/2
        amplification = 1 - f[4]/f[3]*(psi - eta)
        
        # VOLATILITY COMPUTATION
        sigma_eta_eta = sigma*(psi - eta)/amplification  # sigma_eta *times* eta
        sigma_q = sigma_eta_eta*f[4]/f[3]
        sigma_theta = sigma_eta_eta*f[2]/f[1]
        risk_premium = - sigma_theta*(sigma + sigma_q)
        
        household_premium = (a_ - a)/f[3] + delta - delta_ + risk_premium
               
        if household_premium > 0 # households want to hold more
            psi_R = psi
        else
            psi_L = psi
        end
    end
    (psi, sigma_eta_eta, sigma_q,
     sigma_theta, risk_premium, household_premium) = (
     float64(psi), float64(sigma_eta_eta), float64(sigma_q),
     float64(sigma_theta), float64(risk_premium), float64(household_premium))

    mu_q = r - (a - iota)/f[3] - Phi + delta - sigma*sigma_q + risk_premium
    mu_eta_eta = (-(psi - eta)*(sigma + sigma_q)*(sigma + sigma_q + sigma_theta)
                  + eta*(a - iota)/f[3] + eta*(1 - psi)*(delta_ - delta))
    qpp = 2*(mu_q*f[3] - f[4]*mu_eta_eta)/sigma_eta_eta^2
    thetapp = 2*((rho - r)*f[1] - f[2]*mu_eta_eta)/sigma_eta_eta^2 

    fp = [f[2], thetapp, f[4], qpp];

    leverage = psi/eta
    rk = r + risk_premium
    r_k = r + household_premium

    dyn = [psi, sigma_eta_eta, sigma_q, mu_eta_eta, mu_q, iota, leverage, rk, r_k]
    return fp, dyn

end

# MAIN PART OF THE CODE:
# solves the system of ODE's starting with boundary conditions F0, i.e.
# theta(0) = 1, theta'(0) = #large negative number# and q(0) = q_; 
# searching for q'(0) such that q'(eta*) and theta'(eta*) reach 0 
# at the same boundary eta*

etaspan = linspace(0,1,50)
#F0 = [1 -1e+10 q_ 0]'  # [theta(0), theta'(0), q(0), q'(0)]
F0 = [1 -1000 q_ 0]'

# note that if theta(eta) satisfies the system of ODE's then so does
# const*theta(eta), for any positive constant const
# hence we set theta(0) = 1, and then normalize

odefun(f, eta) = fnct(eta, f, r, rho, a, a_, delta, delta_, sigma)[1]

QL = 0; QR = 1000#QR = 1e+15;
for iter = 1:50
    F0[4] = (QL + QR)/2  # this is q'(0)
    func = odefun; y0 = F0; xspan = etaspan
    [etaout,fout,TE,YE,IE] = odelay(odefun, F0, etaspan, events)  
    if IE[1] == 3 # if q'(eta) has reached zero, we 
        QL[1] = F0[4]  # increase q'(0)
    else        # if q(eta) reached qmax or theta'(0) reached 0 
        QR[1] = F0[4]  # we reduce q'(0)
    end

end

# here we are basically done... let me just compute all other variables
# from the same function fnct
# N = length(etaout)
# dynout = zeros(N, 9)
# for n = 1:N
#     [fp dynout(n,:)] = fnct(etaout[n], fout[n,:], r, rho, a, a_, delta, delta_, sigma);
# end

# # normalize theta, to make sure that theta(eta*) = 1
# normalization = fout(N,1);
# fout(:,1:2) = fout(:,1:2)/normalization;

































