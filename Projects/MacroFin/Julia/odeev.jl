#=
	author: Jonathan Payne
=#

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
        x0, t0 = rk4step(fdash, t0, x0, h)
        wi[1:neqn, i+1] = x0
    end
    return wi, ti
end

# ------------------------------
# Usage example 1:
# ------------------------------
t0, x0, tf, N = 0, [0,-4], 10, 100
A = [[1 2], [3 2]]
fdash(x, t) = A*x
xtrue(t) = [ 8/5*exp(-t) - 8/5*exp(4*t), -8/5*exp(-t) - 12/5*exp(4*t) ]
tspan = [t0, tf]
wi2, ti2 = rk4uniform(fdash, x0, tspan, N)

using ODE
t0, x0, tf, N = 0, [0,-4], 10, 100
A = [[1 2], [3 2]]
fdash(t, x) = A*x
xtrue(t) = [ 8/5*exp(-t) - 8/5*exp(4*t), -8/5*exp(-t) - 12/5*exp(4*t) ]
tspan = [t0, tf]
span = linspace(t0, tf, N)
ode23s(fdash, x0, span)


# ------------------------------
# odelay (from the pycse package)
# ------------------------------
# Solve an ODE with events.
#
# Inputs:
# func        = is callable, with signature func(Y, x)
# y0          = are the initial conditions
# xspan       = what you want to integrate over
# events      = is a list of callable functions with signature event(Y, x).
# 			    These functions return zero when an event has happened.
# [value, isterminal, direction] = event(Y, x)
# value       = the value of the event function. When value = 0, an event
# 			    is triggered
# isterminal  = True if the integration is to terminate at a zero of
# 			    this event function, otherwise, False.
# direction   = 0 if all zeros are to be located (the default),
#			  = +1 if only zeros where the event function is increasing, and
#             = -1 if only zeros where the event function is decreasing.
# TOLERANCE   = what is used to identify when an event has occurred.
# fsolve_args = a dictionary of options for fsolve
# kwargs  	  = additional options you want to send to odeint.
#
# Returns [x, y, te, ye, ie]
# x           = the independent variable array
# y           = the solution
# te          = an array of independent variable values where events occurred
# ye          = an array of the solution at the points where events occurred
# ie          = an array of indices indicating which event function occurred.

function odelayold(func, y0, xspan,
				 events, TOLERANCE = 1e-6)#, fsolve_args=[], args...)
	#if fsolve_args == None
    #    fsolve_args = []
    #end
    N = length(xspan)
    x0 = xspan[1] # initial point

    sol = [float64(y0)]
    TE, YE, IE = [], [], [] # to store where events occur

    # initial value of events
    e = zeros(length(events), N)
    for (i, event) in enumerate(events)
        e[i, 1], isterminal, direction = event(y0, x0)
    end

    # now we step through the integration
	for (i, x1) in enumerate(xspan[1:end-1])
		println("i", i)
        x2 = xspan[i + 1]
        f1 = sol[i]

        f2, wz = rk4uniform(func, f1', [x1, x2], 1)#, args...)
       	
        #X[i+1] = x2
        #println(X[i+1])
        X = push!(X, x2)
        sol = push!(sol, f2[:,end][1])
        
        # check event functions. At each step we compute the event
        # functions, and check if they have changed sign since the
        # last step. If they changed sign, it implies a zero was
        # crossed.        
        for (j, event) in enumerate(events)
     		e[j, i + 1], isterminal, direction = event(sol[i + 1], X[i + 1])
                
            if ((e[j, i + 1] * e[j, i] < 0) || abs(e[j, i + 1]) < TOLERANCE || abs(e[j, i]) < TOLERANCE)
                xLt = X[end]       # Last point
                fLt = sol[end]
                eLt = e[j, i+1]
        
	            function objective(x)
	            	# evaluate ode from xLT to x
	             	txspan = [xLt, x]
	             	tempsol, wyz = rk4uniform(func, fLt, txspan, length(txspan))#, args...)
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
                tempsol, wyz = rk4uniform(func, fLt, txspan)#, args...)
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
                        sol[end] = fZ[1]
                        return X, sol, TE, YE, IE
                    end
                end

            end
        end
    end
    # at the end, return what we have
    return X, sol, TE, YE, IE

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





# ------------------------------
# Example usage : 1 dimension
# ------------------------------

# myode(f, x) = 3*x^2 + 12*x -4

# function event1(f, x)
#     # An event is when f = 0 and event is decreasing'
#     isterminal = true
#     direction = -1
#     return f, isterminal, direction
# end

# function event2(f, x)
#     # An event is when f = 0 and increasing'
#     isterminal = false
#     direction = 1
#     return f, isterminal, direction
# end

# f0 = -120

# xspan = linspace(-8, 4, 50)
# events=[event1, event2]
# func = myode; y0 = f0;
# X, F, TE, YE, IE = odelay(myode, f0, xspan, events)

# ------------------------------
# Example usage : 2 dimensions
# ------------------------------

# function myode(z, t)
#     # PROJ: ODE for projectile motion with linear air resistance
#     # z = [z0, z1, z2, z3]
#     g = 9.81 # Specify approximate graviatational constant
#     b = 0.28 # Representative value
#     return [z[3], z[4], -g-b*z[3], -b*z[4]]
# end

# function linevent(z, t)
#     # LINEVENT: Contains the event we are looking for
#     # In this event, z(1) = 0 (hitting the ground)
#     lookfor = z[1] # Sets this to 0
#     isterminal = true # Stop when event is located
#     direction = -1 # specify downward direction
#     return lookfor, isterminal, direction
# end

# f0 = [0.18, 0, 5, 10]
# xspan = linspace(0, 10, 50)
# events = [linevent]
# func = myode
# y0 = f0
# Xm, Fm, TEm, YEm, IEm = odelay(myode, f0, xspan, events)















