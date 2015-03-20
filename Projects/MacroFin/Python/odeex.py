"""
Example implementation of odes in python

"""

import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt

# --------------------
# odelay (from the pycse package)
# --------------------

def odelay(func, y0, xspan, events, TOLERANCE = 1e-6, fsolve_args=None, **kwargs):
    '''Solve an ODE with events.
    func is callable, with signature func(Y, x)
    y0 are the initial conditions xspan is what you want to integrate
    over
    events is a list of callable functions with signature event(Y, x).
    These functions return zero when an event has happened.
    TOLERANCE is what is used to identify when an event has occurred.
    
    [value, isterminal, direction] = event(Y, x)
    value is the value of the event function. When value = 0, an event
    is triggered
    isterminal = True if the integration is to terminate at a zero of
    this event function, otherwise, False.
    direction = 0 if all zeros are to be located (the default), +1
    if only zeros where the event function is increasing, and -1 if
    only zeros where the event function is decreasing.
    fsolve_args is a dictionary of options for fsolve
    
    kwargs are any additional options you want to send to odeint.
    Returns [x, y, te, ye, ie]
    x is the independent variable array
    y is the solution
    te is an array of independent variable values where events occurred
    ye is an array of the solution at the points where events occurred
    ie is an array of indices indicating which event function occurred.
    '''
    if 'full_output' in kwargs:
        raise Exception('full_output not supported as an option')

    if fsolve_args is None:
        fsolve_args = {}

    x0 = xspan[0]  # initial point

    X = [x0]
    sol = [y0]
    TE, YE, IE = [], [], [] # to store where events occur
    
    # initial value of events
    e = np.zeros((len(events), len(xspan)))
    for i,event in enumerate(events):
        e[i,0], isterminal, direction = event(y0, x0)

    # now we step through the integration
    for i, x1 in enumerate(xspan[0:-1]):
        x2 = xspan[i + 1]
        f1 = sol[i]

        f2 = odeint(func, f1, [x1, x2])#, **kwargs)
        
        X += [x2]
        sol += [f2[-1,:]]

        # check event functions. At each step we compute the event
        # functions, and check if they have changed sign since the
        # last step. If they changed sign, it implies a zero was
        # crossed.        
        for j, event in enumerate(events):
            e[j, i + 1], isterminal, direction = event(sol[i + 1], X[i + 1])
                
            if ((e[j, i + 1] * e[j, i] < 0)        # sign change in
                                                   # event means zero
                                                   # crossing
                or np.abs(e[j, i + 1]) < TOLERANCE # this point is
                                                   # practically 0
                or np.abs(e[j, i]) < TOLERANCE):

                xLt = X[-1]       # Last point
                fLt = sol[-1]
                eLt = e[j, i+1]

                # we need to find a value of x that makes the event zero
                def objective(x):
                    # evaluate ode from xLT to x
                    txspan = [xLt, x]
                    tempsol = odeint(func, fLt, txspan)#, **kwargs)
                    sol = tempsol[-1, :]
                    val, isterminal, direction = event(sol, x)
                    return val

                from scipy.optimize import fsolve
                xZ, = fsolve(objective, xLt, **fsolve_args)  # this should be the
                                              # value of x that makes
                                              # the event zero

                # now evaluate solution at this point, so we can
                # record the function values here.
                txspan = [xLt, xZ]
                tempsol = odeint(func, fLt, txspan, **kwargs)
                fZ = tempsol[-1,:]

                vZ, isterminal, direction = event(fZ, xZ)

                COLLECTEVENT = False
                if direction == 0:
                    COLLECTEVENT = True
                elif (e[j, i + 1] > e[j, i] ) and direction == 1:
                    COLLECTEVENT = True
                elif (e[j, i + 1] < e[j, i] ) and direction == -1:
                    COLLECTEVENT = True
                
                if COLLECTEVENT:
                    TE.append(xZ)
                    YE.append(fZ)
                    IE.append(j)

                    if isterminal:
                        X[-1] = xZ
                        sol[-1] = fZ
                        return (np.array(X), 
                                np.array(sol), 
                                np.array(TE), 
                                np.array(YE), 
                                np.array(IE))

    # at the end, return what we have
    return (np.array(X), 
            np.array(sol), 
            np.array(TE), 
            np.array(YE), 
            np.array(IE))

# ----------------------------------------
# Example usage : Univariate
# ----------------------------------------
# def myode(f, x):
#     return 3*x**2 + 12*x -4

# def event1(f, x):
#     'an event is when f = 0 and event is decreasing'
#     isterminal = True
#     direction = -1
#     return f, isterminal, direction

# def event2(f, x):
#     'an event is when f = 0 and increasing'
#     isterminal = False
#     direction = 1
#     return f, isterminal, direction

# f0 = -120

# xspan = np.linspace(-8, 4)
# events=[event1, event2]
# func = myode
# y0 = f0
# X, F, TE, YE, IE = odelay(myode, f0, xspan, events)

#import matplotlib.pyplot as plt
#plt.plot(X, F, '.-')

# plot the event locations.use a different color for each event
#colors = 'rg'

#for x,y,i in zip(TE, YE, IE):
#    plt.plot([x], [y], 'o', color=colors[i])

#plt.savefig('images/event-ode-2.png')
#plt.show()
#print(TE, YE, IE)


# ----------------------------------------
# Example usage: Multivariate
# ----------------------------------------
# System of differential equations:
# z1' = z3
# z2' = z4
# z3' = -g - b*z3
# z4' = -bz4

def myode(z, t):
    # PROJ: ODE for projectile motion with linear air resistance
    # z = [z0, z1, z2, z3]
    g = 9.81 # Specify approximate graviatational constant
    b = 0.28 # Representative value
    return [z[2], z[3], -g-b*z[2], -b*z[3]]

def linevent(z, t):
    # LINEVENT: Contains the event we are looking for
    # In this event, z(1) = 0 (hitting the ground)
    lookfor = z[0] # Sets this to 0
    stop = 1 # Stop when event is located
    direction = -1 # specify downward direction
    return lookfor, stop, direction

f0 = [0.18, 0, 5, 10]
xspan = np.linspace(0, 10)
events = [linevent]
func = myode
y0 = f0
X, F, TE, YE, IE = odelay(myode, f0, xspan, events)





# --------------------
# odelaym (adapated from the pycse package to account for multivariate case)
# --------------------
# I think there are mistakes in the above code that only show up in the
# multivariate case

def odelaym(func, y0, xspan, events, TOLERANCE = 1e-6, fsolve_args=None, **kwargs):
    '''Solve an ODE with events.
    func is callable, with signature func(Y, x)
    y0 are the initial conditions xspan is what you want to integrate
    over
    events is a list of callable functions with signature event(Y, x).
    These functions return zero when an event has happened.
    TOLERANCE is what is used to identify when an event has occurred.
    
    [value, isterminal, direction] = event(Y, x)
    value is the value of the event function. When value = 0, an event
    is triggered
    isterminal = True if the integration is to terminate at a zero of
    this event function, otherwise, False.
    direction = 0 if all zeros are to be located (the default), +1
    if only zeros where the event function is increasing, and -1 if
    only zeros where the event function is decreasing.
    fsolve_args is a dictionary of options for fsolve
    
    kwargs are any additional options you want to send to odeint.
    Returns [x, y, te, ye, ie]
    x is the independent variable array
    y is the solution
    te is an array of independent variable values where events occurred
    ye is an array of the solution at the points where events occurred
    ie is an array of indices indicating which event function occurred.
    '''
    if 'full_output' in kwargs:
        raise Exception('full_output not supported as an option')

    if fsolve_args is None:
        fsolve_args = {}

    x0 = xspan[0]  # initial point

    X = [x0]
    sol = [y0]
    TE, YE, IE = [], [], [] # to store where events occur
    
    # initial value of events
    e = np.zeros((len(events), len(xspan)))
    for i,event in enumerate(events):
        e[i,0], isterminal, direction = event(y0, x0)

    # now we step through the integration
    for i, x1 in enumerate(xspan[0:-1]):
        x2 = xspan[i + 1]
        f1 = sol[i]

        f2 = odeint(func, f1, [x1, x2])#, **kwargs)
        
        X += [x2]
        sol += [f2[-1, :]]

        # check event functions. At each step we compute the event
        # functions, and check if they have changed sign since the
        # last step. If they changed sign, it implies a zero was
        # crossed.        
        for j, event in enumerate(events):
            e[j, i + 1], isterminal, direction = event(sol[i + 1], X[i + 1])
                
            if ((e[j, i + 1] * e[j, i] < 0)        # sign change in
                                                   # event means zero
                                                   # crossing
                or np.abs(e[j, i + 1]) < TOLERANCE # this point is
                                                   # practically 0
                or np.abs(e[j, i]) < TOLERANCE):

                xLt = X[i+1]       # Last point
                fLt = sol[i+1]
                eLt = e[j, i+1]

                # we need to find a value of x that makes the event zero
                def objective(x):
                    # evaluate ode from xLT to x
                    txspan = [xLt, x]
                    tempsol = odeint(func, fLt, txspan)#, **kwargs)
                    sol = tempsol[-1, :]
                    val, isterminal, direction = event(sol, x)
                    return val

                from scipy.optimize import fsolve
                xZ, = fsolve(objective, xLt, **fsolve_args)  # this should be the
                                              # value of x that makes
                                              # the event zero

                # now evaluate solution at this point, so we can
                # record the function values here.
                txspan = [xLt, xZ]
                tempsol = odeint(func, fLt, txspan, **kwargs)
                fZ = tempsol[-1,:]

                vZ, isterminal, direction = event(fZ, xZ)

                COLLECTEVENT = False
                if direction == 0:
                    COLLECTEVENT = True
                elif (e[j, i + 1] > e[j, i] ) and direction == 1:
                    COLLECTEVENT = True
                elif (e[j, i + 1] < e[j, i] ) and direction == -1:
                    COLLECTEVENT = True
                
                if COLLECTEVENT:
                    TE.append(xZ)
                    YE.append(fZ)
                    IE.append(j)

                    if isterminal:
                        X[i+1] = xZ
                        sol[i+1] = fZ
                        return (np.array(X), 
                                np.array(sol), 
                                np.array(TE), 
                                np.array(YE), 
                                np.array(IE))

    # at the end, return what we have
    return (np.array(X), 
            np.array(sol), 
            np.array(TE), 
            np.array(YE), 
            np.array(IE))


# ----------------------------------------
# Example usage: Multivariate
# ----------------------------------------
# System of differential equations:
# z1' = z3
# z2' = z4
# z3' = -g - b*z3
# z4' = -bz4

def myode(z, t):
    # PROJ: ODE for projectile motion with linear air resistance
    # z = [z0, z1, z2, z3]
    g = 9.81 # Specify approximate graviatational constant
    b = 0.28 # Representative value
    return [z[2], z[3], -g-b*z[2], -b*z[3]]

def linevent(z, t):
    # LINEVENT: Contains the event we are looking for
    # In this event, z(1) = 0 (hitting the ground)
    lookfor = z[0] # Sets this to 0
    stop = True # Stop when event is located
    direction = -1 # specify downward direction
    return lookfor, stop, direction

f0 = [0.18, 0, 5, 10]
xspan = np.linspace(0, 10)
events = [linevent]
func = myode
y0 = f0
Xm, Fm, TEm, YEm, IEm = odelaym(myode, f0, xspan, events)


































