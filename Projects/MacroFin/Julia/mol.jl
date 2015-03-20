#=
    Author: Jonathan Payne
    This is simple illustration of the method lines technique for solving
    differential equations.
    cd("/Users/Juddy/Dropbox/Studies/Classes/2014A/Financial Economics/6_Replication/2_Code")
=#

# Runge-Kutta Method
function rk4mono(a, b, h, yinit, fdash)
    # This code calculates ODE using Runge-Kutta 4th order method
    # Solves the problem:
    # dy(x)/dx = fdash(x,y)
    x = a:h:b
    y = zeros(1, length(x))
    y[1] = yinit;
    for i = 1:(length(x)-1)
        k_1 = fdash(x[i], y[i]);
        k_2 = fdash(x[i] + 0.5*h, y[i] + 0.5*h*k_1)
        k_3 = fdash(x[i] + 0.5*h, y[i] + 0.5*h*k_2)
        k_4 = fdash(x[i] + h, y[i] + k_3*h)
        y[i+1] = y[i] + (1/6)*(k_1 + 2*k_2 + 2*k_3 + k_4)*h
    end
    return y
end 
a, b, h, yinit = 0, 10, 0.1, 48
fdash(x,y) = 9.8 - 0.196*y
ymono = rk4mono(a, b, h, yinit, fdash)


# Multi-dimensional Runge-Kutta Method
function rk4old(y, dydx, n, x, h, yout, derivs)
    # Given values for n variables y[1,...,n] and their derviatives
    # dydx[1,...,n] known at x, this algorithm uses the fourth-order
    # Runge-Kutta method to advance the solution over an interval h
    # and return the incremented variables as yout[1,...n], which need
    # not be a distinct array from y. The user supplies the routine
    # derivs(x, y, dydx) which returns derivatives dydx at x

    dym = zeros(n)
    dyt = zeros(n)
    yt = zeros(n)
    hh = h*0.5
    h6 = h/6.0
    xh = x+hm
    
end

function rk4(fdash, t0, x0, tf, N)
    # Approximate the solution of the initial value problem:
    # x'(t) = fdash(t, x),
    # x(t0) = x0
    # using the classical fourth order Runge-Kutta method.
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
    neqn = length(x0)
    ti = linspace(t0, tf, N+1)
    wi = zeros(neqn, N+1)
    wi[1:neqn, 1] = x0
    h = (tf - t0) / N
    for i = 1:N
        k1 = h*fdash(t0, x0)
        k2 = h*fdash(t0 + h/2, x0 + k1/2)
        k3 = h*fdash(t0 + h/2, x0 + k2/2)
        k4 = h*fdash(t0 + h, x0 + k3)
        x0 = x0 + (k1 + 2*k2 + 2k3 + k4) / 6
        t0 = t0 + h
        wi[1:neqn, i+1] = x0
    end
    return wi, ti
end
t0, x0, tf, N = 0, 48, 10, 100
fdash(t, x) = 9.8 - 0.196*x
wi, ti = rk4(fdash, t0, x0, tf, N)

t0, x0, tf, N = 0, [0,-4], 10, 100
A = [[1 2], [3 2]]
fdash(t, x) = A*x
xtrue(t) = [ 8/5*exp(-t) - 8/5*exp(4*t), -8/5*exp(-t) - 12/5*exp(4*t) ]
wi, ti = rk4(fdash, t0, x0, tf, N)

function rk4step(fdash, t, x, h)
    # Calculates the the nth step in the classical fourth order Runga-Kutta
    # procedure for solving a differential equation of the form:
    #   x'(t) = fdash(t, x),
    # Inputs:
    #   * fdash: user defined function
    #   * t: previous value of independent variable
    #   * x: previous value of dependent variable
    #   * h: increment
    # Outputs:
    #   * xnew: next step
    #   * tnew: new time
    k1 = h*fdash(t, x)
    k2 = h*fdash(t + h/2, x + k1/2)
    k3 = h*fdash(t + h/2, x + k2/2)
    k4 = h*fdash(t + h, x + k3)
    xnew = x + (k1 + 2*k2 + 2k3 + k4) / 6
    tnew = t + h
    return xnew, tnew
end

function rk4uniform(fdash, t0, x0, tf, N)
    # Approximate the solution of the initial value problem using the
    # classical fourth order Runge-Kutta method with uniform step sizes:
    #   * x'(t) = fdash(t, x),
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

# Usage example 1:
t0, x0, tf, N = 0, [0,-4], 10, 100
A = [[1 2], [3 2]]
fdash(t, x) = A*x
xtrue(t) = [ 8/5*exp(-t) - 8/5*exp(4*t), -8/5*exp(-t) - 12/5*exp(4*t) ]
wi2, ti2 = rk4uniform(fdash, t0, x0, tf, N)

# Usage example 2:
function Fdash(t, x)
    xp = zeros(2,1)
    xp[1] = x[2]
    xp[2] = -t*x[1] + exp(t)*x[2] + 3*sin(2*t)
    return xp
end
t0, x0, tf, N = 0, [2,8], 4, 100
xi5, ti5 = rk4uniform(Fdash, t0, x0, tf, N)



# ----------------------------------------
# METHOD OF LINES
# ----------------------------------------
# Solves the heat (or diffusion equation) using the method of lines technique:
#   * u_t = kappa u_xx
#   * u(x, 0) = n(x)
#   * u(0, t) = g_0(t), t>0
#   * u(1, t) = g_1(t), t>0

xs, xf = 0, 5
M = 11 # number of mesh points
h = (xf - xs)/M # size of the space mesh
k = 0.2 # size of the time mesh
xgrid = linspace(xs, xf, M)
uxs(t) = 1
uxf(t) = 1
ut0(x) = x

# Set up the discretisation in the space dimension
# Need to discretise into the form:
#   * U'(t) = A U(t) + g(t)
A = (1/h^2)*(diagm(ones(M-1),-1)+diagm(-2*ones(M))+diagm(ones(M-1),1))
function gfunc(t)
    gvec = zeros(M);
    gvec[1] = uxs(t)
    gvec[end] = uxf(t)
    return gvec
end

# Run the ode solver
fdash(t, Uall) = A*Uall + gfunc(t)
t0, Unit, tf, N = 0, ut0(xgrid), 10, 100
wi2, ti2 = rk4uniform(fdash, t0, Unit, tf, N)


########################################################################
#                                                                      %
#    Example of ADI Method for 2D heat equation                        %
#                                                                      %
#          u_t = u_{xx} + u_{yy} + f(x,t)                              %
#                                                                      %
#    Test problme:                                                     %
#      Exact solution: u(t,x,y) = exp(-t) sin(pi*x) sin(pi*y)          %
#      Source term:  f(t,x,y) = exp(-t) sin(pi*x) sin(pi*y) (2pi^2-1)  %
#                                                                      %
#    Files needed for the test:                                        %
#                                                                      %
#     adi.m:      This file, the main calling code.                    %
#     f.m:        The file defines the f(t,x,y)                        %
#     uexact.m:    The exact solution.                                 %
#                                                                      %
#     Results:         n              e            ratio               %
#                     10           0.0041                              %
#     t_final=0.5     20           0.0010           4.1                %
#                     40           2.5192e-04       3.97               %       
#                     80           6.3069e-05       3.9944             %
########################################################################

uexact(t,x,y) = exp(-t)*sin(pi*x)*sin(pi*y)
f(t,x,y) = exp(-t)*sin(pi*x)*sin(pi*y)*(2*pi^2-1)

a = 0; b = 1; c = 0; d = 1; n = 40; tfinal = 0.5
m = n
h = (b-a)/n; dt = h; h1 = h*h
x = a:h:b; y = c:h:d

# ----- Initial condition
t = 0
u1 = zeros(m+1, m+1)
for i=1:m+1
    for j=1:m+1
        u1[i,j] = uexact(t, x[i], y[j])
    end
end

# ----- Loop for time t
k_t = floor(tfinal/dt)
for k = 1:k_t
    t1 = t + dt
    t2 = t + dt/2
    # A. Sweep in the x-direction
    # -- Set the boundary condition
    u2 = zeros(m+1, n+1)
    for i = 1:m+1
        u2[i, 1] = uexact(t2, x[i], y[1])
        u2[i, n+1] = uexact(t2, x[i], y[n+1])
        u2[1, i] = uexact(t2, x[1], y[i])
        u2[m+1, i] = uexact(t2, x[m+1], y[i])
    end
    # -- Iterate
    for j = 2:n
        A = spzeros(m-1, m-1)
        b = zeros(m-1,1)
        for i = 2:m
            b[i-1] = ( (u1[i,j-1] - 2*u1[i,j] + u1[i,j+1])/h1
                      + f(t2, x[i], y[j]) + 2*u1[i,j]/dt )
            if i == 2
                b[i-1] = b[i-1] + uexact(t2, x[i-1], y[j])/h1
                A[i-1, i] = -1/h1
            else
                if i == m
                    b[i-1] = b[i-1] + uexact(t2, x[i+1], y[j])/h1
                    A[i-1, i-2] = -1/h1
                else
                    A[i-1, i] = -1/h1
                    A[i-1, i-2] = -1/h1
                end
            end
            A[i-1, i-1] = 2/dt + 2/h1
        end
        ut = A\b
        for i = 1:m-1
            u1[i+1, j] = ut[i]
        end
    end
    # B. Sweep in the y-direction
    # -- Set the boundary condition
    u1 = zeros(m+1, n+1)
    for i = 1:m+1
        u1[i,1] = uexact(t1, x[i], y[1])
        u1[i,n+1] = uexact(t1, x[i], y[m+1])
        u1[1,i] = uexact(t1, x[1], y[i])
        u1[m+1,i] = uexact(t1, x[m+1], y[i])
    end
    # -- Iterate
    for i = 2:m
        A = spzeros(m-1, m-1)
        b = zeros(m-1,1)
        for j = 2:n
            b[j-1] = ( (u2[i-1,j] - 2*u2[i,j] + u2[i+1,j])/h1
                      + f(t2, x[i], y[j]) + 2*u2[i,j]/dt )
            if j == 2
                b[j-1] = b[j-1] + uexact(t1, x[i], y[j-1])/h1
                A[j-1,j] = -1/h1
            else
                if j == n
                    b[j-1] = b[j-1] + uexact(t1, x[i], y[j+1])/h1
                    A[j-1,j-2] = -1/h1
                else
                    A[j-1,j] = -1/h1
                    A[j-1,j-2] = -1/h1
                end
            end
            A[j-1,j-1] = 2/dt + 2/h1
        end
        # -- Solve the system
        ut = A\b
        for j = 1:n-1
            u1[i,j+1] = ut[j]
        end
    end
t = t + dt
end

# ----- Data analysis
ue = zeros(m+1, m+1)
for i = 1:m+1
    for j = 1:n+1
        ue[i,j] = uexact(tfinal, x[i], y[i])
    end
end
err = maximum(abs(u1-ue))
