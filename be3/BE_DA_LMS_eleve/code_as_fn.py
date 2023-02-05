
import numpy as np
from numpy import zeros, eye, exp
from numpy.linalg import \
    (inv,# To invert a matrix
    norm,# To compute the Euclidean norm
    cholesky) # To compute Cholesky factorization
from numpy.random import randn # To generate samples from a normalized Gaussian
import matplotlib.pyplot as plt # To plot a grapH

from numpy.core.numeric import zeros_like
from numpy import\
(round, shape,
copy,
dot,# Matrix-matrix or matrix-vector product
eye,# To generate an identity matrix
ones, # To generate an array full of ones
zeros, # To generate an array full of zeros
linspace)# To get space and time position indices for observations 
from numpy.linalg import \
(inv,# To invert a matrix
norm) # To compute the Euclidean norm
from numpy.random import randn # To generate samples from a normalized Gaussian
import matplotlib.pyplot as plt # To plot a graph
from .VariationalMethods.operators import Hessian3dVar, obs_operator, Rmatrix, Bmatrix
from .VariationalMethods.operators import Precond, Hessian4dVar
from .Model.models import lorenz95 
from .VariationalMethods.solvers import pcg
from numpy.lib.function_base import append
import math

def threedvar_analysis(b_kind='diagonal', sigma_b=0.8, length_scale=10, smoothing=4, total_space_obs=20, sigmaR=0.2, verbose=True, plot=True):
    np.random.seed(seed=42)

    def funcval(x):
        eo = y - obs.hop(x)
        eb = x-xb
        J = eb.dot(B.invdot(eb)) + eo.dot(R.invdot(eo))
        return J
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                          (1) Initialization
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n  = 40           # state space dimension 

    # Model class initialization
    dt = 0.025        # time step for 4th order Runge-Kutta
    F = 8             # Forcing term
    model = lorenz95(F,dt)

    # Observation class initialization
    sigmaR = sigmaR # observation error std
    total_space_obs = total_space_obs # total number of observations at fixed time 
    space_inds_obs = round(linspace(0, n, total_space_obs, endpoint = False)).astype(int) # observation locations in the space
    obs = obs_operator(sigmaR, space_inds_obs, n)
    R = Rmatrix(sigmaR)

    # Background class initialization
    sigmaB = sigma_b # background error std
    if(b_kind == 'diagonal'):
        B = Bmatrix(sigmaB,'diagonal')
    else:
        B = Bmatrix(sigmaB,'diffusion', D=length_scale, M=smoothing)

    # Minimization initialization
    max_outer = 5 # number of maximum outer loops
    max_inner = 20 # number of maximum inner loops
    tol = 1e-6  # tolerance for the inner loop
    tol_grad = 1e-6 # tolerance for the outer loop
    In = eye(n)
    F = Precond(B) # Define the preconditioner

    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                          (2) Generate the truth
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xt = 3.0*ones(n)+randn(n) # the true state
    xt = model.traj(xt,5000)  # spin-up


    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                          (3) Generate the background
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #xb = xt + sigmaB*randn(n)
    #The square root of B can be used to create correlated errors
    xb = xt + B.sqrtdot(randn(n))

    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                          (4) Generate the observations
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y = obs.hop(xt) + sigmaR*randn(len(space_inds_obs)) # TO DO

    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                          (5) Variational data assimilation
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iter_outer = 0 
    xa = copy(xb)   # Choose the initial vector
    if(verbose):
        print('')
        print('iter', '         f(x)', '          ||grad(x)||')
    while iter_outer < max_outer: # Gauss-Newton loop
        A = Hessian3dVar(obs,R,B)    # Define the Hessian matrix (Binv + HtRinvH)
        # TO DO                      # Complete Hessian3dVar in operators.py

        d = y - obs.hop(xa) # obs.misfit(y, xt) # MISFIT
        b = B.invdot(xb - xa) + obs.adj_hop(R.invdot(d))  # Right hand side vector (Binv(xb - xa) + Ht*Rinv*d)
        if verbose :
            print('{:<9d}{:<20.2f}{:<9.2f}'.format(iter_outer, funcval(xa), norm(b)))
        if norm(b) < tol_grad:
            break
        # Calculate the increment dx such that 
        # (Binv + HtRinvH) dx = Binv(xb - x) + Ht Rinv d
        # Solve the linear system by using an iterative solver PCG
        dxs, error, i, flag = pcg(A, zeros_like(xa), b, F, max_inner, tol) 
        xa += dxs[-n:]
        iter_outer += 1

    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                          (5) Diagnostics
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    err_b = norm(xt - xb)/norm(xt)
    err_a = norm(xt - xa)/norm(xt)
    if(verbose):
        print('')
        print('||xt - xb||_2 / ||xt||_2 = ', err_b)
        print('||xt - xa||_2 / ||xt||_2 = ', err_a)

    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                          (6) Plots
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(plot):
        xrange = range(0,n)
        fig, ax = plt.subplots()
        ax.plot(xrange,xt, '-k', label='Truth')
        ax.plot(xrange, xb, '-.b', label='Background')
        ax.plot(space_inds_obs, y, 'og', label='Observations')
        ax.plot(xrange, xa, '-r', label='Analysis')
        leg = ax.legend()
        plt.xlabel('x-coordinate')
        plt.ylabel('Temperature')
        plt.show()

    return {
        'relative_err_xb': err_b,
        'relative_err_xa': err_a,
        'xt': xt,
        'xb': xb,
        'y': y,
        'inds_obs': space_inds_obs,
        'xa': xa,
        'b_kind': b_kind,
        'length_scale': length_scale,
        'smoothing': smoothing,
        'b_sigma': sigma_b,
        'total_space_obs': total_space_obs
    }


def statistical_analysis(sigma0=0.4, locations=(3,), btype='diagonal', sigmaB=0.8, correlation_length=1.0, plot=True):
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                          (1) Initialization
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    np.random.seed(seed=42)
    n  = 4  # state space dimension 

    # Observation operator
    I = eye(n)
    inds_obs = list(locations) # [3] # Location of the observations (1 dimensional grid)
    H = I[inds_obs] # H as a selection operator - H is size m x n
    m = len(inds_obs) # number of observations

    # Observation errors
    # R is a diagonal matrix
    sigmaR = sigma0 # 0.4 # observation error std
    R = zeros((m,m))
    for ii in range(m):
        R[ii,ii] = sigmaR*sigmaR 

    # Background errors
    sigmaB = sigmaB # 0.8 # background error std
    L = correlation_length # 1.0 # correlation length scale
    btype = btype #'diagonal'
    B = zeros((n,n))
    if btype == 'diagonal':
        # no correlation between the grid points
        for ii in range(n):
            B[ii,ii] = sigmaB*sigmaB  
    if btype == 'soar':
        # correlation between the grid points
        for ii in range(n):
            for jj in range(n):
                rij = abs(jj-ii)
                rho = (1 + rij/L)*exp(-rij/L)
                B[ii,jj] = sigmaB*sigmaB*rho          
    # B = B12 * B12^T
    B12 = cholesky(B)
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                          (2) Generate the truth
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xt = randn(n) # the true state

    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                          (3) Generate the background
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xb = xt + np.einsum('ij, j', B12, np.random.randn(n))
    # size of xb is same as xt
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                          (4) Generate the observations
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y = np.einsum('mn, n-> m', H, xt) + np.sqrt(sigmaR) * randn(m) # as long as R is diagonal
    # y should be a vector of size m
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                          (5) Obtain analysis from BLUE
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #Kalman Gain matrix
    BGt = np.matmul(B, H.T) # nn x nm = nm
    K = np.matmul(BGt, inv(np.matmul(H, BGt) + R)) # nm x mm = nm
    #BLUE analysis
    second_term = np.einsum('nm,m->n', K, (y - np.einsum('mn, n -> m', H, xb)))
    xa = xb + second_term

    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                          (5) Diagnostics
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    relative_err_xb = norm(xt - xb)/norm(xt)
    relative_err_xa = norm(xt - xa)/norm(xt)
    if(plot):
        print('')

        print('')
        print('||xt - xb||_2 / ||xt||_2 = ', relative_err_xb)
        print('||xt - xa||_2 / ||xt||_2 = ', relative_err_xa)
        print('\n')

    # Analysis covariance matrix
    Sinv = inv(B) + (H.T).dot(inv(R).dot(H))
    S = inv(Sinv)
    if(plot):
        print('Analysis covariance matrix: \n')
    trace = []
    for ii in range(n):
        if(plot):
            print('S[{}, {}]: {}'.format(ii, ii, S[ii,ii]))
        trace.append(S[ii, ii])

    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                          (6) Plots
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(plot):
        xrange = range(0,n)
        fig, ax = plt.subplots()
        ax.plot(xrange, xt, '+k', label='Truth')
        ax.plot(xrange, xb, 'db', label='Background')
        ax.plot(inds_obs, y, 'og', label='Observations')
        ax.plot(xrange, xa, 'xr', label='Analysis')
        leg = ax.legend()
        plt.xlabel('x-coordinate')
        plt.ylabel('Temperature')
        plt.show()
    return {
        'relative_err_xb': relative_err_xb,
        'relative_err_xa': relative_err_xa,
        'covariance_mat_trace': trace,
        'xt': xt,
        'xb': xb,
        'y': y,
        'xa': xa,
        'inds_obs': inds_obs
    }


def fourdvar_analysis(n=100, nt=10, b_kind='diffusive', D=5, M=4, f_obs_time=30,
                      f_obs_space=10, sigma_obs=1e-2, max_outer=5,
                      max_inner=200, precond='B', verbose=True, plot=True):
    def nonlinear_funcval(x):
        eo = y - obs.gop(x)
        eb = x-xb
        J = eb.dot(B.invdot(eb)) + eo.dot(R.invdot(eo))
        return J

    def quadratic_funcval(x, dx):
        eo = obs.gop(x) - y + obs.tlm_gop(x, dx)
        eb = x-xb+dx
        J = eb.dot(B.invdot(eb)) + eo.dot(R.invdot(eo))
        return J    

    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                          (1) Initialization
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    np.random.seed(42)
    n  = n             # state space dimension 
    nt = nt            # number of time steps

    # Model class initialization
    dt = 0.025        # time step for 4th order Runge-Kutta
    F = 8             # Forcing term
    model = lorenz95(F,dt)

    # Observation class initialization
    sigmaR = sigma_obs # 1e-2# observation error std
    total_space_obs = f_obs_space # total number of observations at fixed time 
    total_time_obs = f_obs_time # 5 # total number of observations at fixed location 
    space_inds_obs = round(linspace(0, n, total_space_obs, endpoint = False)).astype(int) # observation locations in the space
    time_inds_obs = round(linspace(0, nt, total_time_obs, endpoint = False)).astype(int) # observation locations along the time 
    m = total_space_obs*total_time_obs
    obs = obs_operator(sigmaR, space_inds_obs, n, time_inds_obs, nt, model)
    R = Rmatrix(sigmaR)

    # Background class initialization
    sigmaB = 0.8 # background error std
    if(b_kind == 'diagonal'):
        B = Bmatrix(sigmaB,'diagonal')
    else:
        B = Bmatrix(sigmaB,'diffusion', D=D, M=M) # 5, 4

    # Minimization initialization
    max_outer = max_outer # 10 # number of maximum outer loops
    max_inner = max_inner # 500 # number of maximum inner loops
    tol = 1e-6  # tolerance for the inner loop
    tol_grad = 1e-6 # tolerance for the outer loop
    In = eye(n)
    if(precond=='B'):
        F = Precond(B) # Define the preconditioner
    else:
        F = Precond(In)
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                          (2) Generate the truth
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xt = 3.0*ones(n)+randn(n) # the true state
    xt = model.traj(xt,5000)  # spin-up
    #xt = Bmatrix(sigmaB,'diffusion', D=5, M=4).dot(xt)

    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                          (3) Generate the background
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xb = xt + B.sqrtdot(randn(n))
    #The square root of B can be used to create correlated errors
    #xb = xt + B.sqrtdot(randn(n)) 

    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                          (4) Generate the observations
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y = obs.generate_obs(xt)

    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                          (5) Variational data assimilation
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iter_outer = 0 
    dxs = []
    xa = copy(xb)   # Choose the initial vector
    quadcost = math.log10(quadratic_funcval(xa, zeros_like(xa)))
    if(verbose):
        print('')
        print('iter', '  CGiter', '        f(x)', '          ||grad(x)||')
    iter = 0
    total_number_of_inner_iterations = 0
    while iter_outer < max_outer: # Gauss-Newton loop
        A = Hessian4dVar(obs,R,B, xa)    # Define the Hessian matrix (Binv + HtRinvH)
        d = obs.misfit(y, xa) # misfit calculation (y - G(xa))
        b = B.invdot(xb-xa) + obs.adj_gop(xa, R.invdot(d))
        if(verbose):
            print('{:<9d}{:<9d}{:<20.2f}{:<9.2f}'.format(iter_outer, iter, nonlinear_funcval(xa), norm(b)))
        if norm(b) < tol_grad:
            break
        # Calculate the increment dx such that 
        # (Binv + HtRinvH) dx = Binv(xb - x) + Ht Rinv d
        dxs, error, iter, flag = pcg(A, zeros_like(xa), b, F, max_inner, tol)
        total_number_of_inner_iterations += iter
        dx = dxs[-n:]
        for i in range(iter):
            qval = math.log10(quadratic_funcval(xa, dxs[i*n:(i+1)*n]))
            quadcost = append(quadcost, qval)
        xa += dx
        iter_outer += 1

    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                          (5) Diagnostics
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(verbose):
        print('')
        print('||xt - xb||_2 / ||xt||_2 = ', norm(xt - xb)/norm(xt))
        print('||xt - xa||_2 / ||xt||_2 = ', norm(xt - xa)/norm(xt))

    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                          (6) Plots
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(plot):
        xrange = range(0,n)
        fig, ax = plt.subplots()
        ax.plot(xrange,xt, '-k', label='Truth')
        ax.plot(xrange, xb, '-.b', label='Background')
        ax.plot(space_inds_obs, y[:len(space_inds_obs)], 'og', label='Observations')
        ax.plot(xrange, xa, '-r', label='Analysis')
        leg = ax.legend()
        plt.xlabel('x-coordinate')
        plt.ylabel('Temperature')
        #plt.show()

        #plt.figure()
        #plt.plot(quadcost,'r*')
        #plt.show()

    return {
        'xt': xt,
        'xb': xb,
        'y': y,
        'inds_obs': space_inds_obs,
        'xa': xa,
        'total_space_obs': total_space_obs,
        'n ': n,
        'nt ': nt,
        'b_kind ': b_kind,
        'D ': D,
        'M ': M,
        'f_obs_time ': f_obs_time,
        'f_obs_space ': f_obs_space,
        'sigma_obs ': sigma_obs,
        'max_outer ': max_outer,
        'max_inner ': max_inner,
        'total_number_of_inner_iterations': total_number_of_inner_iterations
    }










