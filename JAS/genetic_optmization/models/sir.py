
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Aug 20 2020
@author: Juliane Oliveira julianlanzin@gmail.com 
@author: Moreno rodrigues rodriguesmsb@gmail.com

"""


from scipy.integrate import solve_ivp
from scipy import optimize
from scipy.optimize import least_squares
import numpy as np
import warnings
warnings.filterwarnings('ignore')


class sir:
    def __init__(self, pop):
        self.pop = pop
    
    def get_parmeters(self):
        try:
            return self.beta, self.gamma, self.i0
        except:
            print("Model must be fitted to get parameters")
    
    def __sir_model(self, t, X, beta, gamma):
    
        S, I, R, Nw = X 
    
        Nw_prime = beta * S * I
    
        S_prime = -beta * S * I
    
        I_prime = beta * S * I - gamma * I
    
        R_prime = gamma * I 
    
        return S_prime, I_prime, R_prime, Nw_prime
    
    #Define ode solver using ivp
    def sir_ode_solver(self, time, beta, gamma, i0):
    
        y0 = 1-i0, i0, 0, i0
    
        time_span = [min(time), max(time)]
        time_eval = np.arange(time_span[0],time_span[1] + 1, 1)
    
        solution_ode = solve_ivp(
            fun = lambda t, y: self.__sir_model(t = t, X = y, beta = beta, gamma = gamma),
            t_span = time_span,
            y0 = y0,
            t_eval = time_eval,
            method = "LSODA")
        return solution_ode
    
    def sir_least_squares_error_ode(self, par, time_exp, f_exp, fitting_model):
        args = par
        time = (time_exp.min(), time_exp.max())
    
        y_model = fitting_model(time, *args)
        simulated_time = y_model.t
        simulated_ode_solution = y_model.y
        sf, i_f, rf, nwf = simulated_ode_solution
    
        residual = (f_exp - (nwf * self.pop))**2

        return np.sum(residual/np.sqrt(self.pop * nwf + 1))
    
    ## Run sir
    def fit(self, x, y, ncores = -1, maxiter = 100, seed = 12, bounds = [(0,2), (1/21,1/4.99)]):
        bounds.append((1,1/self.pop))
        self.x = x.values.astype(np.float)
        y = y.values.astype(np.float)/self.pop
        result_sir = optimize.differential_evolution(
            self.sir_least_squares_error_ode,
            bounds = bounds,
            args = (x, y, self.sir_ode_solver),
            popsize = 300,
            strategy = "best1bin",
            tol = 1e-5,
            recombination = 0.5,
            maxiter = maxiter,
            disp = False,
            seed = seed,
            workers = ncores)
        self.beta, self.gamma, self.i0 = result_sir.x[0], result_sir.x[1], result_sir.x[2]
    
    def predict(self, time):
        time_pred = np.arange(min(self.x), max(self.x) + 1 + time, 1)

        solution_ode = sir_ode_solver(time_pred, self.beta, self.gamma, self.i0)

        self.S, self.I, self.R, self.Nw = solution_ode.y
        return self.Nw * pop
    
    

        

    
