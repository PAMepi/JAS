
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Aug 11 2020
@author: Moreno rodrigues rodriguesmsb@gmail.com

"""


from scipy.integrate import solve_ivp, odeint
from scipy import optimize
import numpy as np
import warnings
warnings.filterwarnings('ignore')




class compartimental_models:
    def __init__(self):
        pass

    #define least square errors
    def sir_least_squares_error_ode(self,par, time_exp, f_exp, fitting_model, initial_conditions):
        args = par
        time = (time_exp.min(), time_exp.max())
    
        y_model = fitting_model(initial_conditions, time, *args)
        simulated_time = y_model.t
        simulated_ode_solution = y_model.y
        _, simulated_qoi, _ = simulated_ode_solution
    
        residual = f_exp - simulated_qoi

        return np.sum(residual ** 2.0)



class SIR(compartimental_models):
    def __init__(self, pop):
        self.pop = pop


        #Define initial conditions for compartiments
        self.__S0 = self.pop/self.pop
        self.__I0 = 1/self.pop
        self.__R0 = 0.
        self.__y0 = [self.__S0, self.__I0, self.__R0]

  
    def __sir_model(self, time, init, beta, gamma):
    
        S, I, R = init

        S_prime = -beta * S * I
        I_prime = beta * S * I - gamma * I
        R_prime = gamma * I 
    
        return S_prime, I_prime, R_prime


    #Define ode solver using ivp
    def sir_ode_solver(self, y0, time, beta, gamma):
 
        #Create a time span
        time_span = [min(time), max(time)]

        #Create a list with points to evaluate the function 
        time_eval = np.arange(time_span[0],time_span[1] + 1, 1) 
    
        solution_ode = solve_ivp(
            fun = lambda time, init: self.__sir_model(time, init, beta = beta, gamma = gamma),
            t_span = time_span,
            y0 = y0,
            t_eval = time_eval,
            method = "LSODA")

        return solution_ode
    
    
    def fit(self, x, y, workers = -1, bounds = [(0,2), (1/21,1/4.99)]):
        x = x.values.astype(np.float)
        y = y.values.astype(np.float)/self.pop
        result_sir = optimize.differential_evolution(
            super().sir_least_squares_error_ode,
            bounds = bounds, #bounds are given in the following order: beta, gamma
            args = (x, y, self.sir_ode_solver, self.__y0),
            popsize = 300,
            strategy = "best1bin",
            tol = 1e-5,
            recombination = 0.5,
            maxiter = 100,
            disp = False,
            seed = 12,
            workers = workers)
        self.beta, self.gamma = result_sir.x[0], result_sir.x[1]
    
    def predict(self, time):
        time = time
        solution_ode_sir = self.sir_ode_solver(self.__y0, time, self.beta, self.gamma)
        self.S, self.I, self.R = solution_ode_sir.y
        return self.I * self.pop

    def RO(self):
        return self.beta/self.gamma

    
        


   

