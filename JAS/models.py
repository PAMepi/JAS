
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Aug 20 2020
@author: Juliane Oliveira julianlanzin@gmail.com 
@author: Moreno rodrigues rodriguesmsb@gmail.com

"""


from scipy.integrate import solve_ivp, odeint
from scipy import optimize
from scipy.integrate import odeint
from scipy.optimize import least_squares
import numpy as np
import warnings
warnings.filterwarnings('ignore')




class compartimental_models:
    def __init__(self):
        pass

    def r0(self):
        return self.beta/self.gamma

#     #define least square errors
#     def sir_least_squares_error_ode(self,par, time_exp, f_exp, fitting_model, initial_conditions):
#         args = par
#         time = (time_exp.min(), time_exp.max())
    
#         y_model = fitting_model(initial_conditions, time, *args)
#         simulated_time = y_model.t
#         simulated_ode_solution = y_model.y
#         _, simulated_qoi, _ = simulated_ode_solution
    
#         residual = f_exp - simulated_qoi

#         return np.sum(residual ** 2.0)



# class SIR(compartimental_models):
#     def __init__(self, pop):
#         self.pop = pop


#         #Define initial conditions for compartiments
#         self.__S0 = self.pop/self.pop
#         self.__I0 = 1/self.pop
#         self.__R0 = 0.
#         self.__y0 = [self.__S0, self.__I0, self.__R0]

  
#     def __sir_model(self, time, init, beta, gamma):
    
#         S, I, R = init

#         S_prime = -beta * S * I
#         I_prime = beta * S * I - gamma * I
#         R_prime = gamma * I 
    
#         return S_prime, I_prime, R_prime


#     #Define ode solver using ivp
#     def sir_ode_solver(self, y0, time, beta, gamma):
 
#         #Create a time span
#         time_span = [min(time), max(time)]

#         #Create a list with points to evaluate the function 
#         time_eval = np.arange(time_span[0],time_span[1] + 1, 1) 
    
#         solution_ode = solve_ivp(
#             fun = lambda time, init: self.__sir_model(time, init, beta = beta, gamma = gamma),
#             t_span = time_span,
#             y0 = y0,
#             t_eval = time_eval,
#             method = "LSODA")

#         return solution_ode
    
    
#     def fit(self, x, y, workers = -1, bounds = [(0,2), (1/21,1/4.99)]):
#         x = x.values.astype(np.float)
#         y = y.values.astype(np.float)/self.pop
#         result_sir = optimize.differential_evolution(
#             super().sir_least_squares_error_ode,
#             bounds = bounds, #bounds are given in the following order: beta, gamma
#             args = (x, y, self.sir_ode_solver, self.__y0),
#             popsize = 300,
#             strategy = "best1bin",
#             tol = 1e-5,
#             recombination = 0.5,
#             maxiter = 100,
#             disp = False,
#             seed = 12,
#             workers = workers)
#         self.beta, self.gamma = result_sir.x[0], result_sir.x[1]
    
#     def predict(self, time):
#         time = time
#         solution_ode_sir = self.sir_ode_solver(self.__y0, time, self.beta, self.gamma)
#         self.S, self.I, self.R = solution_ode_sir.y
#         return self.I * self.pop


class sir(compartimental_models):
    def __init__(self, pop):
        self.pop = pop
    
    #Defining the Model
    def __sir(self, f ,t, parametros):
    
        #parameters
        b, gamma = parametros
    
        #variables
        S, I, R, Nw = f
        #S = f[0]
        #I = f[1]
        #R = f[2]
        #Nw = f[3]
    

        #equations
        dNw_dt = b * S * I
        dS_dt = - b * S * I
        dI_dt = b * S * I - gamma * I
        dR_dt = gamma * I 
    
        #Returning the derivatives
        return [dS_dt, dI_dt, dR_dt, dNw_dt]

    
    #Calculating the Error
    def __least_square_error(self,pars, ts0, y):

        y = self.y
    
        #assaing the Least square given parameters to the numerical integration for the error calculation
        b, gamma, i0 = pars
    
        #initial conditions
        q0 = [1-i0,i0,0,i0]
    
        #parameters to feed the E.D.Os
        parode = b, gamma
    
        #Integrating
        qs = odeint(self.__sir, q0, ts0, args = (parode,),mxstep = 1000000)

        #sinf the epidemic curve
        #sdth the death curve

        sinf = qs[:,-1]

        #define the standardized residuals
        erri = (self.pop * sinf - y) / np.sqrt(self.pop * sinf + 1.0)
    
        return np.r_[erri]

    
    def fit(self, x, y, n_tries, bounds = [[0., 2.0], [1/14, 1/5]]):
        self.x = x
        self.y = y
        #
        ts0 = np.arange(1, len(self.x) + 1)
    
        #Convert bounds to array
        bounds.append([0,50/self.pop])
        bounds = np.array(bounds)
      
       

        #Add infinity to try to avoid local
        self.best_res = None
        best_cost = np.inf
        for i in range(n_tries):
        
            #create a set of ramdom parameters
            par0 = np.random.rand(len(bounds[0]))
        
            #Limit those parameters to the interval defined
            par0 = bounds[:,0] + par0 * (bounds[:,1] - bounds[:,0])
        
            res = optimize.least_squares(lambda pars: self.__least_square_error(pars, ts0, self.y), par0, bounds = (bounds[:,0],bounds[:,1]))
        
            if res.cost < best_cost:
                best_cost = res.cost
                self.best_res = res
        
        self.beta, self.gamma, self.i0 = self.best_res.x
    

    def predict(self,time):
        
        q0 = [1 - self.i0, self.i0,0, self.i0]

        #pick best parameters
        parode = self.beta, self.gamma
        predicted = odeint(self.__sir, q0, np.arange(1, len(time) + 1), args = (parode,), mxstep = 1000000)

        self.S = predicted[:,0]
        self.I = predicted[:,1]
        self.R = predicted[:,2]
        self.Nw = predicted[:,3]

        #Compute cases
        return self.Nw * self.pop

class sir_bv(compartimental_models):
    def __init__(self, pop):
        self.pop = pop
    
    #Defining the Steap Function 

    def __H(self,t):
        h = 1.0/(1.0+ np.exp(-2.0 * 50 * t))
        return h


    #Defining the Beta(t)
    def __beta(self, t, t1, b, b1):
        beta = b * self.__H(t1 - t) + b1 * self.__H(t - t1) 
        return beta
    

    #Defining the Model
    def __sir(self, f, t, parametros):
    
        #parameters
        b, b1, gamma, t1 = parametros
        
        
        #variables
        S, I, R, Nw = f
    
        #equations
        dNw_dt = self.__beta(t, t1, b, b1) * S * I
        
        dS_dt = - self.__beta(t, t1, b, b1) * S * I
        dI_dt = self.__beta(t, t1, b, b1) * S * I - gamma * I
        dR_dt = gamma * I 
    
        #Returning the derivatives
        return [dS_dt, dI_dt, dR_dt, dNw_dt]

    #Calculating the Error
    def __least_square_error_bv(self, pars, ts0, y):

        y = self.y
    
        #assaing the Least square given parameters to the numerical integration for the error calculation
        b, b1, gamma, t1, i0 = pars
    
        #initial conditions
        q0 = [1-i0, i0, 0, i0]
    
        #parameters to feed the E.D.Os
        parode = b, b1, gamma, t1
    
        #Integrating
        qs = odeint(self.__sir, q0, ts0, args = (parode,),mxstep=1000000)
        

        #sinf the epidemic curve
        #sdth the death curve

        sinf = qs[:,-1]

        #define the standardized residuals
        erri = (self.pop * sinf - y) / np.sqrt(self.pop * sinf + 1.0)
    
        return np.r_[erri]
    
    def fit(self, x, y, n_tries, bounds = [[0., 2.0], [0., 2.0], [1/14, 1/5], [15,45]]):
        self.x = x
        self.y = y
        #
        ts0 = np.arange(1, len(self.x) + 1)
    
        #Convert bounds to array
        bounds.append([0,50/self.pop])
        bounds = np.array(bounds)
        
      
    
        #Add infinity to try to avoid local
        self.best_res = None
        best_cost = np.inf
        for i in range(n_tries):
        
            #create a set of ramdom parameters
            par0 = np.random.rand(len(bounds))
        
            #Limit those parameters to the interval defined
            par0 = bounds[:,0] + par0 * (bounds[:,1] - bounds[:,0])
            
        
            res = optimize.least_squares(lambda pars: self.__least_square_error_bv(pars, ts0, self.y), par0, bounds = (bounds[:,0],bounds[:,1]))
        
            if res.cost < best_cost:
                best_cost = res.cost
                self.best_res = res
        
        self.beta, self.beta1, self.gamma, self.t1, self.i0 = self.best_res.x
        #note to transform t into the correct date use int(t1 - 1)
    
    def predict(self, time):
        
        q0 = [1 - self.i0, self.i0,0, self.i0]

        #pick best parameters
        parode = self.beta, self.beta1, self.gamma, self.t1
        predicted = odeint(self.__sir, q0, np.arange(1, len(time) + 1), args = (parode,), mxstep = 1000000)

        self.S = predicted[:,0]
        self.I = predicted[:,1]
        self.R = predicted[:,2]
        self.Nw = predicted[:,3]

        #Compute cases
        return self.Nw * self.pop

