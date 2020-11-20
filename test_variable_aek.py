# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 20:00:25 2020

@author: user
"""

import pandas as pd 
import cplex

c = cplex.Cplex()
data = pd.ExcelFile(r'C:\Users\user\Desktop\Research_Project_Shkurte_Esati\Code\FT_problem_modular_capacities\example_inputdata_modular_capacities.xlsx')


set_IPLinks = pd.read_excel(data, sheet_name='IPlinks' , usecols= ['IP links E'])
set_modules = pd.read_excel(data, sheet_name='Modules', usecols= ['capacity per module uk'])
costs_per_module = pd.read_excel(data, sheet_name='CostsPerModule', usecols=['Set of costs C of modules K'])
# Calc column length
N_IPlinks = len(set_IPLinks)
N_Modules = len(set_modules)


cost = costs_per_module["Set of costs C of modules K"].tolist()

#cost = [float(i) for i in cost]

varnames_typeInteger = ["ae" + str(e + 1) + "k" + str(k + 1) for e in range(N_IPlinks) for k in range(N_Modules)]

obj1 = cost * N_IPlinks
lb1 = [0.0] * N_IPlinks * N_Modules

var_typeInteger = list(c.variables.add(obj = obj1, 
                                lb= lb1,
                                ub=[cplex.infinity] * N_IPlinks * N_Modules,
                                types= [c.variables.type.integer]* N_IPlinks * N_Modules,
                                names = varnames_typeInteger
                                ))