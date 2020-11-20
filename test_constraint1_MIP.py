# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 20:26:46 2020

@author: user
"""

import pandas as pd 
import cplex

c = cplex.Cplex()
data = pd.ExcelFile(r'C:\Users\user\Desktop\Research_Project_Shkurte_Esati\Code\FT_problem_modular_capacities\example_inputdata_modular_capacities.xlsx')


set_IPLinks = pd.read_excel(data, sheet_name='IPlinks' , usecols= ['IP links E'])
set_modules = pd.read_excel(data, sheet_name='Modules', usecols= ['capacity per module uk'])
costs_per_module = pd.read_excel(data, sheet_name='CostsPerModule', usecols=['Set of costs C of modules K'])
set_IPLinks = pd.read_excel(data, sheet_name='IPlinks' , usecols= ['IP links E'])
set_Demands = pd.read_excel(data, sheet_name='Demands' , usecols= ['D Demands'])
set_TrafficVolume = pd.read_excel(data, sheet_name='TrafficVolume' , usecols= ['h(d)'])
set_States = pd.read_excel(data, sheet_name='States', usecols=['States S'])    
set_Paths = pd.read_excel(data, sheet_name='Paths', usecols=['Paths P(d)'])
set_paths_per_demand = pd.read_excel(data, sheet_name='Paths', usecols=['paths per demand'])
set_Paths_Total = pd.read_excel(data, sheet_name = "Paths")
subset_states = pd.read_excel(data, sheet_name='subset states', usecols= ['Subset of States S€'])
set_paths_per_link = pd.read_excel(data, sheet_name='Paths R(d,e)', usecols= ['Paths per link'])

set_modules = pd.read_excel(data, sheet_name='Modules', usecols= ['capacity per module uk'])
costs_per_module = pd.read_excel(data, sheet_name='CostsPerModule', usecols=['Set of costs C of modules K'])
#modules = pd.read_excel(data, sheet_name='Modules')
# Calc column length
N_IPlinks = len(set_IPLinks)
N_Demands = len(set_Demands)
N_States = len(set_States)
N_Paths = len(set_Paths)

N_IPlinks = len(set_IPLinks)
N_Modules = len(set_modules)


cost = costs_per_module["Set of costs C of modules K"].tolist()
modules = set_modules["capacity per module uk"].tolist()
modules = [float(i) for i in modules]

# extract column from pandas dataframe 
demandList = set_Paths_Total["D Demands"].tolist()


# loop through all entries
counter = 0
numPaths_of_d = []
loop_counter = 0

for d in demandList:
    
        loop_counter = loop_counter + 1;  
        if isinstance(d,str) == 0:
                counter = counter + 1
                if (loop_counter == len(demandList)): 
                        numPaths_of_d.append(counter)
            
            
        else:
                if loop_counter == 1:
                        counter = 1
                    
                else:
                    numPaths_of_d.append(counter)
                    counter = 1


varnames_typeInteger = ["ae" + str(e + 1) + "k" + str(k + 1) for e in range(N_IPlinks) for k in range(N_Modules)]

obj1 = cost * N_IPlinks
lb1 = [0.0] * N_IPlinks * N_Modules

var_typeInteger = list(c.variables.add(obj = obj1, 
                                lb= lb1,
                                ub=[cplex.infinity] * N_IPlinks * N_Modules,
                                types= [c.variables.type.integer]* N_IPlinks * N_Modules,
                                names = varnames_typeInteger
                                ))

varnames_type2 = ["xd" + str(d + 1) + "p" + str(p + 1)  for d in range(N_Demands) for p in range(numPaths_of_d[d])]    
var_type2 = list(c.variables.add(lb=[0.0] * N_Paths,
                            ub=[cplex.infinity] * N_Paths,
                            types= [c.variables.type.continuous] * N_Paths,
                            names = varnames_type2
                            ))


for e in range(N_IPlinks):
        
        # names of the constraints
        constnames_type1 = ["c1_" + str(e)]
        
        # right hand side is a list of zeros since the variables of type one we can bring them to the left hand side
        rhs1 = [0.0]
        
        # In set_of_p_perlink we save the elements of the set of Paths per link - for example: p0,p3,p4,p7,p9,p11
        set_of_p_perlink = set_paths_per_link.at[e,'Paths per link']
        
        # In newstr we save the numbers of the elements of the set of Paths per link - for example: 0,3,4,7,9,11
        newstr = ''.join((ch if ch in '0123456789.-e' else ' ') for ch in set_of_p_perlink)
        
        # in listOfNumbers we save the elements of newstr as a list - for example: [0,3,4,7,9,11]
        listOfNumbers = [int(i) for i in newstr.split()]
        
        #ind = [var_type3[p + N_Paths * s]] + [var_type2[p]] + var_AFT[p]
                
                
        
        # indicies get the variables 
        # we use the list listOfNumbers to loop through the elements of variable of type 2
        ind = [var_type2[p] for p in listOfNumbers] + [var_typeInteger[k + e*N_Modules] for k in range(N_Modules)]
        
        modules1 = [element * (-1) for element in modules]
        
        # val gives the coefficients of the indicies, length of the list is the same as the length of listOfNumbers
        val = [1.0] * len(listOfNumbers) + modules1
    
        
        rows = [[ind, val]]
        
        # lin_expr is a matrix in list-of-lists format.
        # senses specifies the senses of the linear constraint -L means less-than
        c.linear_constraints.add(lin_expr = rows,
                                 senses="L", 
                                 rhs=rhs1,
                                 names = constnames_type1)

b = c.variables.get_num_integer()   
b1 = c.variables.get_num()  
b2 = c.variables.get_names()  