# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 21:47:42 2020

@author: user
"""

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
alpha = pd.read_excel(data, sheet_name='Alpha')
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
capacity_permodule = [float(i) for i in modules]

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

varnames_type3 = ["xd" + str(d + 1) + "p" + str(p + 1) + "s" + str(s + 1) for s in range(N_States) for d in range(N_Demands) for p in range(numPaths_of_d[d])]
var_type3 = list(c.variables.add(lb=[0.0] * N_Paths * N_States,
                                ub=[cplex.infinity] * N_Paths * N_States,
                                types= [c.variables.type.continuous]* N_Paths * N_States,
                                names = varnames_type3
                                ))
        
#counter needed to name the constraints
cnt_link_states = 0
for e in range(N_IPlinks):
    
    # In subset_of_states we save the elements of the subset of states per link - for example: s1,s2
    subset_of_states = subset_states.at[e,'Subset of States S€']
    
    # In newstr we save the numbers of the elements of the subset of states per link - for example: 1,2 
    newstr = ''.join((ch if ch in '0123456789.-e' else ' ') for ch in subset_of_states)
    
    # in listOfNumbers we save the elements of newstr as a list - for example: [1,2]
    listOfNumbers = [int(i) for i in newstr.split()]
    
    #counter needed to name the constraints
    ind_naming = 0
    
    for s in (listOfNumbers):
        
        # names of the constraints
        constnames_type3 = ["c3_" + str(ind_naming + cnt_link_states)]
        
        #counter increased by 1
        ind_naming = ind_naming + 1
        
        # right hand side is a list of zeros since the variables of type one we can bring them to the left hand side
        rhs1 = [0.0]
        
        # In set_of_p_perlink we save the elements of the set of Paths per link - for example: p0,p3,p4,p7,p9,p11
        set_of_p_perlink = set_paths_per_link.at[e,'Paths per link']
        
        # In newstr1 we save the numbers of the elements of the set of Paths per demand - for example: 0,3,4,7,9,11
        newstr1 = ''.join((ch if ch in '0123456789.-e' else ' ') for ch in set_of_p_perlink)
        
        # in listOfNumbers we save the elements of newstr as a list - for example: [0,3,4,7,9,11]
        listOfNumbers1 = [int(i) for i in newstr1.split()]
        
        # indicies get the variables
        # we use the list listOfNumbers1 to loop through the elements of variable of type 3
        ind = [var_type3[k + (s - 1) * N_Paths] for k in listOfNumbers1] + [var_typeInteger[k + e*N_Modules] for k in range(N_Modules)]
        
        # val gives the coefficients of the indicies, length of the list is the same as the length of listOfNumbers1
        # alpha.iloc[e,s] gets elements as a matrix from excel sheet alpha
        result = [element * (- float(alpha.iloc[e,s])) for element in capacity_permodule] 
        
        val = [1.0] * len(listOfNumbers1) + result
        
        rows = [[ind, val]]
        
        # lin_expr is a matrix in list-of-lists format.
        # senses specifies the senses of the linear constraint -L means less-than
        c.linear_constraints.add(lin_expr = rows,
                                 senses="L", 
                                 rhs=rhs1,
                                 names = constnames_type3) 
    
    # the counter is increased by the length of the list listOfNumbers
    cnt_link_states = cnt_link_states + len(listOfNumbers)

b = c.variables.get_num_integer()   
b1 = c.variables.get_num()  
b2 = c.variables.get_names()  
