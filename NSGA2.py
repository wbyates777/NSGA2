# NSGA2 07/06/2021
# Program Name: NSGA2.py
# Description: This is a python 3.7 implementation of the NSGA-II algorithm 
# Author: W.B. Yates - University of Exeter

# Employs integer solution vectors and assumes positive, real valued, objective functions to be minimised

# K. Deb, A. Pratap, S. Agarwal and T. Meyarivan, "A fast and elitist multiobjective genetic algorithm: NSGA-II," 
# in IEEE Transactions on Evolutionary Computation, 
# vol. 6, no. 2, pp. 182-197, April 2002, 
# doi: 10.1109/4235.996017

# Hypervolume code see https://ls11-www.cs.tu-dortmund.de/rudolph/hypervolume/start
# Note for Python 3 compatability you must replace xrange with range, and 
# on line 166, replace 'decorated.sort()' with 'sorted(decorated, key=lambda n: n[0])' 

import math
import random
import matplotlib.pyplot as plt
import json
import csv
import sys
from  hv3 import HyperVolume
from typing import NamedTuple


# 
# File reading/writting
# 

def readJSONFile(filename):
    with open(filename, 'r') as myfile:
        return json.loads(myfile.read())


def saveJSONFile(filename, json_str):
    with open(filename,'w') as myfile:
        json.dump(json_str, myfile)

# 
# NSGA2 configuration and output
# 

SolutionSpec = NamedTuple('SolutionSpec', [
    ("sol_num", int),           # length of this part of solution vector
    ("min_val", int),           # minimum value for element of this part of solution vector 
    ("max_val", int),           # maximum value for element of this part of solution vector 
])

Config = NamedTuple('Config', [
    ("seed", int),              # random number seed
    ("pop_size", int),          # the number of solutions to work with; typically less than 100
    ("gen_num", int),           # the number of generations or iterations
    ("sol_spec", list),         # specify the solution vector [(sol_spec1), (sol_spec2), ...]
    ("mutate_prob", float),     # probability of mutating element of solution vector
    ("obj_func", object),       # the multi-objective function to be minimised
    ("obj_num", int),           # the number of objectives to be minimised
    ("ref_point", list),        # reference point for hypervolume calculation
    ("trace", int)              # trace level - 0 implies 'silent'
])


Results = NamedTuple('Results', [
    ("solutions", list),        # the solution set
    ("obj_values", list),       # their objective function values 
    ("pareto_fronts", list),    # the solution indexes sorted into Pareto fronts (rank 0,1,2...) 
])

class NSGA2:
  
    # constructor
    def __init__(self, config:Config):
        """
        A python 3.7 implementation of NSGA2 by W.B. Yates - University of Exeter 
        Assumes positive objective values to be minimised
        K. Deb, A. Pratap, S. Agarwal and T. Meyarivan, "A fast and elitist multiobjective genetic algorithm: NSGA-II," 
        in IEEE Transactions on Evolutionary Computation, 
        vol. 6, no. 2, pp. 182-197, April 2002, 
        doi: 10.1109/4235.996017
        """
        
        random.seed(config.seed)
        
        self.pop_size    = config.pop_size
        self.gen_num     = config.gen_num
        
        self.sol_nums   = []
        self.min_vals   = []
        self.max_vals   = []
        self.sol_length = 0
        for a in config.sol_spec:
            self.sol_nums.append(a.sol_num)
            self.min_vals.extend([a.min_val for i in range(a.sol_num)])
            self.max_vals.extend([a.max_val for i in range(a.sol_num)])
            self.sol_length = self.sol_length + a.sol_num
        
        self.mutate_prob = config.mutate_prob
        self.trace       = config.trace
        self.evalFunc    = config.obj_func
        self.obj_num     = config.obj_num
        self.hv          = HyperVolume(config.ref_point);
        
        if self.trace > 0:
            print('\nConstructed population of ', self.pop_size, ' solutions of length ', self.sol_length) 
            print('Running for ', self.gen_num, ' generations')


    # begin output methods
    def print_solution(self, idx, solutions, obj_values):
        print("%d)" % idx, end = " ")
        for i in range(self.sol_length):
            print(solutions[idx][i], end = " ")
        print("[", end = "")
        for i in range(self.obj_num):
            print(round(obj_values[idx][i], 4), end = "")
            if i < self.obj_num - 1:
                print(", ", end = "" )
        print("]")
        
    # save the results to file as json
    def save(self, fileName, data:Results):
        results = []
        solutions = []
        for i in range(len(data.solutions)):
            solutions.append({
                'solution_id': i,
                'solution': data.solutions[i],
                'objective_values': data.obj_values[i],
            }) 
        results.append({'solutions':solutions})
        pareto_fronts = []
        for i in range(len(data.pareto_fronts)):
            pareto_fronts.append({
                'rank': i,
                'solution_ids': data.pareto_fronts[i]
            })
            
        results.append({'pareto_fronts':pareto_fronts})
        saveJSONFile(fileName, results) 
        
    # plot a Pareto front 
    def plot_front(self, data:Results):

        def onpick(event):
            ind = event.ind
            print('Selected solutions: ')
            for i in range(len(ind)):
                self.print_solution(ind[i], solutions, obj_values)

        fig, ax = plt.subplots()
        plt.xlabel('function 1', fontsize = 14)
        plt.ylabel('function 2', fontsize = 14)
        plt.title('Pareto Front')
        for i in range(len(data.pareto_fronts)):
            # obj_slice[obj_func_idx][front_idx] = obj_values[front_idx][obj_func_idx]
            obj_slices = [[(data.obj_values[i][j]) for i in data.pareto_fronts[i]] for j in range(self.obj_num)] 
            ax.plot(obj_slices[0], obj_slices[1], marker='o', linestyle='', markersize=8, label= 'Rank ' + str(i), picker = True)
        plt.legend()
        fig.canvas.mpl_connect('pick_event', onpick)
      
        plt.show()
    # end output methods
    
    
    # begin genetic operators
    def crossover(self, solution1, solution2):
        cutpoint1 = random.randint(0, self.sol_length)
        cutpoint2 = random.randint(0, self.sol_length)
        if cutpoint1 > cutpoint2:
            tmp = cutpoint1
            cutpoint1 = cutpoint2
            cutpoint2 = tmp
        
        child = []
        if random.random() > 0.5:
            child = solution1[:]            # you *must* copy by value here
            for i in range(cutpoint1, cutpoint2):
                child[i] = solution2[i]
        else:
            child = solution2[:]            # you *must* copy by value here
            for i in range(cutpoint1, cutpoint2):
                child[i] = solution1[i]

        return self.mutation(child)


    def mutation(self, solution):
        for i in range(self.sol_length):
            if random.random() < self.mutate_prob:
                solution[i] = random.randint(self.min_vals[i], self.max_vals[i])
        return solution
    # end of genetic operators
    
    
    # do objective values 'a' dominate objective values 'b'? 
    # assume positive values and minimisation -- add constraint checks here
    def dominates(self, a, b): 
        any = False
        for i in range(self.obj_num):
             if a[i] > b[i]:
                 return False
             if a[i] < b[i]:
                 any = True
        return any
        

    # fast non dominated sort
    def fast_nd_sort(self, obj_values):
        pop_ids = range(len(obj_values))
        S = [[] for i in pop_ids]
        n = [0 for i in pop_ids]
        rank = [0 for i in pop_ids]
        front = [[]]

        for p in pop_ids:
            S[p] = []           # S[p] is the set of points dominated by p          
            n[p] = 0            # n[p] the number of points that dominate p            
            for q in pop_ids:
                if self.dominates(obj_values[p], obj_values[q]):
                    if q not in S[p]:
                        S[p].append(q)
                elif self.dominates(obj_values[q], obj_values[p]):
                    n[p] = n[p] + 1
            # if p is not dominated; add it to the Pareto front
            if n[p] == 0:       
                rank[p] = 0
                if p not in front[0]:
                    front[0].append(p)

        i = 0
        while front[i] != []:
            Q = []
            for p in front[i]:
                for q in S[p]:
                    n[q] = n[q] - 1
                    if  n[q] == 0:
                        rank[q] = i + 1
                        if q not in Q:
                            Q.append(q)
            i = i + 1
            front.append(Q)

        # delete empty set
        del front[len(front)-1]
        return front


    # return the order of values - smallest to largest
    def ordering(self, values):
        value_ids = range(len(values))
        sorted_values = [(values[i], i) for i in value_ids]
        sorted_values.sort(key = lambda a:a[0])
        return [sorted_values[i][1] for i in value_ids]
        
        
    # calculate crowding distance of the solutions in front
    def crowding_distance(self, obj_values, front):
        dist = [0 for i in range(len(front))]
        dist[0]              = 999999       # arbitrary large int
        dist[len(front) - 1] = 999999
        obj_ids = range(self.obj_num)
        
        # obj_slices[obj_func_idx][front_idx] = obj_values[front_idx][obj_func_idx]
        obj_slices = [[obj_values[i][j] for i in front] for j in obj_ids] 

        scale_factor = [max(obj_slices[i]) - min(obj_slices[i]) for i in obj_ids] 
        
        for i in obj_ids:
            if scale_factor[i] == 0:
                scale_factor[i] = 1
 
        order = [self.ordering(obj_slices[i]) for i in obj_ids] 
        
        for i in obj_ids:
            for j in range(1, len(front)-1):
                dist[j] = dist[j] + ((obj_slices[i][order[i][j+1]] - obj_slices[i][order[i][j-1]]) / scale_factor[i])

        return dist
                

    # tournament selection - favour dominant solutions, or 'isolated' solutions with large crowding scores
    def tournament(self, a, b, solutions, obj_values, crowding_values): 
        if self.dominates(obj_values[a], obj_values[b]):
            return solutions[a]
        if self.dominates(obj_values[b], obj_values[a]):
            return solutions[b]
            
        if crowding_values[a] > crowding_values[b]:
            return solutions[a]
        if crowding_values[b] > crowding_values[a]:
            return solutions[b]
            
        if random.random() > 0.5:
            return solutions[a]
        return solutions[b]
        
        
    # construct a matrix of obj_values[pop_size][obj_num]
    def evaluate(self, solutions):
        return  [self.evalFunc(solutions[j]) for j in range(len(solutions))]
        
    def evolve(self):
        # create and evaluate the solutions in the population using the multiple objective functions
        solutions = [[random.randint(self.min_vals[i], self.max_vals[i]) for i in range(self.sol_length)] for j in range(self.pop_size)]
        obj_values = self.evaluate(solutions)

        if self.trace > 0:
            print("\nInitial HV is ", round(self.hv.compute(obj_values),4))
            
        gen_no = 0
        while gen_no < self.gen_num:

            # calculate crowding numbers for the population as a whole
            crowding_values = self.crowding_distance(obj_values, range(self.pop_size))
            
            # construct new generation of solutions using tournament selection and genetic operators
            new_solutions = []
            while len(new_solutions) != self.pop_size:
                idxs = [random.randint(0, self.pop_size-1) for i in range(4)]
                parent1 = self.tournament(idxs[0], idxs[1], solutions, obj_values, crowding_values)
                parent2 = self.tournament(idxs[2], idxs[3], solutions, obj_values, crowding_values)
                new_solutions.append(self.crossover(parent1, parent2))   

            # evaluate them 
            new_obj_values = self.evaluate(new_solutions)

            # elitist - add old solutions and objective values
            new_solutions.extend(solutions)
            new_obj_values.extend(obj_values)  
   
            # sort solution indexes - best have rank = 0, worst have rank = len(pareto_fronts) - 1
            pareto_fronts = self.fast_nd_sort(new_obj_values)

            # get the indexes of the best solutions
            best_solution_ids = []
            rank = 0
            while (len(best_solution_ids) + len(pareto_fronts[rank])) <= self.pop_size:
                best_solution_ids.extend(pareto_fronts[rank])
                rank = rank + 1

            if len(best_solution_ids) < self.pop_size:
                finished = False
                while finished == False:
                    crowding_values = self.crowding_distance(new_obj_values, pareto_fronts[rank])
                    # order pareto_fronts[rank] by crowding values - largest to smallest
                    order = self.ordering(crowding_values)
                    order.reverse()
                    for i in order:
                        best_solution_ids.append(pareto_fronts[rank][i])
                        if len(best_solution_ids) == self.pop_size:
                            finished = True
                            break
                    rank = rank + 1
        
            # construct the next generation (no objective function evaluations required)
            solutions = [new_solutions[i] for i in best_solution_ids]
            obj_values = [new_obj_values[i] for i in best_solution_ids] 
            gen_no = gen_no + 1
                        
            if self.trace > 1:
                front_values  = [(new_obj_values[i]) for i in pareto_fronts[0]]
                volume = self.hv.compute(front_values)
                print(gen_no, ") HV is ", round(volume,4))
        
        pareto_fronts = self.fast_nd_sort(obj_values)
   
        if self.trace > 0:
            front_vals = [(obj_values[i]) for i in pareto_fronts[0]]
            print("Final HV is ", round(self.hv.compute(front_vals), 4))
            print("\nThe Pareto front on generation ", gen_no, " is")
            print(pareto_fronts[0])
        
        if self.trace > 1:
            print("\nThe best solutions on generation ", gen_no, " are")
            for i in pareto_fronts[0]:
                self.print_solution(i, solutions, obj_values)
            print("\n") 
        
        return Results(solutions, obj_values, pareto_fronts)


# 
# Dummy multiple objective function - multi_obj_func
# 

# first function to minimise
def function1(solution):
    value = 0
    for i in range(len(solution)-2):
        tmp = abs((solution[i] - i))
        value = value + (tmp * tmp)
    return value

# second function to minimise
def function2(solution):
    value = 0
    for i in range(2, len(solution)):
        tmp = abs(solution[i] - (len(solution) - i))
        value = value + (tmp * tmp) 
    return value

def multi_obj_func(solution):
    return function1(solution), function2(solution)


#
#
#


def main():
    
    fileName = "NSGA2.config"
    
    args = sys.argv[1:]
    if len(args) == 1:
        fileName = args[0]
        
    jsondata = readJSONFile(fileName)

    sol_spec = [SolutionSpec(a['sol_num'], a['min_val'], a['max_val']) for a in jsondata['sol_spec']]
   
    config = Config (
        jsondata['seed'],           # random number seed 
        jsondata['pop_size'],       # the number of solutions to work with; typically between 20 and 100
        jsondata['gen_num'],        # the number of generations or iterations; typically between 100 and 10000
        sol_spec,                   # specify the solution vector [(sol_spec1), (sol_spec2), ...]
        jsondata['mutate_prob'],    # probability of mutating element of solution vector
        multi_obj_func,             # the multi-objective function to be minimised
        2,                          # the number of objectives to be minimised
        [1000.0, 1000.0],           # reference point for hypervolume calculation
        jsondata['trace']           # trace level - 0 implies 'silent'   
    )

    optimiser = NSGA2(config)
    
    results = optimiser.evolve()

    if config.trace > 1:
        optimiser.plot_front(results)
        
    optimiser.save('results.json', results)
      

if __name__ == '__main__':
    main()

