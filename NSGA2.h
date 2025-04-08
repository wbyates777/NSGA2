/* NSGA2 19/03/2025

 $$$$$$$$$$$$$$$
 $   NSGA2.h   $
 $$$$$$$$$$$$$$$

 by W.B. Yates
 Copyright (c) W.B. Yates. All rights reserved.
 History:

 A fast and elitist multiobjective genetic algorithm: NSGA-II
 
 Assumes positive objective values to be minimised.
 
 K. Deb, A. Pratap, S. Agarwal and T. Meyarivan, "A fast and elitist multiobjective genetic algorithm: NSGA-II,"
 in IEEE Transactions on Evolutionary Computation,
 vol. 6, no. 2, pp. 182-197, April 2002,
 doi: 10.1109/4235.996017
 
 Example 1
 
 
 Config config;
 
 config.seed = 21;
 config.pop_size  = 100;
 config.gen_num = 50;
 
 config.trace = 3;
 config.mutate_prob = 0.07;
 
 SolSpec solspec1;
 solspec1.sol_num =  15;
 solspec1.min_val =  -10;
 solspec1.max_val =  10;
 
 SolSpec solspec2;
 solspec2.sol_num =  15;
 solspec2.min_val =  0;
 solspec2.max_val =  20;
 
 config.obj_num = 2; // 
 config.con_num = 2; // constraints
 config.evalFunc = multi_obj_func;
 
 config.sol_spec = { solspec1, solspec2 };
 config.refPoint = { 10000.0,10000.0};
 
 NSGA2 nsga(config);
 
 ///
 
 // first function to minimise
 double 
 function1(const Chrom &solution)
 {
     double value = 0;
     for (int i = 0; i < solution.size()-2; ++i)
     {
         double tmp = std::fabs(solution[i] - i);
         value += (tmp * tmp);
     }
     return value;
 }

 // second function to minimise
 double 
 function2(const Chrom &solution)
 {
     double value = 0;
     int N = (int)solution.size();
     for (int i = 2; i < solution.size(); ++i)
     {
         double tmp = std::fabs(solution[i] - (N - i) );
         value += (tmp * tmp);
     }
     return value;
 }


 std::pair<ObjValue, Constraint>
 multi_obj_func(const Chrom &solution)
 {
     std::vector<double> obj_values = { function1(solution), function2(solution) };
     std::vector<int> constraints = {0, 0};
     if (solution[0] == 3) 
         constraints[0]++;
     return { obj_values, constraints };
 }

*/


#ifndef __NSGA2_H__
#define __NSGA2_H__


#include <vector>

#ifndef __HYPERVOLUME_H__
#include "HyperVolume.h"
#endif

#ifndef __URAND_H__
#include "URand.h"
#endif

typedef std::vector<int> Chrom; // chromosome of length n
typedef std::vector<Chrom> Solutions;

typedef std::vector<int> Constraint; // number of violations for each constraint
typedef std::vector<Constraint> Constraints;

typedef std::vector<double> ObjValue; // value of each objective to be optimised
typedef std::vector<ObjValue> ObjValues;

typedef std::vector<std::vector<int>> ParetoFronts;

// a function EvalFunction that takes a Chrom and returns an ObjValue and a Constraint
// must be provided by the user
typedef std::function<std::pair<ObjValue, Constraint>(const Chrom&)> EvalFunction;

struct SolSpec
{
    int sol_num;
    int min_val;
    int max_val;
};

struct Config
{
    int    seed;
    int    pop_size;  // population size (usually 50-500)
    int    gen_num;   // number of iterations or generations
    int    obj_num;   // number of doubles in an ObjValue
    int    con_num;   // number of ints in a Constraint
    int    trace;
    double mutate_prob;
    EvalFunction evalFunc; // takes a Chrom and returns an ObjValue and a Constraint
    std::vector<SolSpec> sol_spec;
    Point refPoint;   // see HyperVolume.h
};

struct Results 
{
    Results( void )=default;
    Results(const Solutions &s, const ObjValues &ov, const Constraints &cs, const ParetoFronts &pf): solutions(s), obj_values(ov), constraints(cs), pareto_fronts(pf) {}
    ~Results( void )=default;
    
    Solutions    solutions; 
    Constraints  constraints;
    ObjValues    obj_values; 
    ParetoFronts pareto_fronts; // best solution idx are in pareto_fronts[0]
};



class NSGA2
{
public:

    NSGA2( const Config& config );
    ~NSGA2( void )=default;

    Results 
    evolve( void );
    
private:

    bool
    is_in( int x, const std::vector<int> &v ) const { return (std::find(v.begin(), v.end(), x) != v.end()); }
    
    void
    print_solution(int idx, const Solutions &solutions, const ObjValues &obj_values, const Constraints &constraints) const;
    
    Chrom
    crossover( const Chrom &solution1, const Chrom &solution2);
    
    Chrom
    mutation(Chrom solution);
    
    bool
    dominates(const ObjValue &a, const ObjValue &b, const Constraint &ca, const Constraint &cb) const;
    
    ParetoFronts
    fast_nd_sort(const ObjValues &obj_values, const Constraints &constraints) const;
    
    std::vector<int>
    ordering( const std::vector<double> &values ) const;
    
    std::vector<double>
    crowding_distance(const ObjValues &obj_values, const std::vector<int> &front) const;
    
    Chrom
    tournament(int a, int b, const Solutions &solutions, 
               const ObjValues &obj_values, const Constraints &constraints, 
               std::vector<double> &crowding_values);
    
    std::pair<ObjValues, Constraints>
    evaluate(const Solutions &solutions) const;

    int m_pop_size;
    int m_gen_num;
    int m_sol_length;
    int m_trace;
    int m_obj_num;
    int m_con_num;
    
    double m_mutate_prob;
    
    std::vector<int> m_sol_nums;
    std::vector<int> m_min_vals;
    std::vector<int> m_max_vals;
    
    EvalFunction m_evalFunc;
    HyperVolume m_hv;
    URand m_ran;
    
};

#endif


