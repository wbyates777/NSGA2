/* NSGA2 19/03/2025

 $$$$$$$$$$$$$$$$$
 $   NSGA2.cpp   $
 $$$$$$$$$$$$$$$$$

 by W.B. Yates
 Copyright (c) W.B. Yates. All rights reserved.
 History:

 A fast and elitist multiobjective genetic algorithm: NSGA-II
 
 Assumes positive objective values to be minimised.
 
 K. Deb, A. Pratap, S. Agarwal and T. Meyarivan, "A fast and elitist multiobjective genetic algorithm: NSGA-II,"
 in IEEE Transactions on Evolutionary Computation,
 vol. 6, no. 2, pp. 182-197, April 2002,
 doi: 10.1109/4235.996017
 
*/


#ifndef __NSGA2_H__
#include "NSGA2.h"
#endif

#include <iostream>

NSGA2::NSGA2( const Config& config ) 
{
    m_ran.seed(config.seed); 
    
    m_pop_size    = config.pop_size;
    m_gen_num     = config.gen_num; 
    m_trace       = config.trace;
    m_obj_num     = config.obj_num;
    m_con_num     = config.con_num;
    m_evalFunc    = config.evalFunc;
    
    m_sol_length  = 0;
    for (const SolSpec &cs : config.sol_spec)
    {
        m_sol_nums.push_back(cs.sol_num);

        m_min_vals.insert(m_min_vals.end(), cs.sol_num, cs.min_val);
        m_max_vals.insert(m_max_vals.end(), cs.sol_num, cs.max_val);

        m_sol_length += cs.sol_num;
    }
    
    m_hv.setRefPoint(config.refPoint);
    
    if (m_trace > 0)
    {
        std::cout << "\nConstructed population of " << m_pop_size << " solutions of length " << m_sol_length<< std::endl;
        std::cout << "Running for " << m_gen_num << " generations" << std::endl;
    }
}

void
NSGA2::print_solution(int idx, const Solutions &solutions, const ObjValues &obj_values, const Constraints &constraints) const
{
    std::cout << idx << ") ";
    for (int j = 0; j < m_sol_length; ++j)
        std::cout << solutions[idx][j] << " "; 
    
    // objective function values
    std::cout << "[";
    for (int i = 0; i < m_obj_num; ++i) 
    {
        std::cout << std::format("{:.2f}", obj_values[idx][i]);
        if (i < m_obj_num - 1)
            std::cout << ", ";
    } 
    std::cout << "]";
    
    // constraint violations
    std::cout << "[";
    for (int i = 0; i < m_con_num; ++i) 
    {
        std::cout << constraints[idx][i]; 
        if (i < m_con_num - 1)
            std::cout << ", ";
    } 
    std::cout << "]" << std::endl;
}

//
// begin genetic operators
//

Chrom
NSGA2::crossover( const Chrom &solution1, const Chrom &solution2)
{
    int cutpoint1 = m_ran.rndInt(m_sol_length);
    int cutpoint2 = m_ran.rndInt(m_sol_length);
    
    if (cutpoint1 > cutpoint2)
        std::swap(cutpoint1, cutpoint2);
 
    std::vector<int> child;
    if (m_ran.rndFloat() > 0.5)
    {
        child = solution1;           
        for (int i = cutpoint1; i < cutpoint2; ++i)
            child[i] = solution2[i];
    }
    else
    {
        child = solution2;          
        for (int i = cutpoint1; i < cutpoint2; ++i)
            child[i] = solution1[i];
    }
    
    return mutation(child);
}

Chrom
NSGA2::mutation(Chrom solution)
{
    for (int i = 0; i < m_sol_length; ++i)
    {
        if (m_ran.rndFloat() < m_mutate_prob)
            solution[i] = m_ran.rndInt(m_min_vals[i], m_max_vals[i] + 1); // inclusive of max_vals[j]
    }
    return solution;
}

//
// end of genetic operators
//

bool
NSGA2::dominates(const ObjValue &a, const ObjValue &b, const Constraint &ca, const Constraint &cb ) const
// do objective values 'a' dominate objective values 'b'?
// assume positive values and minimisation
{
    // constraints
    int count_a = 0, count_b = 0, count_a_beats_b = 0;
    
    for (int i = 0; i < m_con_num; ++i)
    {
        if (ca[i] > 0)
            count_a++;
        if (cb[i] > 0)
            count_b++;
        if (ca[i] < cb[i])
            count_a_beats_b++;
    }
     
    //  solution a is feasible and solution b is not        
    if (count_a == 0 && count_b > 0)
        return true;
     
    // solutions a and b are both infeasible, but solution a has less overall constraint violations.
    if (count_a > 0 && count_b > 0)
    {
        if (count_a <= count_b && count_a_beats_b > 0)
            return true;
        else return false;
    }
 
    // objective values
    bool any = false;
    for (int i = 0; i < m_obj_num; ++i)
    {
        if (a[i] > b[i])
            return false;
        
        if (a[i] < b[i])
            any = true;
    }
    return any;
}

ParetoFronts
NSGA2::fast_nd_sort(const ObjValues &obj_values, const Constraints &constraints) const
// fast non-dominated sort
{
    int N = (int) obj_values.size();
    
    std::vector<int>  pop_ids(N);
    std::iota(pop_ids.begin(), pop_ids.end(), 0);
    
    std::vector<std::vector<int>> S(N);
    std::vector<int> n(N, 0), rank(N, 0);

    ParetoFronts fronts(1);

    for (int p : pop_ids)
    {   
        S[p].clear();           // S[p] is the set of points dominated by p
        n[p] = 0;               // n[p] the number of points that dominate p
        for (int q : pop_ids)
        {
            if (dominates(obj_values[p], obj_values[q], constraints[p], constraints[q]))
            {    
                if (!is_in(q, S[p])) //  if (q not in S[p])   
                    S[p].push_back(q);
            } 
            else if (dominates(obj_values[q], obj_values[p], constraints[q], constraints[p]))
                n[p] += 1;
            
        }
        // if p is not dominated; add it to the Pareto front
        if (n[p] == 0)
        {
            rank[p] = 0;
            if (!is_in(p, fronts[0])) // if (p not in front[0])
                fronts[0].push_back(p);
        }
    }
    
    int i = 0;
    while (!fronts[i].empty())
    {    
        std::vector<int> Q;
        for (int p : fronts[i])
        {
            for (int q : S[p])
            {
                n[q] -= 1;
                if  (n[q] == 0)
                {
                    rank[q] = i + 1;
                    if (!is_in(q, Q))
                        Q.push_back(q);
                }
            }
        }
        i++;
        fronts.push_back(Q);
    }
    
    // delete empty set
    fronts.pop_back();
    return fronts;
}

std::vector<int>
NSGA2::ordering( const std::vector<double> &values ) const
{
    std::vector<int> indexes(values.size());
    std::iota(indexes.begin(), indexes.end(), 0);
    auto op = [&values](int i, int j)  { return values[i] < values[j]; };
    std::sort( indexes.begin(), indexes.end(), op );
    
    return indexes;
}


std::vector<double>
NSGA2::crowding_distance(const ObjValues &obj_values, const std::vector<int> &front) const
// calculate crowding distance of the solutions in front
{
    std::vector<double> dist(front.size(), 0.0); 
    dist[0]                = std::numeric_limits<double>::max();  
    dist[front.size() - 1] = std::numeric_limits<double>::max();
    
    ObjValues obj_slices(m_obj_num, std::vector<double>(front.size()));
    for (int i = 0; i < m_obj_num; ++i)
    {
        for (int j = 0; j < front.size(); ++j)
        {
            obj_slices[i][j] = obj_values[front[j]][i];
        }
    }

    for (int i = 0; i < m_obj_num; ++i)
    {
        std::vector<int> order = ordering(obj_slices[i]);
 
        // val = maxval - minval
        double val = obj_slices[i][order.back()] - obj_slices[i][order[0]];
        double scale_factor = (val == 0.0) ? 1.0 : val;
        
        for (int j = 1; j <  front.size()-1; ++j)
        {
            dist[j] += ((obj_slices[i][order[j+1]] - obj_slices[i][order[j-1]]) / scale_factor);
        }
    }
    return dist;
}

Chrom
NSGA2::tournament( int a, int b, const Solutions &solutions, 
                   const ObjValues &obj_values, const Constraints &constraints, 
                   std::vector<double> &crowding_values )
{
    if (dominates(obj_values[a], obj_values[b], constraints[a], constraints[b]))
        return solutions[a];
    if (dominates(obj_values[b], obj_values[a], constraints[b], constraints[a]))
        return solutions[b];

    if (crowding_values[a] > crowding_values[b])
        return solutions[a];
    if (crowding_values[b] > crowding_values[a])
        return solutions[b];

    if (m_ran.rndFloat() > 0.5)
        return solutions[a];
    return solutions[b];
}


std::pair<ObjValues, Constraints>
NSGA2::evaluate(const Solutions &solutions) const
// construct a matrix of obj_values v[pop_num][obj_num] and constraint violations  c[pop_num][con_num]
{
    std::pair<ObjValues, Constraints> retVal(solutions.size(), solutions.size());

    for (int i = 0; i <  solutions.size(); ++i)
    {
        std::pair<ObjValue, Constraint> res = m_evalFunc(solutions[i]);
        retVal.first[i] = res.first; 
        retVal.second[i] = res.second; 
    }
    return retVal;
}


Results 
NSGA2::evolve( void )
{
    // create and evaluate the solutions in the population using the multiple objective functions
    Solutions solutions(m_pop_size, std::vector<int>(m_sol_length));
    for (int i = 0; i < m_pop_size; ++i)
    {
        for (int j = 0; j < m_sol_length; ++j)
        {
            solutions[i][j] = m_ran.rndInt(m_min_vals[j], m_max_vals[j] + 1); // inclusive of max_vals[j]
        }
    }
    
    std::pair<ObjValues, Constraints> initial_res = evaluate(solutions);
    ObjValues &obj_values     = initial_res.first;
    Constraints &constraints  = initial_res.second; 
    
    if (m_trace > 0) 
        std::cout << "\nInitial HV is " << std::format("{:.4f}", m_hv.compute(obj_values)) << std::endl;
    
    int idxs[4];
    int gen_no = 0;
    ParetoFronts pareto_fronts;
    
    while (gen_no < m_gen_num)
    {
        // calculate crowding numbers for the population as a whole
        std::vector<int> pop_front(m_pop_size);
        std::iota(pop_front.begin(), pop_front.end(), 0); 
        std::vector<double> crowding_values = crowding_distance(obj_values, pop_front);
        
        // construct new generation of new solutions using tournament selection and genetic operators
        Solutions new_solutions;
        new_solutions.reserve(m_pop_size);
        while (new_solutions.size() != m_pop_size)
        {
            for (int i = 0; i < 4; ++i)
                idxs[i] = m_ran.rndInt(m_pop_size);
            
            Chrom parent1 = tournament(idxs[0], idxs[1], solutions, obj_values, constraints, crowding_values);
            Chrom parent2 = tournament(idxs[2], idxs[3], solutions, obj_values, constraints, crowding_values);
            new_solutions.push_back(crossover(parent1, parent2));
        }
        
        // evaluate them
        std::pair<ObjValues, Constraints> new_res = evaluate(new_solutions);
        ObjValues   &new_obj_values  = new_res.first;
        Constraints &new_constraints = new_res.second; 
        
        // elitist - add old solutions, objective, and constraint values
        new_solutions.insert(new_solutions.end(), solutions.begin(), solutions.end());
        new_obj_values.insert(new_obj_values.end(), obj_values.begin(), obj_values.end());
        new_constraints.insert(new_constraints.end(), constraints.begin(), constraints.end());
        
        // sort solution indexes - 'best' have rank = 0, 'worst' have rank = pareto_fronts.size() - 1
        pareto_fronts = fast_nd_sort(new_obj_values, new_constraints);
        
        // get the indexes of the best solutions
        std::vector<int> best_solution_ids;
        
        int rank = 0;
        while ((best_solution_ids.size() + pareto_fronts[rank].size()) <= m_pop_size)
        {   
            best_solution_ids.insert(best_solution_ids.end(), pareto_fronts[rank].begin(), pareto_fronts[rank].end());
            rank++;
        }
        
        if (best_solution_ids.size() < m_pop_size)
        {   
            bool finished = false;
            while (finished == false)
            {    
                crowding_values = crowding_distance(new_obj_values, pareto_fronts[rank]);
                // order pareto_fronts[rank] by crowding values - largest to smallest
                std::vector<int> order = ordering(crowding_values);
                order = std::vector<int>(order.rbegin(), order.rend()); // reverse order
                for (int i : order)
                {   
                    best_solution_ids.push_back(pareto_fronts[rank][i]);
                    if (best_solution_ids.size() == m_pop_size)
                    {   
                        finished = true;
                        break;
                    }
                }
                rank++;
            }
        }
        
        // construct the next generation (no objective function evaluations required)
        int count = 0;
        for (int i : best_solution_ids)
        {
            solutions[count]   = new_solutions[i];
            obj_values[count]  = new_obj_values[i];
            constraints[count] = new_constraints[i];
            count++;
        }
        
        gen_no++;
        
        if (m_trace > 1)
        {
            std::vector<std::vector<double>> front_values(pareto_fronts[0].size());
            int count = 0;
            for (int i : pareto_fronts[0])
                front_values[count++] = new_obj_values[i];
            double volume = m_hv.compute(front_values);
            std::cout << gen_no <<  ") HV is " << std::format("{:.4f}", volume) << std::endl;
        }
    }
    
    // finished evolution 
    pareto_fronts = fast_nd_sort(obj_values, constraints);
    
    if (m_trace > 0)
    {
        std::vector<std::vector<double>> front_values(pareto_fronts[0].size());
        for (int i : pareto_fronts[0])
            front_values[i] = obj_values[i];
        
        std::cout << "Final HV is " << m_hv.compute(front_values)  << std::endl;
        std::cout << "\nThe Pareto front on generation " << gen_no << " is (" <<  pareto_fronts[0].size() << ") " << std::endl;
        
        for (int i : pareto_fronts[0])
            std::cout << i << " ";
        
        std::cout << std::endl;
    }
    
    if (m_trace > 1)
    {
        std::cout << "\nThe best solutions on generation " << gen_no <<  " are:" << std::endl;
        for (int i : pareto_fronts[0])
        {
            print_solution(i, solutions, obj_values, constraints);
        }
    }
    return Results(solutions, obj_values, constraints, pareto_fronts);
}

