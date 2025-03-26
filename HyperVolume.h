/* HyperVolume 21/03/2025

 $$$$$$$$$$$$$$$$$$$$$
 $   HyperVolume.h   $
 $$$$$$$$$$$$$$$$$$$$$

 by W.B. Yates
 Copyright (c) W.B. Yates. All rights reserved.
 History:

 Translated from python to c++. Tested against original.
 
 """
 Hypervolume computation based on variant 3 of the algorithm in the paper:
 C. M. Fonseca, L. Paquete, and M. Lopez-Ibanez. An improved dimension-sweep
 algorithm for the hypervolume indicator. In IEEE Congress on Evolutionary
 Computation, pages 1157-1163, Vancouver, Canada, July 2006.

 Minimization is implicitly assumed here!

 """
 
 Original code
 
 #    Copyright (C) 2010 Simon Wessing
 #    TU Dortmund University
 #
 #    This program is free software: you can redistribute it and/or modify
 #    it under the terms of the GNU General Public License as published by
 #    the Free Software Foundation, either version 3 of the License, or
 #    (at your option) any later version.
 #
 #    This program is distributed in the hope that it will be useful,
 #    but WITHOUT ANY WARRANTY; without even the implied warranty of
 #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 #    GNU General Public License for more details.
 #
 #    You should have received a copy of the GNU General Public License
 #    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 __author__ = "Simon Wessing"
 
  
*/


#ifndef __HYPERVOLUME_H__
#define __HYPERVOLUME_H__


#include <vector>

typedef  std::vector<double> Point;


class HyperVolume
{

public:

    HyperVolume( void );
    HyperVolume( const Point &refPoint );
    ~HyperVolume( void )=default;

    void
    setRefPoint( const Point &refPoint );
    
    double 
    compute(const std::vector<Point> &front);

    
private:

    typedef struct Node
    /* 
       A special data structure needed by FonsecaHyperVolume. 
       It consists of several doubly linked lists that share common nodes. So, 
       every node has multiple predecessors and successors, one in every list.
    */
    {
        Node( void ): ignore(0) {}
        Node(const Point& p): ignore(0), point(p)
        {
            next.resize(point.size(), nullptr);  
            prev.resize(point.size(), nullptr);
            volume.resize(point.size(), 0.0); 
            area.resize(point.size(), 0.0);
        }
        ~Node( void ) { ignore = 0; }


        int    ignore;     
        
        Point  point;
        
        std::vector<Node*>  next;
        std::vector<Node*>  prev;
        std::vector<double> area;    
        std::vector<double> volume;  
    } Node;
    
    void 
    append(Node* node, int index); 

    void
    remove(Node *node, int index, Point &bounds);
        
    void
    reinsert(Node *node, int index, Point &bounds);
    // 
    
    template<typename T>
    bool
    any( const std::vector<T>& vec ) { return (*std::min_element(vec.begin(), vec.end()) > 0); }

    void
    setSentinal( void );
    
    bool 
    weaklyDominates(const Point &point, const Point &other) const;
    
    void
    preProcess(const std::vector<Point> &front);
    
    double 
    hvRecursive(int dimIndex, int length, Point &bounds);
   

    Point m_refPoint;
    Node m_sentinel; 
    std::vector<Node> m_nodes;
     
};

#endif


