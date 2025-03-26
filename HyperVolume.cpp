/* HyperVolume 21/03/2025

 $$$$$$$$$$$$$$$$$$$$$$$
 $   HyperVolume.cpp   $
 $$$$$$$$$$$$$$$$$$$$$$$

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
 
 Test example point
 
  HV  = 7348493500.00000000 
 
  when 
  
  m_refPoint = {2000.0, 2000.0, 2000.0};
 
  std::vector<Point> front = {
    { 495.0, -417.0, 0.0 }, 
    { 658.0, 366.0, 1.0 }, 
    { 471.0, 733.0, 0.5 }, 
    { 697.0, 258.0, 10.0 }, 
    { 1111.0, 214.0, 11.0  }, 
    { 876.0, 253.0, 12.0  }, 
    { 476.0, 713.0, 13.0  }, 
    { 908.0, 237.0, 10.0  }, 
    { 1133.0, 213.0, 10.0  }, 
    { 672.0, 306.0, 4.0  }, 
    { 467.0, 815.0, 3.0  }, 
    { 1321.0, 200.0, -1.0  }, 
    { 657.0, 374.0, -1.0  }
  }; 

*/


#ifndef __HYPERVOLUME_H__
#include "HyperVolume.h"
#endif

#include <iostream>
#include <cassert>
#include <numeric>

HyperVolume::HyperVolume( void ) {}

HyperVolume::HyperVolume( const Point &refPoint )
{
    m_refPoint = refPoint;
    setSentinal();
}


void 
HyperVolume::append(Node* node, int index)
// Appends a node to the end of the list at the given index.
{
    Node *lastButOne = m_sentinel.prev[index];
    node->next[index] = &m_sentinel;
    node->prev[index] = lastButOne;
    // set the last element as the new one
    m_sentinel.prev[index] = node;
    lastButOne->next[index] = node;
}

void
HyperVolume::remove(Node *node, int index, Point &bounds)
{
// Removes and returns 'node' from all lists in [0, 'index'[. 
    for (int i = 0; i < index; ++i)
    {
        Node* predecessor = node->prev[i];
        Node* successor   = node->next[i];
        predecessor->next[i] = successor;
        successor->prev[i]   = predecessor;  
        
        if (bounds[i] > node->point[i])
            bounds[i] = node->point[i];
    }
}   
    
void
HyperVolume::reinsert(Node *node, int index, Point &bounds)
/*  
    Inserts 'node' at the position it had in all lists in [0, 'index'[
    before it was removed. This method assumes that the next and previous 
    nodes of the node that is reinserted are in the list.
  */
{
    for (int i = 0; i < index; ++i)
    {
        node->prev[i]->next[i] = node;
        node->next[i]->prev[i] = node;
        if (bounds[i] > node->point[i])
            bounds[i] = node->point[i];
    }
}


bool 
HyperVolume::weaklyDominates(const Point &point, const Point &other) const
{
    for (int i = 0; i < point.size(); ++i)
    {
        if (point[i] > other[i])
            return false;
    }
    return true;
}


double 
HyperVolume::hvRecursive(int dimIndex, int length, Point &bounds)
/*  
      Recursive call to hypervolume calculation.
      In contrast to the paper, the code assumes that the reference point
      is [0, ..., 0]. This allows the avoidance of a few operations.
  */
{
    double hvol = 0.0;
  
    if (length == 0)
        return hvol;

    if (dimIndex == 0)
    {
        // special case: only one dimension - why using hypervolume at all?
        return -m_sentinel.next[0]->point[0];
    }
    else if (dimIndex == 1)
    {   
        // special case: two dimensions, end recursion
        Node *q  = m_sentinel.next[1];
        double h = q->point[0];
        Node *p  = q->next[1];
        while (p != &m_sentinel)
        {
            const Point &pCargo = p->point;
            hvol += h * (q->point[1] - pCargo[1]);
            if (pCargo[0] < h)
            {
                h = pCargo[0];
            }
            q = p;
            p = q->next[1];
        }
        hvol += h * q->point[1];
        return hvol;
    }
    else
    {
        Node *p = &m_sentinel;
        Node *q = p->prev[dimIndex];
        while (q != &m_sentinel)
        {   
            if (q->ignore < dimIndex)
                q->ignore = 0;
            q = q->prev[dimIndex];
        }
        
        q = p->prev[dimIndex];
        while (length > 1 && (q->point[dimIndex] > bounds[dimIndex] || q->prev[dimIndex]->point[dimIndex] >= bounds[dimIndex]))
        {
            p = q;
            remove(p, dimIndex, bounds);
            q = p->prev[dimIndex];
            length -= 1;
        }
        
        std::vector<double> &qArea = q->area;
        const Point &qCargo = q->point;
        Node *qPrevDimIndex = q->prev[dimIndex];
        if (length > 1)
        {
            hvol = qPrevDimIndex->volume[dimIndex] + qPrevDimIndex->area[dimIndex] * (qCargo[dimIndex] - qPrevDimIndex->point[dimIndex]);
        }
        else
        {   
            qArea[0] = 1;
            for (int j = 0; j < dimIndex; ++j) 
                qArea[j+1] = qArea[j] * -qCargo[j];
        }
        q->volume[dimIndex] = hvol;
        if (q->ignore >= dimIndex) 
        {
            qArea[dimIndex] = qPrevDimIndex->area[dimIndex];
        }
        else 
        {
            qArea[dimIndex] = hvRecursive(dimIndex - 1, length, bounds);
            if (qArea[dimIndex] <= qPrevDimIndex->area[dimIndex])
                q->ignore = dimIndex;
        }
        
        while (p != &m_sentinel)
        {
            double pCargoDimIndex = p->point[dimIndex];
            hvol += q->area[dimIndex] * (pCargoDimIndex - q->point[dimIndex]);
            bounds[dimIndex] = pCargoDimIndex;
            reinsert(p, dimIndex, bounds);
            length++;
            q = p;
            p = p->next[dimIndex];
            q->volume[dimIndex] = hvol;
            
            if (q->ignore >= dimIndex)
                q->area[dimIndex] = q->prev[dimIndex]->area[dimIndex];
            else q->area[dimIndex] = hvRecursive(dimIndex - 1, length, bounds);
            
            if (q->area[dimIndex] <= q->prev[dimIndex]->area[dimIndex])
                q->ignore = dimIndex;
        }
        hvol -= q->area[dimIndex] * q->point[dimIndex];
        
        return hvol;
    }
}

void
HyperVolume::setRefPoint( const Point &refPoint )
{
    m_refPoint = refPoint;
    setSentinal();
}

void
HyperVolume::setSentinal( void )
{
    int dims = (int) m_refPoint.size();
    
    m_sentinel.point = Point(dims, 0.0);
    
    m_sentinel.ignore = 0;
    m_sentinel.next.clear();
    m_sentinel.prev.clear();
    m_sentinel.area.clear();
    m_sentinel.volume.clear();
    
    m_sentinel.next.resize(dims, nullptr);
    m_sentinel.prev.resize(dims, nullptr);
    m_sentinel.area.resize(dims, 0.0);
    m_sentinel.volume.resize(dims, 0.0);

    for (int i = 0; i < dims; ++i)
    {
        m_sentinel.next[i] = &m_sentinel;
        m_sentinel.prev[i] = &m_sentinel; 
    }
}


void
HyperVolume::preProcess(const std::vector<Point> &front)
// Sets up the list data structure needed for calculation.
{
    setSentinal();
    
    int dims = (int) m_refPoint.size();

    m_nodes.resize(front.size());
    for (int i = 0; i < front.size(); ++i)
    {
        m_nodes[i] = Node(front[i]);  
    }
    
    std::vector<int> order(m_nodes.size());
    std::iota(order.begin(), order.end(), 0);
    
    for (int i = 0; i < dims; ++i)
    {
        // order the nodes by dimention
        auto op1 = [&](int idx1, int idx2)  { return m_nodes[idx1].point[i] < m_nodes[idx2].point[i]; };
        std::sort(order.begin(), order.end(), op1 );
                               
        for (int j : order)
        {
            Node *lastButOne = m_sentinel.prev[i];
            m_nodes[j].next[i] = &m_sentinel;
            m_nodes[j].prev[i] = lastButOne;
            // set the last element as the new one
            m_sentinel.prev[i] = &m_nodes[j];
            lastButOne->next[i] = &m_nodes[j]; 
        } 
    }
}



double
HyperVolume::compute(const std::vector<Point> &front)
{
    if (front.empty()) 
        return 0.0;
    
    assert(m_refPoint.size() == front[0].size());
    
    std::vector<Point> relevantPoints;
    
    for (const Point &point : front)
    {
        // only consider points that dominate the reference point
        if (weaklyDominates(point, m_refPoint))
            relevantPoints.push_back(point);
    }

    int dims = (int) m_refPoint.size();
    if (any(m_refPoint))
    {
        // shift points so that refPoint == [0, ..., 0]
        // this way the reference point doesn't have to be explicitly used
        // in the HV computation
        for (int i = 0; i < relevantPoints.size(); ++i)
        {
            for (int j = 0; j < dims; ++j)
            {
                relevantPoints[i][j] -= m_refPoint[j];
            }
        }
    }
 
    preProcess(relevantPoints);
    
    Point bounds(dims, -std::numeric_limits<double>::max()); 
    
    return hvRecursive(dims - 1, (int) relevantPoints.size(), bounds);
}
