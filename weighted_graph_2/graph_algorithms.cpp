#ifndef GRAPH_ALGS
#define GRAPH_ALGS

#include <map>
#include <vector>
#include <queue>
#include <stack>
#include <list>
#include <deque>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <utility>
#include <algorithm>
#include <limits>
#include "weighted_graph.hpp"
#include "easy_weighted_graph_algorithms.cpp"

//Returns true if the graph is connected, false otherwise.
template <typename vertex>
bool is_connected(const weighted_graph<vertex>& g){
	// check for empty graph
	if (g.num_vertices() == 0) {
		return true;
	}
	
	// do traversal (dft is faster), if it visits all vertices, its connected
	std::vector<vertex> vertices = depth_first(g, *g.begin());
	if (vertices.size() == g.num_vertices()) {
		return true;
	}
	return false;	
}

//Returns a vector of weighted graphs, where each weighted graph is a connected
//component of the input graph.
template <typename vertex>
std::vector<weighted_graph<vertex>> connected_components(const weighted_graph<vertex>& g){	
	// traverse, mark visited vertex -> create/form a new graph
	// for every unvisited vertex, traverse again, and form new graph
	std::vector<weighted_graph<vertex> > con_components;
	// check for empty graph
	if(g.num_vertices() <= 0) return con_components;
	// set of visited vertices
	std::unordered_set<vertex> visited;
	// give graph iterator start as start vertex
	vertex starting_vertex = *g.begin();
	bool endLoop = false;
	// whilst graph isn't broken loop through vertices
	if (!is_connected(g)) {
		while(visited.size() < g.num_vertices()){
			std::vector<vertex> conn_vertices = depth_first(g, starting_vertex);
			visited.insert(starting_vertex);
			weighted_graph<vertex> cloned_graph = weighted_graph<vertex>(g);

			// check if all vertex in graph are not in vertices

			for( auto i = g.begin(); i != g.end(); i++){
				bool vertexFound = false;
				for( auto j = conn_vertices.begin(); j!= conn_vertices.end(); j ++){
					
					if( *i == *j) {
						visited.insert(  *i);
						vertexFound = true;
					}
				}
				
				if( !vertexFound) {
					cloned_graph.remove_vertex(*i);				
				}
			}	
			//look for the next starting vertex
			for( auto i = g.begin(); i != g.end(); i++){
				if (visited.count(*i) == 0) {
					starting_vertex = *i;
					break;
					}
			}					
			con_components.push_back(cloned_graph);
		}
	} else {
		con_components.push_back(g);
	}
	return con_components;
}

//Returns a map of the vertices of the weighted graph g and their distances from
//the given starting vertex v.
template <typename vertex> 
std::map<vertex, int> dijkstras(const weighted_graph<vertex>& g, const vertex& v){
	// load the vertices from the graph
	std::unordered_set<vertex> unvisited(g.begin(), g.end());
	// map data structure to store vertex and its distance
	std::map<vertex, int> dist;
	
	// loop through, d[v] = 0, all others infinity (or large number)
    for( auto it = g.begin(); it != g.end(); it++ ){
        auto s = *it;
        if( s == v)
            dist[s] = 0;
        else
            dist[s] = std::numeric_limits<int>::max();
            
    }
	
	// set current vertex to v;
	vertex current = v;
	
    
	// while there are still vertices to visit
	while (!unvisited.empty()) {
		for (auto n_it = g.cneighbours_begin(current); n_it != g.cneighbours_end(current); ++n_it){
			if (unvisited.count(n_it->first) > 0 && dist[n_it->first] > dist[current] + g.get_edge_weight(current, n_it->first)) {
					// update distance if better path found
					dist[n_it->first] = dist[current] + g.get_edge_weight(current, n_it->first);
				}
            
        }
		// mark vertex as visited
		unvisited.erase(current);
        
		// min pattern to find current nearest unvisited vertex;
        int min_dist = std::numeric_limits<int>::max();
		for (auto i: unvisited) {
            if ( dist[i] < min_dist ){
                min_dist = dist[i];
                current = i;
                }
            }
		}
	// return distances map	
	return dist;
}


//Returns a vector containing all the articulation points of the
//input weighted graph g.
template <typename vertex>
std::vector<vertex> articulation_points(const weighted_graph<vertex>& g){
	// clone graph, iterate over g and call is_connected
	std::vector<vertex> articulation_points;
	for(auto i = g.begin(); i != g.end(); i++){
		weighted_graph<vertex> cloned_graph = weighted_graph<vertex>(g);
		cloned_graph.remove_vertex(*i);
		// if removing ith vertex breaks graph, add to articulation points
		if (!is_connected(cloned_graph)) {
			articulation_points.push_back(*i);
			}
		}
	return articulation_points;
}

#endif

