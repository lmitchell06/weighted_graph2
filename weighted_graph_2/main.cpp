#include <iostream>
#include "weighted_graph.hpp"
#include "graph_algorithms.cpp"
//#include "tests/test_helper.cpp"

int main(){
	/*weighted_graph<int> g;
		
		auto r = (std::rand()%20) + 5;
		
		for (auto i = 0; i < r; ++i){
			g.add_vertex(i);
		}
		
		std::vector<int> vertices(g.begin(), g.end());
		
		auto min_edges = random_tree(vertices);
		
		for (auto e : min_edges){
			g.add_edge(e.first, e.second, (std::rand()%10) + 1);
		}
		
		auto shortest_path_tree_weight = g.total_weight();
		
		auto start_vertex = vertices[std::rand()%vertices.size()];
		
		auto actual_distances = compute_tree_distances(g, start_vertex);
				
		auto extra_edges = random_tree(vertices);
		
		for (auto e : extra_edges){
			if (!g.are_adjacent(e.first, e.second)){
				g.add_edge(e.first, e.second, (std::rand()%10) + shortest_path_tree_weight + 1);
			}
		}
		
		auto computed_shortest_paths = dijkstras(g, start_vertex);
		
		std::cout << (actual_distances, computed_shortest_paths) << std::endl;
		
		for (auto p : actual_distances){
			std::cout << (computed_shortest_paths[p.first], p.second) << std::endl;
		}*/
		
}
