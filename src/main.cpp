#include <cassert>
#include <iostream>
#include <stdint.h>
#include <unordered_map>
#include <vector>

using NodeID = size_t;

#define TODO(msg) assert(false && msg)

template<class Type>
class Node {
 public:
  explicit Node(Type value) : value_(value) {}
  Type get_value() const { return value_; }
  void set_value(Type value) { value_ = value; }

 private:
  Type value_;
};

struct Edge
{
  NodeID start;
  NodeID end;
};

class AdjacencyList {
  using NeighbourList = std::vector<NodeID>;

 public:
  const NeighbourList& get_neighbours_of(NodeID node_id) {
	return neighbour_list_.at(node_id);
  }

 private:
  std::unordered_map<NodeID, NeighbourList> neighbour_list_;
};

double random_double() { return static_cast<double>(rand()) / RAND_MAX; }

template<class NodeType>
class Graph {
  bool directed = false;

 public:
  static Graph* make_erdos_renyi(size_t n, double p) {
	auto* graph = new Graph<double>(n, std::ceil(n * p));
	// OPTIMIZE: Its O(n+n^2)
	for (size_t i = 0; i < n; i++) { graph->add_node(i); }
	for (size_t i = 0; i < n; i++) {
	  for (size_t j = 0; i < n; i++) {
		if (random_double() < p) {
		  graph->add_edge(i, j);
		}
	  }
	}
	return graph;
  }
  static Graph* make_watts_strogatz(size_t n, size_t k) {
	auto* graph = new Graph<double>(n, (n * k) / 2);
	// OPTIMIZE: Its O(n+n^2)
	for (size_t i = 0; i < n; i++) { graph->add_node(i); }
	for (int64_t i = 0; i < n; i++) {
	  for (int64_t j = 0; i < n; i++) {
		int64_t coef = std::abs(i - j) % (n - 1 - k / 2);
		if (coef > 0 && coef < k / 2) {
		  graph->add_edge(i, j);
		}
	  }
	}
	TODO("There is a step 2 missing");
	for (size_t i = 0; i < n; i++) { graph->add_node(i); }
	return graph;
  }

  static Graph* make_barabasi_albert(size_t n) {
	auto* graph = new Graph<double>(n, n);
	for (size_t i = 0; i < n; i++) {
	  size_t sum_k = 0;
	  for (size_t j = 0; j < n; j++) { sum_k += graph->degree(j); }
	  double p_i = static_cast<double>(graph->degree(i)) / sum_k;
	  size_t current = graph->add_node(i);
	  if (random_double() < p_i) {
		graph->add_edge(current, i);
	  }
	}
	return graph;
  }

  explicit Graph(size_t future_nodes_count = 0, size_t future_edges_count = 0) {
	nodes_.reserve(future_nodes_count);
	edges_.reserve(future_edges_count);
  }
  void add_edge(NodeID n_a, NodeID n_b) {
	TODO("Check if this edge is already in the graph and if permutation is "
		 "already in graph");
	edges_.push_back(Edge{n_a, n_b});
  }

  NodeID add_node(Node<NodeType> node) {
	nodes_.push_back(node);
	return nodes_.size() - 1;
  }
  NodeID add_node(NodeType value) { return add_node(Node<NodeType>(value)); }
  void remove_node(NodeID node_id) { TODO("Not implemented"); }
  void remove_edge(NodeID from, NodeID to) { TODO("Not implemented"); }

  std::vector<NodeID> get_neighbours(NodeID node_id) {
	std::vector<NodeID> neighbours;

	// TODO: That could and should be cached
	for (auto& edge : edges_) {
	  if (edge.start == node_id) {
		neighbours.push_back(edge.end);
	  }
	  if (!directed && edge.start == node_id) {
		neighbours.push_back(edge.start);
	  }
	}
	return neighbours;
  }

  size_t neighbour_count(NodeID node_id) {
	return get_neighbours(node_id).size();
  }
  size_t degree(NodeID node_id) { return neighbour_count(node_id); }

  double global_clustering_coefficient();
  double node_clustering_coefficient(NodeID node_id);
  double average_clustering_coefficient();

  size_t radius();
  size_t diameter();

 private:
  std::vector<Node<NodeType>> nodes_;
  std::vector<Edge> edges_;
};

int main() {
  size_t n = 10'000;

  std::clog << "Making barabasi albert:";
  auto* ba = Graph<double>::make_barabasi_albert(n);
  std::clog << " done" << std::endl;
  std::clog << "Making watts srogatz:";
  auto* ws = Graph<double>::make_watts_strogatz(n, 10);
  std::clog << " done" << std::endl;
  std::clog << "Making erdos renyi:";
  auto* er = Graph<double>::make_erdos_renyi(n, 0.1);
  std::clog << " done" << std::endl;

  return 0;
}
