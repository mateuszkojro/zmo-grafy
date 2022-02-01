#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <unordered_map>
#include <utility>
#include <vector>

#define TIME_SCOPE(msg, content)                                               \
  do {                                                                         \
	auto start = std::chrono::high_resolution_clock::now();                    \
	content auto stop = std::chrono::high_resolution_clock::now();             \
	std::clog << msg << " took: "                                              \
			  << std::chrono::duration_cast<std::chrono::microseconds>(        \
					 stop - start)                                             \
					 .count()                                                  \
			  << "us" << std::endl;                                            \
  } while (0)

using NodeID = size_t;
using NeighbourList = std::vector<NodeID>;

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

 public:
  const NeighbourList& get_neighbours_of(NodeID node_id) {
	return neighbour_list_[node_id];
  }

  bool are_connected(NodeID n_a, NodeID n_b) {
	auto& node_a_neighbours = neighbour_list_.at(n_a);
	auto& node_b_neighbours = neighbour_list_.at(n_b);

	// FIXME: We are sorting so we could use that this and dont need to look two
	// times

	if (std::find(node_a_neighbours.begin(), node_a_neighbours.end(), n_b)
		!= node_a_neighbours.end()) {
	  return true;
	}

	if (std::find(node_b_neighbours.begin(), node_b_neighbours.end(), n_a)
		!= node_b_neighbours.end()) {
	  return true;
	}

	return false;
  }

  void reserve(size_t node_count = 0, size_t edges_per_node = 0) {
	neighbour_list_.reserve(node_count);
	nodes_.reserve(node_count);
  }

  NodeID add_node(double value, NeighbourList initial_connections = {}) {
	NodeID node_id = ++last_node_;
	nodes_.push_back(node_id);
	neighbour_list_[node_id] = std::move(initial_connections);
	return node_id;
  }

  const std::vector<NodeID>& get_nodes() { return nodes_; }

  void add_edge(NodeID node_a, NodeID node_b) {
	size_t idx = 0;
	NodeID target_node = 0;
	if (node_a < node_b) {
	  idx = node_a;
	  target_node = node_b;
	} else {
	  idx = node_b;
	  target_node = node_a;
	}
	auto& neighbours = neighbour_list_[idx];
	if (std::find(neighbours.begin(), neighbours.end(), target_node)
		== neighbours.end()) {
	  neighbours.push_back(target_node);
	}
  }

  void remove_edge(NodeID node_a, NodeID node_b) {
	size_t idx = 0;
	NodeID target_node = 0;
	if (node_a < node_b) {
	  idx = node_a;
	  target_node = node_b;
	} else {
	  idx = node_b;
	  target_node = node_a;
	}
	auto& neighbours = neighbour_list_[idx];
	std::remove(neighbours.begin(), neighbours.end(), target_node);
  }

 private:
  size_t edges_per_node_ = 0;// TODO: Figure out if that caching would be
							 // working
  NodeID last_node_ = 0;
  std::vector<NodeID> nodes_;
  std::unordered_map<NodeID, NeighbourList> neighbour_list_;
};

double random_double() { return static_cast<double>(rand()) / RAND_MAX; }

template<class NodeType>
class Graph {
  bool directed_ = false;

 public:
  //  static std::random_device rd_;
  //  static std::mt19937 random_engine_;

  static Graph* make_erdos_renyi(size_t n, double p) {
	auto* graph = new Graph<double>(n, std::ceil(n * p));
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
  static Graph* make_watts_strogatz(size_t n, size_t k, double beta = 0.1) {
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

	std::uniform_int_distribution<> int_distribution(0, n);
	for (size_t i = 0; i < n; i++) {
	  for (size_t j = i + 1; j < i + k / 2; j++) {
		if (random_double() < beta) {
		  // rewire with uniform distribution and avoid self loops
		  size_t rewire_index = rand();
		  if (rewire_index != j) {
			graph->remove_edge(i, j);
			graph->add_edge(i, rewire_index);
		  }
		}
	  }
	}
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
	content_.reserve(future_nodes_count, future_edges_count);
  }
  void add_edge(NodeID n_a, NodeID n_b) { content_.add_edge(n_a, n_b); }

  NodeID add_node(Node<NodeType> node) { return content_.add_node(0); }
  NodeID add_node(NodeType value) { return add_node(Node<NodeType>(value)); }

  void remove_node(NodeID node_id) { TODO("Not implemented"); }
  void remove_edge(NodeID from, NodeID to) { content_.remove_edge(from, to); }

  std::vector<NodeID> get_neighbours(NodeID node_id) {
	return content_.get_neighbours_of(node_id);
  }

  size_t neighbour_count(NodeID node_id) {
	return get_neighbours(node_id).size();
  }
  size_t degree(NodeID node_id) { return neighbour_count(node_id); }

  size_t distance_between(NodeID from, NodeID to) { return 0; }

  size_t acentricity(NodeID node_id) {
	size_t max_distance = 0;
	for (const auto& node : content_.get_nodes()) {
	  if (size_t d = distance_between(node_id, node); d > max_distance) {
		max_distance = d;
	  }
	}
	return max_distance;
  }

  double global_clustering_coefficient() {
	uint64_t number_of_closed_triplets = 0;
	uint64_t number_of_open_triplets = 0;
	for (NodeID a : content_.get_nodes()) {
	  for (NodeID b : content_.get_nodes()) {
		for (NodeID c : content_.get_nodes()) {
		  if (is_triplet(a, b, c)) {
			if (is_triangle(a, b, c)) {
			  number_of_closed_triplets++;
			}
			number_of_open_triplets++;
		  }
		}
	  }
	}
	return static_cast<double>(number_of_closed_triplets)
		   / number_of_open_triplets;
  }

  size_t count_connections(const std::vector<NodeID>& nodes) {
	size_t number_of_connections = 0;
	for (const auto& a : nodes) {
	  for (const auto& b : nodes) {
		const auto& neighbours_of_a = content_.get_neighbours_of(a);
		if (std::find(neighbours_of_a.begin(), neighbours_of_a.end(), b)
			!= neighbours_of_a.end()) {
		  number_of_connections++;
		}
	  }
	}
	return number_of_connections;
  }

  double node_clustering_coefficient(NodeID node_id) {
	const NeighbourList& neighbour_list = content_.get_neighbours_of(node_id);
	size_t k = neighbour_list.size();

	if (k == 0 || k == 1)
	  return 0;

	size_t number_of_links_between_neighbours =
		count_connections(neighbour_list);
	return static_cast<double>(number_of_links_between_neighbours)
		   / (k * (k - 1));
  }

  double average_clustering_coefficient() {
	double sum_of_coeficients = 0;
	for (const NodeID& node : content_.get_nodes())
	  sum_of_coeficients += node_clustering_coefficient(node);

	return sum_of_coeficients / content_.get_nodes().size();
  }

  bool is_triangle(NodeID n_a, NodeID n_b, NodeID n_c) {
	assert(n_a != n_b && n_b != n_c && n_c != n_a
		   && "Triangle from the same nodes does not make sense");

	bool a_conected_b = content_.are_connected(n_a, n_b);
	bool b_conected_c = content_.are_connected(n_b, n_c);
	bool c_conected_a = content_.are_connected(n_c, n_a);

	return a_conected_b || b_conected_c || c_conected_a;
  }

  bool is_triplet(NodeID n_a, NodeID n_b, NodeID n_c) {
	assert(n_a != n_b && n_b != n_c && n_c != n_a
		   && "Triplet from the same nodes does not make sense");
	bool a_conected_b = content_.are_connected(n_a, n_b);
	bool b_conected_c = content_.are_connected(n_b, n_c);
	bool c_conected_a = content_.are_connected(n_c, n_a);

	if (a_conected_b || c_conected_a) {
	  return true;
	}

	if (b_conected_c || a_conected_b) {
	  return true;
	}

	if (c_conected_a || b_conected_c) {
	  return true;
	}

	return false;
  }

  size_t radius() {
	size_t min_accentricity = SIZE_MAX;
	for (const auto& node : content_.get_nodes()) {
	  size_t a = acentricity(node);
	  if (a < min_accentricity) {
		min_accentricity = a;
	  }
	}
	return min_accentricity;
  }

  size_t diameter() {
	size_t max_distance = 0;
	for (const auto& node_a : content_.get_nodes()) {
	  for (const auto& node_b : content_.get_nodes()) {
		size_t d = distance_between(node_a, node_b);
		if (d > max_distance) {
		  max_distance = d;
		}
	  }
	}
	return max_distance;
  }

 private:
  AdjacencyList content_;
};

int main() {
  size_t n = 25'000;

  TIME_SCOPE("barabasi albert", {
	std::clog << "Making barabasi albert:";
	auto* ba = Graph<double>::make_barabasi_albert(n);
	std::clog << " done" << std::endl;
  });

  TIME_SCOPE("erdos renyi", {
	std::clog << "Making erdos renyi:";
	auto* er = Graph<double>::make_erdos_renyi(n, 0.1);
	std::clog << " done" << std::endl;
  });

  Graph<double>* ws;
  TIME_SCOPE("watts srogatz", {
	std::clog << "Making watts srogatz:";
	ws = Graph<double>::make_watts_strogatz(n, 10);
	std::clog << " done" << std::endl;
  });

  double clustering_coef = 0;
  TIME_SCOPE("Clustering coef",
			 { clustering_coef = ws->average_clustering_coefficient(); });
  std::cout << "Average clustering coef: " << clustering_coef << std::endl;
  return 0;
}
