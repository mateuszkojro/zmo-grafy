#include <cassert>
#include <iostream>
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

template<class NodeType>
class Graph {
 public:
  static Graph* make_erdos_renyi();
  static Graph* make_watts_strogatz();
  static Graph* make_barabasi_albert();

 public:
  Graph();
  void add_edge(NodeID n_a, NodeID n_b) {
	TODO("Check if this edge is already in the graph");
	edges_.push_back(Edge{n_a, n_b});
  }
  NodeID add_node(Node<NodeType> node) { nodes_.push_back(node); }
  NodeID add_node(NodeType value) { add_node(Node<NodeType>(value)); }

  std::vector<Node<NodeType>> get_neighbours(NodeID node_id) {}
  size_t neighbour_count(NodeID node_id) {
	return get_neighbours(node_id).size();
  }

  double clustering_coefficient();

  size_t radius();
  size_t diameter();

 private:
  std::vector<Node<NodeType>> nodes_;
  std::vector<Edge> edges_;
};

int main() {
  std::cout << "Hello world!" << std::endl;
  return 0;
}
