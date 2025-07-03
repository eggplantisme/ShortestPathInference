#define max_nodes 2000001
#define max_edges 14000000

#ifndef NETWORK_H
#define NETWORK_H

#include <vector>

//**********************************************************************************************
//**********************************************************************************************
//*************************class for networks*dynamic*definitions*******************************
//**********************************************************************************************
//**********************************************************************************************
class node;
class edge;
class network;

//**********************************************************************************************
//*************************node*****************************************************************
//**********************************************************************************************
class node 
{
public:
   int degree;
   int ass_degree; // assigned degree by SF generation constructor
   std::vector<int> edges;      
   std::vector<int> nodes;
   std::vector<int> sim;
   std::vector<int> ip;	
   int flag;
   int box; // contains the number of box this node belongs to
   double distance; // distance form the origin for optimal path
   double centrality;
   node();

};

//***********************************************************************************************
//*************************edge*****************************************************************
//**********************************************************************************************
class edge
{
public:
   int flag;
   int nodes[2];
   double w; //
   edge();};


//***********************************************************************************************
//*************************network**************************************************************
//**********************************************************************************************
class network
{
public:
	int num_nodes, num_edges;
	node * nod;
	edge * edg;

	network();              // empty nw

	void grab(char ch[]);
	void random_link_prediction_split(double q, network& training, network& testing, long * idum);
	void fast_get_cls(int node, char temp[]); // gets cluster with using a temporary file; the result is the network with many zero degree nodes
	void write(char filename[]); // outputs a network in a usual format
	void erase();// erases network 
	void s1_finite_size_param_infer(double gamma, double T, double tpr, double& nnn, double& kp_max, double& avg_k, double& mu, double& Rh);
	int isedge(int, int); //returns the edge number if 2 nodes are by a single edge otherwise returns -1;
	int add_node(); // returns the number of nodes, if fail returns -1;
	int add_edge(int,int); // returns -1 if can't add an edge, number of edge otherwise
	void refresh(); // clears all the flags on nodes and edges;
	void write_w(char filename []); // outputs a nw with weights on the links
    int distance_no_refresh_out_num_nodes(int nn);

 };

#endif