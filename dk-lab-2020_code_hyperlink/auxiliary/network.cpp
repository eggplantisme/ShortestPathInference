#define PI 3.14159265

#include <assert.h>
#include <cstdio>
#include <cmath>
#include <iostream>
#include "network.h"
#include "global.h"

//*****************node*constructor********************************************************************
node :: node()
{
    int i;
    ip.resize(4);
    for (i = 0; i < 4; i++)
    {
        ip[i] = -1;
    }
   distance = -1;
//   flag = 0;
   box = -1;
   degree = 0;
   ass_degree = 0;

};

//*****************edge*constructor********************************************************************
edge :: edge()
{
     flag = 0;
     nodes[0] = nodes[1] = -1;
     w = 1;
};


network :: network()
{
    num_nodes = 0;
    num_edges = 0;
    nod = new node[max_nodes];
    edg = new edge[max_edges];
};

int network :: add_node()
{
    if (num_nodes < max_nodes)
    {
        nod[num_nodes].flag = 0;
        nod[num_nodes].degree = 0;
        num_nodes++;
        return num_nodes - 1;
    }
    else
    {
        return -1;  
    }
}

int network :: add_edge(int n1, int n2)
{
    int num;
    int d1, d2;
    assert(n1 >= 0);
    assert(n1 < num_nodes);
    assert(n2 >= 0);
    assert(n2 < num_nodes);
    assert(n1 != n2);
    if ( num_edges < max_edges )
    {
        num = isedge(n1,n2);
        if (num != -1)
        {
            return num;
        }
        else 
        {
            edg[num_edges].flag = 0;
            edg[num_edges].w = 1;
            edg[num_edges].nodes[0] = n1;
            edg[num_edges].nodes[1] = n2;
            num_edges++;
            
            d1 = nod[n1].nodes.size();
            d2 = nod[n2].nodes.size();
            nod[n1].nodes.resize(d1 + 1);
            assert(nod[n1].nodes.size() == d1 + 1);
            
            nod[n1].edges.resize(d1 + 1);
            assert(nod[n1].edges.size() == d1 + 1);
                
            nod[n2].nodes.resize(d2 + 1);
            assert(nod[n2].nodes.size() == d2 + 1);
            
            nod[n2].edges.resize(d2 + 1);
            assert(nod[n2].edges.size() == d2 + 1);
                        
            nod[n1].edges[nod[n1].degree] = num_edges - 1;
            nod[n1].nodes[nod[n1].degree] = n2;
            nod[n2].edges[nod[n2].degree] = num_edges - 1;
            nod[n2].nodes[nod[n2].degree] = n1;
            nod[n1].degree++;
            nod[n2].degree++;
        }
    }
    else
    {
        return -1;
    }
    return num_edges - 1;
}

void network :: grab(char filename[])
{
    
    FILE * ff;
    char str[1001];
    int i;
    int node1, node2, edge;
    std::cout << "grab network from file " << filename << std::endl;
    ff = fopen(filename,"r");
    if (ff == NULL) {
        perror("Failed to open file");
        return;
    }
    assert(fgets(str,1000,ff) != NULL);
    // std::cout << "start network grab." << std::endl;
    erase();
    std::cout << "grab" << std::endl;
    i = 0;
    do
    {
        if ( (i%10000) == 0)
        {
//          std::cout << "String read : "<< i << std::endl;
        }
        i++;
        node1 = get_string_number(1,str);
        node2 = get_string_number(2,str);
       
        assert(node1 >= 0);
        assert(node2 >= 0);
        if ((node1 < max_nodes)&&(node2 < max_nodes))
        {
            if (node1 >= num_nodes)
                do
                {
                    add_node();
                }
                while(node1 >= num_nodes);
                if (node2 >= num_nodes)
                    do
                    {
                        add_node();
                    }
                    while(node2 >= num_nodes);
        
            if (node1 != node2)
            {
                edge = isedge(node1, node2);
                if (edge == -1)
                {
                    edge = add_edge(node1, node2);
                    edg[edge].w = 1;
                }
                else
                {
                    edg[edge].w += 1;
                }
            }
        }
    }
    while(fgets(str,1000,ff) != NULL);
    for (i = 0; i < num_edges; i++)
    {
        edg[i].w = 1/edg[i].w;
    }
    fclose(ff);
    return;
}

void network :: random_link_prediction_split(double q, network& training, network& testing, long * idum)
{
    int i, node1, node2;
    for (i = 0; i < num_nodes; i++)
    {
        training.add_node();
        testing.add_node();
    }
    for (i = 0; i < num_edges; i++)
    {
        node1 = edg[i].nodes[0];
        node2 = edg[i].nodes[1];
        if (ran1(idum)  > q) //output edge
        {
            testing.add_edge(node1, node2);
        }
        else
        {
            training.add_edge(node1, node2);
        }
    }
    return;
}

void network :: fast_get_cls(int node, char temp[])
{
    int j, n_nodes, ctr = 0;
    int node1, node2;
    FILE * tmp;
    tmp = fopen(temp, "w");
    refresh();
//  distance(node);
    distance_no_refresh_out_num_nodes(node);
    std::cout << "save temp " << std::endl;
    tmp = fopen(temp, "w");
    for (j = 0; j < num_nodes; j++)
    {
        if (nod[j].distance != -1)
        ctr++;
    }
    std::cout << " ctr = " << ctr << std::endl;
    //  cin >> j;
    for (j = 0; j < num_edges; j++)
    {
        node1 = edg[j].nodes[0];
        node2 = edg[j].nodes[1];
        if (( nod[node1].distance != -1) && ( nod[node2].distance != -1))
        {
            fprintf(tmp, "%d %d \n", node1, node2);
        }
    }
    fclose(tmp);
    n_nodes = num_nodes;
    grab(temp);
}

void network :: write(char filename[]) // outputs a network in a usual format
{
    int i, j;
    FILE * ff;
    int node1, node2, degree;
    ff = fopen(filename, "w");
    std::cout << "num_nodes" << std::endl;
    for (i = 0; i < num_nodes; i++)
    {
        node1 = i;
        degree = nod[node1].degree;
        for (j = 0; j < degree; j++)
        {
            node2 = nod[node1].nodes[j];
            if (node2 > node1)
            {
                //cout << node1 << " " << node2 << endl;
                fprintf(ff, "%d %d \n", node1, node2);  
            }
        }
    }
    fclose(ff);
}   


void network :: write_w(char filename[]) // outputs a network in a usual format
{
    int i, j;
    FILE * ff;
    double weight;
    int edge;
    int node1, node2, degree;
    ff = fopen(filename, "w");
    std::cout << "num_nodes" << std::endl;
    for (i = 0; i < num_nodes; i++)
    {
        node1 = i;
        degree = nod[node1].degree;
        for (j = 0; j < degree; j++)
        {
            node2 = nod[node1].nodes[j];
            edge = isedge(node1, node2);
            assert(edge >= 0);
            assert(edge < num_edges);
            weight = edg[edge].w;
            if (node2 > node1)
            {
//              std::cout << node1 << " " << node2 << " " << weight << std::endl;
                fprintf(ff, "%d %d  %f\n", node1, node2, weight);   
            }
        }
    }
    fclose(ff);
}   

void network :: erase()
{
    for (int i = 0; i < num_nodes; i++)
    {
        nod[i].nodes.resize(0);
        nod[i].edges.resize(0);
    }
    num_nodes = 0;
    num_edges = 0;
    refresh();
    return;
}

void network :: s1_finite_size_param_infer(double gamma, double T, double tpr, double& nnn, double& kp_max, double& avg_k, double& mu, double& Rh)
{
    
    int i;
    int N_obs = 0, E;
    double k_obs, avg_kp, k_max = 0, pk0, beta;
    double Integral;
    double error;
    //std::cout << "Program estimates parameters of the given network" << std::endl;
    //std::cout << "Program takes four parameters" << std::endl;
    //std::cout << "1) network filename" << std::endl;
    //std::cout << "2) gamma, between 2 and 3" << std::endl;
    //std::cout << "3) 0 < T < 1" << std::endl;
    //std::cout << "4) Link True Positive Rate" << std::endl;
    //std::cout << "gamma = " << gamma << std::endl;
    assert(gamma <=3.0);
    assert(gamma > 2.0);
    beta = 1./T;
    //std::cout << "T = " << T << std::endl;
    assert(beta >= 1.0);;
    //std::cout << "TPR = " << tpr << std::endl;
    assert(tpr > 0);
    assert(tpr <= 1.0);
    
    //std::cout << "Basic parameter estimation" << std::endl;
    
    for (i = 0; i < num_nodes; i++)
    {
        if (nod[i].degree > 0)
        {
            N_obs++;
        }
        if (nod[i].degree > k_max)
        {
            k_max = nod[i].degree;
        }
    }
    E = num_edges;
    std::cout << "N observable = " << N_obs << std::endl;
    std::cout << "E = " << E << std::endl;
    std::cout << "max k  = " << k_max << std::endl;
    k_obs = 2.*E /N_obs;
    std::cout << "<k>_obs = " << k_obs << std::endl;
    avg_kp = (gamma - 1)/(gamma - 2);
    std::cout << "<kp> = " << avg_kp << std::endl;
    Integral = (PI/beta)/(std::sin(PI/beta));
    //std::cout << "I(beta) = " << Integral << std::endl;
    std::cout << "<k> = " << avg_k << std::endl;
    //std::cout << "kp_max estimation" << std::endl;
    double xleft = 10, xright = 100000;
    double ff1, ff2, ff3;
    do
    {
        ff1 = incomplete_negative_gamma_function(xleft, avg_k, avg_kp, tpr, gamma, k_obs, k_max);
        ff2 = incomplete_negative_gamma_function(( xleft + xright)/2., avg_k, avg_kp, tpr, gamma, k_obs, k_max);
        ff3 = incomplete_negative_gamma_function( xright, avg_k, avg_kp, tpr, gamma, k_obs, k_max);

        if (ff1 * ff3 >= 0){
            std::cout << "********** ERROR **********" << std::endl;
            std::cout << "Model parameter inference failed: bisection solver failed to find a solution for kp_max parameter!" << std::endl;
            std::cout << "Try tweaking 'xleft' and 'xright' parameters in 'network.cpp'." << std::endl;
            std::cout << "NOTE: the underlying model assumes the network has a heavy-tailed distribution with the tail exponent gamma. Solver may fail to converge if it is not the case." << std::endl;
            std::cout << "***************************" << std::endl;
        }

        assert(ff1 * ff3 < 0);
        if (ff2 * ff3 < 0)
        {
            xleft = (xleft + xright)/2.;
        }
        else
        {
            xright = (xleft + xright)/2.;
        }
        //std::cout << "xleft = " << xleft << std::endl;
        //std::cout << "xright = " << xright << std::endl;
        //std::cout << "ff = " <<  ff2 << std::endl;
    }
    while (fabs(ff2) > 1e-6);
    kp_max = xleft;
    //std::cout << "kp max = " << kp_max << std::endl;
    //std::cout << "incomplete_negative_gamma_function check = " << incomplete_negative_gamma_function(kp_max, avg_k, avg_kp, tpr, gamma, k_obs, k_max) << std::endl;
    avg_k = finite_avg_degree(k_max, kp_max, avg_kp, tpr, gamma);
    //std::cout << "avg k = " << avg_k << std::endl;
    //std::cout << "finite_alpha = " << finite_alpha(kp_max, gamma) << std::endl;
    pk0 = (gamma - 1) * std::exp((gamma - 1)  * std::log(    (avg_k/ avg_kp)  * tpr * finite_alpha(kp_max, gamma)  ) ) * incomplete_gamma(1-gamma, (avg_k/ avg_kp)  * tpr * finite_alpha(kp_max, gamma)  );
    //std::cout << "Estimated P(0) = " << pk0 << std::endl;
    assert(pk0 >= 0);
    assert(pk0 <= 1.0);
    nnn = N_obs / (1-pk0);
    //std::cout << "true N = " << nnn << std::endl;
    mu = (1./2.)*(avg_k / (avg_kp * avg_kp * Integral));
    //std::cout << "mu = " << mu << std::endl;
    Rh = 2 * std::log (nnn / (PI* mu));
    //std::cout << "Rh = " << Rh << std::endl;
    
    return;
}

int network :: isedge(int n1, int n2)
{
    int i,j;
    int k1, k2; 
    int edge;
    k1 = nod[n1].degree;
    k2 = nod[n2].degree;
    if (k1 < k2)
    {
        for(i = 0; i < k1 ; i++)
        {
            if (nod[n1].nodes[i] == n2)
            {
                edge = nod[n1].edges[i];
                return edge;    
            }
    
        }   
    }
    else
    {
        for(i = 0; i < k2 ; i++)
        {
            if (nod[n2].nodes[i] == n1)
            {
                edge = nod[n2].edges[i];
                return edge;    
            }

        }
    }
    return -1;
}

void network :: refresh()
{
    int i;
    for (i = 0; i < num_nodes; i++)
    {
        nod[i].flag = 0;
//      nod[i].box = -1;
        nod[i].distance = -1;
    };
    for (i = 0; i < num_edges; i++)
    {
        edg[i].flag = 0;
    };
    
    return;
}

  

int network :: distance_no_refresh_out_num_nodes(int nn) // returns the total number of nodes travelled
{ 
    int i,j,k,l;
    int size, degree;
    int node, node1;
    int node_ctr = 1;
    int ctr = 0; // counter of the number of iterations
    //cout << "distance" << endl;
    std::vector <int> list[2]; // keeps the list of current nodes
    // initialization
    list[0].resize(1);
    list[0][0] = nn;
    list[1].resize(0);
    nod[nn].distance = 0;
    //cout << "node[" << nn << "] : " << nod[nn].distance << endl;
    // find neighbors of the current list 
    do
    {
        ctr++;
        k = (ctr + 1) % 2;
        l = (ctr) % 2;
        //  cout << "k = " << k << endl;
        //  cout << "l = " << l << endl;
        list[l].resize(0);
        size = list[k].size();
        for (i = 0; i < size; i++) // goes over curr nodes
        {   
            //      cout << "i = " << i << endl;
            //      cout << "nod[0].distance = " << nod[0].distance << endl;
            node = list[k][i];
            degree = nod[node].degree;
            for (j = 0; j < degree; j++) // goes over connected nodes
            {
                //          cout << "j = " << j << endl;
                node1 = nod[node].nodes[j];
                // if the distance is smaller than one already assigned update distance
                if ((nod[node1].distance > ctr) || (nod[node1].distance == -1))
                {
                    //              cout << "node[" << node1 << "] : " << nod[node1].distance << endl;
                    //              cout << "ctr =" << ctr << endl; 
                    nod[node1].distance = ctr;
                    node_ctr++;
                    list[l].resize(list[l].size() + 1);
                    list[l][ list[l].size() - 1 ] = node1;
                    // update current list  
                    //  cout << "node[" << node1 << "] : " << nod[node1].distance << endl;
                }
            }
        }
        // output current list
        //cout << "current list = ";
        //  cout << endl;
        // add newly acquired nodes to thel list
    }
    while(list[l].size() >= 1);
    return node_ctr;
}
