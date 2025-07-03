#define PI 3.14159265

#include "mle.h"
#include "global.h"
#include <vector>
#include <cmath>
#include <assert.h>
#include <cstdio>
#include <iostream>

double global_h2_find_node_coords(int node1, std::vector<int>& id, network& nw, std::vector<double>& angles, std::vector<double>& real_radii, int nnn, double Rh, double beta, double tpr, int grid_size)
{
    // NOT SUITABLE FOR PARALLEL CODE
    double max_theta, xhi, angle, fine_angle;
    double R = nnn / (2*PI), L, L_max, p;
    int i, j, t, node2;
    
    assert(grid_size > 0);
    
    max_theta = -100;
    L_max = -1e100;
    
    for (j = 0; j < id.size(); j++)
    {
        node2 = id[j];
        if (nw.isedge(node1, node2) != -1)
        {
            nw.nod[node2].box = 1; // pre-mark nodes in the club that are connected to node1
        }
    }
    
    
    for (j = 0; j < grid_size; j++)    // local optimization
    {
        angle = (2*PI*j)/ (grid_size) - PI;
        L = 0;
        for (i = 0; i < id.size(); i++)
        {
            node2 = id[i];
            assert(nw.nod[node2].flag == 1);
            xhi = distH2approx(real_radii[node1], angle, real_radii[node2],angles[node2], 1.0);
            p = tpr / (1. + std::exp(beta * (xhi - Rh)/2.) );
            if (nw.nod[node2].box ==1)
            {
                L += std::log(p);
            }
            else
            {
                L += std::log(1. -p);
            }
        }
        if (L > L_max)
        {
            L_max = L;
            max_theta = angle;
        }
    }
    
    for (j = 0; j < id.size(); j++)
    {
        node2 = id[j];
        nw.nod[node2].box = 0; // release labellings
    }
    return max_theta;
}

double h2_find_node_coords(int node1, std::vector<int>& id, network& nw, std::vector<double>& angles, std::vector<double>& real_radii, int nnn, double Rh, double beta, double tpr, int grid_size)
{
    // NOT SUITABLE FOR PARALLEL CODE
    double max_theta, xhi, angle, fine_angle;
    double R = nnn / (2*PI), L, L_max, p;
    int i, j, t, node2;
    
    assert(grid_size > 0);
    
    max_theta = -100;
    L_max = -1e100;
    
    for (j = 0; j < grid_size; j++)    // local optimization
    {
        angle = (2*PI*j)/ (grid_size) - PI;
        L = 0;
        for (i = 0; i < nw.nod[node1].degree; i++)
        {
            node2 = nw.nod[node1].nodes[i];
            if (nw.nod[node2].flag == 1) // flag=1 means node is in the club
            {
                xhi = distH2approx(real_radii[node1], angle, real_radii[node2],angles[node2], 1.0);
                p = tpr / (1. + std::exp(beta * (xhi - Rh)/2.) );
                L += std::log(p);
            }
        }
        if (L > L_max)
        {
            L_max = L;
            max_theta = angle;
        }
    }
    
    
    for (j = 0; j < id.size(); j++)
    {
        node2 = id[j];
        if (nw.isedge(node1, node2) != -1)
        {
            nw.nod[node2].box = 1; // pre-mark nodes in the club that are connected to node1
        }
    }
    fine_angle = max_theta;
    L_max = -1e100;
    for (t = - grid_size*300./nnn; t < grid_size * 300./nnn; t++)
    {
        angle = max_theta + ((2*PI*t)/(grid_size));
        if (angle > PI)
        {
            do
            {
                angle -= 2 * PI;
            }
            while(angle > PI);
        }
        if (angle < -PI)
        {
            do
            {
                angle += 2* PI;
            }
            while(angle < -PI);
        }
        assert(angle <= PI);
        assert(angle >= -PI);
        L =    0;
        for (j = 0; j < id.size(); j++) // global optmization
        {
            node2 = id[j];
            xhi = distH2approx(real_radii[node1], angle, real_radii[node2],angles[node2], 1.0);
            p = tpr / (1. + std::exp(beta * (xhi - Rh)/2.) );
            if(nw.nod[node2].box == 1) // node is connected to the node in the club
            {
                L += std::log(p);
            }
            else
            {
                L += std::log(1. - p);
            }
        }
        if (L_max < L)
        {
            L_max = L;
            fine_angle = angle;
        }
    }
    for (j = 0; j < id.size(); j++)
    {
        node2 = id[j];
        nw.nod[node2].box = 0; // release labellings
    }
    
    
    return fine_angle;
}

double global_h2_find_core_coords_single_round(std::vector<int>& id, network& nw, std::vector<double>& angles, std::vector<double>& radii, std::vector<double>& real_angles, int nnn, double Rh, double beta, double error_value, double noise, long * idum, double tpr, int grid_size) // returns Likelihood
{
    double max_error, error;
    double L;
    double angle;
    int ctr, i, j, k, node1, node2;
    max_error = 0;
    ctr = 0;
    for (i = 0; i < id.size(); i++)
    {
        node1 = id[i];
        angles[node1] = noise*(ran1(idum) - 0.5) + angles[node1];
        if (angles[node1] > PI)
        {
            angles[node1] -= 2*PI;
        }
        if (angles[node1] < -PI)
        {
            angles[node1] += 2*PI;
        }
    }
    ctr = 0;
    do
    {
        ctr++;
        max_error = 0;
        for (i = 0; i < id.size(); i++)
        {
            
            node1 = id[i];
            angle = global_h2_find_node_coords(node1, id, nw, angles, radii, nnn, Rh, beta, tpr, grid_size);
            error = distS1(angle,angles[node1]);
            if (error > max_error)
            {
                max_error = error;
            }
            angles[node1] = angle;
        }
        
    }
    while((max_error > error_value)&&(ctr < 10));
    L = h2loglikelihood(id, nw, angles, radii, nnn,  Rh, beta, tpr);
    return L;
}

double h2_find_core_coords_single_round(std::vector<int>& id, network& nw, std::vector<double>& angles, std::vector<double>& radii, std::vector<double>& real_angles, int nnn, double Rh, double beta, double error_value, double noise, long * idum, double tpr, int grid_size) // returns Likelihood
{
    double max_error, error;
    double L;
    double angle;
    int ctr, i, j, k, node1, node2;
    max_error = 0;
    ctr = 0;
    for (i = 0; i < id.size(); i++)
    {
        node1 = id[i];
        angles[node1] = noise*(ran1(idum) - 0.5) + angles[node1];
        if (angles[node1] > PI)
        {
            angles[node1] -= 2*PI;
        }
        if (angles[node1] < -PI)
        {
            angles[node1] += 2*PI;
        }
    }
    ctr = 0;
    do
    {
        ctr++;
        max_error = 0;
        for (i = 0; i < id.size(); i++)
        {
            
            node1 = id[i];
            angle = h2_find_node_coords(node1, id, nw, angles, radii, nnn, Rh, beta, tpr, grid_size);
            error = distS1(angle,angles[node1]);
            if (error > max_error)
            {
                max_error = error;
            }
            angles[node1] = angle;
        }
        
    }
    while((max_error > error_value)&&(ctr < 10));
    L = h2loglikelihood(id, nw, angles, radii, nnn,  Rh, beta, tpr);
    return L;
}

double h2loglikelihood(std::vector <int>& id, network& nw, std::vector <double>& angles, std::vector <double>& real_radii, int nnn,  double Rh, double beta, double tpr)
{
    int i, j, node1, node2;
    
    double R = nnn / (2*PI), p, xhi, L = 0;
    for (i = 0; i < id.size(); i++)
    {
        node1 = id[i];
        for (j = i + 1; j < id.size(); j++)
        {
            node2 = id[j];
            xhi = distH2(real_radii[node1], angles[node1], real_radii[node2], angles[node2], 1.0);
            assert(xhi >= 0);
            p = tpr / (1. + std::exp(beta * (xhi - Rh)/2.) );
            if (nw.isedge(node1,node2) == -1)
            {
                L += std::log(1-p);
            }
            else
            {
                L += std::log (p);
            }
        }
    }
    return L;
}

void h2_find_remaining_coords(std::vector<int>& id, network& nw_old, std::vector<double>& a_list, std::vector<double>& r_list, std::vector<double>& real_a_list, int nnn, double Rh, double T, double error_value, int num_iter, double max_noise, long * idum, double tpr, char * out_filename) // returns Likelihood
{
    int size, ctr, i, j, max_k, max_node, max_index;
    int node1, node2;
    std::vector<int> embed_list;
    FILE * ff;
    
    ff = fopen(out_filename, "w");
    // doublecheck
    for (i = 0; i < nw_old.num_nodes; i++)
    {
        if (nw_old.nod[i].degree >= 1)
        {
            fprintf(ff,"%d\t%d\t%f\t%f\t%f\n", i, nw_old.nod[i].degree,  r_list[i], a_list[i], real_a_list[i]);
        }
    }
    fclose(ff);
    std::cout << "embed remaining nodes" << std::endl;
    
    
    ctr = 0;
    for (i = 0; i < nw_old.num_nodes; i++)
    {
        if ((nw_old.nod[i].degree >= 1.0) && (nw_old.nod[i].flag == 0))
        {
            embed_list.resize(ctr + 1);
            embed_list[ctr] = i;
            ctr++;
        }
    }
    
    std::cout << "embed list compiled" << std::endl;
    std::cout << "embed list size = " << ctr << std::endl;
    
    // make sure we have some k = 1 nodes
    if (embed_list.size() != 0){
        do
        {
            max_k = 0;
            max_node = embed_list[0];
            for (i = 0; i < embed_list.size(); i++)
            {
                node1 = embed_list[i];
                ctr = 0;
                for (j = 0; j < nw_old.nod[node1].degree; j++)
                {
                    node2 = nw_old.nod[node1].nodes[j];
                    if (nw_old.nod[node2].flag == 1)
                    {
                        ctr++; // k to the core
                    }
                }
                if (max_k < ctr)
                {
                    max_k = ctr;
                    max_node = node1;
                    max_index =  i;
                }
            }
            
            a_list[max_node] = h2_find_node_coords(max_node, id, nw_old, a_list, r_list, nnn, Rh, 1./T, tpr, nnn);
            
            
            size = embed_list.size();
            embed_list[max_index] = embed_list[size - 1];
            embed_list.resize(size - 1);
            // add the node to the core/
            
            nw_old.nod[max_node].flag = 1;
            size = id.size();
            id.resize(size  + 1);
            id[size] = max_node;
            if (embed_list.size() % 100 == 0)
            {
                std::cout << "remaining num to embed = " << embed_list.size() << std::endl;
            }
        }
        while(embed_list.size() >= 1);
    }
    
    
    
    
    // output resulting nodes
    std::cout << "output resulting angles" << std::endl;
    std::cout << "output filename = " << out_filename << std::endl;
    ff = fopen(out_filename, "w");
    // doublecheck
    for (i = 0; i < nw_old.num_nodes; i++)
    {
        if (nw_old.nod[i].degree >= 1)
        {
            fprintf(ff,"%d\t%d\t%f\t%f\t%f\n", i, nw_old.nod[i].degree,  r_list[i], a_list[i], real_a_list[i]);
        }
    }
    fclose(ff);

    return;
}

void h2_massage_all(std::vector<int>& id, network& nw_old, std::vector<double>& a_list, std::vector<double>& r_list, std::vector<double>& real_a_list, int nnn, double Rh, double T, double error_value, int num_iter, long * idum, double tpr, int grid_size, char * out_filename) // returns Likelihood
{
    int i, j, node1;
    double noise, L;
    FILE * ff;
    std::cout << "error value = " << error_value << std::endl;
    std::cout << "core size = " << id.size() << std::endl;
    std::cout << "start core massaging..." << std::endl;
    for (i = 0; i < num_iter; i++)
    {
        noise = 0;
        shuffle_vector(id, idum);
        L = h2_find_core_coords_single_round_progress(id, nw_old, a_list, r_list, real_a_list, nnn, Rh, 1./T, error_value, noise, idum, tpr, grid_size);
        std::cout << "iter = " << i << "; L = " << L << std::endl;
        ff = fopen(out_filename, "w");
        for (j = 0; j < id.size(); j++)
        {
            node1 = id[j];
            fprintf(ff, "%d\t%f\t%f\t%f\t%d\n", node1, a_list[node1], real_a_list[node1], r_list[node1], nw_old.nod[node1].degree);
        }
        fclose(ff);
    }
    
    return;
}

double  h2_find_core_coords_single_round_progress(std::vector<int>& id, network& nw, std::vector<double>& angles, std::vector<double>& radii, std::vector<double>& real_angles, int nnn, double Rh, double beta, double error_value, double noise, long * idum, double tpr, int grid_size) // returns Likelihood
{
    double max_error, error;
    double L;
    double angle;
    int ctr, i, j, k, node1, node2;
    max_error = 0;
    ctr = 0;
    for (i = 0; i < id.size(); i++)
    {
        node1 = id[i];
        angles[node1] = noise*(ran1(idum) - 0.5) + angles[node1];
        if (angles[node1] > PI)
        {
            angles[node1] -= 2*PI;
        }
        if (angles[node1] < -PI)
        {
            angles[node1] += 2*PI;
        }
    }
    ctr = 0;
    {
        ctr++;
        max_error = 0;
        for (i = 0; i < id.size(); i++)
        {
            
            node1 = id[i];
            angle = h2_find_node_coords(node1, id, nw, angles, radii, nnn, Rh, beta, tpr, grid_size);
            error = distS1(angle,angles[node1]);
            if (error > max_error)
            {
                max_error = error;
            }
            angles[node1] = angle;
        }
        std::cout << "max_error = " << max_error << std::endl;
    }
    
    L = h2loglikelihood(id, nw, angles, radii, nnn,  Rh, beta, tpr);
    return L;
}