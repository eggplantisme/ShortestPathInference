#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <assert.h>
#include <iterator>
#include <string.h>
#include <sstream>
#include <chrono>


#define max_nodes 2000001
#define max_edges 14000000
#define max_degree 100000
#define max_cent 100000000
#define max_distance 1e48
#define max_shells 200
#define PI 3.14159265
#define ee 2.71828

#include "./auxiliary/network.h"
#include "./auxiliary/global.h"
#include "./auxiliary/mle.h"

inline double conn_prob(double dist, double mu, double T)
{
    return 1./(1.+ std::exp( (dist - mu) / (2*T) ));
    
}


double dist(double x1, double x2)
{
    return PI- std::fabs(PI - fabs(x1 - x2));
}




int main(int argc, char *argv[])
{
    
    int iter, num_iter = 100, mid_avg;
    double max_noise = PI/4;
    int core_size1 = 200;
    int core_size2 = 1000;
    
    int grid_mult = 1;
    int i, j, k, ctr, node_ctr, N, best_iter;
    double R0, Rh, gamma, avg, T, pp, nnn, mu, avg_k, kp_max, avg_kp, avg_r2, k_cut1, error_value, noise, L, Lmax, x, min, max, mult;
    double a, a1, a2, a3;  // angles
    int max_k, node1, node2, max_node, size, max_index;
    FILE * ff;
    
    long idum;
    int num_avg = 20;
    std::vector<int> id, embed_list, k_list, knode_list;
    
    std::vector <double> best_coords, real_r_list, real_a_list, r_list, a_list, kp_list, score_list, true_list, invcosh;
    network nw_old;
    clock_t t1, t2;

    char* edgelist_filename;
    
    //************************************************************************************************
    std::cout << "Initialization ..." << std::endl;;
    if (argc > 0)
    {
        if (argc == 9)
        {
            
            
            edgelist_filename = basename(argv[1]);
            std::cout << "Here is the filename: " << edgelist_filename << std::endl;

            nw_old.grab(argv[1]);
            N = nw_old.num_nodes;
                        std::cout << "N = " << N << std::endl;
            gamma = atof(argv[2]);
                        std::cout << "gamma = " << gamma << std::endl;
            avg_kp = (gamma - 1)/(gamma - 2);
                        std::cout << "avg_kp = " << avg_kp << std::endl;
            T = atof(argv[3]);
                        std::cout << "T = " << T << std::endl;
            idum = atoi(argv[4]);
                        std::cout << "idum = " << idum << std::endl;
            pp = 1. - atof(argv[6]);
                        std::cout << "q = " << pp << std::endl;
            assert(pp >= 0);
            assert(pp <= 1.0);
            num_avg = atoi(argv[7]);
                        std::cout << "num avg = " << num_avg << std::endl;
            mid_avg = roundx(0.5*num_avg);
            
            grid_mult = atoi(argv[8]);
                        std::cout << "grid size multiplier coefficient = " << grid_mult << std::endl;
            assert(num_avg >= 1);
        }
        else
        {
            std::cout << std::endl;
            std::cout << "Program embeds a given network with 1-q randomly removed links" << std::endl;
            std::cout << "Program consists of the following steps:" << std::endl;
            std::cout << "A. Extract GCC of the network" << std::endl;
            std::cout << "B. Infer all model parameters except gamma, T and q" << std::endl;
            std::cout << "C. Infer radial coordinates" << std::endl;
            std::cout << "D. Infer angular core coordinates" << std::endl;
            std::cout << "E. Infer remaining node coordinates" << std::endl;
            std::cout << "**********" << std::endl;
            std::cout << "Program takes following parameters:" << std::endl;
            
            for (i = 0; i < argc; i++)
            {
                std::cout << argv[i] << std::endl;
            }
            
            std::cout << "1) network filename" << std::endl;
            std::cout << "2) gamma" << std::endl;
            std::cout << "3) T temperature" << std::endl;
            std::cout << "4) random seed" << std::endl;
            std::cout << "5) coords filename " << std::endl;
            std::cout << "6) fraction of removed links " << std::endl;
            std::cout << "7) number of layers" << std::endl;
            std::cout << "8) grid size multiplier coefficient" << std::endl;
            std::cout << "Exiting program" << std::endl;
            
            
            return -1;
        }
        
        //************************************************************************************************
        std::cout << "A. Extract GCC of the network" << std::endl;
        max_k = 0;
        node1 = -1;
        for (i = 0; i < nw_old.num_nodes; i++)
        {
            if(max_k < nw_old.nod[i].degree)
            {
                max_k = nw_old.nod[i].degree;
                node1 = i;
            }
        }
        assert(node1 >= 0);
        
        
        // get gcc and save it to a uniquely named tmp file
        const auto current_time = std::chrono::system_clock::now();
        auto seconds_since_epoch = std::chrono::duration_cast<std::chrono::seconds>(current_time.time_since_epoch()).count();
        std::stringstream ss;
        ss << "./" << seconds_since_epoch << "_" << edgelist_filename << "_" << "tmp_gcc.net";
        char* temp_file_path = new char[ss.str().length() + 1];
        ss >> temp_file_path;
        std::cout << "Saving temporary largest connected component to: " << temp_file_path << std::endl;
        nw_old.fast_get_cls(node1, temp_file_path);
        nw_old.write(temp_file_path);
        nw_old.grab(temp_file_path);

        std::cout << "B. Infer all model parameters except gamma, T and q" << std::endl;
    
        nw_old.s1_finite_size_param_infer(gamma, T, pp, nnn, kp_max, avg_k, mu, Rh);
        std::cout << "inferred N = " << nnn << std::endl;
        std::cout << "inferred kp_max = " << kp_max << std::endl;
        std::cout << "inferred avg_k = " << avg_k << std::endl;
        std::cout << "inferred mu = " << mu << std::endl;
        std::cout << "inferred Rh = " << Rh << std::endl;
        //************************************************************************************************
        std::cout << "C. Infer radial coordinates" << std::endl;
        kp_list.resize(N);
        r_list.resize(N);
        real_r_list.resize(N);
        real_a_list.resize(N);
        
        for (i = 0; i < N; i++)
        {
            real_r_list[i] = -1.;
            real_a_list[i] = -100.;
        }
        
        avg_r2 = avg_kp * std::exp(-Rh/2.);  // <exp(-r/2)>
        
        
        for (i = 0; i < N; i++)
        {
            r_list[i] = Rh;
            kp_list[i] = 1.;
        }
        for (i = 0; i < nw_old.num_nodes; i++)
        {
            if (nw_old.nod[i].degree > T * gamma)
            {
                kp_list[i] = nw_old.nod[i].degree - T* gamma;
                kp_list[i] /= (pp * finite_alpha(kp_max, gamma) * avg_k);
                kp_list[i] *= avg_kp;
                r_list[i] = - 2. * std::log( avg_r2 *(nw_old.nod[i].degree + (gamma - 1) * T ) / ( pp * finite_alpha(kp_max, gamma) *avg_k)  );
            }
        }
        
        //************************************************************************************************
        std::cout << "D. Infer angular core coordinates" << std::endl;
        
        
        k_list.resize(0);
        knode_list.resize(0);
        ctr = 0;
        error_value = 0.0001;
        for (i = 0; i < nw_old.num_nodes; i++)
        {
            if (nw_old.nod[i].degree >= 2)
            {
                knode_list.resize(ctr + 1);
                k_list.resize(ctr + 1);
                knode_list[ctr] = i;
                k_list[ctr] = nw_old.nod[i].degree;
                ctr++;
            }
        }
        std::cout << "total nodes for the core = " << ctr << std::endl;
        quickSort(k_list, knode_list);
        
        std::cout << "initializing angles" << std::endl;
        a_list.resize(N);
        
        Lmax = -1e100;
        
        for (iter = 0; iter < 20; iter ++)
        {
            std::cout << "Replica " << iter+1 << std::endl;
            std::cout << "**********************" << std::endl;
            for(i = 0; i < N; i++)
            {
                a_list[i] = -100;
            }
            id.resize(0);
            min = 0;
            max = 20;
            mult = std::exp( (1./double(num_avg)) * std::log(double(knode_list.size() / 20. ) ));
            ctr = 0;
            node_ctr = 0;
            for (i = 1; i < mid_avg; i++)  // initial initializations repeat 10 times choose the highest likelihood
            {
                for (j = roundx(min); j < roundx(max); j++)
                {
                    ctr++;
                    node1  = knode_list[j];
                    if (i == 1)
                    {
                        a_list[node1] = 2*PI * (ran1(&idum) - 0.5);
                    }
                    else
                    {
                        if ((id.size() < 500)||(nnn < 1500))
                        {
                            a_list[node1] = global_h2_find_node_coords(node1, id, nw_old, a_list, r_list, nnn, Rh, 1./T, pp, grid_mult*id.size());
                        }
                        else
                        {
                            a_list[node1] = h2_find_node_coords(node1, id, nw_old, a_list, r_list, nnn, Rh, 1./T, pp, grid_mult*id.size());
                        }
                    }
                    id.resize(node_ctr + 1);
                    id[node_ctr] = node1;
                    node_ctr++;
                    nw_old.nod[node1].flag = 1;
                }
                // embed nodes single round
                noise = (i)*(0.01-max_noise)/num_avg + max_noise;
                shuffle_vector(id, &idum);
                if ((id.size() < 500)||(nnn < 1500))
                {
                    L = global_h2_find_core_coords_single_round(id, nw_old, a_list, r_list, real_a_list, nnn, Rh, 1./T, error_value, noise, &idum, pp, grid_mult*id.size()); // returns Likelihood
                }
                else
                {
                    L = h2_find_core_coords_single_round(id, nw_old, a_list, r_list, real_a_list, nnn, Rh, 1./T, error_value, noise, &idum, pp, grid_mult*id.size()); // returns Likelihood
                }
                    
                ff = fopen(argv[5], "w");
                for (k = 0; k < id.size(); k++)
                {
                    node1 = id[k];
                    fprintf(ff, "%d\t%f\t%f\t%f\t%d\n", node1, a_list[node1], real_a_list[node1], r_list[node1], nw_old.nod[node1].degree);
                }
                fclose(ff);
                min = max;
                max = 20 * std::exp(i*std::log(mult));
            }
            L = h2loglikelihood(id, nw_old, a_list, r_list, nnn,  Rh, 1./T, pp);
            if (Lmax < L)
            {
                Lmax = L;
                best_iter = iter;
                best_coords = a_list;
            }
            std::cout << "size = " << id.size() << std::endl;
            std::cout << "iter = " << iter << std::endl;
            std::cout << "best_iter = " << best_iter << std::endl;
            std::cout << "L = " << L << std::endl;
            std::cout << "Lmax = " << Lmax << std::endl;
            
            
            
        }
        
        a_list = best_coords;
        L = h2loglikelihood(id, nw_old, a_list, r_list, nnn,  Rh, 1./T, pp);
        std::cout << "Loptimal = " << L << std::endl;
        
        std::cout << "resuming the rest of the calculations" << std::endl;
        for (i = mid_avg; i <= num_avg + 1; i++)
        {
            std::cout << "iter = " << i << std::endl;
            for (j = roundx(min); j < roundx(max); j++)
            {
                ctr++;
                node1  = knode_list[j];
                if (i == 1)
                {
                    a_list[node1] = 2*PI * (ran1(&idum) - 0.5);
                }
                else
                {
                    if ((id.size() < 500)||(nnn < 1500))
                    {
                        a_list[node1] = global_h2_find_node_coords(node1, id, nw_old, a_list, r_list, nnn, Rh, 1./T, pp, grid_mult*id.size());
                    }
                    else
                    {
                        a_list[node1] = h2_find_node_coords(node1, id, nw_old, a_list, r_list, nnn, Rh, 1./T, pp, grid_mult*id.size());
                    }
                }
                id.resize(node_ctr + 1);
                id[node_ctr] = node1;
                node_ctr++;
                nw_old.nod[node1].flag = 1;
            }
            // embed nodes single round
            noise = (i)*(0.01-max_noise)/num_avg + max_noise;
            shuffle_vector(id, &idum);
            if ((id.size() < 500)||(nnn < 1500))
            {
                L = global_h2_find_core_coords_single_round(id, nw_old, a_list, r_list, real_a_list, nnn, Rh, 1./T, error_value, noise, &idum, pp, grid_mult*id.size()); // returns Likelihood
            }
            else
            {
                L = h2_find_core_coords_single_round(id, nw_old, a_list, r_list, real_a_list, nnn, Rh, 1./T, error_value, noise, &idum, pp, grid_mult*id.size()); // returns Likelihood
            }
            std::cout << "noise = " << noise << ";    iter= " << i << "/" << num_avg << "; L = " << L << std::endl;
            ff = fopen(argv[5], "w");
            for (k = 0; k < id.size(); k++)
            {
                node1 = id[k];
                fprintf(ff, "%d\t%f\t%f\t%f\t%d\n", node1, a_list[node1], real_a_list[node1], r_list[node1], nw_old.nod[node1].degree);
            }
            fclose(ff);
            min = max;
            max = 20 * std::exp(i*std::log(mult));
            std::cout << "i = " << i << std::endl;
            std::cout << "mult = " << mult << std::endl;
            std::cout << "min = " << min << std::endl;
            std::cout << "max = " << max << std::endl;
        }
        std::cout << "Control sum = " << ctr << std::endl;    
        
        
        
        //************************************************************************************************
        std::cout << "E. Infer remaining node coordinates" << std::endl;
        h2_find_remaining_coords(id, nw_old, a_list, r_list, real_a_list, nnn, Rh, T, error_value, num_iter, max_noise, &idum, pp, argv[5]);
        
        std::cout << "E2. Massage All" << std::endl;
        h2_massage_all(id, nw_old, a_list, r_list, real_a_list, nnn, Rh, T, error_value, 20, &idum, pp, grid_mult*nnn, argv[5]);
        
        std::cout << "Add small noise to angular coordinates" << std::endl;
        for (i = 0; i < nw_old.num_nodes; i++)
        {
            if (nw_old.nod[i].degree > 0)
            {
                a_list[i] = 1e-4*(ran1(&idum) - 0.5) + a_list[i];
                if (a_list[i] > PI)
                {
                    a_list[i] -= 2*PI;
                }
                if (a_list[i] < -PI)
                {
                    a_list[i] += 2*PI;
                }
            }    
        }

        //************************************************************************************************
        std::cout << "output coordinates" << std::endl;
        ff = fopen(argv[5], "w");
        for (i = 0; i < nw_old.num_nodes; i++)
        {
            if (nw_old.nod[i].degree > 0)
            {
                fprintf(ff, "%d\t%d\t%f\t%f\n", i, nw_old.nod[i].degree, r_list[i], a_list[i]);
            }
        }
        fclose(ff);
        std::remove(temp_file_path); // @Jiaze remove temporary file
        return 0;
        //************************************************************************************************
    }
    
}

