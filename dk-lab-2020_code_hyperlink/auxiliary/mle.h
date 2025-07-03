#ifndef MLE_H
#define MLE_H

#include <vector>
#include "network.h"

double global_h2_find_node_coords(int node1, std::vector<int>& id, network& nw, std::vector<double>& angles, std::vector<double>& real_radii, int nnn, double Rh, double beta, double tpr, int grid_size);
double h2_find_node_coords(int node1, std::vector<int>& id, network& nw, std::vector<double>& angles, std::vector<double>& real_radii, int nnn, double Rh, double beta, double tpr, int grid_size);
double  global_h2_find_core_coords_single_round(std::vector<int>& id, network& nw, std::vector<double>& angles, std::vector<double>& radii, std::vector<double>& real_angles, int nnn, double Rh, double beta, double error_value, double noise, long * idum, double tpr, int grid_size);
double h2_find_core_coords_single_round(std::vector<int>& id, network& nw, std::vector<double>& angles, std::vector<double>& radii, std::vector<double>& real_angles, int nnn, double Rh, double beta, double error_value, double noise, long * idum, double tpr, int grid_size);
double h2loglikelihood(std::vector <int>& id, network& nw, std::vector <double>& angles, std::vector <double>& real_radii, int nnn,  double Rh, double beta, double tpr);
void h2_find_remaining_coords(std::vector<int>& id, network& nw_old, std::vector<double>& a_list, std::vector<double>& r_list, std::vector<double>& real_a_list, int nnn, double Rh, double T, double error_value, int num_iter, double max_noise, long * idum, double tpr, char * out_filename);
void h2_massage_all(std::vector<int>& id, network& nw_old, std::vector<double>& a_list, std::vector<double>& r_list, std::vector<double>& real_a_list, int nnn, double Rh, double T, double error_value, int num_iter, long * idum, double tpr, int grid_size, char * out_filename);
double  h2_find_core_coords_single_round_progress(std::vector<int>& id, network& nw, std::vector<double>& angles, std::vector<double>& radii, std::vector<double>& real_angles, int nnn, double Rh, double beta, double error_value, double noise, long * idum, double tpr, int grid_size);

#endif