#ifndef GLOBAL_H
#define GLOBAL_H

#include <vector>

int roundx(double x);
double finite_alpha(double kp_max, double gamma);
void quickSort(std::vector<int>& numbers, std::vector<int>& addition);
float ran1(long *idum);
void shuffle_vector(std::vector<int>& v1, long *idum); //shuffles the entries of the vector 
double distH2(double r1, double angle1, double r2, double angle2, double zeta);
double distH2approx(double r1, double angle1, double r2, double angle2, double zeta);
int get_string_number(int num, char str[]);
double distS1(double x1, double x2);
int get_word_beginning(int num, char str[]);
int get_word_ending(int ctr1, char str[]);
int get_number(int ctr1, int ctr2, char str[]);
void q_sort(std::vector<int>& numbers, std::vector<int>& addition, int left, int right);
double incomplete_gamma(double a, double z);
double finite_avg_degree(int kmax, double kp_max, double avg_kp, double tpr, double gamma);
double incomplete_negative_gamma_function(double x, double avg_k, double avg_kp, double tpr, double gamma, double k_obs, double kmax);
void reverse(std::vector<int>& numbers);

#endif