#pragma once
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
//#include <cstdio>
//#include <cstdlib>
#include <iomanip>
using namespace std;
void get_norm_info(double* x, int num, double* x_min, double* x_max);

vector<int> select_data(double* t, int num, double threshold);

void normalize(double* x, int num, double* x_norm, double* x_min, double* x_max);

double get_real_val(double x0, double* x_max, double* x_min);

void select(double* x_norm, int num_select, vector<double>& x_norm_select, vector<int>& id_select);

double sum_2matrix_multiply(const vector<double>& a, const vector<double>& b, int num_select);

double sum_3matrix_multiply(const vector<double>& a, const vector<double>& b, const vector<double>& c, int num_select);

double sum_4matrix_multiply(const vector<double>& a, const vector<double>& b, const vector<double>& c, const vector<double>& d, int num_select);

void matrix_inverse(vector<vector<double>>& matrix_ori, vector<vector<double>>& matrix_inverse, int N);

void matrix_multiply(vector<vector<double>>& A_matrix_inverse, vector<vector<double>>& B_matrix, vector<vector<double>>& X_matrix, int rows, int cols);

void polynomial_fitting_7parameters_caculate_max_val_point(const vector<double>& x_norm_select, const vector<double>& y_norm_select, const vector<double>& z_norm_select, const vector<double>& t_select, int num_select, int num_parameter, double* x0_ptr, double* y0_ptr, double* z0_ptr);

void polynomial_fitting_10parameters_caculate_max_val_point(const vector<double>& x_norm_select, const vector<double>& y_norm_select, const vector<double>& z_norm_select, const vector<double>& t_select, int num_select, int num_parameter, double* x0_ptr, double* y0_ptr, double* z0_ptr);

bool get_max_val_point(double* x, double* y, double* z, double* t,int num, double threshold, int num_parameter,double* x0_real_ptr,double* y0_real_ptr,double* z0_real_ptr);