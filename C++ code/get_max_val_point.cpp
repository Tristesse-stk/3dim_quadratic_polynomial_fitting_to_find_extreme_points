#include "get_max_val_point.h"

void get_norm_info(double* x, int num, double* x_min, double* x_max)
{
	float min_temp = x[0];
	float max_temp = x[0];
	for (int i = 1; i < num; i++)
	{
		if (x[i] < min_temp)
		{
			min_temp = x[i];
		}

		if (x[i] > max_temp)
		{
			max_temp = x[i];
		}
	}
	x_min[0] = min_temp;
	x_max[0] = max_temp;
}

vector<int> select_data(double* t, int num, double threshold)
{
	vector<int> id_select;
	for (int i = 0; i < num; i++)
	{
		if (t[i] > threshold)
		{
			id_select.push_back(i);
		}
	}
	return id_select;
}

void normalize(double* x, int num, double* x_norm, double* x_min, double* x_max)
{
	float len_max_min = x_max[0] - x_min[0];
	float center = (x_max[0] + x_min[0]) / 2;
	for (int i = 0; i < num; ++i)
	{
		x_norm[i] = (x[i] - center) / (len_max_min / 2);
	}
}

double get_real_val(double x0, double* x_max, double* x_min)
{
	return x0 * ((x_max[0] - x_min[0]) / 2) + ((x_max[0] + x_min[0]) / 2);
}

void select(double* x_norm, int num_select, vector<double>& x_norm_select, vector<int>& id_select)
{
	for (int i = 0; i < num_select; i++)
	{
		x_norm_select[i] = x_norm[id_select[i]];
	}
}

double sum_2matrix_multiply(const vector<double>& a, const vector<double>& b, int num_select)
{
	double temp_sum = 0.0;
	for (int i = 0; i < num_select; i++)
	{
		temp_sum += a[i] * b[i];
	}
	return temp_sum;
}

double sum_3matrix_multiply(const vector<double>& a, const vector<double>& b, const vector<double>& c, int num_select)
{
	double temp_sum = 0.0;
	for (int i = 0; i < num_select; i++)
	{
		temp_sum += a[i] * b[i] * c[i];
	}
	return temp_sum;
}

double sum_4matrix_multiply(const vector<double>& a, const vector<double>& b, const vector<double>& c, const vector<double>& d, int num_select)
{
	double temp_sum = 0.0;
	for (int i = 0; i < num_select; i++)
	{
		temp_sum += a[i] * b[i] * c[i] * d[i];
	}
	return temp_sum;
}


//void matrix_inverse(double(*a)[N], double(*b)[N])
//求逆矩阵
//matrix_ori:原矩阵，matrix_inverse:逆矩阵
void matrix_inverse(vector<vector<double>>& matrix_ori, vector<vector<double>>& matrix_inverse, int N)
{
	int i, j, k;
	double max, temp;
	// 定义一个临时矩阵t
	//double t[N][N];
	vector<vector<double>> t(N, vector<double>(N));
	// 将a矩阵临时存放在矩阵t[n][n]中
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			t[i][j] = matrix_ori[i][j];
		}
	}
	// 初始化B矩阵为单位矩阵
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			matrix_inverse[i][j] = (i == j) ? (double)1 : 0;
		}
	}
	// 进行列主消元，找到每一列的主元
	for (i = 0; i < N; i++)
	{
		max = t[i][i];
		// 用于记录每一列中的第几个元素为主元
		k = i;
		// 寻找每一列中的主元元素
		for (j = i + 1; j < N; j++)
		{
			if (fabs(t[j][i]) > fabs(max))
			{
				max = t[j][i];
				k = j;
			}
		}
		//cout<<"the max number is "<<max<<endl;
		// 如果主元所在的行不是第i行，则进行行交换
		if (k != i)
		{
			// 进行行交换
			for (j = 0; j < N; j++)
			{
				temp = t[i][j];
				t[i][j] = t[k][j];
				t[k][j] = temp;
				// 伴随矩阵B也要进行行交换
				temp = matrix_inverse[i][j];
				matrix_inverse[i][j] = matrix_inverse[k][j];
				matrix_inverse[k][j] = temp;
			}
		}
		if (t[i][i] == 0)
		{
			//std::cout << "\nthe matrix does not exist inverse matrix\n";
			std::cout << "\n无法求解逆矩阵\n";
			break;
		}
		// 获取列主元素
		temp = t[i][i];
		// 将主元所在的行进行单位化处理
		//cout<<"\nthe temp is "<<temp<<endl;
		for (j = 0; j < N; j++)
		{
			t[i][j] = t[i][j] / temp;
			matrix_inverse[i][j] = matrix_inverse[i][j] / temp;
		}
		for (j = 0; j < N; j++)
		{
			if (j != i)
			{
				temp = t[j][i];
				//消去该列的其他元素
				for (k = 0; k < N; k++)
				{
					t[j][k] = t[j][k] - temp * t[i][k];
					matrix_inverse[j][k] = matrix_inverse[j][k] - temp * matrix_inverse[i][k];
				}
			}

		}

	}
}

void matrix_multiply(vector<vector<double>>& A_matrix_inverse, vector<vector<double>>& B_matrix, vector<vector<double>>& X_matrix, int rows, int cols)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			for (int ii = 0; ii < rows; ii++)
			{
				X_matrix[i][j] += A_matrix_inverse[i][ii] * B_matrix[ii][j];
			}
		}
	}
}

void polynomial_fitting_7parameters_caculate_max_val_point(const vector<double>& x_norm_select, const vector<double>& y_norm_select, const vector<double>& z_norm_select, const vector<double>& t_select,int num_select, int num_parameter, double* x0_ptr, double* y0_ptr, double* z0_ptr)
{
	vector<vector<double>> A_matrix(num_parameter, vector<double>(num_parameter));
	vector<vector<double>> B_matrix(num_parameter, vector<double>(1));
	vector<vector<double>> A_matrix_inverse(num_parameter, vector<double>(num_parameter));
	vector<vector<double>> X_matrix(num_parameter, vector<double>(1, 0));

	A_matrix[0][0] = sum_2matrix_multiply(x_norm_select, x_norm_select, num_select);
	A_matrix[0][1] = accumulate(x_norm_select.begin(), x_norm_select.end(), 0.0);
	A_matrix[0][2] = sum_2matrix_multiply(y_norm_select, y_norm_select, num_select);
	A_matrix[0][3] = accumulate(y_norm_select.begin(), y_norm_select.end(), 0.0);
	A_matrix[0][4] = sum_2matrix_multiply(z_norm_select, z_norm_select, num_select);
	A_matrix[0][5] = accumulate(z_norm_select.begin(), z_norm_select.end(), 0.0);
	A_matrix[0][6] = num_select;
	B_matrix[0][0] = accumulate(t_select.begin(), t_select.end(), 0.0);

	A_matrix[1][0] = sum_3matrix_multiply(x_norm_select, x_norm_select, z_norm_select, num_select);
	A_matrix[1][1] = sum_2matrix_multiply(x_norm_select, z_norm_select, num_select);
	A_matrix[1][2] = sum_3matrix_multiply(y_norm_select, y_norm_select, z_norm_select, num_select);
	A_matrix[1][3] = sum_2matrix_multiply(y_norm_select, z_norm_select, num_select);
	A_matrix[1][4] = sum_3matrix_multiply(z_norm_select, z_norm_select, z_norm_select, num_select);
	A_matrix[1][5] = sum_2matrix_multiply(z_norm_select, z_norm_select, num_select);
	A_matrix[1][6] = accumulate(z_norm_select.begin(), z_norm_select.end(), 0.0);
	B_matrix[1][0] = sum_2matrix_multiply(t_select, z_norm_select, num_select);

	A_matrix[2][0] = sum_4matrix_multiply(x_norm_select, x_norm_select, z_norm_select, z_norm_select, num_select);
	A_matrix[2][1] = sum_3matrix_multiply(x_norm_select, z_norm_select, z_norm_select, num_select);
	A_matrix[2][2] = sum_4matrix_multiply(y_norm_select, y_norm_select, z_norm_select, z_norm_select, num_select);
	A_matrix[2][3] = sum_3matrix_multiply(y_norm_select, z_norm_select, z_norm_select, num_select);
	A_matrix[2][4] = sum_4matrix_multiply(z_norm_select, z_norm_select, z_norm_select, z_norm_select, num_select);
	A_matrix[2][5] = sum_3matrix_multiply(z_norm_select, z_norm_select, z_norm_select, num_select);
	A_matrix[2][6] = sum_2matrix_multiply(z_norm_select, z_norm_select, num_select);
	B_matrix[2][0] = sum_3matrix_multiply(t_select, z_norm_select, z_norm_select, num_select);

	A_matrix[3][0] = sum_3matrix_multiply(x_norm_select, x_norm_select, y_norm_select, num_select);
	A_matrix[3][1] = sum_2matrix_multiply(x_norm_select, y_norm_select, num_select);
	A_matrix[3][2] = sum_3matrix_multiply(y_norm_select, y_norm_select, y_norm_select, num_select);
	A_matrix[3][3] = sum_2matrix_multiply(y_norm_select, y_norm_select, num_select);
	A_matrix[3][4] = sum_3matrix_multiply(z_norm_select, z_norm_select, y_norm_select, num_select);
	A_matrix[3][5] = sum_2matrix_multiply(z_norm_select, y_norm_select, num_select);
	A_matrix[3][6] = accumulate(y_norm_select.begin(), y_norm_select.end(), 0.0);
	B_matrix[3][0] = sum_2matrix_multiply(t_select, y_norm_select, num_select);

	A_matrix[4][0] = sum_4matrix_multiply(x_norm_select, x_norm_select, y_norm_select, y_norm_select, num_select);
	A_matrix[4][1] = sum_3matrix_multiply(x_norm_select, y_norm_select, y_norm_select, num_select);
	A_matrix[4][2] = sum_4matrix_multiply(y_norm_select, y_norm_select, y_norm_select, y_norm_select, num_select);
	A_matrix[4][3] = sum_3matrix_multiply(y_norm_select, y_norm_select, y_norm_select, num_select);
	A_matrix[4][4] = sum_4matrix_multiply(z_norm_select, z_norm_select, y_norm_select, y_norm_select, num_select);
	A_matrix[4][5] = sum_3matrix_multiply(z_norm_select, y_norm_select, y_norm_select, num_select);
	A_matrix[4][6] = sum_2matrix_multiply(y_norm_select, y_norm_select, num_select);
	B_matrix[4][0] = sum_3matrix_multiply(t_select, y_norm_select, y_norm_select, num_select);

	A_matrix[5][0] = sum_3matrix_multiply(x_norm_select, x_norm_select, x_norm_select, num_select);
	A_matrix[5][1] = sum_2matrix_multiply(x_norm_select, x_norm_select, num_select);
	A_matrix[5][2] = sum_3matrix_multiply(y_norm_select, y_norm_select, x_norm_select, num_select);
	A_matrix[5][3] = sum_2matrix_multiply(y_norm_select, x_norm_select, num_select);
	A_matrix[5][4] = sum_3matrix_multiply(z_norm_select, z_norm_select, x_norm_select, num_select);
	A_matrix[5][5] = sum_2matrix_multiply(z_norm_select, x_norm_select, num_select);
	A_matrix[5][6] = accumulate(x_norm_select.begin(), x_norm_select.end(), 0.0);
	B_matrix[5][0] = sum_2matrix_multiply(t_select, x_norm_select, num_select);

	A_matrix[6][0] = sum_4matrix_multiply(x_norm_select, x_norm_select, x_norm_select, x_norm_select, num_select);
	A_matrix[6][1] = sum_3matrix_multiply(x_norm_select, x_norm_select, x_norm_select, num_select);
	A_matrix[6][2] = sum_4matrix_multiply(y_norm_select, y_norm_select, x_norm_select, x_norm_select, num_select);
	A_matrix[6][3] = sum_3matrix_multiply(y_norm_select, x_norm_select, x_norm_select, num_select);
	A_matrix[6][4] = sum_4matrix_multiply(z_norm_select, z_norm_select, x_norm_select, x_norm_select, num_select);
	A_matrix[6][5] = sum_3matrix_multiply(z_norm_select, x_norm_select, x_norm_select, num_select);
	A_matrix[6][6] = sum_2matrix_multiply(x_norm_select, x_norm_select, num_select);
	B_matrix[6][0] = sum_3matrix_multiply(t_select, x_norm_select, x_norm_select, num_select);

	for (int i = 0; i < 7; i++)
	{
		for (int j = 0; j < 7; j++)
		{
			//cout<< A_matrix[i][j]<<" ";
			cout << setiosflags(ios::left) << setw(10) << A_matrix[i][j];
		}
		cout<< endl;
	}

	//for (int i = 0; i < 7; i++)
	//{
	//	cout << B_matrix[i][0] << endl;
	//}

	matrix_inverse(A_matrix, A_matrix_inverse, num_parameter);
	//for (int i = 0; i < 7; i++)
	//{
	//	for (int j = 0; j < 7; j++)
	//	{
	//		cout << A_matrix_inverse[i][j] << " ";
	//	}
	//	cout << endl;
	//}
	matrix_multiply(A_matrix_inverse, B_matrix, X_matrix, num_parameter, 1);

	//for (int i = 0; i < 7; i++)
	//{
	//	cout << X_matrix[i][0] << endl;
	//}
	//for (int i = 0; i < 10000; i++)
	//{
	//	cout << endl;
	//}

	double a0 = X_matrix[0][0]; double a1 = X_matrix[1][0]; double a2 = X_matrix[2][0];
	double a3 = X_matrix[3][0]; double a4 = X_matrix[4][0]; double a5 = X_matrix[5][0];
	x0_ptr[0] = -a1 / (2 * a0); y0_ptr[0] = -a3 / (2 * a2); z0_ptr[0] = -a5 / (2 * a4);
}

void polynomial_fitting_10parameters_caculate_max_val_point(const vector<double>& x_norm_select, const vector<double>& y_norm_select, const vector<double>& z_norm_select, const vector<double>& t_select, int num_select, int num_parameter,double* x0_ptr, double* y0_ptr, double* z0_ptr)
{
	//A_matrix，B_matrix:A_matrix大小10*10，B_matrix大小10*1
	//对10个待定参数求偏导令它为0得到的10个方程组，把系数和变量分离成矩阵相乘的形式，左边就是A_matrix，右边就是B_matrix
	vector<vector<double>> A_matrix(num_parameter, vector<double>(num_parameter));
	vector<vector<double>> B_matrix(num_parameter, vector<double>(1));
	//A_matrix_inverse:A_matrix的逆矩阵
	vector<vector<double>> A_matrix_inverse(num_parameter, vector<double>(num_parameter));
	//X_matrix:A_matrix*X_matrix=B_matrix，求解出X_matrix，也就是10个参数
	vector<vector<double>> X_matrix(num_parameter, vector<double>(1, 0));

	//C_matrix,D_matrix:C_matrix大小3*3，D_matrix大小3*1
	//对拟合好的多项式中的x,y,z求偏导令它为0得到的3个方程组，把系数和变量分离成矩阵相乘的形式，左边就是C_matrix，右边就是D_matrix
	vector<vector<double>> C_matrix(3, vector<double>(3));
	//C_matrix_inverse:C_matrix的逆矩阵
	vector<vector<double>> C_matrix_inverse(3, vector<double>(3));
	vector<vector<double>> D_matrix(3, vector<double>(1));

	//max_val_point:C_matrix*max_val_point=D_matrix，求解出max_val_point，也就是极大值点坐标的x0,y0,z0
	vector<vector<double>> max_val_point(3, vector<double>(1));

	//A_matrix[9][0] = sum_2matrix_multiply(x_norm_select, x_norm_select, num_select);
	//A_matrix[9][1] = sum_2matrix_multiply(y_norm_select, y_norm_select, num_select);
	//A_matrix[9][2] = sum_2matrix_multiply(z_norm_select, z_norm_select, num_select);
	//A_matrix[9][3] = sum_2matrix_multiply(x_norm_select, y_norm_select, num_select);
	//A_matrix[9][4] = sum_2matrix_multiply(x_norm_select, z_norm_select, num_select);
	//A_matrix[9][5] = sum_2matrix_multiply(y_norm_select, z_norm_select, num_select);
	//A_matrix[9][6] = accumulate(x_norm_select.begin(), x_norm_select.end(), 0.0);
	//A_matrix[9][7] = accumulate(y_norm_select.begin(), y_norm_select.end(), 0.0);
	//A_matrix[9][8] = accumulate(z_norm_select.begin(), z_norm_select.end(), 0.0);
	//A_matrix[9][9] = num_select;
	//B_matrix[9][0] = accumulate(t_select.begin(), t_select.end(), 0.0);

	//A_matrix[8][0] = sum_3matrix_multiply(x_norm_select, x_norm_select, z_norm_select, num_select);
	//A_matrix[8][1] = sum_3matrix_multiply(y_norm_select, y_norm_select, z_norm_select, num_select);
	//A_matrix[8][2] = sum_3matrix_multiply(z_norm_select, z_norm_select, z_norm_select, num_select);
	//A_matrix[8][3] = sum_3matrix_multiply(x_norm_select, y_norm_select, z_norm_select, num_select);
	//A_matrix[8][4] = sum_3matrix_multiply(x_norm_select, z_norm_select, z_norm_select, num_select);
	//A_matrix[8][5] = sum_3matrix_multiply(y_norm_select, z_norm_select, z_norm_select, num_select);
	//A_matrix[8][6] = sum_2matrix_multiply(x_norm_select, z_norm_select, num_select);//A_matrix[9][4]; //
	//A_matrix[8][7] = sum_2matrix_multiply(y_norm_select, z_norm_select, num_select);//A_matrix[9][5]; //
	//A_matrix[8][8] = sum_2matrix_multiply(z_norm_select, z_norm_select, num_select);//A_matrix[9][2]; //
	//A_matrix[8][9] = accumulate(z_norm_select.begin(), z_norm_select.end(), 0.0);//A_matrix[9][8];// 
	//B_matrix[8][0] = sum_2matrix_multiply(t_select, z_norm_select, num_select);

	//A_matrix[7][0] = sum_3matrix_multiply(x_norm_select, x_norm_select, y_norm_select, num_select);
	//A_matrix[7][1] = sum_3matrix_multiply(y_norm_select, y_norm_select, y_norm_select, num_select);
	//A_matrix[7][2] = sum_3matrix_multiply(z_norm_select, z_norm_select, y_norm_select, num_select);//A_matrix[8][5]; //
	//A_matrix[7][3] = sum_3matrix_multiply(x_norm_select, y_norm_select, y_norm_select, num_select);
	//A_matrix[7][4] = sum_3matrix_multiply(x_norm_select, z_norm_select, y_norm_select, num_select);//A_matrix[8][3]; //
	//A_matrix[7][5] = sum_3matrix_multiply(y_norm_select, z_norm_select, y_norm_select, num_select);//A_matrix[8][1]; //
	//A_matrix[7][6] = sum_2matrix_multiply(x_norm_select, y_norm_select, num_select);//A_matrix[9][3]; //
	//A_matrix[7][7] = sum_2matrix_multiply(y_norm_select, y_norm_select, num_select);//A_matrix[9][1]; //
	//A_matrix[7][8] = sum_2matrix_multiply(z_norm_select, y_norm_select, num_select);//A_matrix[9][5]; //
	//A_matrix[7][9] = accumulate(y_norm_select.begin(), y_norm_select.end(), 0.0);//A_matrix[9][7]; //
	//B_matrix[7][0] = sum_2matrix_multiply(t_select, y_norm_select, num_select);

	//A_matrix[6][0] = sum_3matrix_multiply(x_norm_select, x_norm_select, x_norm_select, num_select);
	//A_matrix[6][1] = sum_3matrix_multiply(y_norm_select, y_norm_select, x_norm_select, num_select);//A_matrix[7][3]; //
	//A_matrix[6][2] = sum_3matrix_multiply(z_norm_select, z_norm_select, x_norm_select, num_select);//A_matrix[8][4]; //
	//A_matrix[6][3] = sum_3matrix_multiply(x_norm_select, y_norm_select, x_norm_select, num_select);//A_matrix[7][0]; //
	//A_matrix[6][4] = sum_3matrix_multiply(x_norm_select, z_norm_select, x_norm_select, num_select);//A_matrix[8][0]; //
	//A_matrix[6][5] = sum_3matrix_multiply(y_norm_select, z_norm_select, x_norm_select, num_select);//A_matrix[8][3]; //
	//A_matrix[6][6] = sum_2matrix_multiply(x_norm_select, x_norm_select, num_select);//A_matrix[9][0]; //
	//A_matrix[6][7] = sum_2matrix_multiply(y_norm_select, x_norm_select, num_select);//A_matrix[9][3]; //
	//A_matrix[6][8] = sum_2matrix_multiply(z_norm_select, x_norm_select, num_select);//A_matrix[9][4]; //
	//A_matrix[6][9] = accumulate(x_norm_select.begin(), x_norm_select.end(), 0.0);//A_matrix[9][6]; //
	//B_matrix[6][0] = sum_2matrix_multiply(t_select, x_norm_select, num_select);

	//A_matrix[5][0] = sum_4matrix_multiply(x_norm_select, x_norm_select, y_norm_select, z_norm_select, num_select);
	//A_matrix[5][1] = sum_4matrix_multiply(y_norm_select, y_norm_select, y_norm_select, z_norm_select, num_select);
	//A_matrix[5][2] = sum_4matrix_multiply(z_norm_select, z_norm_select, y_norm_select, z_norm_select, num_select);
	//A_matrix[5][3] = sum_4matrix_multiply(x_norm_select, y_norm_select, y_norm_select, z_norm_select, num_select);
	//A_matrix[5][4] = sum_4matrix_multiply(x_norm_select, z_norm_select, y_norm_select, z_norm_select, num_select);
	//A_matrix[5][5] = sum_4matrix_multiply(y_norm_select, z_norm_select, y_norm_select, z_norm_select, num_select);
	//A_matrix[5][6] = sum_3matrix_multiply(x_norm_select, y_norm_select, z_norm_select, num_select); //A_matrix[8][3]; //
	//A_matrix[5][7] = sum_3matrix_multiply(y_norm_select, y_norm_select, z_norm_select, num_select);//A_matrix[8][1]; //
	//A_matrix[5][8] = sum_3matrix_multiply(z_norm_select, y_norm_select, z_norm_select, num_select);//A_matrix[8][5]; //
	//A_matrix[5][9] = sum_2matrix_multiply(y_norm_select, z_norm_select, num_select);//A_matrix[9][5]; //
	//B_matrix[5][0] = sum_3matrix_multiply(t_select, y_norm_select, z_norm_select, num_select);

	//A_matrix[4][0] = sum_4matrix_multiply(x_norm_select, x_norm_select, x_norm_select, z_norm_select, num_select);
	//A_matrix[4][1] = sum_4matrix_multiply(y_norm_select, y_norm_select, x_norm_select, z_norm_select, num_select);//A_matrix[5][3]; //
	//A_matrix[4][2] = sum_4matrix_multiply(z_norm_select, z_norm_select, x_norm_select, z_norm_select, num_select);
	//A_matrix[4][3] = sum_4matrix_multiply(x_norm_select, y_norm_select, x_norm_select, z_norm_select, num_select);//A_matrix[5][0]; //
	//A_matrix[4][4] = sum_4matrix_multiply(x_norm_select, z_norm_select, x_norm_select, z_norm_select, num_select);
	//A_matrix[4][5] = sum_4matrix_multiply(y_norm_select, z_norm_select, x_norm_select, z_norm_select, num_select);//A_matrix[5][4]; //
	//A_matrix[4][6] = sum_3matrix_multiply(x_norm_select, x_norm_select, z_norm_select, num_select);//A_matrix[8][0]; //
	//A_matrix[4][7] = sum_3matrix_multiply(y_norm_select, x_norm_select, z_norm_select, num_select);//A_matrix[8][3]; //
	//A_matrix[4][8] = sum_3matrix_multiply(z_norm_select, x_norm_select, z_norm_select, num_select);//A_matrix[8][4]; //
	//A_matrix[4][9] = sum_2matrix_multiply(x_norm_select, z_norm_select, num_select);//A_matrix[9][4]; //
	//B_matrix[4][0] = sum_3matrix_multiply(t_select, x_norm_select, z_norm_select, num_select);

	//A_matrix[3][0] = sum_4matrix_multiply(x_norm_select, x_norm_select, x_norm_select, y_norm_select, num_select);
	//A_matrix[3][1] = sum_4matrix_multiply(y_norm_select, y_norm_select, x_norm_select, y_norm_select, num_select);
	//A_matrix[3][2] = sum_4matrix_multiply(z_norm_select, z_norm_select, x_norm_select, y_norm_select, num_select);//A_matrix[5][4]; //
	//A_matrix[3][3] = sum_4matrix_multiply(x_norm_select, y_norm_select, x_norm_select, y_norm_select, num_select);
	//A_matrix[3][4] = sum_4matrix_multiply(x_norm_select, z_norm_select, x_norm_select, y_norm_select, num_select);//A_matrix[5][0]; //
	//A_matrix[3][5] = sum_4matrix_multiply(y_norm_select, z_norm_select, x_norm_select, y_norm_select, num_select);//A_matrix[5][3]; //
	//A_matrix[3][6] = sum_3matrix_multiply(x_norm_select, x_norm_select, y_norm_select, num_select);//A_matrix[7][0]; //
	//A_matrix[3][7] = sum_3matrix_multiply(y_norm_select, x_norm_select, y_norm_select, num_select);//A_matrix[7][3]; //
	//A_matrix[3][8] = sum_3matrix_multiply(z_norm_select, x_norm_select, y_norm_select, num_select);//A_matrix[8][3]; //
	//A_matrix[3][9] = sum_2matrix_multiply(x_norm_select, y_norm_select, num_select);//A_matrix[7][6]; //
	//B_matrix[3][0] = sum_3matrix_multiply(t_select, x_norm_select, y_norm_select, num_select);

	//A_matrix[2][0] = sum_4matrix_multiply(x_norm_select, x_norm_select, z_norm_select, z_norm_select, num_select);//A_matrix[4][4]; //
	//A_matrix[2][1] = sum_4matrix_multiply(y_norm_select, y_norm_select, z_norm_select, z_norm_select, num_select);//A_matrix[5][5]; //
	//A_matrix[2][2] = sum_4matrix_multiply(z_norm_select, z_norm_select, z_norm_select, z_norm_select, num_select);
	//A_matrix[2][3] = sum_4matrix_multiply(x_norm_select, y_norm_select, z_norm_select, z_norm_select, num_select);//A_matrix[5][4]; //
	//A_matrix[2][4] = sum_4matrix_multiply(x_norm_select, z_norm_select, z_norm_select, z_norm_select, num_select);//A_matrix[4][2]; //
	//A_matrix[2][5] = sum_4matrix_multiply(y_norm_select, z_norm_select, z_norm_select, z_norm_select, num_select);//A_matrix[5][2]; //
	//A_matrix[2][6] = sum_3matrix_multiply(x_norm_select, z_norm_select, z_norm_select, num_select);//A_matrix[8][4]; //
	//A_matrix[2][7] = sum_3matrix_multiply(y_norm_select, z_norm_select, z_norm_select, num_select);//A_matrix[8][5]; //
	//A_matrix[2][8] = sum_3matrix_multiply(z_norm_select, z_norm_select, z_norm_select, num_select);//A_matrix[8][2]; //
	//A_matrix[2][9] = sum_2matrix_multiply(z_norm_select, z_norm_select, num_select);//A_matrix[9][2]; //
	//B_matrix[2][0] = sum_3matrix_multiply(t_select, z_norm_select, z_norm_select, num_select);

	//A_matrix[1][0] = sum_4matrix_multiply(x_norm_select, x_norm_select, y_norm_select, y_norm_select, num_select);//A_matrix[3][3]; //
	//A_matrix[1][1] = sum_4matrix_multiply(y_norm_select, y_norm_select, y_norm_select, y_norm_select, num_select);
	//A_matrix[1][2] = sum_4matrix_multiply(z_norm_select, z_norm_select, y_norm_select, y_norm_select, num_select);//A_matrix[5][5]; //
	//A_matrix[1][3] = sum_4matrix_multiply(x_norm_select, y_norm_select, y_norm_select, y_norm_select, num_select);//A_matrix[3][1]; //
	//A_matrix[1][4] = sum_4matrix_multiply(x_norm_select, z_norm_select, y_norm_select, y_norm_select, num_select);//A_matrix[5][3]; //
	//A_matrix[1][5] = sum_4matrix_multiply(y_norm_select, z_norm_select, y_norm_select, y_norm_select, num_select);//A_matrix[5][1]; //
	//A_matrix[1][6] = sum_3matrix_multiply(x_norm_select, y_norm_select, y_norm_select, num_select);//A_matrix[7][3]; //
	//A_matrix[1][7] = sum_3matrix_multiply(y_norm_select, y_norm_select, y_norm_select, num_select);//A_matrix[7][1]; //
	//A_matrix[1][8] = sum_3matrix_multiply(z_norm_select, y_norm_select, y_norm_select, num_select);//A_matrix[8][1]; //
	//A_matrix[1][9] = sum_2matrix_multiply(y_norm_select, y_norm_select, num_select);//A_matrix[9][1]; //
	//B_matrix[1][0] = sum_3matrix_multiply(t_select, y_norm_select, y_norm_select, num_select);

	//A_matrix[0][0] = sum_4matrix_multiply(x_norm_select, x_norm_select, x_norm_select, x_norm_select, num_select);
	//A_matrix[0][1] = sum_4matrix_multiply(y_norm_select, y_norm_select, x_norm_select, x_norm_select, num_select);//A_matrix[3][3]; //
	//A_matrix[0][2] = sum_4matrix_multiply(z_norm_select, z_norm_select, x_norm_select, x_norm_select, num_select);//A_matrix[4][4]; //
	//A_matrix[0][3] = sum_4matrix_multiply(x_norm_select, y_norm_select, x_norm_select, x_norm_select, num_select);//A_matrix[3][0]; //
	//A_matrix[0][4] = sum_4matrix_multiply(x_norm_select, z_norm_select, x_norm_select, x_norm_select, num_select);//A_matrix[4][0]; //
	//A_matrix[0][5] = sum_4matrix_multiply(y_norm_select, z_norm_select, x_norm_select, x_norm_select, num_select);//A_matrix[5][0]; //
	//A_matrix[0][6] = sum_3matrix_multiply(x_norm_select, x_norm_select, x_norm_select, num_select);//A_matrix[6][0]; //
	//A_matrix[0][7] = sum_3matrix_multiply(y_norm_select, x_norm_select, x_norm_select, num_select);//A_matrix[7][0]; //
	//A_matrix[0][8] = sum_3matrix_multiply(z_norm_select, x_norm_select, x_norm_select, num_select);//A_matrix[8][0]; //
	//A_matrix[0][9] = sum_2matrix_multiply(x_norm_select, x_norm_select, num_select);//A_matrix[9][0]; //
	//B_matrix[0][0] = sum_3matrix_multiply(t_select, x_norm_select, x_norm_select, num_select);

	A_matrix[9][0] = sum_2matrix_multiply(x_norm_select, x_norm_select, num_select);
	A_matrix[9][1] = sum_2matrix_multiply(y_norm_select, y_norm_select, num_select);
	A_matrix[9][2] = sum_2matrix_multiply(z_norm_select, z_norm_select, num_select);
	A_matrix[9][3] = sum_2matrix_multiply(x_norm_select, y_norm_select, num_select);
	A_matrix[9][4] = sum_2matrix_multiply(x_norm_select, z_norm_select, num_select);
	A_matrix[9][5] = sum_2matrix_multiply(y_norm_select, z_norm_select, num_select);
	A_matrix[9][6] = accumulate(x_norm_select.begin(), x_norm_select.end(), 0.0);
	A_matrix[9][7] = accumulate(y_norm_select.begin(), y_norm_select.end(), 0.0);
	A_matrix[9][8] = accumulate(z_norm_select.begin(), z_norm_select.end(), 0.0);
	A_matrix[9][9] = num_select;
	B_matrix[9][0] = accumulate(t_select.begin(), t_select.end(), 0.0);

	A_matrix[8][0] = sum_3matrix_multiply(x_norm_select, x_norm_select, z_norm_select, num_select);
	A_matrix[8][1] = sum_3matrix_multiply(y_norm_select, y_norm_select, z_norm_select, num_select);
	A_matrix[8][2] = sum_3matrix_multiply(z_norm_select, z_norm_select, z_norm_select, num_select);
	A_matrix[8][3] = sum_3matrix_multiply(x_norm_select, y_norm_select, z_norm_select, num_select);
	A_matrix[8][4] = sum_3matrix_multiply(x_norm_select, z_norm_select, z_norm_select, num_select);
	A_matrix[8][5] = sum_3matrix_multiply(y_norm_select, z_norm_select, z_norm_select, num_select);
	A_matrix[8][6] = A_matrix[9][4]; //sum_2matrix_multiply(x_norm_select, z_norm_select, num_select);//
	A_matrix[8][7] = A_matrix[9][5]; //sum_2matrix_multiply(y_norm_select, z_norm_select, num_select);//
	A_matrix[8][8] = A_matrix[9][2]; //sum_2matrix_multiply(z_norm_select, z_norm_select, num_select);//
	A_matrix[8][9] = A_matrix[9][8];// accumulate(z_norm_select.begin(), z_norm_select.end(), 0.0);//
	B_matrix[8][0] = sum_2matrix_multiply(t_select, z_norm_select, num_select);

	A_matrix[7][0] = sum_3matrix_multiply(x_norm_select, x_norm_select, y_norm_select, num_select);
	A_matrix[7][1] = sum_3matrix_multiply(y_norm_select, y_norm_select, y_norm_select, num_select);
	A_matrix[7][2] = A_matrix[8][5]; //sum_3matrix_multiply(z_norm_select, z_norm_select, y_norm_select, num_select);//
	A_matrix[7][3] = sum_3matrix_multiply(x_norm_select, y_norm_select, y_norm_select, num_select);
	A_matrix[7][4] = A_matrix[8][3]; //sum_3matrix_multiply(x_norm_select, z_norm_select, y_norm_select, num_select);//
	A_matrix[7][5] = A_matrix[8][1]; //sum_3matrix_multiply(y_norm_select, z_norm_select, y_norm_select, num_select);//
	A_matrix[7][6] = A_matrix[9][3]; //sum_2matrix_multiply(x_norm_select, y_norm_select, num_select);//
	A_matrix[7][7] = A_matrix[9][1]; //sum_2matrix_multiply(y_norm_select, y_norm_select, num_select);//
	A_matrix[7][8] = A_matrix[9][5]; //sum_2matrix_multiply(z_norm_select, y_norm_select, num_select);//
	A_matrix[7][9] = A_matrix[9][7]; //accumulate(y_norm_select.begin(), y_norm_select.end(), 0.0);//
	B_matrix[7][0] = sum_2matrix_multiply(t_select, y_norm_select, num_select);

	A_matrix[6][0] = sum_3matrix_multiply(x_norm_select, x_norm_select, x_norm_select, num_select);
	A_matrix[6][1] = A_matrix[7][3]; //sum_3matrix_multiply(y_norm_select, y_norm_select, x_norm_select, num_select);//
	A_matrix[6][2] = A_matrix[8][4]; //sum_3matrix_multiply(z_norm_select, z_norm_select, x_norm_select, num_select);//
	A_matrix[6][3] = A_matrix[7][0]; //sum_3matrix_multiply(x_norm_select, y_norm_select, x_norm_select, num_select);//
	A_matrix[6][4] = A_matrix[8][0]; //sum_3matrix_multiply(x_norm_select, z_norm_select, x_norm_select, num_select);//
	A_matrix[6][5] = A_matrix[8][3]; //sum_3matrix_multiply(y_norm_select, z_norm_select, x_norm_select, num_select);//
	A_matrix[6][6] = A_matrix[9][0]; //sum_2matrix_multiply(x_norm_select, x_norm_select, num_select);//
	A_matrix[6][7] = A_matrix[9][3]; //sum_2matrix_multiply(y_norm_select, x_norm_select, num_select);//
	A_matrix[6][8] = A_matrix[9][4]; //sum_2matrix_multiply(z_norm_select, x_norm_select, num_select);//
	A_matrix[6][9] = A_matrix[9][6]; //accumulate(x_norm_select.begin(), x_norm_select.end(), 0.0);//
	B_matrix[6][0] = sum_2matrix_multiply(t_select, x_norm_select, num_select);

	A_matrix[5][0] = sum_4matrix_multiply(x_norm_select, x_norm_select, y_norm_select, z_norm_select, num_select);
	A_matrix[5][1] = sum_4matrix_multiply(y_norm_select, y_norm_select, y_norm_select, z_norm_select, num_select);
	A_matrix[5][2] = sum_4matrix_multiply(z_norm_select, z_norm_select, y_norm_select, z_norm_select, num_select);
	A_matrix[5][3] = sum_4matrix_multiply(x_norm_select, y_norm_select, y_norm_select, z_norm_select, num_select);
	A_matrix[5][4] = sum_4matrix_multiply(x_norm_select, z_norm_select, y_norm_select, z_norm_select, num_select);
	A_matrix[5][5] = sum_4matrix_multiply(y_norm_select, z_norm_select, y_norm_select, z_norm_select, num_select);
	A_matrix[5][6] = A_matrix[8][3]; //sum_3matrix_multiply(x_norm_select, y_norm_select, z_norm_select, num_select); //
	A_matrix[5][7] = A_matrix[8][1]; //sum_3matrix_multiply(y_norm_select, y_norm_select, z_norm_select, num_select);//
	A_matrix[5][8] = A_matrix[8][5]; //sum_3matrix_multiply(z_norm_select, y_norm_select, z_norm_select, num_select);//
	A_matrix[5][9] = A_matrix[9][5]; //sum_2matrix_multiply(y_norm_select, z_norm_select, num_select);//
	B_matrix[5][0] = sum_3matrix_multiply(t_select, y_norm_select, z_norm_select, num_select);

	A_matrix[4][0] = sum_4matrix_multiply(x_norm_select, x_norm_select, x_norm_select, z_norm_select, num_select);
	A_matrix[4][1] = A_matrix[5][3]; //sum_4matrix_multiply(y_norm_select, y_norm_select, x_norm_select, z_norm_select, num_select);//
	A_matrix[4][2] = sum_4matrix_multiply(z_norm_select, z_norm_select, x_norm_select, z_norm_select, num_select);
	A_matrix[4][3] = A_matrix[5][0]; //sum_4matrix_multiply(x_norm_select, y_norm_select, x_norm_select, z_norm_select, num_select);//
	A_matrix[4][4] = sum_4matrix_multiply(x_norm_select, z_norm_select, x_norm_select, z_norm_select, num_select);
	A_matrix[4][5] = A_matrix[5][4]; //sum_4matrix_multiply(y_norm_select, z_norm_select, x_norm_select, z_norm_select, num_select);//
	A_matrix[4][6] = A_matrix[8][0]; //sum_3matrix_multiply(x_norm_select, x_norm_select, z_norm_select, num_select);//
	A_matrix[4][7] = A_matrix[8][3]; //sum_3matrix_multiply(y_norm_select, x_norm_select, z_norm_select, num_select);//
	A_matrix[4][8] = A_matrix[8][4]; //sum_3matrix_multiply(z_norm_select, x_norm_select, z_norm_select, num_select);//
	A_matrix[4][9] = A_matrix[9][4]; //sum_2matrix_multiply(x_norm_select, z_norm_select, num_select);//
	B_matrix[4][0] = sum_3matrix_multiply(t_select, x_norm_select, z_norm_select, num_select);

	A_matrix[3][0] = sum_4matrix_multiply(x_norm_select, x_norm_select, x_norm_select, y_norm_select, num_select);
	A_matrix[3][1] = sum_4matrix_multiply(y_norm_select, y_norm_select, x_norm_select, y_norm_select, num_select);
	A_matrix[3][2] = A_matrix[5][4]; //sum_4matrix_multiply(z_norm_select, z_norm_select, x_norm_select, y_norm_select, num_select);//
	A_matrix[3][3] = sum_4matrix_multiply(x_norm_select, y_norm_select, x_norm_select, y_norm_select, num_select);
	A_matrix[3][4] = A_matrix[5][0]; //sum_4matrix_multiply(x_norm_select, z_norm_select, x_norm_select, y_norm_select, num_select);//
	A_matrix[3][5] = A_matrix[5][3]; //sum_4matrix_multiply(y_norm_select, z_norm_select, x_norm_select, y_norm_select, num_select);//
	A_matrix[3][6] = A_matrix[7][0]; //sum_3matrix_multiply(x_norm_select, x_norm_select, y_norm_select, num_select);//
	A_matrix[3][7] = A_matrix[7][3]; //sum_3matrix_multiply(y_norm_select, x_norm_select, y_norm_select, num_select);//
	A_matrix[3][8] = A_matrix[8][3]; //sum_3matrix_multiply(z_norm_select, x_norm_select, y_norm_select, num_select);//
	A_matrix[3][9] = A_matrix[7][6]; //sum_2matrix_multiply(x_norm_select, y_norm_select, num_select);//
	B_matrix[3][0] = sum_3matrix_multiply(t_select, x_norm_select, y_norm_select, num_select);

	A_matrix[2][0] = A_matrix[4][4]; //sum_4matrix_multiply(x_norm_select, x_norm_select, z_norm_select, z_norm_select, num_select);//
	A_matrix[2][1] = A_matrix[5][5]; //sum_4matrix_multiply(y_norm_select, y_norm_select, z_norm_select, z_norm_select, num_select);//
	A_matrix[2][2] = sum_4matrix_multiply(z_norm_select, z_norm_select, z_norm_select, z_norm_select, num_select);
	A_matrix[2][3] = A_matrix[5][4]; sum_4matrix_multiply(x_norm_select, y_norm_select, z_norm_select, z_norm_select, num_select);//
	A_matrix[2][4] = A_matrix[4][2]; //sum_4matrix_multiply(x_norm_select, z_norm_select, z_norm_select, z_norm_select, num_select);//
	A_matrix[2][5] = A_matrix[5][2]; //sum_4matrix_multiply(y_norm_select, z_norm_select, z_norm_select, z_norm_select, num_select);//
	A_matrix[2][6] = A_matrix[8][4]; //sum_3matrix_multiply(x_norm_select, z_norm_select, z_norm_select, num_select);//
	A_matrix[2][7] = A_matrix[8][5]; //sum_3matrix_multiply(y_norm_select, z_norm_select, z_norm_select, num_select);//
	A_matrix[2][8] = A_matrix[8][2]; //sum_3matrix_multiply(z_norm_select, z_norm_select, z_norm_select, num_select);//
	A_matrix[2][9] = A_matrix[9][2]; //sum_2matrix_multiply(z_norm_select, z_norm_select, num_select);//
	B_matrix[2][0] = sum_3matrix_multiply(t_select, z_norm_select, z_norm_select, num_select);

	A_matrix[1][0] = A_matrix[3][3]; //sum_4matrix_multiply(x_norm_select, x_norm_select, y_norm_select, y_norm_select, num_select);//
	A_matrix[1][1] = sum_4matrix_multiply(y_norm_select, y_norm_select, y_norm_select, y_norm_select, num_select);
	A_matrix[1][2] = A_matrix[5][5]; //sum_4matrix_multiply(z_norm_select, z_norm_select, y_norm_select, y_norm_select, num_select);//
	A_matrix[1][3] = A_matrix[3][1]; //sum_4matrix_multiply(x_norm_select, y_norm_select, y_norm_select, y_norm_select, num_select);//
	A_matrix[1][4] = A_matrix[5][3]; //sum_4matrix_multiply(x_norm_select, z_norm_select, y_norm_select, y_norm_select, num_select);//
	A_matrix[1][5] = A_matrix[5][1]; //sum_4matrix_multiply(y_norm_select, z_norm_select, y_norm_select, y_norm_select, num_select);//
	A_matrix[1][6] = A_matrix[7][3]; //sum_3matrix_multiply(x_norm_select, y_norm_select, y_norm_select, num_select);//
	A_matrix[1][7] = A_matrix[7][1]; //sum_3matrix_multiply(y_norm_select, y_norm_select, y_norm_select, num_select);//
	A_matrix[1][8] = A_matrix[8][1]; //sum_3matrix_multiply(z_norm_select, y_norm_select, y_norm_select, num_select);//
	A_matrix[1][9] = A_matrix[9][1]; //sum_2matrix_multiply(y_norm_select, y_norm_select, num_select);//
	B_matrix[1][0] = sum_3matrix_multiply(t_select, y_norm_select, y_norm_select, num_select);

	A_matrix[0][0] = sum_4matrix_multiply(x_norm_select, x_norm_select, x_norm_select, x_norm_select, num_select);
	A_matrix[0][1] = A_matrix[3][3]; //sum_4matrix_multiply(y_norm_select, y_norm_select, x_norm_select, x_norm_select, num_select);//
	A_matrix[0][2] = A_matrix[4][4]; //sum_4matrix_multiply(z_norm_select, z_norm_select, x_norm_select, x_norm_select, num_select);//
	A_matrix[0][3] = A_matrix[3][0]; //sum_4matrix_multiply(x_norm_select, y_norm_select, x_norm_select, x_norm_select, num_select);//
	A_matrix[0][4] = A_matrix[4][0]; //sum_4matrix_multiply(x_norm_select, z_norm_select, x_norm_select, x_norm_select, num_select);//
	A_matrix[0][5] = A_matrix[5][0]; //sum_4matrix_multiply(y_norm_select, z_norm_select, x_norm_select, x_norm_select, num_select);//
	A_matrix[0][6] = A_matrix[6][0]; //sum_3matrix_multiply(x_norm_select, x_norm_select, x_norm_select, num_select);//
	A_matrix[0][7] = A_matrix[7][0]; //sum_3matrix_multiply(y_norm_select, x_norm_select, x_norm_select, num_select);//
	A_matrix[0][8] = A_matrix[8][0]; //sum_3matrix_multiply(z_norm_select, x_norm_select, x_norm_select, num_select);//
	A_matrix[0][9] = A_matrix[9][0]; //sum_2matrix_multiply(x_norm_select, x_norm_select, num_select);//
	B_matrix[0][0] = sum_3matrix_multiply(t_select, x_norm_select, x_norm_select, num_select);




	for (int i = 0; i < num_parameter; i++)
	{
		for (int j = 0; j < num_parameter; j++)
		{
			cout << setiosflags(ios::left) << setw(10) << A_matrix[i][j];
		}
		cout << endl;
	}


	//for (int i = 0; i < num_parameter; i++)
	//{
	//	cout << B_matrix[i][0] << endl;
	//}

	//对矩阵A求逆
	matrix_inverse(A_matrix, A_matrix_inverse, num_parameter);
	//for (int i = 0; i < num_parameter; i++)
	//{
	//	for (int j = 0; j < num_parameter; j++)
	//	{
	//		cout << setiosflags(ios::left) << setw(10) << A_matrix_inverse[i][j];
	//	}
	//	cout << endl;
	//}

	//求解出X_matrix
	matrix_multiply(A_matrix_inverse, B_matrix, X_matrix, num_parameter, 1);

	//for (int i = 0; i < 7; i++)
	//{
	//	cout << X_matrix[i][0] << endl;
	//}
	//for (int i = 0; i < 10000; i++)
	//{
	//	cout << endl;
	//}

	
	double a0 = X_matrix[0][0]; double a1 = X_matrix[1][0]; double a2 = X_matrix[2][0];
	double a3 = X_matrix[3][0]; double a4 = X_matrix[4][0]; double a5 = X_matrix[5][0];
	double a6 = X_matrix[6][0]; double a7 = X_matrix[7][0]; double a8 = X_matrix[8][0];

	C_matrix[0][0] = 2 * a0; C_matrix[0][1] = a3; C_matrix[0][2] = a4;
	C_matrix[1][0] = a3; C_matrix[1][1] = 2 * a1; C_matrix[1][2] = a5;
	C_matrix[2][0] = a4; C_matrix[2][1] = a5; C_matrix[2][2] = 2 * a2;

	D_matrix[0][0] = -a6; D_matrix[1][0] = -a7; D_matrix[2][0] = -a8;

	//对C求逆矩阵
	matrix_inverse(C_matrix, C_matrix_inverse, 3);
	//求解得max_val_point
	matrix_multiply(C_matrix_inverse, D_matrix, max_val_point, 3, 1);
	x0_ptr[0] = max_val_point[0][0];
	y0_ptr[0] = max_val_point[1][0];
	z0_ptr[0] = max_val_point[2][0];
}

bool get_max_val_point(double* x, double* y, double* z, double* t, int num,double threshold,int num_parameter, double* x0_real_ptr, double* y0_real_ptr, double* z0_real_ptr)
{
	//x_norm，y_norm，z_norm是正则化之后存放的数据
	double* x_norm = new double[num]; double* y_norm = new double[num]; double* z_norm = new double[num];
	//x_min，y_min，z_min是x,y,z数组中的最小值；x_max，y_max，z_max是x,y,z数组中的最大值
	double* x_min = new double; double* x_max = new double;
	double* y_min = new double; double* y_max = new double;
	double* z_min = new double; double* z_max = new double;

	//get_norm_info():获取x,y,z的最大值和最小值，用于正则化,x_min和x_max记录
	get_norm_info(x, num, x_min, x_max); get_norm_info(y, num, y_min, y_max); get_norm_info(z, num, z_min, z_max);
	//select_data():对t进行筛选，只有得分在threshold以上的才能作为有效的数据，返回的是满足条件的序号向量
	vector<int> id_select = select_data(t, num, threshold);
	//num_select:筛选过之后的数据个数
	int num_select = id_select.size();
	cout << "用于拟合的有效数据数量为:" << num_select << endl;
	if (num_select == 0)
	{
		cout << "阈值threshold设置过高导致没有数据满足条件，请调小threshold" << endl;
		return false;
	}
	//x_norm_select,y_norm_select,z_norm_select是经过预处理之后用来拟合的数据
	vector<double> x_norm_select(num_select); vector<double> y_norm_select(num_select); vector<double> z_norm_select(num_select); vector<double> t_select(num_select);
	
	//对原始数据x,y,z正则化操作压缩到-1到1之间，得到x_norm，y_norm，z_norm
	normalize(x, num, x_norm, x_min, x_max); normalize(y, num, y_norm, y_min, y_max); normalize(z, num, z_norm, z_min, z_max);
	
	//利用上述的id_select对正则化之后的数据x_norm，y_norm，z_norm,t进行筛选，得到x_norm_select，y_norm_select，z_norm_select,t_select(t已经在0到1之间无需正则化)
	select(x_norm, num_select, x_norm_select, id_select); select(y_norm, num_select, y_norm_select, id_select); select(z_norm, num_select, z_norm_select, id_select); select(t, num_select, t_select, id_select);
	
	//x0_ptr,y0_ptr,z0_ptr是正则化之后的数据的极大值点坐标的指针
	double* x0_ptr = new double; double* y0_ptr = new double; double* z0_ptr = new double;

	//目待拟合多项式参数有7个或者10个
	if (num_parameter == 7)
	{
		//待拟合多项式:t=a0*x^2+a1*y^2+a2*z^2+a3*x+a4*y+a5*z+a6
		//输入x_norm_select, y_norm_select, z_norm_select, t_select, num_select, num_parameter
		//输出x0_ptr, y0_ptr, z0_ptr
		polynomial_fitting_7parameters_caculate_max_val_point(x_norm_select, y_norm_select, z_norm_select, t_select, num_select, num_parameter, x0_ptr, y0_ptr, z0_ptr);
	}
	else if (num_parameter == 10)
	{
		//待拟合多项式:t=a0*x^2+a1*y^2+a2*z^2+a3*x*y+a4*x*z+a5*y*z+a6*x+a7*y+a8*z+a9
		//输入和输出同上
		polynomial_fitting_10parameters_caculate_max_val_point(x_norm_select, y_norm_select, z_norm_select, t_select, num_select, num_parameter, x0_ptr, y0_ptr, z0_ptr);
	}

	//get_real_val()将正则化数据求得的极大值点返回到真实坐标下的极大值
	double x0_real = get_real_val(x0_ptr[0], x_max, x_min); double y0_real = get_real_val(y0_ptr[0], y_max, y_min); double z0_real = get_real_val(z0_ptr[0], z_max, z_min);
	x0_real_ptr[0] = x0_real; y0_real_ptr[0] = y0_real; z0_real_ptr[0] = z0_real;
	delete x0_ptr; delete y0_ptr; delete z0_ptr;
	delete x_min; delete x_max;delete y_min; delete y_max;delete z_min; delete z_max;
	delete[] x_norm; delete[] y_norm; delete[] z_norm;
	return true;
}

