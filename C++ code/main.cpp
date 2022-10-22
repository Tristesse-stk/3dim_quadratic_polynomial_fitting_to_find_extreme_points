#include <fstream>
#include <windows.h>  
#include <ctime>
#include <io.h>
#include "get_max_val_point.h"
using namespace std;

void read_data(const char* data_txt_path, double* x, double* y, double* z, double* t, int rows)
{
	int num = 0;
	std::ifstream data_txt(data_txt_path);
	std::string data_line;
	while (num <rows && data_txt >> data_line)
	{
		for (int col_id = 0; col_id < 4; col_id++)
		{
			string::size_type pos = data_line.find(",");
			if (pos != data_line.npos)
			{
				std::string data_str_ = data_line.substr(0, pos);
				data_line = data_line.substr(pos + 1, data_line.length() - pos - 1);
				if (col_id == 0)
				{
					x[num]= atof(data_str_.c_str());
				}
				else if (col_id == 1)
				{
					y[num] = atof(data_str_.c_str());
				}
				else if (col_id == 2)
				{
					z[num] = atof(data_str_.c_str());
				}
			}
			else {
				t[num] = atof(data_line.c_str());
			}
		}
		num += 1;
	}

}

int main()
{
	const char* data_txt_path = "D:\\Datasets\\origin\\拟合\\1.txt";
	if ((_access(data_txt_path, 0)) == -1)
	{
		cout << "txt文件路径不存在，请重新检查路径定义" << endl;
		return 0;
	}
	int num_points = 125;
	int num_parameters = 10;
	double threshold = 0.8;
	//x,y,z代表3个输入变量的数组，t代表得分结果标签
	double* x = new double[num_points]; double* y = new double[num_points]; double* z = new double[num_points]; double* t = new double[num_points];
	//x0_real_ptr,y0_real_ptr,z0_real_ptr代表极值点坐标的指针
	double* x0_real_ptr = new double; double* y0_real_ptr = new double; double* z0_real_ptr = new double;
	//把txt里的数据读上来，把数据给到x, y, z, t里去
	read_data(data_txt_path, x, y, z, t, num_points);
	clock_t startTime, endTime;
	startTime = clock();//计时开始
	//start = GetTickCount();
	//Sleep(1);
	//get_max_val_point函数参数解释:
	//x,y,z代表3个输入变量的数组，t代表得分结果标签
	//num_points代表有多少组数据,这里都是125
	//threshold代表得分阈值，只有得分超过此阈值的数据会作为有效数据用来拟合，小于等于它的都会被剔除，我这里设的是0.8
	//x0_real_ptr,y0_real_ptr,z0_real_ptr就是最终得到的极值点的坐标值，用指针传入，执行完后保存的就是最终结果
	
	//result代表有没有计算出结果，计算出来了就是true，中间有环节出错导致没法计算结果就是false
	bool result;
	result = get_max_val_point(x, y, z, t, num_points, threshold, num_parameters, x0_real_ptr, y0_real_ptr, z0_real_ptr);
	//for (int i = 0; i < 1000; i++)
	//{
	//	result = get_max_val_point(x, y, z, t, num_points, threshold, num_parameters, x0_real_ptr, y0_real_ptr, z0_real_ptr);
	//}
	
	endTime = clock();
	cout << "The run time is: " << (double)(endTime - startTime)/1000 / CLOCKS_PER_SEC << "s" << endl;
	if(result)
	{ 
		cout << "极大值点x:" << x0_real_ptr[0] << endl;
		cout << "极大值点y:" << y0_real_ptr[0] << endl;
		cout << "极大值点z:" << z0_real_ptr[0] << endl;
	}
	else
	{
		cout << "求极大值点过程有环节出错导致无法求得极大值点坐标" << endl;
	}
	delete[] x;delete[] y;delete[] z; delete[] t;
	delete x0_real_ptr; delete y0_real_ptr; delete z0_real_ptr;
	system("pause");
	return 0;
}