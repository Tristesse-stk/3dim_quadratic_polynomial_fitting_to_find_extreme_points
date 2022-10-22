此项目利用已知的多组3维数据来做3元2次多项式拟合，进一步求解极值点。

包括了Matlab和C++两份语言的代码。

MATLAB：

使用MATLAB软件打开3dim_quadratic_polynomial_fitting_to_find_extreme_points.m,修改文件路径，然后执行即可，

在命令行窗口x0_real、y0_real、z0_real就是最终求解结果。

C++：

在Visual Studio创建新项目，将main.cpp和get_max_val_point.cpp添加到源文件,将get_max_val_point.h添加到头文件;

然后修改txt文件路径，执行即可。

在命令行会弹出结果的提示。