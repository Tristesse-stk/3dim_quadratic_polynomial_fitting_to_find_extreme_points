clear all;
close all;
clc

%对于类似的txt文件，不含有字符，只有数字
txt_path = 'D:\Datasets\origin\拟合\1.txt';
num_points = 125;
threshold = 0.8;
data=load(txt_path);
x = data(:,1);
y = data(:,2);
z = data(:,3);
t = data(:,4);
x_max = max(x);
y_max = max(y);
z_max = max(z);

x_min = min(x);
y_min = min(y);
z_min = min(z);

x_center = (x_max+x_min)/2;
y_center = (y_max+y_min)/2;
z_center = (z_max+z_min)/2;

x_norm = (x-x_center)/((x_max-x_min)/2);
y_norm = (y-y_center)/((y_max-y_min)/2);
z_norm = (z-z_center)/((z_max-z_min)/2);

[t_sort,t_sort_id] = sort(t);
id_effect = 0;
for id=1:num_points
    if((t_sort(id,1)<threshold) && (t_sort(id+1,1)>=threshold))
        id_effect = id+1;
        break;
    end
end

t_sort_id_select = t_sort_id(id_effect:size(t_sort_id),:);
x_norm_select = x_norm(t_sort_id_select,:);
y_norm_select = y_norm(t_sort_id_select,:);
z_norm_select = z_norm(t_sort_id_select,:);
t_select = t(t_sort_id_select,:);

num_parameters = 10;
A = zeros(num_parameters,num_parameters);
B = zeros(num_parameters,1);
A(10,1)=sum(x_norm_select.*x_norm_select);
A(10,2)=sum(y_norm_select.*y_norm_select);
A(10,3)=sum(z_norm_select.*z_norm_select);
A(10,4)=sum(x_norm_select.*y_norm_select);
A(10,5)=sum(x_norm_select.*z_norm_select);
A(10,6)=sum(y_norm_select.*z_norm_select);
A(10,7)=sum(x_norm_select);
A(10,8)=sum(y_norm_select);
A(10,9)=sum(z_norm_select);
A(10,10)=size(y_norm_select,1);
B(10,1)= sum(t_select);

A(9,1)=sum(x_norm_select.*x_norm_select.*z_norm_select);
A(9,2)=sum(y_norm_select.*y_norm_select.*z_norm_select);
A(9,3)=sum(z_norm_select.*z_norm_select.*z_norm_select);
A(9,4)=sum(x_norm_select.*y_norm_select.*z_norm_select);
A(9,5)=sum(x_norm_select.*z_norm_select.*z_norm_select);
A(9,6)=sum(y_norm_select.*z_norm_select.*z_norm_select);
A(9,7)=sum(x_norm_select.*z_norm_select);
A(9,8)=sum(y_norm_select.*z_norm_select);
A(9,9)=sum(z_norm_select.*z_norm_select);
A(9,10)=sum(z_norm_select);
B(9,1)= sum(t_select.*z_norm_select);

A(8,1)=sum(x_norm_select.*x_norm_select.*y_norm_select);
A(8,2)=sum(y_norm_select.*y_norm_select.*y_norm_select);
A(8,3)=sum(z_norm_select.*z_norm_select.*y_norm_select);
A(8,4)=sum(x_norm_select.*y_norm_select.*y_norm_select);
A(8,5)=sum(x_norm_select.*z_norm_select.*y_norm_select);
A(8,6)=sum(y_norm_select.*z_norm_select.*y_norm_select);
A(8,7)=sum(x_norm_select.*y_norm_select);
A(8,8)=sum(y_norm_select.*y_norm_select);
A(8,9)=sum(z_norm_select.*y_norm_select);
A(8,10)=sum(y_norm_select);
B(8,1)= sum(t_select.*y_norm_select);


A(7,1)=sum(x_norm_select.*x_norm_select.*x_norm_select);
A(7,2)=sum(y_norm_select.*y_norm_select.*x_norm_select);
A(7,3)=sum(z_norm_select.*z_norm_select.*x_norm_select);
A(7,4)=sum(x_norm_select.*y_norm_select.*x_norm_select);
A(7,5)=sum(x_norm_select.*z_norm_select.*x_norm_select);
A(7,6)=sum(y_norm_select.*z_norm_select.*x_norm_select);
A(7,7)=sum(x_norm_select.*x_norm_select);
A(7,8)=sum(y_norm_select.*x_norm_select);
A(7,9)=sum(z_norm_select.*x_norm_select);
A(7,10)=sum(x_norm_select);
B(7,1)= sum(t_select.*x_norm_select);

A(6,1)=sum(x_norm_select.*x_norm_select.*y_norm_select.*z_norm_select);
A(6,2)=sum(y_norm_select.*y_norm_select.*y_norm_select.*z_norm_select);
A(6,3)=sum(z_norm_select.*z_norm_select.*y_norm_select.*z_norm_select);
A(6,4)=sum(x_norm_select.*y_norm_select.*y_norm_select.*z_norm_select);
A(6,5)=sum(x_norm_select.*z_norm_select.*y_norm_select.*z_norm_select);
A(6,6)=sum(y_norm_select.*z_norm_select.*y_norm_select.*z_norm_select);
A(6,7)=sum(x_norm_select.*y_norm_select.*z_norm_select);
A(6,8)=sum(y_norm_select.*y_norm_select.*z_norm_select);
A(6,9)=sum(z_norm_select.*y_norm_select.*z_norm_select);
A(6,10)=sum(y_norm_select.*z_norm_select);
B(6,1)= sum(t_select.*y_norm_select.*z_norm_select);

A(5,1)=sum(x_norm_select.*x_norm_select.*x_norm_select.*z_norm_select);
A(5,2)=sum(y_norm_select.*y_norm_select.*x_norm_select.*z_norm_select);
A(5,3)=sum(z_norm_select.*z_norm_select.*x_norm_select.*z_norm_select);
A(5,4)=sum(x_norm_select.*y_norm_select.*x_norm_select.*z_norm_select);
A(5,5)=sum(x_norm_select.*z_norm_select.*x_norm_select.*z_norm_select);
A(5,6)=sum(y_norm_select.*z_norm_select.*x_norm_select.*z_norm_select);
A(5,7)=sum(x_norm_select.*x_norm_select.*z_norm_select);
A(5,8)=sum(y_norm_select.*x_norm_select.*z_norm_select);
A(5,9)=sum(z_norm_select.*x_norm_select.*z_norm_select);
A(5,10)=sum(x_norm_select.*z_norm_select);
B(5,1)= sum(t_select.*x_norm_select.*z_norm_select);

A(4,1)=sum(x_norm_select.*x_norm_select.*x_norm_select.*y_norm_select);
A(4,2)=sum(y_norm_select.*y_norm_select.*x_norm_select.*y_norm_select);
A(4,3)=sum(z_norm_select.*z_norm_select.*x_norm_select.*y_norm_select);
A(4,4)=sum(x_norm_select.*y_norm_select.*x_norm_select.*y_norm_select);
A(4,5)=sum(x_norm_select.*z_norm_select.*x_norm_select.*y_norm_select);
A(4,6)=sum(y_norm_select.*z_norm_select.*x_norm_select.*y_norm_select);
A(4,7)=sum(x_norm_select.*x_norm_select.*y_norm_select);
A(4,8)=sum(y_norm_select.*x_norm_select.*y_norm_select);
A(4,9)=sum(z_norm_select.*x_norm_select.*y_norm_select);
A(4,10)=sum(x_norm_select.*y_norm_select);
B(4,1)= sum(t_select.*x_norm_select.*y_norm_select);

A(3,1)=sum(x_norm_select.*x_norm_select.*z_norm_select.*z_norm_select);
A(3,2)=sum(y_norm_select.*y_norm_select.*z_norm_select.*z_norm_select);
A(3,3)=sum(z_norm_select.*z_norm_select.*z_norm_select.*z_norm_select);
A(3,4)=sum(x_norm_select.*y_norm_select.*z_norm_select.*z_norm_select);
A(3,5)=sum(x_norm_select.*z_norm_select.*z_norm_select.*z_norm_select);
A(3,6)=sum(y_norm_select.*z_norm_select.*z_norm_select.*z_norm_select);
A(3,7)=sum(x_norm_select.*z_norm_select.*z_norm_select);
A(3,8)=sum(y_norm_select.*z_norm_select.*z_norm_select);
A(3,9)=sum(z_norm_select.*z_norm_select.*z_norm_select);
A(3,10)=sum(z_norm_select.*z_norm_select);
B(3,1)= sum(t_select.*z_norm_select.*z_norm_select);

A(2,1)=sum(x_norm_select.*x_norm_select.*y_norm_select.*y_norm_select);
A(2,2)=sum(y_norm_select.*y_norm_select.*y_norm_select.*y_norm_select);
A(2,3)=sum(z_norm_select.*z_norm_select.*y_norm_select.*y_norm_select);
A(2,4)=sum(x_norm_select.*y_norm_select.*y_norm_select.*y_norm_select);
A(2,5)=sum(x_norm_select.*z_norm_select.*y_norm_select.*y_norm_select);
A(2,6)=sum(y_norm_select.*z_norm_select.*y_norm_select.*y_norm_select);
A(2,7)=sum(x_norm_select.*y_norm_select.*y_norm_select);
A(2,8)=sum(y_norm_select.*y_norm_select.*y_norm_select);
A(2,9)=sum(z_norm_select.*y_norm_select.*y_norm_select);
A(2,10)=sum(y_norm_select.*y_norm_select);
B(2,1)= sum(t_select.*y_norm_select.*y_norm_select);

A(1,1)=sum(x_norm_select.*x_norm_select.*x_norm_select.*x_norm_select);
A(1,2)=sum(y_norm_select.*y_norm_select.*x_norm_select.*x_norm_select);
A(1,3)=sum(z_norm_select.*z_norm_select.*x_norm_select.*x_norm_select);
A(1,4)=sum(x_norm_select.*y_norm_select.*x_norm_select.*x_norm_select);
A(1,5)=sum(x_norm_select.*z_norm_select.*x_norm_select.*x_norm_select);
A(1,6)=sum(y_norm_select.*z_norm_select.*x_norm_select.*x_norm_select);
A(1,7)=sum(x_norm_select.*x_norm_select.*x_norm_select);
A(1,8)=sum(y_norm_select.*x_norm_select.*x_norm_select);
A(1,9)=sum(z_norm_select.*x_norm_select.*x_norm_select);
A(1,10)=sum(x_norm_select.*x_norm_select);
B(1,1)= sum(t_select.*x_norm_select.*x_norm_select);

%parm = inv(A)*b;
parm = A\B;
a0 = parm(1,1);a1 = parm(2,1);a2 = parm(3,1);a3 = parm(4,1);a4 = parm(5,1);a5 = parm(6,1);a6 =parm(7,1);a7 = parm(8,1);a8 =parm(9,1);a9 =parm(10,1);

C = zeros(3,3);
D = zeros(3,1);

C(1,1)=2*a0;C(1,2)=a3;C(1,3)=a4;
C(2,1)=a3;C(2,2)=2*a1;C(2,3)=a5;
C(3,1)=a4;C(3,2)=a5;C(3,3)=2*a2;
D(1,1)=-a6;D(2,1)=-a7;D(3,1)=-a8;
max_val_point = C\D;

error_square = 0;
t_pred = zeros(size(z_norm_select));
for i=1:size(z_norm_select,1)
    t_pred(i,1) = a0*x_norm_select(i,1)^2+a1*y_norm_select(i,1)^2+a2*z_norm_select(i,1)^2+a3*x_norm_select(i,1)*y_norm_select(i,1)+a4*x_norm_select(i,1)*z_norm_select(i,1)+a5*y_norm_select(i,1)*z_norm_select(i,1)+a6*x_norm_select(i,1)+a7*y_norm_select(i,1)+a8*z_norm_select(i,1)+a9;
    error_square = error_square+(t_pred(i,1)-t_select(i,1))^2;
end

all_compare = [x_norm_select,y_norm_select,z_norm_select,t_pred,t_select];

x0_real = max_val_point(1)*((x_max-x_min)/2)+x_center
y0_real = max_val_point(2)*((y_max-y_min)/2)+y_center
z0_real = max_val_point(3)*((z_max-z_min)/2)+z_center


