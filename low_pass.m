function res=low_pass(y,ny,weights,m,s,max_l)
% 本函数计算ny*ny 滤波矩阵
% weights和y都是行向量
% 本函数得到的滤波矩阵res 需右乘函数值行向量 其作用
res=zeros(ny);
A=-40;
p=16;

% 左乘列向量的版本
% f=@(l) exp(A*(l/max_l)^p) * swal(y,l,m,s)...
%         .*( weights .* swal(y',l,m,s) );

% 右乘行向量的版本
f=@(l) exp(A*(l/max_l)^p) * swal(y,l,m,s)...
        .*( weights' .* swal(y',l,m,s) );

for l=max( [abs(s),abs(m)] ):max_l
    res=res+f(l);
end