function res=laplacian_D(y,ny,weights,m,s,max_l)
% 本函数计算ny*ny spin-weighted s Laplace-Beltrami 求导矩阵
% weights和y都是行向量
% 本函数得到的求导矩阵res右乘函数值行向量得到其导数
res=zeros(ny);
f=@(l) -(l-s)*(l+s+1) ...
        * swal(y,l,m,s)...
        .*( weights' .* swal(y',l,m,s) ); %在点乘中，列矢量的指标为点乘后矩阵的行指标
for l=max(abs(s),abs(m)):max_l
    res=res+f(l);
end
