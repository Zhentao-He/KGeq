function res=swal(y,l,m,s)
% computes values for normalized spin-weighted associated Legendre (swaL)
% 调用格式swal(y,l,m,s)
% y 是向量.
% l,m,s are ints.
% 输出结果res与y同size
alpha=abs(m-s);
beta=abs(m+s);
if mod(alpha+beta,2)==1
    msg="The entered m and s are incorrect!";
    error(msg)
end
n=l-(alpha+beta)/2;
if n<0
    res=0;
    return 
end
norm=(-1)^max(m,-s)*sqrt(...
    (2*n + alpha + beta + 1)* 2^(-alpha-beta-1) * factorial(n) ...
    *factorial(n + alpha + beta) ...
    /factorial( n + alpha ) ...
    /factorial( n + beta)...
    );
% 当y是矢量时，jacobiP(n,alpha,beta,y)只接受标量n,alpha,beta
% 或者与y相同size的vec
res=norm * (1-y).^(alpha/2) .* (1+y).^(beta/2) .* jacobiP(n,alpha,beta,y);

