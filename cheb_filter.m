function F=cheb_filter(f)
% 本函数对Chebshev多项式进行滤波

% to coef
re_coef=myFCT1(real(f));
im_coef=myFCT1(imag(f));

% 判断size
typeflag=1;
if isrow(f)
    N=length(f)-1;
    typeflag=0; %代表是行向量
else
    [N,~]=size(f);
    N=N-1;
end

% 指数过滤
A = -40; p = 16;
filter_exp = exp(A*( (0:N)/N ).^p);
if typeflag==1
    filter_exp = filter_exp';
end
re_coef = filter_exp.*re_coef;
im_coef = filter_exp.*im_coef;

% to real 
F=(myFCT1(re_coef) + 1i*myFCT1(im_coef))*2/N;

