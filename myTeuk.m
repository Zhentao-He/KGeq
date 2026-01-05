clear;
%----parameter----
M=1;    % Black hole Mass
a=0.7*M;% Black hole Spin 
L=1;    % compactification_length
% spin weight
s = 0; 
% 磁量子数m
m = 1;
 
% Horizon
rH = M + sqrt(M^2-a^2); 
RH = L^2/rH;

% ----grid----
nR=176;
ny=38;
max_l=23;

[y,weights] = gauss(ny);
y=y.';   %行向量y=-cos(\theta)

[D1R,R]=cheb(nR);
R=(R+1)*RH/2; D1R=D1R/(RH/2); D2R = D1R^2;

[RR,yy] = ndgrid(R,y);
syy = sqrt(1-yy.^2);      % sin(\theta)

% 角方向的求导、滤波矩阵（与s,m有关）
DLaplace = laplacian_D(y,ny,weights,m,s,max_l);
swal_filter = low_pass(y,ny,weights,m,s,max_l);

%--------------------------------------------------------------------------
% 计算方程系数
cTT = 8*M*(2*M - a^2/L^2*R) .* (1+2*M*R/L^2) - a^2*(1-y.^2); % R-y 自动扩充为RR-yy

cTR = -2*(L^2 - (8*M^2-a^2)*RR.^2/L^2 + 4*a^2*M*RR.^3/L^4 );

cRR = -(L^2 - 2*M*RR + a^2*RR.^2/L^2).*RR.^2/L^2;

% 与m,s有关的方程系数
cT_fun = @(s,m) 1i*m*(2*a*(1+4*M*RR/L^2)) ...
    + 2*( 2*M*(-s+ (2+s)*2*M*RR/L^2 -3*a^2*RR.^2/L^4) - a^2*RR/L^2 + 1i*s*a*(-y)); % RR-y 自动扩充为RR-yy

cR_fun = @(s,m) 1i*m*2*a*RR.^2/L^2 ...
    + 2*RR.*(-(1+s)+(s+3)*M*RR/L^2-2*a^2*RR.^2/L^4);

c_fun = @(s,m) 1i*m*2*a*RR/L^2 ...
    +2*( (1+s)*M*RR/L^2 - a^2*RR.^2/L^4);

cT = cT_fun(s,m);

cR = cR_fun(s,m);

c = c_fun(s,m);

%--------------------------------------------------------------------------
clear RR yy
%--------------------------------------------------------------------------
% linear initial data

% bump fuction 
% ----radial----
w = nR/2-20 : nR/2; % 控制bump所在的区域
r = L^2./R( w );
rl=r(1); ru=r(end);
width = ru-rl;
A=0*R;
A( w ) = 1 * ( (r-rl)/width ).^2 .* ((ru-r)/width).^2 ...
    .* exp(-width./(r-rl)-2*width./(ru-r));
% ----angular----
l = 2;
psi_init = A.*swal(y,l,m,s); % initial data of psi4 引力扰动

% ingoing
% components of vec n
n_T = 2 + 4*M*R/L^2;
n_R = R.^2/L^2;

% 辅助演化变量
Q_init = D1R*psi_init;
P_init =  cTT.*( -n_R./n_T .* Q_init ) + cTR .* Q_init + cT .* psi_init;
clear Q_init
%--------------------------------------------------------------------------
% time step for RK4
dt = 9/((nR+1)^2);
tn = floor(0.1*M/dt);
dT = tn*dt;
TL = 60*M; % 总时长
TN = ceil(TL/dT);

T = (0:TN)*dT;

% 开始演化---------------------
eq = @(P,psi)KGeq(P,psi,cTT,cTR,cRR,cT,cR,c,D1R,D2R,DLaplace);
P_sol = zeros(nR+1,ny,TN);
psi_sol = zeros(nR+1,ny,TN);
P_sol(:,:,1) = P_init;
psi_sol(:,:,1) = psi_init;
tic 
P = P_sol(:,:,1);
psi = psi_sol(:,:,1);
for ii = 2:TN+1
    for jj = 1:tn
          [P,psi] = teuk_time_step(eq,dt,P,psi,swal_filter);
    end
    psi_sol(:,:,ii) = psi;
    P_sol(:,:,ii) = P;
end
toc

figure
plot(T,reshape(abs(psi_sol(end,end,:)),[length(T),1]))

%%
function [kP,kphi] = KGeq(P,psi,cTT,cTR,cRR,cT,cR,c,D1R,D2R,DLaplace)
% 本函数根据线性阶的Teuk方程计算演化变量P,psi的时间导数
% 方程系数cT,cR,c和求导矩阵DLaplace与s,m有关，必须在调用此函数时确保一致
    kphi = (P - ( cTR.*(D1R*psi) + cT.*psi ) )./cTT;
    kP = -(cRR.*(D2R*psi) + cR.*(D1R*psi) )  - c.*psi + psi*DLaplace;
end

function [P_next,psi_next] = teuk_time_step(eq,dT,P,psi,swal_filter)
    % teuk time step = RK4 + filter
    % 方程系数cT,cR,c和求导矩阵DLaplace与s,m有关，必须在调用此函数时确保一致
    
    % RK4
    [kP1,kpsi1] = eq(P,psi);
    [kP2,kpsi2] = eq(P+kP1*dT/2,psi+kpsi1*dT/2);
    [kP3,kpsi3] = eq(P+kP2*dT/2,psi+kpsi2*dT/2);
    [kP4,kpsi4] = eq(P+kP3*dT  ,psi+kpsi3*dT);
    P_next = P + dT/6*(kP1+2*kP2+2*kP3+kP4);
    psi_next = psi + dT/6*(kpsi1+2*kpsi2+2*kpsi3+kpsi4);
    
    % filter
    % cheb_filter
    P_next = cheb_filter(P_next);
    psi_next = cheb_filter(psi_next);
    % swal_filter
    P_next = P_next*swal_filter;
    psi_next = psi_next*swal_filter;
end

