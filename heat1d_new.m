%------------------------------------------------------------------------
%--- Heat Equation in two dimensions-------------------------------------
%--- Solves Ut=alpha*(Uxx+Uyy)-------------------------------------------
%------------------------------------------------------------------------
clc;
close all;
clear all;
%--dimensions...........................................................
N = 1000;  
DX=0.1; % step size
%DY=0.1;
Nx=1000; 
%Ny=5;
X=0:DX:Nx;
%Y=0:DY:Ny;
alpha=5; % arbitrary thermal diffusivity 
%--boundary conditions----------------------------------------------------
U(1:N) = 0 ;
%U(1,1:N) = 0; 
%U(N,1:N) = 0;  
%U(1:N,1) = 0;  
%U(1:N,N) = 0;  
%--initial condition------------------------------------------------------
U(150)=1000;% a heated patch at the center
U(500)=1000;
U(750)=1000;
%U=U+10*rand(1,1000);
U=abs(U);

Umax=max(max(U)); 
%-------------------------------------------------------------------------
DT = DX^2/(2*alpha); % time step 
s=1/(2*alpha);
M=2000; % maximum number of allowed iteration
%---finite difference scheme----------------------------------------------
fram=0;
Ncount=M;
loop=1;
t=[(1-2*s) s zeros(1,N-3) s];
A = gallery('circul',t);
A(1,:)=zeros(1,N);
A(N,:)=zeros(1,N);
A(:,1)=zeros(1,N);
A(:,N)=zeros(1,N);
Unit=U;
%%
Afin=A^Ncount;
% m=N;
% t=randperm(N,m);
% S=zeros(m,N);
% for(i=1:m)
%      S(i,t(i))=1;
% end
M=Afin; 
meas=M*Unit';
%%
% estimate u0 and N0 together
[N,upred] = Npredict_new(meas,500,A);
figure;plot([upred,Unit']);
legend('predicted u0','true u0')
