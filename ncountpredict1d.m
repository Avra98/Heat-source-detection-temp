Dx=1; % step size
%DY=0.1;
Nx=500;
%Ny=5;
X=1:Dx:Nx;
%Y=0:DY:Ny;=127; % arbitrary thermal diffusivity (units in mm^2/sec)
alpha_silver=165.6;
alpha_al=97;
alpha_si=88;
alpha_fe=23;
alpha_air=19;
alpha=25;
%--boundary conditions----------------------------------------------------
U(1:Nx) = 0 ;
%U(1,1:N) = 0; 
%U(N,1:N) = 0;  
%U(1:N,1) = 0;  
%U(1:N,N) = 0;  
%--initial condition------------------------------------------------------
%U(150)=100;% a heated patch at the center
%U(500)=100;
%U(750)=100;
%U=U+10*rand(1,1000);
k=5;    
p1 = randperm(Nx-15,k);  %%Dont take sources at the boundary
    %p2 = randperm(n-2,k);

    for(j=1:k)
        U(p1(1,j)+10)=100;
    end
U=abs(U);
%U=U';
Unit=U;

%-------------------------------------------------------------------------
%DT = DX^2/(2*alpha); % time step 
%s=1/(2*alpha);
%sg=1/(2*alpha_gold);
sf=1/(2*alpha_fe);
ss=1/(2*alpha_silver);

s(1:1000)=1/(2*alpha);
%s(301:1000)=ss;
A=makeA1d(alpha,Nx);
%---finite difference scheme----------------------------------------------
fram=0;
Ncount=500;

%%Consider Nx interior points, two boundary points at the extreme. 

M=A^Ncount;
meas=M*Unit'+2*rand(Nx,1);
Unit=Unit';
plot(meas)
%upred=l1optconst(M,meas,eps);
%%
[Npred,beta] = Npredict_new(meas,1000,A);
