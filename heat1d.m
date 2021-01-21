Dx=1; % step size
%DY=0.1;
Nx=1000;
%Ny=5;
X=1:Dx:Nx;
%Y=0:DY:Ny;=127; % arbitrary thermal diffusivity (units in mm^2/sec)
alpha_silver=165.6;
alpha_al=97;
alpha_si=88;
alpha_fe=23;
alpha_air=19;
alpha=35;
%--boundary conditions----------------------------------------------------
U(1:Nx) = 0 ;
s(1:1000)=1/(2*alpha);

Ncount=1500;

%%Consider Nx interior points, two boundary points at the extreme. 
A=zeros(Nx,Nx);

%Filling up intereior pointa
for i=2:Nx-1
        A(i,i)=1-2*s(i);
        A(i,i-1)=s(i);
        A(i,i+1)=s(i);
end

%%Dirchlet boundary conditions 
A(1,:)=zeros(1,Nx);
A(Nx,:)=zeros(1,Nx);
A(:,1)=zeros(1,Nx);
A(:,Nx)=zeros(1,Nx);

sigma=0.5;
k=3;
trails=10;
Npred=zeros(1,35,trails);
for (i=1:trails)
   U(1:Nx) = 0 ;
   U(randperm(Nx-100,k)+100)=100;
   U=abs(U);
   Unit=U;
   meas=(A^Ncount)*Unit';
  % meas=abs(meas+sigma*randn(size(meas)));
   
   [Npred(:,:,i),~] = Npredict1D(meas,500,A);
end