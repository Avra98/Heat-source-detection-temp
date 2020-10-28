clc;
close all;
clear all;
%--dimensions...........................................................
N = 51;  
Ncount=500;
DX=0.1; % step size
DY=0.1;
Nx=5; 
Ny=5;
X=0:DX:Nx;
Y=0:DY:Ny;
alpha=5; % arbitrary thermal diffusivity 
%--boundary conditions----------------------------------------------------
U(1:N,1:N) = 0 ;
U(1,1:N) = 0; 
U(N,1:N) = 0;  
U(1:N,1) = 0;  
U(1:N,N) = 0;  
%--initial condition------------------------------------------------------
U(13,15)=1000;
U(39,39)=1000;
Unit=U;
% a heated patch at the center
Umax=max(max(U));
%-------------------------------------------------------------------------
DT = DX^2/(2*alpha); % time step 
sx=1/(2*alpha);
sy=1/(2*alpha);

%%Reshaping initial 2d heat map into 1d
Unitv= reshape(Unit,[N^2,1]);

%Make the circulant A matrix
 t=[(1-2*sx-2*sy) sx zeros(1,N-2) sy zeros(1,N^2-2*N-1) sy zeros(1,N-2) sx];
 A = gallery('circul',t);
 
 %Set boundary conditions along x-axis
 A([1:N],:)=0; A([N^2-N+1:N^2],:)=0; A(:,[1:N])=0; A(:,[N^2-N+1:N^2])=0;
 
 %Set boundary conditions along y-axis
 for (i=1:N)
     A(i*N,:)=0;
     A(i*N-1,:)=0;
 end
 A=sparse(A);
 
 %Forward operator 
 ut=(A^(Ncount))*Unitv;
 %Reshaping
 utm=reshape(ut,[N,N]);
 utm=utm';
 
 %%Recovering initial U0 from measurements by l1 minimisation. 
 cvx_begin 
    variable beta1(N^2,1)
   minimize norm(beta1,1)
    subject to
      norm(ut-(A^(Ncount))*beta1,2)<=1
cvx_end


 upred=reshape(beta1,[N,N]);
 upred=upred';
 subplot(3,1,1)
imshow(utm)
title('Measurement, heat map after diffusion')
subplot(3,1,2)
imshow(U)
title('Original heat source location')
subplot(3,1,3)
imshow(upred)
title('Reconstrucetd initial condition from measurements')