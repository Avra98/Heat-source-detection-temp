D=D%------------------------------------------------------------------------
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
Ncount=0;
loop=1;
t=[(1-2*s) s zeros(1,N-3) s];
A = gallery('circul',t);
A(1,:)=zeros(1,N);
A(N,:)=zeros(1,N);
A(:,1)=zeros(1,N);
A(:,N)=zeros(1,N);
Unit=U;
while loop==1;
   ERR=0; 
   U_old = U;
for i = 2:N-1

   Residue=s*(U_old(i+1)-2*U_old(i)+U_old(i-1))+ U_old(i)-U(i);
   ERR=ERR+abs(Residue);
  U(i)=U(i)+Residue;
end
if(ERR>=0.000001*Umax)  % allowed error limit is 1% of maximum temperature
    Ncount=Ncount+1;
         if (mod(Ncount,10)==0) % displays movie frame every 10 time steps
              fram=fram+1;
              plot(U);
              %axis([1 N])
              h=gca; 
              get(h,'FontSize') 
              set(h,'FontSize',12)
              colorbar('location','eastoutside','fontsize',12);
              xlabel('X','fontSize',12);
              %ylabel('Y','fontSize',12);
              title('Heat Diffusion','fontsize',12);
              fh = figure(1);
             set(fh, 'color', 'white'); 
            F=getframe;
         end
 
 %--if solution do not converge in 2000 time steps------------------------
 
    if(Ncount>M)
        loop=0;
        disp(['solution do not reach steady state in ',num2str(M),...
            'time steps'])
    end
    
 %--if solution converges within 2000 time steps..........................   
    
else
    loop=0;
    disp(['solution reaches steady state in ',num2str(Ncount) ,'time steps'])
end
end
Afin=A^Ncount;
Afin2=A^((Ncount-1)/2);
Ufin=Afin*Unit';
%------------------------------------------------------------------------
%--display a movie of heat diffusion------------------------------------
subplot(2,1,1)
 movie(F,fram,1)
 subplot(2,1,2)
 plot(Ufin)
 %%%%Measurements%%
 m=N;
 t=randperm(N,m);
 S=zeros(m,N);
 for(i=1:m)
     S(i,t(i))=1;
 end
 M=S*Afin;
% M2=S*Afin2;
 meas=M*Unit';
% upred= OMPerr(M,meas,0.01); 
% upred2=OMPerr(M2,meas,0.01); 
 %upred=abs(upred);
 %upred2=abs(upred2);
 %uon=M\meas;
 cvx_begin
  variable upred(1000,1)
  minimize( norm( upred, 1 ) )
  subject to
    meas==M*upred;
cvx_end
 

%%

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
u0=Unit';
Ufin = (A^1000)*u0;

M1=S*(A^1000);
penalty = 'MCP';           % set penalty function to lasso
penparam = 500;
%upred2 =  lsq_sparsereg(M,meas,penparam,'penalty',penalty);
[upred5] = lsq_sparsereg(M1,meas,1,'penalty',penalty,'penparam',500);
 subplot(2,2,1)
 plot(Unit)
 subplot(2,2,2)
 plot(Ufin)
  subplot(2,2,3)
 plot(upred)
 subplot(2,2,4)
 plot(upred5)
 %%
     
%------END---------------------------------------------------------------
