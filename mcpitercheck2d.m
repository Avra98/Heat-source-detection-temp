N = 51;  
Ncount=400;

%--boundary conditions----------------------------------------------------
U(1:N,1:N) = 0;
U(1,1:N) = 0;  
U(N,1:N) = 0;   
U(1:N,1) = 0;  
U(1:N,N) = 0;  
%--initial condition------------------------------------------------------

%-------------------------------------------------------------------------
alpha_gold=127; % arbitrary thermal diffusivity 
alpha_silver=165.6;
alpha_al=97;
alpha_si=88;
alpha_fe=23;
alpha_air=19;
%DT = DX^2/(2*alpha); % time step 
%%Specifying regions
s=zeros(N,N);
s(1:N,1:N)=1/(2*alpha_fe);
%s(10:20,10:20)= 1/(2*alpha_al);
sx= reshape(s,[N^2,1]);
sy=sx;

Unit=U;
%%Reshaping initial 2d heat map into 1d
Unitv= reshape(Unit,[N^2,1]);

%Make the circulant A matrix
 %t=[(1-2*sx-2*sy) sx zeros(1,N-2) sy zeros(1,N^2-2*N-1) sy zeros(1,N-2) sx];
 %A = gallery('circul',t);
 A=zeros(N^2,N^2);
 for(i=2:(N^2)-1)
     A(i,i)=(1-2*sx(i)-2*sy(i));
     A(i,i-1)=sx(i);
     A(i,i+1)=sx(i);
     A(i,mod(i+N-1,N^2)+1)=sy(i);
     A(i,mod(i-N-1,N^2)+1)=sy(i);
 end
 
     
 
 %Set boundary conditions along x-axis
 A([1:N],:)=0; A([N^2-N+1:N^2],:)=0; A(:,[1:N])=0; A(:,[N^2-N+1:N^2])=0;
 
 %Set boundary conditions along y-axis
 for (i=1:N)
     A(i*N,:)=0;
     A(i*N-1,:)=0;
 end
 %A=sparse(A);
 
 
 %%
 
 
 k=2;
trails=10;
Npred=zeros(1,25,trails);
for (i=1:trails)
   U(1:N,1:N)=0; 
   p1=randperm(N-2,k);
   p2=randperm(N-2,k);
   
   for (j=1:k)
   U(p1(1,j)+1,p2(1,j)+1)=1000;
   end
   
   U=abs(U);
   Unit=U;
   Unitv= reshape(Unit,[N^2,1]);
   meas=(A^Ncount)*Unitv;
   [Npred(:,:,i),~] = Npredict_new2D(meas,70,A);
end