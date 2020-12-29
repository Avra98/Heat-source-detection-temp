function [A]=makeA(alpha,n)
s=zeros(n,n);
s(1:n,1:n)=1/(2*alpha);
%s(10:20,10:20)= 1/(2*alpha_al);
sx= reshape(s,[n^2,1]);
sy=sx;
A=zeros(n^2,n^2);

 for(i=2:(n^2)-1)
     A(i,i)=(1-2*sx(i)-2*sy(i));
     A(i,i-1)=sx(i);
     A(i,i+1)=sx(i);
     A(i,mod(i+n-1,n^2)+1)=sy(i);
     A(i,mod(i-n-1,n^2)+1)=sy(i);
 end
 
%Set boundary conditions along x-axis
 A([1:n],:)=0; A([n^2-n+1:n^2],:)=0; A(:,[1:n])=0; A(:,[n^2-n+1:n^2])=0;
 
 %Set boundary conditions along y-axis
 for (i=1:n)
     A(i*n,:)=0;
     A(i*n-1,:)=0;
 end
 
 
 