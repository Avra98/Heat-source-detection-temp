function [err,k]=l2min(u,u0,A)
n=size(u0,1);
B=eye(n);
u0 = A*u0;
for(i=1:5000)
     e(i)=norm(u-u0);
     u0=A*u0;
end
[err,k]=min(e);
