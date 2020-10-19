function k=l2min(u,u0,A)
n=size(u0,1);
B=eye(n);
for(i=1:5000)
     e(i)=norm(u-B*u0);
     B=B*A;
end
[err,k]=min(e);
