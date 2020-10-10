function k=diagsolve(p,u,u0,c)
%a=p'*u*u0'*p;
%b=p'*u0*u0'*p;%Solving least sqaures instead of direct equality because u may not be in R(A^N).
a=p'*u;
b=p'*u0; % Direct inversion
n=size(u);
for i=1:n
    d(i)=a(i)/b(i);
end
c = c';

[k,~]=polyfit(log(abs(c(1:30))),log(abs(d(1:30))),1);
k=k(1);
end