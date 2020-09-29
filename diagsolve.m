function[D]=diagsolve(p,u,u0)
%a=p'*u*u0'*p;
%b=p'*u0*u0'*p;%Solving least sqaures instead of direct equality because u may not be in R(A^N).
a=p'*u;
b=p'*u0; % Direct inversion
n=size(u);
for(i=1:n)
    d(i)=a(i)/b(i);
end
D=diag(d);
end
[a,~]=polyfit(log(abs(d(1:20))),log(abs(D(1:20))),1);