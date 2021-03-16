function [A]=makeA(alpha,n)
A=zeros(n,n);
s(1:n)=1/(2*alpha);

%Filling up intereior points
for i=2:n-1
        A(i,i)=1-2*s(i);
        A(i,i-1)=s(i);
        A(i,i+1)=s(i);
end

%%Dirchlet boundary conditions 
A(1,:)=zeros(1,n);
A(n,:)=zeros(1,n);
A(:,1)=zeros(1,n);
A(:,n)=zeros(1,n);
end
