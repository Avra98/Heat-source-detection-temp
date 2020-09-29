function [psi,p]= MCP(A,lambda,x)
n=size(x,1);
for(i=2:n-1)
 
        psi(i)= lambda-(abs(x(i)<=(sqrt(2*lambda)/norm(A(:,i))))).*(((norm(A(:,i))^2)/2)*((abs(x(i))-(sqrt(2*lambda)/norm(A(:,i))))^2));
        
end

p=sum(psi);

