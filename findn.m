function [n,q]=findn(D,d)

a=diag(D);
b=diag(d);
j=find(b<0.999 & b>0 );
m=size(j,1);
for (l=1:m)
    q(l)=log(a(j(l)))/log(b(j(l)));
end
n=mean(q(1:10));
end
