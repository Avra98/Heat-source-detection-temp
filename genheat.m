function [source,meas]= genheat(Nimage,k,n,Ncount,alpha)
A=makeA(alpha,n);
source=zeros(n,n,Nimage);
meas=source;
for (i=1:Nimage)
    p1 = randperm(n,k);
    p2 = randperm(n,k);
    for(j=1:k)
        source(p1(1,j),p2(1,j),i)=100;
    end
    meas(:,:,i)=diffuse2d(A,Ncount,source(:,:,i));
end
save('Heatdata.m','source','meas')
