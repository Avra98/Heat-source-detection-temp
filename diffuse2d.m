function [meas]= diffuse2d(A,Ncount,source)
n=size(source,1);
Unitv= reshape(source,[n^2,1]);
ut=(A^(Ncount))*Unitv;
 utm=reshape(ut,[n,n]);
 meas=utm;
 