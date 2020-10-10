function [N0,upred]= Npredict(meas,N0,A)
% N0 is the intial guess of N.
[p,d,~]=svd(A);
for i=1:20
M=(A^(N0));
penalty = 'MCP';
if i<10
   upred = lsq_sparsereg(M,meas,2,'penalty',penalty,'penparam',500);
else
    upred = lsq_sparsereg(M,meas,2,'penalty',penalty,'penparam',2000);
end
N0 =diagsolve_new(p,meas,upred,diag(d));
disp(['N0:',num2str(N0)])
figure(1);plot(upred); pause(.2);
end
%[n1,~]=findn(D,d);
% M1=S*(A^(n1));
% upred1 = lsq_sparsereg(M1,meas,2,'penalty',penalty,'penparam',500);
% D1=diagsolve(p,upred1,u0); 
% [n2,~]=findn(D1,d);
% M2=S*(A^(n2));
% upred2 = lsq_sparsereg(M2,meas,2,'penalty',penalty,'penparam',500);
% D2=diagsolve(p,upred2,u0); 
% [n3,~]=findn(D2,d);
% N_opt=n3; 
end
