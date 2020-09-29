function [N_opt]= Npredict(meas,N0,A,S,u0)
% N0 is the intial guess of N.
[p,d,~]=svd(A);
M=S*(A^(N0));
penalty = 'MCP';
upred = lsq_sparsereg(M,meas,2,'penalty',penalty,'penparam',500);
D=diagsolve(p,upred,u0);
[n1,~]=findn(D,d);
M1=S*(A^(n1));
upred1 = lsq_sparsereg(M1,meas,2,'penalty',penalty,'penparam',500);
D1=diagsolve(p,upred1,u0); 
[n2,~]=findn(D1,d);
M2=S*(A^(n2));
upred2 = lsq_sparsereg(M2,meas,2,'penalty',penalty,'penparam',500);
D2=diagsolve(p,upred2,u0); 
[n3,~]=findn(D2,d);
N_opt=n3; 
