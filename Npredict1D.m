function [N0,upred]= Npredict1D(meas,N0,A)
% N0 is the intial guess of N.
%[p,d,~]=svd(A);
%eps=0.1; %l2 ball constraint 
iter=35;
Npred=zeros(1,iter);
for i=1:iter
M=(A^(N0));
penalty = 'MCP';
%penalty= 'POWER'; 
if i<20
 upred = lsq_sparsereg(M,meas,1,'penalty',penalty,'penparam',100);
else
 upred = lsq_sparsereg(M,meas,.1,'penalty','POWER','penparam',0.01);
end
    %upred= l1opt(M,meas,0.1);
ind = find(abs(upred)>1e-3);
M1 = M(:,ind);

upred(upred~=0) = lsqnonneg(M1,meas); 
%%Either use diagsolve or iterative l2min. For noise measurements, use
%%iterative l2min. Comment out the other one. 
%N0 =diagsolve_new(p,meas,upred,diag(d));  %%Use for diagsolve using polyfit
 [~,N0]=l2min(meas,upred,A);      %% Use for iterative l2 min of(A^n.u(0)-u(n))
 Npred(i)=N0;
%up=reshape(upred,[sqrt(size(meas,1)),sqrt(size(meas,1))]);
% up=up';
disp(['N0:',num2str(N0)]);
figure(1);plot(upred); pause(.2);
end

end