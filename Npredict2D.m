function [N0,upred]= Npredict2D(meas,N0,A)
% N0 is the intial guess of N.
%[p,d,~]=svd(A);
%eps=0.1; %l2 ball constraint 
iter=105;
%Npred=zeros(1,iter);
N1 =N0-0.1;
visited_number= [];
for i=1:iter
M=(A^(N0));
penalty = 'MCP';
%penalty= 'POWER';
i
if ~ismember(N1,visited_number)
   upred = lsq_sparsereg(M,meas,1,'penalty',penalty,'penparam',400);
   visited_number= [visited_number;N1];
   N1 = N0;
else
   upred = lsq_sparsereg(M,meas,10,'penalty','POWER','penparam',0.01);
   %upred = lsq_sparsereg(M,meas,1,'penalty',penalty,'penparam',10);
end
ind = find(upred~=0);
M1 = M(:,ind);
upred(upred~=0) = (M1'*M1)\(M1'*meas); 
%upred(abs(upred)<500)=0;
    %upred= l1opt(M,meas,0.1);

%%Either use diagsolve or iterative l2min. For noise measurements, use
%%iterative l2min. Comment out the other one. 
%N0 =diagsolve_new(p,meas,upred,diag(d));  %%Use for diagsolve using polyfit

 [~,N0]=l2min(meas,upred,A);      %% Use for iterative l2 min of(A^n.u(0)-u(n))
 Npred(i)=N0
 
 end
 
up=reshape(upred,[sqrt(size(meas,1)),sqrt(size(meas,1))]);
% up=up';
disp(['N0:',num2str(N0)]);
figure(1);imagesc(up); pause(.2);


end