function [Npred,upred]= Npredict(meas,N0,A)
% N0 is the intial guess of N.
%[p,d,~]=svd(A);
%eps=0.1; %l2 ball constraint 
iter=40;
Npred=zeros(1,iter);
for i=1:iter
M=(A^(N0));
penalty = 'MCP';
%penalty= 'POWER';

   upred = lsq_sparsereg(M,meas,1,'penalty',penalty,'penparam',100);
   %upred = lsq_sparsereg(M,meas,1,'penalty',penalty,'penparam',1);
    %upred= l1opt(M,meas,0.5);

%%Either use diagsolve or iterative l2min. For noise measurements, use
%%iterative l2min. Comment out the other one. 
%N0 =diagsolve_new(p,meas,upred,diag(d));  %%Use for diagsolve using polyfit
 [~,N0]=l2min(meas,upred,A);      %% Use for iterative l2 min of(A^n.u(0)-u(n))
 Npred(i)=N0;
%up=reshape(upred,[71,71]);
%up=up';
disp(['N0:',num2str(N0)]);
figure(1);plot(upred); pause(.2);
end
%%
 p = plot([1:iter],1000*ones(1,iter),[1:iter],Npred);
 p(1).LineWidth = 2;
 p(2).LineWidth = 4;
 legend('Original Ncount ','Predicted Ncount')
 xlabel('Iterations','fontweight','bold','fontsize',16)
 ylabel('Diffusion time (Ncount)','fontweight','bold','fontsize',16)
 ylim([100 5500])
 title('Ncount','fontweight','bold','fontsize',16)
 ax = gca;
 ax.YAxis.FontWeight = 'bold';
 ax.XAxis.FontWeight = 'bold';
 ax.XAxis.FontSize = 18;
 ax.YAxis.FontSize = 18;
end
