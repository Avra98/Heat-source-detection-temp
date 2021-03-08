function [source,wsource]=wrongsource(n,alpha,k,N)
A=makeA1d(alpha,n);
Unit=zeros(n,1);
source=zeros(n,N);
eps=0.1;
for (i=1:N)
   Unit(1:n,1)=0;
   
   Unit(randperm(n-2,k)+1,1)=100*rand(k,1);
   source(:,i)=Unit(:,1);
   offset=randi([1,100]);
   wsource(:,i)=(A^(offset))*Unit;
end
save('wrongsource.mat','source','wsource');
%%
subplot(2,1,1)
plot(source(:,8))
subplot(2,1,2)
plot(wsource(:,8))
