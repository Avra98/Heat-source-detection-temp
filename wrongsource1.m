function [source,wsource]=wrongsource1(n,alpha,k,N)
A=makeA1d(alpha,n);
Unit=zeros(n,1);
source=zeros(n,N);
eps=0.1;
for i=1:N
   Unit(1:n,1)=0;
   
   Unit(randperm(n-2,k)+1,1)=20*randn(k,1);
   source(:,i)=cumsum(Unit(:,1));
   if sum(source(:,i)<0)<sum(source(:,i)>0)
   source((source(:,i)<0),i)=0;
   else
        source((source(:,i)>0),i)=0;
        source(:,i) = -source(:,i);
   end
   source(:,i) = source(:,i)/norm(source(:,i));
   offset=randi([1,3000]);
   wsource(:,i)=(A^(offset))*(source(:,i));
    wsource(:,i)=(wsource(:,i))/norm(wsource(:,i));
end
save('wrongsource1.mat','source','wsource');
%%
subplot(2,1,1)
plot(source(:,8))
subplot(2,1,2)
plot(wsource(:,8))
