function [beta1]=l1opt(M,meas,eps)
cvx_begin 
    variable beta1(1000,1)
   minimize(norm(beta1,1))
    subject to
      norm(meas-M*beta1,2)<=eps
cvx_end
