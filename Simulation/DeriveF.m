function [f_total,r_total] = DeriveF(a_total,l_total,K,L)
%% initialization 
f_total = ones(127,1);  
r_total = zeros(127,K); 
%%
for j = 1:127 %total 127 possible {r_0,r_1,...,r_6} 
   r_subset = de2bi(j,K,'left-msb');
   r_total(j,:) = r_subset;
   for i = 1:L
       tmp1 = r_subset(1)*a_total(i,1)+ r_subset(2)*a_total(i,2)+r_subset(3)*a_total(i,3)+r_subset(4)*a_total(i,4)+r_subset(5)*a_total(i,5)+r_subset(6)*a_total(i,6)+r_subset(7)*a_total(i,7);
       tmp2 = mod(tmp1,2);
       if tmp2 == 0
           f_i = exp(l_total(i))/(1+exp(l_total(i))); % probability of P(y_i =0);
       else
           f_i = 1/(1+exp(l_total(i))); % probability of P(y_i = 1);
       end
       f_total(j) = f_total(j)*f_i;
   end
end