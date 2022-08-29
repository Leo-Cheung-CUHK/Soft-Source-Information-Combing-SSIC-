function psdu_llr = Soft_NonHTDescramer_llr(PSDULength,a_total,f_total,r_total,yllr,L)
%% Initialization
yllr = yllr(1:(16+8*PSDULength)); % we cut off the meaningful payload
xllr = zeros(size(yllr));

for i=(1+L):size(yllr,1)
    Pdot_y0_i = exp(yllr(i))/(1+exp(yllr(i)));
    Pdot_y1_i = 1/(1+exp(yllr(i)));
    
    Pdot_z0_i = 0;
    Pdot_z1_i = 0;
    
    a_total_index = mod((i-1),127)+1;
    a_set_i = a_total(a_total_index,:);
    
    for j = 1:size(r_total,1)
        r_set_j = r_total(j,:);
        f_j     = f_total(j);
        tmp = mod(r_set_j*a_set_i.',2);
        if tmp == 0
            Pdot_z0_i = Pdot_z0_i + f_j;
        else
            Pdot_z1_i = Pdot_z1_i + f_j;
        end
    end
    
    Pdot_x0_i = Pdot_y0_i*Pdot_z0_i + Pdot_y1_i*Pdot_z1_i;
    Pdot_x1_i = Pdot_y1_i*Pdot_z0_i + Pdot_y0_i*Pdot_z1_i;
    
    xllr(i) = log(Pdot_x0_i/Pdot_x1_i);
end

psdu_llr = xllr(17:end);

end