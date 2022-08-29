function [a_total,S_total] = DetMatRep(L,K)
a_total = zeros(L,K);

%% initialization 
S_total = cell(L+7,1);
S_total(1) = {[1]};
S_total(2) = {[2]};
S_total(3) = {[3]};
S_total(4) = {[4]};
S_total(5) = {[5]};
S_total(6) = {[6]};
S_total(7) = {[7]};

for i = 8:L+7
    % the combined set
    tmp = union(cell2mat(S_total(i-7)),cell2mat(S_total(i-4)));
    tmp2 = intersect(cell2mat(S_total(i-7)),cell2mat(S_total(i-4)));
    [~,nums] = size(tmp2);
    if (nums~=0)
       for j= 1:nums
            tmp(tmp==tmp2(j)) = [];
       end
    end
    
    for k = 1:K
        index = tmp(tmp==k);
        if index>=1 
            a_total(i-7,k) = 1;
        else
            a_total(i-7,k) = 0;
        end
    end
    S_total(i) = mat2cell(tmp,1,length(tmp));
end
S_total = S_total(8:end);

end