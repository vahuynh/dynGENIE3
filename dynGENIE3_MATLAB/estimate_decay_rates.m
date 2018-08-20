function decay_rates = estimate_decay_rates(data,obsTimes)
%function decay_rates = estimate_decay_rates(data,obsTimes)
%
% Estimate the value of the decay rate alpha_i for each gene i.
%
% We suppose that the expression of gene i decreases according to
% A*exp(-alpha_i*t) + C_min 
% between the highest and lowest observed expression levels.
% C_min is set to the minimum expression value over all genes and all
% samples.

nTS = length(data);
ngenes = size(data{1},2);

C_min = min(min(data{1}));
if nTS > 1
    for k=2:nTS
        C_min = min(C_min,min(min(data{k})));
    end
end

decay_rates = zeros(nTS,ngenes);

for k=1:nTS

    t = obsTimes{k};
    
    for i=1:ngenes
    
        x = data{k}(:,i);

        x_min = min(x);
        x_max = max(x);

        idx_min = find(x==x_min,1,'first');
        idx_max = find(x==x_max,1,'first');
        
        x_min = max(x_min-C_min,1e-6);
        x_max = max(x_max-C_min,1e-6);

        decay_rates(k,i) = (log(x_max)-log(x_min)) / abs(t(idx_min)-t(idx_max));
        
    end

   
end

decay_rates = max(decay_rates,[],1);

