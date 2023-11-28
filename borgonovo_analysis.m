%% Borgonovo indices
BorgonovoA_ohc = cell([1,n_regions]);
BorgonovoA_osc = cell([1,n_regions]);
BorgonovoOpts.Type = 'Sensitivity';
BorgonovoOpts.Method = 'Borgonovo';

for i_reg=1:n_regions
    BorgonovoOpts.Borgonovo.Sample.X = Xeval(:,i_reg,:);
    BorgonovoOpts.Borgonovo.Sample.Y = Yeval_ohc{i_reg}; % likely needs some .something
    BorgonovoA_ohc{i_reg} = uq_createAnalysis(BorgonovoOpts);
    BorgonovoOpts.Borgonovo.Sample.Y = Yeval_osc{i_reg}; % likely needs some .something
    BorgonovoA_osc{i_reg} = uq_createAnalysis(BorgonovoOpts);
end


