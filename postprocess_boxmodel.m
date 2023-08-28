function outputs = postprocess_boxmodel(fjord)

    run_length = min(length(fjord.s.t),size(fjord.s.H,2));

    % output variables of interest
    z_profiles=fjord.f.zs;
    temp_profiles=NaN(size(fjord.f.Ts));
    salt_profiles=NaN(size(fjord.f.Ss));
    t_export=NaN(size(fjord.t));
    s_export=NaN(size(fjord.t));
    v_export=NaN(size(fjord.t));
    z_bnds_export=NaN([2,length(fjord.t)]);

    % iterates over all time steps
    for i=1:run_length

        % copies profile values for all depths corresponding to that box
        box_depths=-[0, cumsum(fjord.s.H(:,i)')];
        for k=1:length(box_depths)-1
            inds_box=z_profiles < box_depths(k) & z_profiles > box_depths(k+1);
            temp_profiles(inds_box,i)=fjord.s.T(k,i);
            salt_profiles(inds_box,i)=fjord.s.S(k,i);
        end

        % computes outflow of freshwater
        i_outflow = find(fjord.s.QVs(:,i) > 0);       % finds where flow is going towards the shelf
        [~,i_export] = max(fjord.s.QVs(i_outflow,i)); % finds the strongest flow

        % if something is being exported
        if ~isempty(i_export) 
            z_bnds_export(:,i)=box_depths(i_export:i_export+1); % gets the depth interval
            t_export(i)=fjord.s.Te(i_export,i);                 % copies over the total heat
            s_export(i)=fjord.s.Se(i_export,i);                 % and salt exported,
            v_export(i)=fjord.s.QVs(i_export,i);                % as well as volume fluxes
        end

    end

    % computes total fjord salt (kg) and heat content (J)    
    outputs.sc = fjord.s.S  .*               (fjord.p.betaS*fjord.s.S - fjord.p.betaT*fjord.s.T)  .* fjord.p.L.*fjord.p.W.*fjord.p.H;
    outputs.hc = (fjord.s.T .* fjord.p.cw .* (fjord.p.betaS*fjord.s.S - fjord.p.betaT*fjord.s.T)) .* fjord.p.L.*fjord.p.W.*fjord.p.H;

    % copies all into output structure
    outputs.Tf = temp_profiles;
    outputs.Sf = salt_profiles;
    outputs.Ze = smoothdata(z_bnds_export,2);
    outputs.Te = smoothdata(t_export);
    outputs.Se = smoothdata(s_export);
    outputs.Ve = smoothdata(v_export);

end