function [params,iOpts,probs,fjords_processed] = define_model_param_distrib(datasets,fjords_compilation,glaciers_compilation,i_reg,experiments_taxis,experiments_time_dt,n_years)
%% A more statistically driven approach to sensitivity tests
% derives distributions for all parameters on which our model simulations
% will depend on: L, H, W, Zg, Zs, Ta, Sa, Qa, Da, P0, and C0

%% Get compiled fjord data in pre-processed structure
% also selects the desired period in time
verbose.plot=1;  % change to 1 to produce plots of all parameter distributions
verbose.print=0; % change 1 to print basic fjord statistics

if nargin < 5 || isempty(experiments_taxis)
    datasets.opts.time_start = datetime(2010,01,15);
    datasets.opts.time_end   = datetime(2018,12,15);
else
    datasets.opts.time_start = experiments_taxis(1);
    datasets.opts.time_end   = experiments_taxis(end);
end
datasets.opts.time_interval = [datasets.opts.time_start,datasets.opts.time_end]; 
datasets.opts.dt            = 30.0; % time step (in days) for creating the forcings
if nargin < 6
    experiments_time_dt       = 0.10; % time step (in days) for the actual experiments
end
if nargin < 7
    n_years=1;
end
fjords_processed(size(fjords_compilation)) = struct("p",[],"a",[],"f",[],"t",[],"m",[]);
for i=1:length(fjords_compilation)
    fjords_processed(i) = prepare_boxmodel_input(datasets,fjords_compilation(i));
end

%% Gets distribution of fjord geometry parameters
% since the geometry does not change with time, we can stick to the
% original fjords structure
fjord_stats = print_fjord_statistics(fjords_compilation,verbose);

probs(1) = fjord_stats.L.pd;
iOpts.Marginals(1) = uq_KernelMarginals(fjord_stats.L.total', [min(fjord_stats.L.total) max(fjord_stats.L.total)]);

probs(end+1) = fjord_stats.H.pd;
iOpts.Marginals(end+1) = uq_KernelMarginals(fjord_stats.H.total',[min(fjord_stats.H.total), max(fjord_stats.H.total)]);

% probs(end+1) = fjord_stats.W.pd;
% iOpts.Marginals(end+1) = uq_KernelMarginals(fjord_stats.W.total',[min(fjord_stats.W.total), max(fjord_stats.W.total)]);
probs(end+1) = fjord_stats.a.pd;
iOpts.Marginals(end+1) = uq_KernelMarginals(fjord_stats.a.total',[1, max(fjord_stats.a.total)]);

% probs(end+1) = fjord_stats.Zs.pd;
% iOpts.Marginals(end+1) = uq_KernelMarginals(fjord_stats.Zs.total', [min(fjord_stats.Zs.total) -50]);
probs(end+1) = fjord_stats.Zsr.pd;
iOpts.Marginals(end+1) = uq_KernelMarginals(fjord_stats.Zsr.total', [min(abs(fjord_stats.Zs.total)/max(fjord_stats.H.total)) max(abs(fjord_stats.Zs.total)/max(fjord_stats.H.total))]);
% iOpts.Marginals(end+1) = uq_KernelMarginals(fjord_stats.Zsr.total', [0.05 0.95]);

% probs(end+1) = fjord_stats.Zg.pd;
% iOpts.Marginals(end+1) = uq_KernelMarginals(fjord_stats.Zg.total', [min(fjord_stats.Zg.total) -50]);
probs(end+1) = fjord_stats.Zgr.pd;
iOpts.Marginals(end+1) = uq_KernelMarginals(fjord_stats.Zgr.total', [min(abs(fjord_stats.Zg.total)/max(fjord_stats.H.total)) max(abs(fjord_stats.Zg.total)/max(fjord_stats.H.total))]);
% iOpts.Marginals(end+1) = uq_KernelMarginals(fjord_stats.Zgr.total', [0.05 0.95]);


%% Gets the glacier forcing

runtime_axis = [datetime(datasets.opts.time_start):datasets.opts.dt:datetime(datasets.opts.time_end)]';
% Subglacial discharge - we want the variation in amplitude rather than
% anomalies; otherwise the signal will be dominated by small variations in non-summer months
str_reg = {'SW','SE','CW','CE','NW','NE','NO'};
q_amp_max=NaN(size(glaciers_compilation));
q_amp_min=NaN(size(glaciers_compilation));


if datasets.opts.restrict_to_fjords
    if ~isempty(i_reg)
        [q_forcing,~,~,~] = get_var_forcing_by_region(fjords_processed,'Qsg');
        q_forcing = q_forcing(:,i_reg);
    else
        [q_clim,~,~,~] = get_var_clim_by_region(fjords_processed,'Qsg');
        q_clim = mean(q_clim,2); % average over all regions if none is specified
        q_forcing = repmat(q_clim,n_years); 
    end
    
    for k=1:size(q_forcing,2)
        for peak_prominence=10:-0.5:0.5
            [v_peaks,~] = findpeaks(q_forcing(:,k),'MinPeakProminence',peak_prominence);
            if length(v_peaks) > 3 % at least 3 distinct peaks need to be picked in the series
                break
            end
        end
        try
        q_amp_max(k) = max(v_peaks)./mean(v_peaks);
        q_amp_min(k) = min(v_peaks)./mean(v_peaks);
        catch
            disp('glacial discharge is really low (<0.5 m3/s). Skipping glacier...')
        end
    end
else
    q_reg=NaN([length(glaciers_compilation),length(runtime_axis)]);
    for k=1:length(glaciers_compilation)
        if isempty(i_reg) || strcmp(glaciers_compilation(k).region,str_reg{i_reg}) % will also be true if we do not provide i_reg, hence all data is accounted for
            q_reg(k,:) = get_total_glacier_runoff(glaciers_compilation(k),runtime_axis);
            for peak_prominence=10:-0.5:0.5
                [v_peaks,~] = findpeaks(q_reg(k,:),'MinPeakProminence',peak_prominence);
                if length(v_peaks) > 3 % at least 3 distinct peaks need to be picked in the series
                    break
                end
            end
            try
            q_amp_max(k) = max(v_peaks)./mean(v_peaks);
            q_amp_min(k) = min(v_peaks)./mean(v_peaks);
            catch
                disp('glacial discharge is really low (<0.5 m3/s). Skipping glacier...')
            end
        end
    end
    q_forcing = mean(q_reg,1,'omitnan');
end

% It makes no sense to have Qsg in terms of anomalies, because in the
% frequency distribution, the very small variations in non-summer months would
% dominate the signal. So we prescribe a uniform distribution of
% "amplifying factors" instead
q_amp_max = mean(q_amp_max,'omitnan');
q_amp_min = mean(q_amp_min,'omitnan');
q_pd      = makedist('Uniform','lower',q_amp_min,'upper',q_amp_max);

iOpts.Marginals(end+1).Type     = 'uniform';
iOpts.Marginals(end).Parameters = [q_amp_min q_amp_max];
iOpts.Marginals(end).Bounds     = [0 1.5*q_amp_max];
probs(end+1)                    = q_pd;

%Plot to check parameter space
% figure('Name','Qsg parameter space'); histogram(random(q_pd,1000));

% Solid-ice discharge (same procedure as for T and S)
if datasets.opts.restrict_to_fjords
    if ~isempty(i_reg)
        [d_forcing,d_anom,~,~] = get_var_forcing_by_region(fjords_processed,'D');
        d_forcing=d_forcing(:,i_reg);
        d_anom=d_anom{i_reg};
    else
        [q_clim,~,~,~] = get_var_clim_by_region(fjords_processed,'Qsg');
        q_clim = mean(q_clim,2); % average over all regions if none is specified    
        q_forcing = repmat(q_clim,n_years); 
    end
else
    d_reg=NaN([length(glaciers_compilation),length(runtime_axis)]);
    for k=1:length(glaciers_compilation)
        if isempty(i_reg) || strcmp(glaciers_compilation(k).region,str_reg{i_reg}) % will also be true if we do not provide i_reg, hence all data is accounted for
            time_discharge = glaciers_compilation(k).iceberg.time_axis;
            % we need to convert from total discharge in the month to discharge per second
            seconds_in_year = eomday(time_discharge.Year,time_discharge.Month)*86400*12;
            Gt_to_m3_ice = 1e9*1e3/916.7; % Giga * kg / (kg/m^3) assuming glacier ice density of 916.7 kg/m^3
            D_s = glaciers_compilation(k).iceberg.discharge./seconds_in_year .* Gt_to_m3_ice;
            [~] = data_overlap_check(runtime_axis,time_discharge,'iceberg'); % simple check that the data is available for the simulation period
            d_reg(k,:) = interp1(time_discharge,D_s,runtime_axis,'linear','extrap'); 
        end
    end
    d_forcing = squeeze(median(d_reg,1,'omitnan'));
    d_anom = d_reg - d_forcing;
    d_anom = d_anom(~isnan(d_anom));
end
d_pd                   = fitdist(d_anom,'kernel');
probs(end+1) = d_pd;
iOpts.Marginals(end+1) = uq_KernelMarginals(d_anom,[min(d_anom), max(d_anom)]);

% considering we use an entrainment coefficient of 0.1, P0=[5,30] is
% equivalent to a plume width of 50-300 m
p_pd              = makedist('Uniform','lower',5,'upper',30); 
iOpts.Marginals(end+1).Type     = 'uniform';
iOpts.Marginals(end).Parameters = [5 30];
iOpts.Marginals(end).Bounds     = [5 30];
probs(end+1)                    = p_pd;

% Considering we use a shelf-exchange coefficient between 1e3 and 1e5
p_c0              = makedist('Uniform','lower',1e3,'upper',1e5); 
iOpts.Marginals(end+1).Type     = 'uniform';
iOpts.Marginals(end).Parameters = [1e3 1e5];
iOpts.Marginals(end).Bounds     = [1e3 1e5];
probs(end+1)                    = p_c0;

%% Gets the ocean forcing

% If we are using a gridded ocean dataset, it will take the fields from the fjords structure
if datasets.opts.load_ocean
    % [tocn_forcing, tocn_anom, tocn_decay, depths] = get_var_clim_by_region(fjords_processed,'Ts');
    % [socn_forcing, socn_anom, socn_decay, ~]      = get_var_clim_by_region(fjords_processed,'Ss');
    [tocn_forcing, tocn_anom, tocn_decay, depths] = get_var_forcing_by_region(fjords_processed,'Ts');
    [socn_forcing, socn_anom, socn_decay, ~]      = get_var_forcing_by_region(fjords_processed,'Ss');
    
    % reduce anomaly to a standardised/normalised profile and a "dT factor"
    % Anomaly based on surface, multiplied by an "atenuation profile/decay function"
    tocn_pd = fitdist(tocn_anom{i_reg},'kernel');    
    socn_pd = fitdist(socn_anom{i_reg},'kernel');
    probs(end+1) = tocn_pd;
    iOpts.Marginals(end+1) = uq_KernelMarginals(tocn_anom{i_reg},[min(tocn_anom{i_reg}), max(tocn_anom{i_reg})]);
    
    probs(end+1) = socn_pd;
    iOpts.Marginals(end+1) = uq_KernelMarginals(socn_anom{i_reg},[min(socn_anom{i_reg}), max(socn_anom{i_reg})]);

elseif isfield(datasets.obs,'ctd') % We make no distinction between casts on shelf
    % restrict our probability distribution to select only the casts within the desired period
    date_min = experiments_taxis(1);
    date_max = experiments_taxis(end);
    iobs_min=1;
    iobs_max=length(datasets.obs.ctd_shelf);
    for i=2:length(datasets.obs.ctd_shelf) % we can iterate over the array like this because we sorted it by time before
        if datasets.obs.ctd_shelf(i).time < date_min
            iobs_min=i;
        end
        if datasets.obs.ctd_shelf(i).time < date_max
            iobs_max=i;
        end
    end
    q_cast      = makedist('Uniform','lower',iobs_min,'upper',iobs_max);
    iOpts.Marginals(end+1).Type     = 'uniform';
    iOpts.Marginals(end).Parameters = [iobs_min iobs_max];
    iOpts.Marginals(end).Bounds     = [iobs_min iobs_max];
    probs(end+1)                    = q_cast;

else
    disp('A gridded ocean dataset or CTD casts must be loaded into the datasets structure!')
end

% omeg_pd = makedist('Normal','mu',0.1429,'sigma',0.1167);
% only applies omega to eastern Greenland
% if ismember(i_reg,[2,4,6]) 
%     iOpts.Marginals(end+1).Type     = 'Gaussian';
%     iOpts.Marginals(end).Parameters = [0.1429 0.1167]; % [1/7 days, 3.5/30 days] as per Harden et al. (2011; JClimate)
%     iOpts.Marginals(end).Bounds     = [0 1];
% else
%     omeg_pd = makedist('Normal','mu',0.,'sigma',0.);
%     iOpts.Marginals(end+1).Type = 'Constant' ;
%     iOpts.Marginals(end).Parameters = 0;
% end
% 
% probs(end+1) = omeg_pd;

%% interpolates time series variables to the actual time steps used by the model
datasets.opts.dt            = experiments_time_dt; % time step in days
fjord_dummy = prepare_boxmodel_input(datasets,fjords_compilation(1));
time_axis = fjord_dummy.t;

q_forcing    = interp1(fjords_processed(1).t,q_forcing,time_axis,'linear','extrap');
d_forcing    = interp1(fjords_processed(1).t,d_forcing,time_axis,'linear','extrap');

if datasets.opts.load_ocean
    tocn_forcing = interp1(fjords_processed(1).t,tocn_forcing,time_axis,'linear','extrap');
    tocn_decay = interp1(fjords_processed(1).t,tocn_decay,time_axis,'linear','extrap');
    socn_forcing = interp1(fjords_processed(1).t,socn_forcing,time_axis,'linear','extrap');
    socn_decay = interp1(fjords_processed(1).t,socn_decay,time_axis,'linear','extrap');
end
%% Compiles the parameters into the structure to be used by the wrapper function

params.regID = i_reg;
params.H     = max(fjord_stats.H.total);
params.Qglc  = q_forcing;%(:,i_reg);
params.Dglc  = d_forcing;%(:,i_reg);
params.t     = time_axis;
params.zi    = fjord_dummy.f.zi;
params.xi    = fjord_dummy.f.xi;
params.I0    = fjord_dummy.a.I0;

if datasets.opts.load_ocean
    params.zs    = -depths;
    params.Tocn  = tocn_forcing(:,:,i_reg);
    params.Tdec  = tocn_decay(:,:,i_reg);
    params.Socn  = socn_forcing(:,:,i_reg);
    params.Sdec  = socn_decay(:,:,i_reg);
elseif isfield(datasets.obs,'ctd')
    params.ocn = datasets.obs.ctd_shelf;
end

% Add the names to the input distributions we created
iOpts.Marginals(1).Name = 'L';
% iOpts.Marginals(2).Name = 'W';
iOpts.Marginals(2).Name = 'H';
iOpts.Marginals(3).Name = 'L/W'; % iOpts.Marginals(2).Name = 'W';
iOpts.Marginals(4).Name = 'Zs/H';
iOpts.Marginals(5).Name = 'Zg/H';
% iOpts.Marginals(3).Name = 'Zs';
% iOpts.Marginals(4).Name = 'Zg';
iOpts.Marginals(6).Name = 'Qa';
iOpts.Marginals(7).Name = 'Da';
iOpts.Marginals(8).Name = 'P0';
iOpts.Marginals(9).Name = 'C0';
if datasets.opts.load_ocean
    iOpts.Marginals(10).Name = 'Ta';
    iOpts.Marginals(11).Name = 'Sa';
elseif isfield(datasets.obs,'ctd')
    iOpts.Marginals(10).Name = 'Cast No.';
else
    disp('A gridded ocean dataset or CTD casts must be loaded into the datasets structure!')
end

end


