%% Driver file for simulating the fjords simulated by Tom using MITgcm

%% Configuring paths
run setup_paths

%% Initialise all needed variables

i_yr = 4; % we want only 2019 data for now, because we only have geometry data for glaciers with data for 2019
time_start = datetime(2019,01,01);
time_end   = datetime(2020,02,15);
input_dt   = 30;
model_dt   = 0.1;
tgt_day = 300; % which day of the 400-day run we want to compare

letters = {'B','F'};

fjord_names = {'Nuussuup','Kakiffaat'} ;


%% Compiling data for all fjords around Greenland (requires data-compilation repository)
run compile_process_fjords

% These are the IDs of the corresponding fjords above in the "fjords_processed" data structure
fjord_ids = [9,23];

%% Reading the data from Cowton et al. (2023)
run load_cowton2023_data

tgt_period = 150:250; %time_glaciers > time_start & time_glaciers < time_end;
qsg_all = qsg_glaciers(tgt_period,:);
taxis_qsg = time_glaciers(tgt_period);
clear qsg_glaciers time_glaciers tgt_period % bit of tidying up

qsg_all = [qsg_all(:,2),qsg_all(:,6)];
meta_table = [meta_table(2,:);meta_table(6,:)];
ids_cowton_fjords = [2,6];

%% Reading the MITgcm outputs

mitgcm_kak   = load([inputs_path,'/MITgcm/KAK_10']);
icebergs_kak = load([inputs_path,'/MITgcm/icebergArea_KAK_10']);
mitgcm_nuu   = load([inputs_path,'/MITgcm/NUU_80']);
icebergs_nuu = load([inputs_path,'/MITgcm/icebergArea_NUU_80']);
x_mitgcm=1:400:size(mitgcm_kak.T,1)*400;
y_mitgcm=1:250:size(mitgcm_kak.T,2)*250;
z_mitgcm=1:10:size(mitgcm_kak.T,3)*10;
t_mitgcm=0:10:400;
j_fjord = 1:144;

% Using zero as _fillValue is not a good practice, but should not be
mitgcm_kak.T(mitgcm_kak.S == 0) = NaN;
mitgcm_nuu.T(mitgcm_nuu.S == 0) = NaN;

temp_kak = squeeze(mean(mitgcm_kak.T(j_fjord,:,:,:),[1,2],'omitnan'));
temp_nuu = squeeze(mean(mitgcm_nuu.T(j_fjord,:,:,:),[1,2],'omitnan'));

% Fjord "B"
mitgcm(1).x = x_mitgcm;
mitgcm(1).y = y_mitgcm;
mitgcm(1).z = z_mitgcm;
mitgcm(1).t = t_mitgcm;
mitgcm(1).Tprofile = temp_nuu(:,tgt_day/10);
mitgcm(1).Tupper = mean(temp_nuu(1:5,:),1,'omitnan');
mitgcm(1).Tinter = mean(temp_nuu(5:25,:),1,'omitnan');
mitgcm(1).Tlower = mean(temp_nuu(25:50,:),1,'omitnan');
mitgcm(1).Iprofile = squeeze(mean(icebergs_nuu.icebergArea,[1,2],'omitnan'));
mitgcm(1).Ivolume = icebergs_nuu.icebergArea.*10;

% Fjord "F"
mitgcm(2).x = x_mitgcm;
mitgcm(2).y = y_mitgcm;
mitgcm(2).z = z_mitgcm;
mitgcm(2).t = t_mitgcm;
mitgcm(2).Tprofile = temp_kak(:,tgt_day/10);
mitgcm(2).Tupper = mean(temp_kak(1:5,:),1,'omitnan');
mitgcm(2).Tinter = mean(temp_kak(5:25,:),1,'omitnan');
mitgcm(2).Tlower = mean(temp_kak(25:50,:),1,'omitnan');
mitgcm(2).Iprofile = squeeze(mean(icebergs_kak.icebergArea,[1,2],'omitnan'));
mitgcm(2).Ivolume = icebergs_kak.icebergArea.*10;
clear i_fjord temp_kak temp_nuu t_mitgcm z_mitgcm y_mitgcm x_mitgcm % bit of tidying up

%% Preparing the model runs
n_fjords = size(meta_table,1);

if exist('fjord_model',"var"),       clear fjord_model; end
fjord_model(n_fjords) = struct("p",[],"a",[],"f",[],"t",[],"m",[],"c",[]);
for i_fjord=1:n_fjords

    i_fjord_ctd = ids_cowton_fjords(i_fjord);
    % time axis for the run
    t = 0:model_dt:400; % 400 days (as in the MITgcm simulations)

    % compile metadata and default parameters
    m.ID   = meta_table.Var1(i_fjord);
    m.name = meta_table.GreenlandicName(i_fjord);
    m.lon  = meta_table.Longitude__E_(i_fjord);
    m.lat  = meta_table.Latitude__N_(i_fjord);
    [p,~]  = get_model_default_parameters();
    
    % obtain CTD data
    ctd.temp  = ctd_data.Tdata{i_yr,i_fjord_ctd};
    ctd.salt  = ctd_data.Sdata{i_yr,i_fjord_ctd};
    ctd.depth = ctd_data.zdata{i_yr,i_fjord_ctd};
    

    % fjord geometry
    p.L         = 1e3.*meta_table.Length_km_(i_fjord);
    p.W         = 1e3.*meta_table.Width_km_(i_fjord);
    p.silldepth = -1.* meta_table.SillDepth_m_(i_fjord);
    p.zgl       = -1.* meta_table.GroundingLineDepth_m_(i_fjord);
    if p.silldepth < p.zgl % if the GL sits above the sill
        p.H = abs(p.silldepth-p.Hmin); % we consider the fjord depth to be just Hmin below the sill, 
                                  % so we can avoid fiddling with the number of layers
    else
        p.H = abs(p.zgl); % if the GL is below the sill, then we treat it as the same as the fjord depth
    end

    % ocean forcing
    ts = ctd.temp.s;
    ss = ctd.salt.s;
    zs = -ctd.depth.s;
    if max(abs(zs)) < p.H % if our profile does not extend to the bottom of the fjord
        h_missing = p.H - max(abs(zs));
        depth_range_missing = max(abs(zs))+1:1:p.H;
        new_zs = [zs -depth_range_missing];
        ts = interp1(zs,ts,new_zs,'nearest','extrap');
        ss = interp1(zs,ss,new_zs,'nearest','extrap');
        zs = new_zs;
        clear new_zs h_missing depth_range_missing % tidying up
    end
    f.zs = zs;
    f.Ts = repmat(ts,length(t),1);
    f.Ss = repmat(ss,length(t),1);


    % initial conditions
    H_clim = double(get_fjord_boxes_from_density(f.Ts',f.Ss',f.zs,p)); % negative values, top to bottom
    [T_clim,S_clim] = bin_ocean_profiles(mean(f.Ts,1),mean(f.Ss,1),f.zs,H_clim,p);
    a.H0 = H_clim;
    a.T0 = T_clim;
    a.S0 = S_clim;

    % flipping the ocean forcing like the model wants it
    f.zs = flip(f.zs);
    f.Ts = flip(f.Ts,2)'; 
    f.Ss = flip(f.Ss,2)';
    
    % obtain Qsg
    % interpolate from daily Qsg data to the model time axes
    taxis_qsg_julian = juliandate(taxis_qsg);
    taxis_qsg_num = taxis_qsg_julian-taxis_qsg_julian(1);
    taxis_qsg_num_dtmodel = taxis_qsg_num(1):t(2)-t(1):taxis_qsg_num(end);
    qsg_dtmodel = interp1(taxis_qsg_num,qsg_all(:,i_fjord),taxis_qsg_num_dtmodel);
    f.Qsg = [zeros([1,200/model_dt]), qsg_dtmodel, zeros([1,100/model_dt])];
    clear taxis_qsg_julian taxis_qsg_num taxis_qsg_num_dtmodel qsg_dtmodel % bit of tidying up

    % obtain iceberg data

    % no icebergs
    % p.M0=0; 
    % f.zi = -800:10:0;
    % f.zi = f.zi';
    % f.D  = zeros(1,length(t));
    % f.xi = (p.nu0/sum(a.H0))*exp(p.nu0*f.zi/sum(a.H0))/(1-exp(-p.nu0));
    % a.I0 = 0*f.zi;

    % obs-based
    % runtime_axis = datasets.opts.time_start:model_dt:datasets.opts.time_end;
    [f.zi,ice_d,f.xi,a.I0] = get_total_iceberg_discharge(datasets,fjords_compilation(fjord_ids(i_fjord)).glaciers,runtime_axis',p,a);
    f.D = ice_d(1:length(t)); % we need to interp/extrapolate the 1-year discharge for 400 days
    
    % MITgcm-based (hybrid: using f.xi and f.zi from here and f.D from obs)
    f.xi = flip(mitgcm(i_fjord).Iprofile./trapz(mitgcm(i_fjord).z,mitgcm(i_fjord).Iprofile)); % get avg. iceberg concentration in each layer and make its integral equal 1
    f.zi = -flip(mitgcm(i_fjord).z);
    % f.D = zeros(size(t)); % only consider the initial iceberg concentration in the fjord, and no discharge afterwards
    % f.D(1) = mean(mitgcm(i_fjord).Ivolume(:),'omitnan'); % should this be "mean" or "sum"??
    
    a.I0 = (f.D(1)/(p.W*p.L)).*f.xi;

    fjord_model(i_fjord).t = t;
    fjord_model(i_fjord).m = m;
    fjord_model(i_fjord).p = p;
    fjord_model(i_fjord).a = a;
    fjord_model(i_fjord).f = f;
    fjord_model(i_fjord).c.tf = ctd.temp.f;
    fjord_model(i_fjord).c.sf = ctd.salt.f;
    fjord_model(i_fjord).c.zf = ctd.depth.f;
    fjord_model(i_fjord).c.ts = ctd.temp.s;
    fjord_model(i_fjord).c.ss = ctd.salt.s;
    fjord_model(i_fjord).c.zs = ctd.depth.s;

    % Sanity-check plot to make sure the boxes, T,S profiles, Qsg and D all make sense
    % figure('Position',[100 100 800 500]);
    % ints=[0,-cumsum(a.H0)];
    % subplot(1,3,1); hold on; box on;
    % plot(f.Ts(:,1),f.zs);
    % scatter(a.T0,(ints(1:end-1)+ints(2:end))/2,'filled');
    % for j=1:size(ints,2)
    %     yline(ints(j));
    % end
    % yline(p.zgl,'--k','linewidth',0.5)
    % yline(p.silldepth,'--k','linewidth',0.5)
    % ylim([-p.H 0])
    % subplot(1,3,2); hold on; box on;
    % plot(f.Ss(:,1),f.zs);
    % scatter(a.S0,(ints(1:end-1)+ints(2:end))/2,'filled');
    % for j=1:size(ints,2)
    %     yline(ints(j));
    % end
    % yline(p.zgl,'--k','linewidth',0.5)
    % yline(p.silldepth,'--k','linewidth',0.5)
    % ylim([-p.H 0])
    % subplot(1,3,3); hold on; box on; grid on;
    % plot(t,f.Qsg);
    % plot(t,f.D);
    % legend('Subglacial discharge','Solid-ice discharge');

    clear t m p a f % bit of tidying up
end

%% Run the model


n_fjord_runs = length(fjord_model);
range_C0 = [1e3,5e3,1e4,5e4,1e5];
range_P0 = 5:5:30;
range_M0 = [5e-9,1e-8,2e-8,5e-8];
range_gamma = 0.4:0.05:0.6;
n_C0 = length(range_C0);
n_P0 = length(range_P0);
n_M0 = length(range_M0);
n_gamma = length(range_gamma);
n_combinations=n_gamma*n_M0*n_P0*n_C0;
if exist('ensemble',"var"),clear ensemble; end
% ensemble(n_fjord_runs, length(range_C0), length(range_P0)) = struct("p",[],"a",[],"f",[],"t",[],"m",[],"c",[],"s",[]);
ensemble(n_fjord_runs, length(range_C0), length(range_P0),length(range_M0),length(range_gamma)) = struct("p",[],"t",[],"m",[],"s",[]);
run_counter=0;
for i_fjord=1:n_fjord_runs
tic
for i_C0=1:n_C0
for i_P0=1:n_P0
for i_M0=1:n_M0
for i_gamma=1:n_gamma
    run_counter = run_counter+1;
    cur_fjord = fjord_model(i_fjord);
    cur_fjord.p.C0 = range_C0(i_C0);
    cur_fjord.p.P0 = range_P0(i_P0);
    cur_fjord.p.M0 = range_M0(i_M0);
    cur_fjord.p.gamma = range_gamma(i_gamma);
    % cur_fjord.p.plot_runtime=1;
    try
        [cur_fjord.s,cur_fjord.f] = boxmodel(cur_fjord.p, cur_fjord.t, cur_fjord.f, cur_fjord.a);
        fprintf('run %d complete. ',run_counter)

        ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).p = cur_fjord.p;
        ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).t = cur_fjord.t;
        ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).m = cur_fjord.m;
        % ensemble(i_fjord,i_C0,i_P0).s.Tfinal = mean(cur_fjord.s.T(:,end-365:end),2); % last 30 days(?)
        % ensemble(i_fjord,i_C0,i_P0).s.Sfinal = mean(cur_fjord.s.S(:,end-365:end),2);
        % ensemble(i_fjord,i_C0,i_P0).s.Hfinal = mean(cur_fjord.s.H(:,end-365:end),2);
        ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).s.Tfinal = mean(cur_fjord.s.T(:,(tgt_day-5:tgt_day+5)/model_dt),2); % 10-day avg centered at the target day (since MITgcm outputs a 10-day average)
        ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).s.Sfinal = mean(cur_fjord.s.S(:,(tgt_day-5:tgt_day+5)/model_dt),2);
        ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).s.Hfinal = mean(cur_fjord.s.H(:,(tgt_day-5:tgt_day+5)/model_dt),2);

        Tregular=NaN([length(mitgcm(i_fjord).z),length(cur_fjord.s.t)]);
        for k=1:length(cur_fjord.s.t)
            ints=[0; cumsum(cur_fjord.s.H(:,k))];
            z_box = (ints(1:end-1)+ints(2:end))/2;
            Tregular(:,k) = interp1(z_box,cur_fjord.s.T(:,k),mitgcm(i_fjord).z,'nearest','extrap');
        end
        
        ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).s.Tupper = mean(Tregular(1:5,:),1,'omitnan');   % avg temp of the upper 50 m
        ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).s.Tinter = mean(Tregular(5:25,:),1,'omitnan');  % avg temp of the 50 - 250 m layer
        ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).s.Tlower = mean(Tregular(25:50,:),1,'omitnan'); % avg temp of the 250 - 500 m layer
        fprintf('Output interpolation complete. ')
    catch ME
        fprintf('run %d failed: %s. ',run_counter,ME.message)
    end
    fprintf('\n')
end % for i_gamma
end % for i_M0
end % for i_P0
end % for i_C0
fprintf('Done with fjord %d. ',i_fjord)
toc
fprintf('\n')
end % for i_fjord
save([outs_path,'boxmodel_runs_MITgcm_hybice_comp_n',num2str(n_combinations),'_day',num2str(tgt_day)],'-v7.3','ensemble','fjord_model');
disp('Outputs saved.')
% load([outs_path,'boxmodel_runs_MITgcm_comp_n',num2str(n_combinations),'_day',num2str(tgt_day)]);
%% post-processing to make things 1:1 comparable


if exist('res_obs',"var"),       clear res_obs; end
res_obs(size(ensemble)) = struct("tf",[],"sf",[],"zf",[]);
if exist('res_box',"var"),       clear res_box; end
res_box(size(ensemble)) = struct("tf",[],"sf",[],"zf",[],"ID",[],"name",[]);
for i=1:size(ensemble,1)
    zf_obs = fjord_model(i).c.zf;
    tf_obs = fjord_model(i).c.tf;
    sf_obs = fjord_model(i).c.sf;
    n_completed = 0; % completed runs for that fjord

    % n_layers = ensemble(i,1,1).p.N+ensemble(i,1,1).p.sill;
    tf_box_comp = NaN([length(zf_obs),n_fjord_runs,n_C0,n_P0,n_M0,n_gamma]);
    sf_box_comp = NaN(size(tf_box_comp));
    tupper_box_comp = NaN([length(fjord_model(i).t)-1,n_fjord_runs,n_C0,n_P0,n_M0,n_gamma]);
    tinter_box_comp = NaN(size(tupper_box_comp));
    tlower_box_comp = NaN(size(tupper_box_comp));
    for j=1:size(ensemble,2)
    for k=1:size(ensemble,3)
    for l=1:size(ensemble,4)
    for n=1:size(ensemble,5) % we skip 'm' because it is being used for the metadata structure elsewhere in the code
        if ~isempty(ensemble(i,j,k,l,n).s) & ~isnan(ensemble(i,j,k,l,n).s.Hfinal)
            ints=[0; cumsum(ensemble(i,j,k,l,n).s.Hfinal)];
            z_box = (ints(1:end-1)+ints(2:end))/2;
            tf_box = ensemble(i,j,k,l,n).s.Tfinal; 
            sf_box = ensemble(i,j,k,l,n).s.Sfinal; 
            
            tf_box_comp(:,i,j,k,l,n) = interp1(z_box,tf_box,zf_obs,'nearest','extrap');
            sf_box_comp(:,i,j,k,l,n) = interp1(z_box,sf_box,zf_obs,'nearest','extrap');

            tupper_box_comp(:,i,j,k,l,n) = ensemble(i,j,k,l,n).s.Tupper;
            tinter_box_comp(:,i,j,k,l,n) = ensemble(i,j,k,l,n).s.Tinter;
            tlower_box_comp(:,i,j,k,l,n) = ensemble(i,j,k,l,n).s.Tlower;

            n_completed = n_completed+1;
        end
    end % n
    end % l
    end % k
    end % j
    res_obs(i).zf = zf_obs;
    res_obs(i).tf = tf_obs;
    res_obs(i).sf = sf_obs;
    
    res_obs(i).zs = fjord_model(i).c.zs;
    res_obs(i).ts = fjord_model(i).c.ts;
    res_obs(i).ss = fjord_model(i).c.ss;

    res_box(i).t  = ensemble(i,1,1,1,1).t(2:end);
    res_box(i).zf = z_box;
    res_box(i).tf = median(tf_box_comp(:,i,:,:,:),[3,4,5,6],'omitnan');
    res_box(i).tfmin = min(tf_box_comp(:,i,:,:,:),[],[3,4,5,6],'omitnan');
    res_box(i).tfmax = max(tf_box_comp(:,i,:,:,:),[],[3,4,5,6],'omitnan');
    res_box(i).sf = median(sf_box_comp(:,i,:,:,:),[3,4,5,6],'omitnan');


    res_box(i).Tupper = median(tupper_box_comp(:,i,:,:,:,:),[3,4,5,6],'omitnan');
    res_box(i).Tupper_min = min(tupper_box_comp(:,i,:,:,:,:),[],[3,4,5,6],'omitnan');
    res_box(i).Tupper_max = max(tupper_box_comp(:,i,:,:,:,:),[],[3,4,5,6],'omitnan');

    res_box(i).Tinter = median(tinter_box_comp(:,i,:,:,:,:),[3,4,5,6],'omitnan');
    res_box(i).Tinter_min = min(tinter_box_comp(:,i,:,:,:,:),[],[3,4,5,6],'omitnan');
    res_box(i).Tinter_max = max(tinter_box_comp(:,i,:,:,:,:),[],[3,4,5,6],'omitnan');

    res_box(i).Tlower = median(tlower_box_comp(:,i,:,:,:,:),[3,4,5,6],'omitnan');
    res_box(i).Tlower_min = min(tlower_box_comp(:,i,:,:,:,:),[],[3,4,5,6],'omitnan');
    res_box(i).Tlower_max = max(tlower_box_comp(:,i,:,:,:,:),[],[3,4,5,6],'omitnan');

    res_box(i).id = fjord_model(i).m.ID{1};
    res_box(i).name = fjord_model(i).m.name{1};
    res_box(i).n = n_completed;
end

%% Plotting results
lcolor = lines(3);
figure('Name','Temperature comparison','Position',[40 40 1200 800]);
tiledlayout(2,2);
for i=1:n_fjord_runs
    
    nexttile; hold on; box on; grid on
    text(0.02,1.05,sprintf("(%s) %s",res_box(i).id,res_box(i).name),'units','normalized','fontsize',14)
    text(0.02,0.05,sprintf("n=%d",res_box(i).n),'Units','normalized','FontSize',14)

    % figure; hold on
    % y2 = [res_obs(i).zf; flip(res_obs(i).zf)]';
    % inBetween = [res_box(i).tfmin, flip(res_box(i).tfmax)];
    % hp = fill(flip(inBetween,1), flip(-y2,2), lcolor(3,:),'edgecolor','none','facealpha',0.2);

    hs = plot(res_obs(i).ts,-res_obs(i).zs,'linewidth',1.5,'color',lcolor(1,:));
    hf = plot(res_obs(i).tf,-res_obs(i).zf,'linewidth',1.5,'color',lcolor(2,:));
    hb = plot(res_box(i).tf,-res_obs(i).zf,'linewidth',1.5,'color',lcolor(3,:));
    hm = plot(mitgcm(i).Tprofile,-mitgcm(i).z,'linewidth',1.5,'color','k');
    plot(res_box(i).tfmin,-res_obs(i).zf,'linewidth',1.5,'color',lcolor(3,:),'LineStyle','--');
    plot(res_box(i).tfmax,-res_obs(i).zf,'linewidth',1.5,'color',lcolor(3,:),'LineStyle','--');
    xlabel('Temperature (^oC)'); ylabel('Depth (m)');
    set(gca,'fontsize',14)
    xlim([-2 6])

    nexttile; hold on; box on; grid on
    if length(res_box(i).t) == length(res_box(i).Tupper)
        hbl = plot(res_box(i).t,res_box(i).Tupper,'linewidth',1.5,'color','k'); % the first two lines are dummy for the legend entries
        hml = plot(res_box(i).t,res_box(i).Tupper,'linewidth',1.5,'color','k','LineStyle','--');
    
        hu = plot(res_box(i).t,res_box(i).Tupper,'linewidth',1.5,'color',lcolor(1,:));
        hi = plot(res_box(i).t,res_box(i).Tinter,'linewidth',1.5,'color',lcolor(2,:));
        hl = plot(res_box(i).t,res_box(i).Tlower,'linewidth',1.5,'color',lcolor(3,:));
    
        plot(mitgcm(i).t,mitgcm(i).Tupper,'linewidth',1.5,'color',lcolor(1,:),'LineStyle','--');
        plot(mitgcm(i).t,mitgcm(i).Tinter,'linewidth',1.5,'color',lcolor(2,:),'LineStyle','--');
        plot(mitgcm(i).t,mitgcm(i).Tlower,'linewidth',1.5,'color',lcolor(3,:),'LineStyle','--');
        % plot(res_box(i).t,res_box(i).Tupper_min,'linewidth',1.5,'color',lcolor(1,:),'LineStyle',':');
        % plot(res_box(i).t,res_box(i).Tupper_max,'linewidth',1.5,'color',lcolor(1,:),'LineStyle',':');
        % plot(res_box(i).t,res_box(i).Tinter_min,'linewidth',1.5,'color',lcolor(2,:),'LineStyle',':');
        % plot(res_box(i).t,res_box(i).Tinter_max,'linewidth',1.5,'color',lcolor(2,:),'LineStyle',':');
        % plot(res_box(i).t,res_box(i).Tlower_min,'linewidth',1.5,'color',lcolor(3,:),'LineStyle',':');
        % plot(res_box(i).t,res_box(i).Tlower_max,'linewidth',1.5,'color',lcolor(3,:),'LineStyle',':');
       
    end
    ylabel('Box model Temperature (^oC)'); xlabel('Model time (days)');
    set(gca,'fontsize',14)
    if i==1 
        hl = legend([hs, hf, hb, hm], {"Shelf","Fjord","Box model","MITgcm"},'fontsize',14,'Location','Southeast'); 
        title(hl,sprintf('Profiles at day %d\n(10-day avg.)',tgt_day))
    end
    if i==2 
        hl = legend([hbl,hml,hu,hi,hl],{"Box model","MITgcm","0-50 m","50-250 m","250-500 m"},'fontsize',14,'Location','Northeast'); 
        title(hl,'Time series')
    end
end


% exportgraphics(gcf,[figs_path,'temp_profiles_MITgcm_hybice_comp_n',num2str(n_combinations),'_day',num2str(tgt_day),'.png'],'Resolution',300)