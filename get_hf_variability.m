%% Configuring paths
run load_local_paths.m % sets data_path, import_path, collation_path, model_path, and project_path
addpath(genpath(import_path))
addpath(genpath(model_path))
addpath(genpath(collation_path))
addpath(genpath(uqlab_path))

[datasets,fjords_compilation,~,~] = compile_datasets(data_path);

%% Compile obs data for each region

sw_data = process_gem_ctd(datasets.obs.folder_gem,'Nuuk'); % Nuuk represents SW data
se_data = process_adc_ctd(datasets.obs.folder_adc,'SF');        % Sermilik represents SE

cw_data = process_gem_ctd(datasets.obs.folder_gem,'Disko'); % Disko represents CW data
% TODO: CE (i.e.., Kangerlussuaq) is contained in the SF (Sermilik) files

nw_data = process_adc_ctd(datasets.obs.folder_adc,'Upernavik'); % Upernarvik represents NW
ne_data = process_gem_ctd(datasets.obs.folder_gem,'Zackenberg'); % Zackenberg represents NE data

% plt_obs_handles = plot_obs(datasets,ne_data); % Plot to check that they are in the right place

%% Rearrange data for a monotonic time axis and select target depths
std_depths = [1,10:10:30,50:25:150,200,250,300:100:1500,1750,2000:500:9000]'; % using 1 instead of 0 to avoid NaN at surface
desired_depths = std_depths(1:7);
data_sel = cell([7,1]);

data_sel{1} = get_casts_table(sw_data,std_depths,'linear',desired_depths,1);
data_sel{2} = get_casts_table(se_data,std_depths,'linear',desired_depths,1);

data_sel{3} = get_casts_table(cw_data,std_depths,'linear',desired_depths,1);

data_sel{5} = get_casts_table(nw_data,std_depths,'linear',desired_depths,1);
data_sel{6} = get_casts_table(ne_data,std_depths,'linear',desired_depths,1);

% some memory management, because things get really nasty here...
clear sw_data se_data cw_data nw_data ne_data

data_daily = cell([7,1]);
for i=1:7
    if ~isempty(data_sel{i})
        data_daily{i} = retime(data_sel{i},'daily','linear');
    end
end

for i=1:7
    if ~isempty(data_sel{i})
        figure; 
        subplot(2,1,1); box on; hold on
        scatter(data_sel{i}.time,data_sel{i}.temp,'.k'); ylabel('Temperature'); 
        plot(data_daily{i}.time,data_daily{i}.temp,'linewidth',0.5);

        subplot(2,1,2); box on; hold on
        scatter(data_sel{i}.time,data_sel{i}.salt,'.k'); ylabel('Salinity');
        plot(data_daily{i}.time,data_daily{i}.salt,'linewidth',0.5);
    end
end

% SW (Nuuk) and NE (Zackenberg) seem to be the only ones with good enough coverage?
[p_sw,f_sw] = pspectrum(data_daily{1});
figure; 
subplot(1,2,1), plot(f_sw/pi,p_sw(:,4)); xlim([0 2]*1e-8)
subplot(1,2,2), plot(f_sw/pi,p_sw(:,5)); xlim([0 2]*1e-8)

[p_ne,f_ne] = pspectrum(data_daily{6});
figure; 
subplot(1,2,1), plot(f_ne/pi,p_ne(:,4)); xlim([0 2]*1e-8)
subplot(1,2,2), plot(f_ne/pi,p_ne(:,5)); xlim([0 2]*1e-8)

Fs = 1/86400; % daily frequency in Hz
sw_temp_fft = fft(data_daily{1}.temp);
Lsw = length(data_daily{1}.time);
P2 = abs(sw_temp_fft./Lsw);
P1 = P2(1:Lsw/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(Lsw/2))/Lsw;
figure; plot(f,P1); xlim([0 20]*1e-8)

%% Alternative: Using T data from a longer mooring in Kangerlussuaq fjord

kf_data = process_noaa_ctd('~/data_common/greenland/obs/NOAA/0127320/1.1/data/1-data/','GP09_GP3_TidbiT');
Fs = 1/86400; % daily frequency in Hz
P_total=zeros([1,ceil(length(kf_data.time)/2)]);
for k=1:length(sf_data.depth)
    kf_temp_fft = fft(kf_data.temp(k,:));
    Lsw = length(kf_data.time);
    P2 = abs(kf_temp_fft./Lsw);
    P1 = P2(1:Lsw/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(Lsw/2))/Lsw;
    P_total = P_total+P1;
end
figure; plot((1./f)./86400,P_total./length(sf_data.depth)); xlim([0 60])
ylabel('Power'); xlabel('Period (days)')% x axis in period(days) instead of frequency
