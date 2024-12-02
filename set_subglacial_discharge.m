% Gets the subglacial discharge for OMG-compiled data
% to be run within `driver_grl_fjrods.m`

folder_qsg = dir([data_path,'/greenland_common/runoff/Karlsson2023/']);
        
% for each existing gateID:
if ~isempty(fjord_matrix.qsg_id1(i_fjord))
    % look for the matching files in folder_qsg (filename contains gateID)
    for i=1:length(folder_qsg)
        if contains(folder_qsg(i).name,num2str(fjord_matrix.qsg_id1(i_fjord)))
            qsg_gate1 = readtable([folder_qsg(i).folder,'/',folder_qsg(i).name]); % read in table
            if ~isempty(qsg_gate1)
                Qsg_sum = qsg_gate1.SurfaceMelt + qsg_gate1.BasalMelt;
                Qsg_time = qsg_gate1.Date_YYYY_MM;
            else
                Qsg_sum = 0;
            end
        end
    end
end
if ~isempty(fjord_matrix.qsg_id2(i_fjord))
    % look for the matching files in folder_qsg (filename contains gateID)
    for i=1:length(folder_qsg)
        if contains(folder_qsg(i).name,num2str(fjord_matrix.qsg_id2(i_fjord)))
            qsg_gate2 = readtable([folder_qsg(i).folder,'/',folder_qsg(i).name]); % read in table
            if ~isempty(qsg_gate2)
                Qsg_sum = Qsg_sum + qsg_gate2.SurfaceMelt + qsg_gate2.BasalMelt;
                if ~exist('Qsg_time','var')
                    Qsg_time = qsg_gate2.Date_YYYY_MM;
                end
            end
        end
    end
end
if ~isempty(fjord_matrix.qsg_id3(i_fjord))
    % look for the matching files in folder_qsg (filename contains gateID)
    for i=1:length(folder_qsg)
        if contains(folder_qsg(i).name,num2str(fjord_matrix.qsg_id3(i_fjord)))
            qsg_gate3 = readtable([folder_qsg(i).folder,'/',folder_qsg(i).name]); % read in table
            if ~isempty(qsg_gate3)
                Qsg_sum = Qsg_sum + qsg_gate3.SurfaceMelt + qsg_gate3.BasalMelt;
                if ~exist('Qsg_time','var')
                    Qsg_time = qsg_gate3.Date_YYYY_MM;
                end
            end
        end
    end
end
clear qsg_gate1 qsg_gate2 qsg_gate3


% trim by Date_YYYY-MM based on our time_start and time_end
taxis_tsg = datetime(Qsg_time,'InputFormat','uuuu-MM') + days(15);
which_dates = taxis_tsg >= time_start & taxis_tsg <= time_end;
taxis_tsg = taxis_tsg(which_dates);
Qsg_sum   = Qsg_sum(which_dates);
n_days_in_month = eomday(taxis_tsg.Year,taxis_tsg.Month);
seconds_in_month = n_days_in_month*86400;
taxis_tsg_daily = NaT([sum(n_days_in_month),1]);
i_doy = 1;
for i_month = 1:length(taxis_tsg)
    for i_day = 1:n_days_in_month(i_month)
        taxis_tsg_daily(i_doy) = datetime(taxis_tsg.Year(i_month),taxis_tsg.Month(i_month),i_day,0,0,0);
        i_doy=i_doy+1;
    end
end

% set f.Qsg and f.tsg accordingly
if exist('n_years','var') % we treat things slightly differently in case we do a repeat of one year, or several years
    tsg = NaT([length(taxis_tsg_daily)*n_years,1]);
    i_beg_yr=1;
    i_end_yr=i_beg_yr+sum(n_days_in_month)-1;
    for i_yr=1:n_years
        tsg(i_beg_yr:i_end_yr) = taxis_tsg_daily + years(i_yr-1);
        i_beg_yr=i_end_yr+1;
        i_end_yr=i_beg_yr+sum(n_days_in_month)-1;
    end
    f.tsg = juliandate(tsg) - juliandate(tsg(1));
else
    f.tsg = juliandate(taxis_tsg) - juliandate(taxis_tsg(1));
end
Qsg = Qsg_sum./seconds_in_month;

% Making the monthly values at daily resolution (FjordRPM will linearly
% interpolate if we do not do that, and we do not want linear interp)
% tt_runoff = timetable(taxis_tsg,Qsg);
% tt_runoff.Properties.VariableContinuity = {'step'};
% tt_runoff_daily = retime(tt_runoff,taxis_tsg_daily);

Qsg_daily = interp1(taxis_tsg,Qsg,taxis_tsg_daily,'nearest','extrap');
% tt_runoff_daily = timetable(taxis_tsg_daily,Qsg_daily);

Qsg = Qsg_daily;
f.tsg = juliandate(tsg) - juliandate(tsg(1));

% makes sure they are in the correct dimensions
if size(Qsg,1) > size(Qsg,2)
    Qsg = Qsg';
end
if size(f.tsg,1) > size(f.tsg,2)
    f.tsg = f.tsg';
end

if exist('n_years','var')
    f.Qsg = repmat(Qsg,1,n_years);
else
    f.Qsg = Qsg;
end