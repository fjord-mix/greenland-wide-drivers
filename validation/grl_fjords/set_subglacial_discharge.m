% Gets the subglacial discharge for OMG-compiled data
% to be run within `driver_grl_fjrods.m`

folder_qsg = dir([data_path,'/greenland/runoff/Karlsson2023/']);
        
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
taxis_tsg = datetime(Qsg_time,'InputFormat','uuuu-MM') + days(14);
which_dates = taxis_tsg > time_start & taxis_tsg < time_end;
taxis_tsg = taxis_tsg(which_dates);
Qsg_sum   = Qsg_sum(which_dates);
seconds_in_month = eomday(taxis_tsg.Year,taxis_tsg.Month)*86400; 

% set f.Qsg and f.tsg accordingly
tsg = NaT([1,length(taxis_tsg)*n_years]);
i_yr_beg=1;
i_yr_end=12;
% tsg(i_yr_beg:i_yr_end) = juliandate(taxis_tsg) - juliandate(taxis_tsg(1));
for i_yr=1:n_years    
    tsg(i_yr_beg:i_yr_end) = taxis_tsg + years(i_yr-1);
    i_yr_beg=i_yr_end+1;
    i_yr_end=i_yr_beg+11;
end
f.tsg = juliandate(tsg) - juliandate(tsg(1));
Qsg = Qsg_sum./seconds_in_month;

% makes sure they are in the correct dimensions
if size(Qsg,1) > size(Qsg,2)
    Qsg = Qsg';
end

f.Qsg = repmat(Qsg,1,n_years);