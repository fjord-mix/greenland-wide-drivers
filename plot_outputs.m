function hf = plot_outputs(hf,fjord_model1,fjord_model2)

model_runtime1 = fjord_model1.s.t(1:size(fjord_model1.s.H,2));
runtime_axis = fjord_model1.m.time_axis;
t0 = convertTo(runtime_axis(1),'datenum');
taxis1 = NaT([size(fjord_model1.s.H,2),1]);
for i_time=1:length(taxis1)
    taxis1(i_time) = datetime(t0+model_runtime1(i_time),'ConvertFrom','datenum');
end


m=3; n=1;
if isempty(hf)
    hf=figure('Position',[20 20 800 600],'Name',fjord_model1.m.name); 
else
    set(gcf,hf);
end
subplot(m,n,1); plot(taxis1,fjord_model1.s.H,'linewidth',1.5); hold on;
ylabel('Thickness (m)'); legend('1','2','3','4','Location','east');
subplot(m,n,2); plot(taxis1,fjord_model1.s.T,'linewidth',1.5); hold on;
ylabel('Temperature (^oC)'); 
subplot(m,n,3); plot(taxis1,fjord_model1.s.S,'linewidth',1.5); hold on;
ylabel('Salinity'); xlabel('Time')

if nargin > 2
    model_runtime2 = fjord_model2.s.t(1:size(fjord_model1.s.H,2));
    runtime_axis = fjord_model2.m.time_axis;
    t0 = convertTo(runtime_axis(1),'datenum');
    taxis2 = NaT([size(fjord_model2.s.H,2),1]);
    for i_time=1:length(taxis2)
        taxis2(i_time) = datetime(t0+model_runtime2(i_time),'ConvertFrom','datenum');
    end

    subplot(m,n,1); plot(taxis2,fjord_model2.s.H,'linewidth',1.5,'LineStyle','--'); 
    subplot(m,n,2); plot(taxis2,fjord_model2.s.T,'linewidth',1.5,'LineStyle','--');
    subplot(m,n,3); plot(taxis2,fjord_model2.s.S,'linewidth',1.5,'LineStyle','--');
end

end