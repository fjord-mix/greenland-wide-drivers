function plot_density(fjord_model1)

if ~isfield(fjord_model1.m,'name'), fjord_model1.m.name='Unamed'; end
model_runtime1 = fjord_model1.s.t(1:size(fjord_model1.s.H,2));
if isfield(fjord_model1.m,'time_axis')
    runtime_axis = fjord_model1.m.time_axis;
    t0 = convertTo(runtime_axis(1),'datenum');
    taxis1 = NaT([size(fjord_model1.s.H,2),1]);
    isdays='';
    for i_time=1:length(taxis1)
        taxis1(i_time) = datetime(t0+model_runtime1(i_time),'ConvertFrom','datenum');
    end
else
    taxis1=fjord_model1.s.t;
    isdays = ' (days)';
end

layer_lbls =cell(size(fjord_model1.a.H0));
for i=1:fjord_model1.p.N+fjord_model1.p.sill
    layer_lbls{i}=num2str(i);
end

dens = gsw_sigma0(fjord_model1.s.S,fjord_model1.s.T);

figure('Position',[20 20 800 300],'Name',fjord_model1.m.name); 
plot(taxis1,dens,'linewidth',1.5); hold on;
ylabel('Density $\sigma_0$ ($kgm^{-3}$)','interpreter','latex'); 
xlabel(['Model time',isdays])
hl = legend(layer_lbls,'Location','southeast'); title(hl,'Layer'); hl.NumColumns=2;
if nargin > 1
    model_runtime2 = fjord_model2.s.t(1:size(fjord_model1.s.H,2));
    runtime_axis = fjord_model2.m.time_axis;
    t0 = convertTo(runtime_axis(1),'datenum');
    taxis2 = NaT([size(fjord_model2.s.H,2),1]);
    for i_time=1:length(taxis2)
        taxis2(i_time) = datetime(t0+model_runtime2(i_time),'ConvertFrom','datenum');
    end
    dens2 = gsw_sigma0(fjord_model2.s.S,fjord_model2.s.T);

    plot(taxis2,dens2,'linewidth',1.5,'LineStyle','--'); 

end

end