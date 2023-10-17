% load('/Users/mmeb1/data_common/greenland/FjordMIX/boxmodel/example_fjord_crash_unknown.mat','fjord_run')
load('/Users/mmeb1/data_common/greenland/FjordMIX/boxmodel/example_adjacent_collapse_unclear.mat','fjord_run')

plot_forcings_summary(fjord_run)

[fjord_run.s,fjord_run.f] = boxmodel(fjord_run.p, fjord_run.t, fjord_run.f, fjord_run.a);