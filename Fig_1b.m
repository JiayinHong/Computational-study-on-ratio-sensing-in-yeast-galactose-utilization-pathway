% this script use experimental data from Escalante et al [ref.12]
% to show ratiometric response and decision front in Fig 1b

% load experimental data
load('expt_trait_table.mat')
expt_96well = trait;
% set up labels
galLabel = {'None','2^{-8}','2^{-7}','2^{-6}','2^{-5}','2^{-4}','2^{-3}','2^{-2}','2^{-1}','1','2','4'};
gluLabel = {'None','2^{-6}','2^{-5}','2^{-4}','2^{-3}','2^{-2}','2^{-1}','1'};
% pre-process expt data
alldata = nan(8,12);
ind1 = find(expt_96well.mask_induction == 0);     % all the rows whose mask_induction == 0
tmp = expt_96well(ind1,:).mask_basal == 0;
% remove tmp from ind1, so that ind1 only contains rows whose mask_induction == 0 while mask_basal ~= 0
ind1(tmp) = [];
% use basal_level to represent induced level in ind1
expt_96well(ind1,:).ind_level = expt_96well(ind1,:).basal_level;
% the GAL1 level under the highest gal/glu ratio corresponds fully ON state
logyfp_to_nm = @(x) exp(x-7.34) * 3000; % exp(7.34) is the average fluorescent intensity while full induction
alldata(:) = logyfp_to_nm(expt_96well{:,'ind_level'});
% normalize to maximum induction level
tmp_exp = log(alldata);
max_exp = max(tmp_exp(:));

figure
set(gcf, 'position', [281 415 383 283])
h1 = imagesc(tmp_exp./max_exp);
h1.Parent.XTick = 1:12;
h1.Parent.YTick = 1:8;
h1.Parent.XTickLabel = galLabel;
h1.Parent.YTickLabel = fliplr(gluLabel);
h2=colorbar;
h2.Limits = [0 1];
h2.Ticks = [0, 0.2, 0.4, 0.6, 0.8, 1];
set(gca,'XTickLabelRotation',45)
set(gca,'FontSize',15,'FontName','Times New Roman','FontWeight','normal')
export_fig('../figures/Fig-1b', '-pdf','-transparent','-c[NaN,NaN,NaN,NaN]')


