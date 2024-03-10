close all, clear, clc

% load file
load('cntrl_diameters.mat');

% Parameters for statistical analysis
alpha = 0.05;
tails = 'both';
varia = 'unequal';
procedure = 'bonferroni';

% experimental groups
doxy_label = {'-D';'+D'};
matr_label = {'-M';'+M'};

diam_dat = [];
diam_lab = [];
count = 1;
for i = 1:numel(doxy_label)
    for j = 1:numel(matr_label)
        label = cellstr([char(doxy_label(i)) '/' char(matr_label(j))]);
        diam_dat = [diam_dat; diameter_final{count}];
        diam_lab = [diam_lab; repmat(label,numel(diameter_final{count}),1)];
        count = count+1;
    end
end

[p,tbl,stats] = anovan(diam_dat, {diam_lab});
[c,m,h,nms] = multcompare(stats,'Estimate','anovan','Alpha',alpha,'Ctype',procedure);
eval(['saveas(gcf,''Anova_' procedure '_contr_diam'',''fig'');']);