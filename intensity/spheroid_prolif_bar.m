close all, clear, clc 

%% Define parameters to read .mat files - YAP

% save_dir = '/Volumes/NikonAX/Codes/Spheroids/yap_prolif_bars';
% target_folder = '/Volumes/NikonAX/Data/2022-01-12/analysis/';
% cd(target_folder)
% cell_line = 'MCF10A YAP 5SA';
% groups = {'group1', 'group2', 'group3'}; % group1 = 0ng/ml; group2 = 100 ng/ml; group 3 = 1000 ng/ml Doxy
% files = {'radii_minusm.mat', 'radii_plusm.mat'};
% n_groups = [14, 16, 17; 11, 17, 16]';
% 
% % looked at the contr diagram and ran next line. 
% % sigstar({[1,3] [1,4] [2,3] [2,4] [3,4]})


%% Define parameters to read .mat files - TAZ

save_dir = '/Volumes/NikonAX/Codes/Spheroids/taz_prolif_bars';
target_folder = '/Volumes/TMR/tmr_data/Nikon A1R/2021-10-10/analysis';
cd(target_folder)
cell_line = 'MCF10A RFPTAZ 4SA';
groups = {'group1', 'group2', 'group3', 'group4'};
files = {'radii'};
n_groups = [11, 11, 10, 10]';

% looked at the contr diagram and ran next line. 
% sigstar({[1,2] [1,3] [1,4] [3,4]})


%% Define parameters to read .mat files - CNTRL

% save_dir = '/Volumes/NikonAX/Codes/Spheroids/cntrl_prolif_bars';
% target_folder = '/Volumes/TMR/tmr_data/Nikon A1R/2021-10-21/analysis';
% cd(target_folder)
% cell_line = 'MCF10A RFPCNTRL';
% groups = {'group1', 'group2', 'group3', 'group4'};
% files = {'radii'};
% % {'-D/-M', '-D/+M', '+D/-M', '+D/+D'}
% n_groups = [12, 15, 13, 13]';
% 
% % looked at the contr diagram and ran next line. 
% % sigstar({[1,2] [1,4] [2,3] [3,4]})




%% Define Variables
% rows are dif dox conc, cols are matrigel
diam_final_avg = zeros(numel(groups), numel(files));
diam_final_sem = zeros(numel(groups), numel(files));
orig_radii = cell(1, numel(files));



%% Read and aggregate data 

for i = 1:numel(files)
    
    load(files{i})
    orig_radii{i} = spheroid_rad_um;

    index = 1;
    for n = 1:size(n_groups, 1)

        findex = index + (n_groups(n, i) - 1);
        diam_final_avg(n, i) = mean((spheroid_rad_um(index:findex) .* 2), 'omitnan');
        diam_final_sem(n, i) = std((spheroid_rad_um(index:findex) .* 2), 'omitnan') / sqrt(n_groups(n, i));
        
        index = findex + 1; 

    end

end







%% Plotting params

N = 4;

% define the size of plotting windows
size_cent = [0.25 0.15 0.50 0.80];
size_oriz = [0    0.15 1    0.80];
size_full = [0    0    1    1   ];

% define the style of the plots
li_type = '-o';
ax_width = 4;
li_width = 3;
bar_width = .8;
fntsiz = 25;
alpha = 0.3;
li_incr = 1;
fnt_incr = 12;
cap_width = 20;

% define colormap
fig0 = figure;
string = '1:10,1:10';
for i = 1:N-1
    string = [string ',1:10,1:10'];
end
h = eval(['plot(' string ')']);
getcolors = get(h,'Color');
for i = 1:N
    colors(i,:) = getcolors(i,1);
end
close(fig0);





%% bar plot of final spheroid diameter

graph_groups = {'-D/-M', '-D/+M', '+D/-M', '+D/+D'};


fig = figure;
set(fig,'Units','Normalized','OuterPosition',size_cent);

% TAZ and CNTRL data are in different shape so dont need this
% diams_for_gr = reshape((diam_final_avg([1 2], :)).', 1, []);
% sems_for_gr = reshape((diam_final_sem([1 2], :)).', 1, []);
diams_for_gr = diam_final_avg;
sems_for_gr = diam_final_sem;


for gr = 1:numel(graph_groups)
    
    b_dia{gr} = bar(gr, diams_for_gr(gr),'EdgeColor', 'k','LineWidth', ax_width,'barwidth', bar_width, 'FaceColor', colors{gr,:});
    hold on
    h_rad{gr} = errorbar(gr, diams_for_gr(gr),sems_for_gr(gr), 'LineStyle', 'none', 'Color', 'k', 'LineWidth', ax_width);
    h_rad{gr}.CapSize = cap_width;
    set(get(get(h_rad{gr}(1), 'Annotation'), 'LegendInformation'), 'IconDisplayStyle','off');

end



box off,
xticks(1:numel(graph_groups)),
xticklabels(graph_groups),
ylim([0 450]),
set(gca,'FontSize',fntsiz,'LineWidth',ax_width),
ylabel('Final Spheroid Diameter (\mum)','FontSize',fntsiz),
set(fig,'PaperPositionMode','auto');
title([cell_line ' Spheroid Diameter'])

print('-dtiff','-r600',[save_dir '/results_otsu_bar']);
% print('-dsvg','-r600',[save_dir '/results_otsu_bar']);
saveas(fig, [save_dir '/results_otsu_bar']);











%% statistical significance

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

orig_diams = cell2mat(orig_radii) .* 2;


% added next line for taz and cntrl
% n_groups = reshape(n_groups.', 2, []);

% move index=1 inside first for loop for yap and before for taz/cntrl


for i = 1:numel(matr_label)
    

index = 1;

    for j = 1:numel(doxy_label)

        findex = index + (n_groups(j, i) - 1);

        
        % first two lines are for yap, next two are for taz and cntrl
        diam_dat = [diam_dat; orig_diams(index:findex, i)];
        label = cellstr([char(doxy_label(j)) '/' char(matr_label(i))]);

        % diam_dat = [diam_dat; orig_diams(index:findex)];
        % label = cellstr([char(doxy_label(i)) '/' char(matr_label(j))]);

        diam_lab = [diam_lab; repmat(label, n_groups(j, i), 1)];
        
        index = findex + 1; 

    end

end

[p,tbl,stats] = anovan(diam_dat, {diam_lab});
figure
[c,m,h,nms] = multcompare(stats,'Estimate','anovan','Alpha',alpha,'Ctype',procedure);
saveas(gcf, [save_dir '/anova_' procedure '_contr_diam']);



























