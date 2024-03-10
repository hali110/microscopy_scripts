close all, clear, clc 

%% Define parameters to read .mat files

origin_folder = pwd;
base_folder = '/Volumes/NikonAX/Data/';
% dates = {'2022-04-16', '2022-04-02', '2022-04-07', '2022-04-19', '2022-04-29', '2022-03-28'};
% dox_groups = {[0,0,1,1] [0,0,1,1,1] [0,0,0,1,1,1] [0,0,0,1,1,1] [0,0,0,1,1] [0,0,0,1,1,1]};
col_con = {'1', '1', '2', '3', '4'};
dates = {'2022-04-16', '2022-04-02', '2022-04-07', '2022-04-29', '2022-03-28'};
dox_groups = {[0,0,1,1] [0,0,1,1,1] [0,0,0,1,1,1] [0,0,0,1,1] [0,0,0,1,1,1]};


%% Read and aggregate data 

orig_radii = cell(1,numel(dates));
for date = 1:numel(dates)
    
    load([base_folder dates{date} '/radii.mat'])
    orig_radii{date} = spheroid_rad_um;
    
    
end

% Rearrange data seperated by Doxycycline
radii_md = cell(4,4);
radii_pd = cell(4,4);

for i = 1:numel(orig_radii)
    
    data = orig_radii{i};
    col = cell2mat(col_con(i));
    groups = logical(cell2mat(dox_groups(i)));
    for day = 1:size(data,1)
        
        radii_md{day, str2double(col)} = [radii_md{day, str2double(col)}; data(day, ~groups)'];
        radii_pd{day, str2double(col)} = [radii_pd{day, str2double(col)}; data(day, groups)'];

    end
end


% Take avgs and stds
nanavg = @(x) (mean(x, 'omitnan'));

radii_md_avg = cell(4,4);
radii_md_sem = cell(4,4);
radii_pd_avg = cell(4,4);
radii_pd_sem = cell(4,4);

for row = 1:size(radii_md, 1)
    for col = 1:size(radii_md, 2)
        
        radii_md_avg{row, col} = cellfun(nanavg, radii_md(row,col));
        radii_pd_avg{row, col} = cellfun(nanavg, radii_pd(row,col));
        
        radii_md_sem{row, col} = std(cell2mat(radii_md(row,col)), 'omitnan') / size(radii_md{row,col}, 1);
        radii_pd_sem{row, col} = std(cell2mat(radii_pd(row,col)), 'omitnan') / size(radii_pd{row,col}, 1);

    end
end


%% Plot

N = size(unique(col_con), 2);

% define the size of plotting windows
size_cent = [0.25 0.15 0.50 0.80];
size_oriz = [0    0.15 1    0.80];
size_full = [0    0    1    1   ];

% define the style of the plots
li_type = '-o';
ax_width = 4;
li_width = 3;
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





%% -D radius vs days (series: conc)
fig = figure;
set(fig,'Units','Normalized','OuterPosition',size_cent);

% 3 because we are discarding 3mgml because it was embedded day3
x = 0:2:6;
conc = [1,2,4];
for i = 1:3
    
    nan_vals = isnan([radii_md_avg{:, conc(i)}]);
    y = [radii_md_avg{:, conc(i)}];
    er_bars = [radii_md_sem{:, conc(i)}];
    
    
    
    plot(x(~nan_vals), y(~nan_vals), li_type, 'DisplayName', [num2str(conc(i)) ' mg/mL'], 'Color', colors{i,1}, 'LineWidth', li_width)
    hold on
    errorbar(x(~nan_vals), y(~nan_vals), er_bars(~nan_vals), '.', 'DisplayName', '', 'Color', colors{i,1}, 'LineWidth', li_width-.75, 'HandleVisibility', 'off')
end

title({'Spheroids invading In Collagen Without TAZ', 'Activation (-Doxycyline)'}, 'FontSize',fntsiz)
legend('FontSize',fntsiz,'Location','northwest');
legend boxoff,
xlabel('Time (Days)','FontSize',fntsiz),
ylabel('Spheroid Radius (\mum)','FontSize',fntsiz),
xlim([0 6]),
ylim([0 800]),
set(gca,'FontSize',fntsiz,'LineWidth',ax_width),
set(fig,'PaperPositionMode','auto');
box off
print('-dtiff','-r600', [origin_folder '/taz_graphs/radius_v_days_minusd']);
saveas(fig,[origin_folder '/taz_graphs/radius_v_days_minusd']);



%% +D radius vs days (series: conc)
fig = figure;
set(fig,'Units','Normalized','OuterPosition',size_cent);

% 3 because we are discarding 3mgml because it was embedded day3
x = 0:2:6;
conc = [1,2,4];
for i = 1:3
    
    nan_vals = isnan([radii_pd_avg{:, conc(i)}]);
    y = [radii_pd_avg{:, conc(i)}];
    er_bars = [radii_pd_sem{:, conc(i)}];
    
    
    
    plot(x(~nan_vals), y(~nan_vals), li_type, 'DisplayName', [num2str(conc(i)) ' mg/mL'], 'Color', colors{i,1}, 'LineWidth', li_width)
    hold on
    errorbar(x(~nan_vals), y(~nan_vals), er_bars(~nan_vals), '.', 'DisplayName', '', 'Color', colors{i,1}, 'LineWidth', li_width-.75, 'HandleVisibility', 'off')
end

title({'Spheroids invading In Collagen With TAZ', 'Activation (+Doxycyline)'}, 'FontSize',fntsiz)
legend('FontSize',fntsiz,'Location','northwest');
legend boxoff,
xlabel('Time (Days)','FontSize',fntsiz),
ylabel('Spheroid Radius (\mum)','FontSize',fntsiz),
xlim([0 6]),
ylim([0 800]),
set(gca,'FontSize',fntsiz,'LineWidth',ax_width),
set(fig,'PaperPositionMode','auto');
box off
print('-dtiff','-r600', [origin_folder '/taz_graphs/radius_v_days_plusd']);
saveas(fig,[origin_folder '/taz_graphs/radius_v_days_plusd']);





%% +D radius vs conc (series: days)
fig = figure;
set(fig,'Units','Normalized','OuterPosition',size_cent);
x = [1, 2, 4];
for i = 1:4
    
    nan_vals = isnan([radii_pd_avg{i, x}]);
    y = [radii_pd_avg{i, x}];
    er_bars = [radii_pd_sem{i, x}];
    
    plot(x(~nan_vals), y(~nan_vals), li_type, 'DisplayName', ['Day ' num2str((2 * (i - 1)))], 'Color', colors{i,1}, 'LineWidth', li_width)
    hold on
    errorbar(x(~nan_vals), y(~nan_vals), er_bars(~nan_vals), '.', 'DisplayName', '', 'Color', colors{i,1}, 'LineWidth', li_width-.75, 'HandleVisibility', 'off')
end
title({'Spheroids invading In Collagen With TAZ', 'Activation (+Doxycyline)'}, 'FontSize',fntsiz)
legend('FontSize',fntsiz,'Location','northwest');
legend boxoff,
xlabel('Collagen Concentration (mg/mL)','FontSize',fntsiz),
ylabel('Spheroid Radius (\mum)','FontSize',fntsiz),
xlim([0 5]),
ylim([0 800]),
set(gca,'FontSize',fntsiz,'LineWidth',ax_width),
set(fig,'PaperPositionMode','auto');
box off
print('-dtiff','-r600', [origin_folder '/taz_graphs/radius_v_conc_plusd']);
saveas(fig,[origin_folder '/taz_graphs/radius_v_conc_plusd']);



%% -D radius vs conc (series: days)
fig = figure;
set(fig,'Units','Normalized','OuterPosition',size_cent);
x = [1, 2, 4];
for i = 1:4
    
    nan_vals = isnan([radii_md_avg{i, x}]);
    y = [radii_md_avg{i, x}];
    er_bars = [radii_md_sem{i, x}];
    
    plot(x(~nan_vals), y(~nan_vals), li_type, 'DisplayName', ['Day ' num2str((2 * (i - 1)))], 'Color', colors{i,1}, 'LineWidth', li_width)
    hold on
    errorbar(x(~nan_vals), y(~nan_vals), er_bars(~nan_vals), '.', 'DisplayName', '', 'Color', colors{i,1}, 'LineWidth', li_width-.75, 'HandleVisibility', 'off')
end
title({'Spheroids invading In Collagen Without TAZ', 'Activation (-Doxycyline)'}, 'FontSize',fntsiz)
legend('FontSize',fntsiz,'Location','northwest');
legend boxoff,
xlabel('Collagen Concentration (mg/mL)','FontSize',fntsiz),
ylabel('Spheroid Radius (\mum)','FontSize',fntsiz),
xlim([0 5]),
ylim([0 800]),
set(gca,'FontSize',fntsiz,'LineWidth',ax_width),
set(fig,'PaperPositionMode','auto');
box off
print('-dtiff','-r600', [origin_folder '/taz_graphs/radius_v_conc_minusd']);
saveas(fig,[origin_folder '/taz_graphs/radius_v_conc_minusd']);





%% Statistical Analysis



% Parameters for statistical analysis
alpha = 0.05;
tails = 'both';
varia = 'unequal';
procedure = 'bonferroni';

% experimental groups
doxy_label = {'-D';'+D'};
conc_label = {'1'; '2'; '4'};

rad_dat = [];
rad_lab = [];


for i = 1:numel(conc_label)
    label = cellstr(char(conc_label(i)));
    rad_dat = [rad_dat; radii_pd{4, str2double(char(conc_label(i)))}];
    rad_lab = [rad_lab; repmat(label, numel(radii_pd{4, str2double(char(conc_label(i)))}), 1)];
end



% comparing collagen concentration on day4 data
[p,tbl,stats] = anova1(rad_dat, rad_lab);
figure
[c,m,h,nms] = multcompare(stats,'Estimate','anova1','Alpha',alpha,'Ctype',procedure);
eval(['saveas(gcf,''Anova_' procedure '_diam'',''fig'');']);







