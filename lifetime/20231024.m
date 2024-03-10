close all, clear, clc

%% PLOT SETTINGS

% Define the size of plotting windows
size_cent = [0.25 0.1 0.5 0.75];
size_full = [0 0 1 1];

% Define the style of the plots
fntsiz = 25;
ax_width = 2;
mar_size = 10;
li_width = 1;

% scale bar
pix_per_um = 1.0057;
scalebar_length = 250;  % scalebar will be 10 micrometer long
unit = sprintf('%sm', '\mu'); % micrometer

%% EXPERIMENTAL SAMPLES AND ANALYSIS PARAMETERS

addpath(genpath('/Users/haider/Documents/tmr/codes/Haider'));

% data folder
folder_flim = '/Volumes/haider/tmr/data/Olympus MPE-RS TWIN/FLIM/Data/20231024-YAPCNTRL_2mg_wbeads/tifs';
figure_name = '20231024_MCF10AYAP_wbeads';

% groups
days = {'day0', 'day2', 'day4'};
days_label = {'Day 0', 'Day 2', 'Day 4'};

stiffness = {'2mgml'};
stiff_label = stiffness;

type = {'+D'};
type_label = {'+Doxy'};

cell_lines = {'CNTRL', 'YAP'};

% number of spheroids per day and cell lines
N_sphe = [3 3 3;            % cntrl - day0, 2, 4
          3 2 3];           % YAP - day0, 2, 4

% excitation wavelengths
wavelength = {'740nm'};
met_factor = {'NADH'};

% folders for ilastik masks and results
thresh_method = 'ilastik';                                                 % {'otsu';'ilastik'}
folder_mask = [folder_flim '/ilastik'];

% imaging parameters
num_pix = 512;
num_tiles = [1 1 4];
overlap = .1;

% analysis parameters to be added here...
filter_size = 5;                                                           % pixel size of the neighborhood used for median filtering
lifetimes_NADH = [0.4,3.0];                                                % reference lifetimes for NADH (ns): {free,bound}
lifetimes_FAD = [0.1,2.6];                                                 % reference lifetimes for FAD (ns): {bound,free}
mask_vals = [2 1 3];                                                       % cytoplasm nucleas background

% plot phasor ?
plot_phasor = 0;

% Define colors for comparison between groups
colors = ['b';'r';'c';'g';'m';'k'];

% result folders
res_folder = [folder_flim '/results'];
mkdir(res_folder)
pdf_folder = [res_folder '/pdfs'];
mkdir(pdf_folder)
heat_folder = [res_folder '/heatmaps'];
mkdir(heat_folder)

%% READ AND MASK INTENSITY IMAGES + CALCULATE alpha1

for i = 1:numel(days)

    for line = 1:numel(cell_lines)

        for u = 1:N_sphe(line, i)

            for w = 1:numel(wavelength)

                folder_data = [folder_flim '/' char(days{i}) '/' char(cell_lines{line}) '/spheroid' num2str(u) '/'];
                folder_resu = [folder_flim '/' char(days{i}) '/' char(cell_lines{line}) '/spheroid' num2str(u) '/results'];

                mkdir(folder_resu)

                % import data
                dir_in = dir([folder_data '/*Intensity.tif']);
                dir_gs = dir([folder_data '/*RawG.tif']);
                dir_ss = dir([folder_data '/*RawS.tif']);

                % order array to match z-stacks
                num_slices = size(dir_in,1);
                files_in = cell(num_slices,1);
                files_gs = cell(num_slices,1);
                files_ss = cell(num_slices,1);
                for n = 1:num_slices
                    files_in(n) = cellstr(dir_in(n).name);
                    files_gs(n) = cellstr(dir_gs(n).name);
                    files_ss(n) = cellstr(dir_ss(n).name);
                end
                files_in = natsort(files_in);
                files_gs = natsort(files_gs);
                files_ss = natsort(files_ss);

                % define spheroid name (remove z slices and wavelengths)
                start_index1 = regexp(files_in{1}, '_z');
                name_mono = files_in{1}(1:start_index1-1);

                clear start_index1 
                
                z_per_slice = num_slices / num_tiles(i);

                for tile = 1:num_tiles(i)

                    for n = 1:z_per_slice

                        index = ((tile-1) * z_per_slice) + n;

                        % read tif seperate for g,s,intensity
                        t = Tiff([folder_data '/' files_gs{index}], 'r');
                        data(tile, n).raw_g = read(t);

                        t = Tiff([folder_data '/' files_ss{index}], 'r');
                        data(tile, n).raw_s = read(t);

                        t = Tiff([folder_data '/' files_in{index}], 'r');
                        data(tile, n).raw_i = read(t);

                        clear t

                        % monolayer name and slice
                        data(tile, n).name = files_in{index}(1:end-14);

                        % define monolayer info
                        data(tile, n).days = days{i};
                        data(tile, n).stif = stiffness{1};
                        data(tile, n).type = type{1};
                        data(tile, n).nsph = u;
                        data(tile, n).wave = wavelength{w};
                        data(tile, n).line = cell_lines{line};

                        % generate mask
                        switch thresh_method
                            case 'otsu'
                                [mask, centroid] = mask_otsu(data(tile, n).raw_i, 2, 5);
                            case 'ilastik'
                                [mask, centroid] = mask_ilastik(folder_mask,data(tile, n).name, [mask_vals(1), mask_vals(2), mask_vals(3)]);
                        end
                        data(tile, n).mask = mask;
                        data(tile, n).mask_centroid = centroid;
                        % subplot(1,2,1), imagesc(data(tile, n).raw_i)
                        % subplot(1,2,2), imagesc(mask)
                        clear mask centroid


                        % median filtering
                        G_mask = medfilt2(data(tile, n).raw_g, [filter_size filter_size]);
                        S_mask = medfilt2(data(tile, n).raw_s, [filter_size filter_size]);

                        % evaluate distance on the phasor plot
                        save_path = [folder_resu '/' data(tile, n).name '_phasor_' char(met_factor(w))];
                        data(tile, n).dist = phasor_distance(G_mask .* data(tile, n).mask, S_mask .* data(tile, n).mask, lifetimes_NADH, char(met_factor(w)), save_path, plot_phasor);

                        % save results
                        if index == num_slices
                            results_NADH = data;
                            save([folder_resu '/' name_mono '_NADH'],'results_NADH');
                            clear results_NADH
                        end
                        clear G_mask S_mask

                    end

                end

            end
            
            clear data

            % plot heatmaps and histograms
            load([folder_resu '/' name_mono '_NADH']);

            % if num_tiles(i) ~= 1
            tile_pixs = floor(2 * num_pix - (overlap * num_pix));
            NADH_alpha1 = NaN(tile_pixs,tile_pixs,z_per_slice);
            % else
                % tile_pixs = num_pix;
                % NADH_alpha1 = NaN(tile_pixs,tile_pixs,num_slices);
            % end


            
            NADH_alpha1_proj_max = NaN(tile_pixs,tile_pixs);
            NADH_alpha1_proj_min = NaN(tile_pixs,tile_pixs);
            NADH_alpha1_hist = [];
            
            NADH_alpha2 = NaN(size(NADH_alpha1));
            NADH_alpha2_proj_max = NaN(tile_pixs,tile_pixs);
            NADH_alpha2_proj_min = NaN(tile_pixs,tile_pixs);
            NADH_alpha2_hist = [];
            
            NADH_alphaR = NaN(size(NADH_alpha1));
            NADH_alphaR_proj_max = NaN(tile_pixs,tile_pixs);
            NADH_alphaR_proj_min = NaN(tile_pixs,tile_pixs);
            NADH_alphaR_hist = [];
            
            mask = NaN(size(NADH_alpha1));
            mask_proj = NaN(tile_pixs,tile_pixs);

            if num_tiles(i) ~= 1

                overlap_pix = floor(num_pix * overlap);
                for z = 1:z_per_slice

                    % calculate parameter stacks
                    NADH_alpha1(num_pix-overlap_pix:end, 1:num_pix, z) = results_NADH(4, z).dist(:, :);
                    NADH_alpha1(num_pix-overlap_pix:end,num_pix-overlap_pix:end,z) = results_NADH(3,z).dist(:, :);
                    NADH_alpha1(1:num_pix,num_pix-overlap_pix:end,z) = results_NADH(2,z).dist(:, :);
                    NADH_alpha1(1:num_pix,1:num_pix,z) = results_NADH(1,z).dist;

                    mask(num_pix-overlap_pix:end, 1:num_pix, z) = results_NADH(4, z).mask(:, :);
                    mask(num_pix-overlap_pix:end,num_pix-overlap_pix:end,z) = results_NADH(3,z).mask(:, :);
                    mask(1:num_pix,num_pix-overlap_pix:end,z) = results_NADH(2,z).mask(:, :);
                    mask(1:num_pix,1:num_pix,z) = results_NADH(1,z).mask;

                end

                NADH_alpha2(:,:,:) = ones(tile_pixs,tile_pixs, z_per_slice) - NADH_alpha1(:,:,:);
                NADH_alphaR(:,:,:) = NADH_alpha1(:,:,:) ./ NADH_alpha2(:,:,:);

                % arrange parameters into 1D arrays for histograms
                NADH_alpha1_hist = NADH_alpha1(isfinite(NADH_alpha1 .* mask));
                NADH_alpha2_hist = NADH_alpha2(isfinite(NADH_alpha2 .* mask));
                NADH_alphaR_hist = NADH_alphaR(isfinite(NADH_alphaR .* mask));

    
            else

                for z = 1:z_per_slice

                    % calculate parameter stacks

                    NADH_alpha1(230:741,230:741,z) = results_NADH(1,z).dist;
                    NADH_alpha2(:,:,z) = ones(tile_pixs,tile_pixs) - NADH_alpha1(:,:,z);
                    NADH_alphaR(:,:,z) = NADH_alpha1(:,:,z) ./ NADH_alpha2(:,:,z);
                    mask(230:741,230:741,z) = results_NADH(1, z).mask;

                    % arrange parameters into 1D arrays for histograms
                    hist_tmp = reshape((NADH_alpha1(:,:,z) .* mask(:,:,z)),tile_pixs*tile_pixs,1);
                    NADH_alpha1_hist = [NADH_alpha1_hist; hist_tmp(isfinite(hist_tmp))];
                    clear hist_tmp
                    hist_tmp = reshape((NADH_alpha2(:,:,z) .* mask(:,:,z)),tile_pixs*tile_pixs,1);
                    NADH_alpha2_hist = [NADH_alpha2_hist; hist_tmp(isfinite(hist_tmp))];
                    clear hist_tmp
                    hist_tmp = reshape((NADH_alphaR(:,:,z) .* mask(:,:,z)),tile_pixs*tile_pixs,1);
                    NADH_alphaR_hist = [NADH_alphaR_hist; hist_tmp(isfinite(hist_tmp))];
                    clear hist_tmp

                end
    
            end

            


            % calculate projections for plotting purposes
            NADH_alpha1_proj_max = max(NADH_alpha1, [], 3);
            NADH_alpha2_proj_max = max(NADH_alpha2, [], 3);
            NADH_alphaR_proj_max = max(NADH_alphaR, [], 3);
            NADH_alpha1_proj_min = min(NADH_alpha1, [], 3);
            NADH_alpha2_proj_min = min(NADH_alpha2, [], 3);
            NADH_alphaR_proj_min = min(NADH_alphaR, [], 3);
            mask_proj = max(mask, [], 3);

            % plot NADH alpha1
            fig = figure;
            set(fig,'Units','Normalized','OuterPosition',size_cent);
            imagesc(NADH_alpha1_proj_max, 'AlphaData', mask_proj);
            daspect([1 1 1])
            c = colorbar;
            c.Label.String = 'NADH \alpha_1 max (-)';
            c.Ticks = [0 0.25 0.5];
            clim([0 .5])
            xticks([]);
            yticks([]);
            set(gca,'FontSize',fntsiz,'LineWidth',0.1),
            box off,
            set(fig,'PaperPositionMode','auto');
            print('-dtiff','-r300',[folder_resu '/' name_mono '_NADH_alpha1_max']);
            saveas(fig,[folder_resu '/' name_mono '_NADH_alpha1_max']);
            close(fig),

            % plot NADH alphaR
            fig = figure;
            ax = axes(fig);
            set(fig,'Units','Normalized','OuterPosition',size_cent);
            imagesc(NADH_alphaR_proj_max, 'AlphaData', mask_proj);
            daspect([1 1 1])
            c = colorbar;
            c.Label.String = 'NADH \alpha_1 / \alpha_2 max (-)';
            c.Ticks = [0 0.35 0.7];
            clim([0 0.7])
            xticks([]);
            yticks([]);
            set(gca,'FontSize',fntsiz,'LineWidth',0.1),
            box off,
            set(fig,'PaperPositionMode','auto');
            hScalebar = scalebar(ax, 'x', scalebar_length, unit, 'Location', 'southeast', 'ConversionFactor', pix_per_um, 'FontSize', 25, 'LineWidth', 5);
            print('-dtiff','-r300',[folder_resu '/' name_mono '_NADH_alphaR_max']);
            print('-dtiff','-r300',[heat_folder '/' name_mono '_NADH_alphaR_max']);
            saveas(fig,[folder_resu '/' name_mono '_NADH_alphaR_max']);
            close(fig),

            fig = figure;
            set(fig,'Units','Normalized','OuterPosition',size_cent);
            imagesc(NADH_alpha1_proj_min, 'AlphaData', mask_proj);
            daspect([1 1 1])
            c = colorbar;
            c.Label.String = 'NADH \alpha_1 min (-)';
            c.Ticks = [0 0.25 0.5];
            clim([0 .5])
            xticks([]);
            yticks([]);
            set(gca,'FontSize',fntsiz,'LineWidth',0.1),
            box off,
            set(fig,'PaperPositionMode','auto');
            print('-dtiff','-r300',[folder_resu '/' name_mono '_NADH_alpha1_min']);
            saveas(fig,[folder_resu '/' name_mono '_NADH_alpha1_min']);
            close(fig),

            % plot NADH alphaR
            fig = figure;
            set(fig,'Units','Normalized','OuterPosition',size_cent);
            imagesc(NADH_alphaR_proj_min, 'AlphaData', mask_proj);
            daspect([1 1 1])
            c = colorbar;
            c.Label.String = 'NADH \alpha_1 / \alpha_2 min (-)';
            c.Ticks = [0 0.25 0.5];
            clim([0 .5])
            xticks([]);
            yticks([]);
            set(gca,'FontSize',fntsiz,'LineWidth',0.1),
            box off,
            set(fig,'PaperPositionMode','auto');
            print('-dtiff','-r300',[folder_resu '/' name_mono '_NADH_alphaR_min']);
            print('-dtiff','-r300',[heat_folder '/' name_mono '_NADH_alphaR_min']);
            saveas(fig,[folder_resu '/' name_mono '_NADH_alphaR_min']);
            close(fig),

          

            clear fig c

            % store data
            FLIM_parameters.NADH_alpha1 = NADH_alpha1;
            FLIM_parameters.NADH_alpha1_proj_max = NADH_alpha1_proj_max;
            FLIM_parameters.NADH_alpha1_proj_min = NADH_alpha1_proj_min;
            FLIM_parameters.NADH_alpha1_hist = NADH_alpha1_hist;
            FLIM_parameters.NADH_alpha2 = NADH_alpha2;
            FLIM_parameters.NADH_alpha2_proj_max = NADH_alpha2_proj_max;
            FLIM_parameters.NADH_alpha2_proj_min = NADH_alpha2_proj_min;
            FLIM_parameters.NADH_alpha2_hist = NADH_alpha2_hist;
            FLIM_parameters.NADH_alphaR = NADH_alphaR;
            FLIM_parameters.NADH_alphaR_proj_max = NADH_alphaR_proj_max;
            FLIM_parameters.NADH_alphaR_proj_min = NADH_alphaR_proj_min;
            FLIM_parameters.NADH_alphaR_hist = NADH_alphaR_hist;
            FLIM_parameters.mask = mask;
            FLIM_parameters.mask_proj = mask_proj;


            save([folder_resu '/' name_mono '_parameters'],'FLIM_parameters', '-v7.3');
            clear results_NADH results_FAD FLIM_parameters
            clear NADH_alpha1 NADH_alpha1_proj NADH_alpha1_hist
            clear NADH_alpha2 NADH_alpha2_proj NADH_alpha2_hist
            clear NADH_alphaR NADH_alphaR_proj NADH_alphaR_hist


        end

    end

end



%% PLOT HISTOGRAMS
close all

% fig3 = figure;
% set(fig3,'Units','Normalized','OuterPosition',size_full);

fig4 = figure;
set(fig4,'Units','Normalized','OuterPosition',size_full);

fig5 = figure;
set(fig5,'Units','Normalized','OuterPosition',size_full);

for i = 1:numel(days)

    fig2 = figure;
    set(fig2,'Units','Normalized','OuterPosition',size_full);

    for line = 1:numel(cell_lines)

        NADH_alpha1_hist_line = [];
        NADH_alpha2_hist_line = [];
        NADH_alphaR_hist_line = [];

        fig1 = figure;
        set(fig1,'Units','Normalized','OuterPosition',size_full);

        for u = 1:N_sphe(line, i)
            
            folder_resu = [folder_flim '/' char(days{i}) '/' char(cell_lines{line}) '/spheroid' num2str(u) '/results'];
            flim_params = dir([folder_resu '/*parameters.mat']);
            load([folder_resu '/' flim_params.name]);
    
            % plot histograms per spheroid
            figure(fig1),
            
            [fx,xi]=ksdensity(FLIM_parameters.NADH_alpha1_hist,[0:.001:1]);
            max_density = max(fx);
            subplot(1,3,1), hold on, 
            plot(xi,fx, 'LineStyle', '-', 'Color', colors(u), 'LineWidth', li_width, 'DisplayName', ['Spheroid ' num2str(u)]);
            % plot([xi(fx == max_density), xi(fx == max_density)], [max_density, 0], 'LineStyle', '--', 'Color', colors(u), 'LineWidth', li_width);
            xlabel('NADH \alpha_1 (-)'),
            xlim([0 1]);
            xticks([0 0.5 1]);
            ylabel('PDF'),
            set(gca,'FontSize',fntsiz/2,'LineWidth',ax_width),
            box off,
            
            [fx,xi]=ksdensity(FLIM_parameters.NADH_alpha2_hist,[0:.001:1]);
            subplot(1,3,2), hold on, plot(xi,fx, 'LineStyle', '-', 'Color', colors(u), 'LineWidth', li_width, 'DisplayName', ['Spheroid ' num2str(u)]);
            xlabel('NADH \alpha_2 (-)'),
            xlim([0 1]);
            xticks([0 0.5 1]);
            ylabel('PDF'),
            set(gca,'FontSize',fntsiz/2,'LineWidth',ax_width),
            box off,
            
            [fx,xi]=ksdensity(FLIM_parameters.NADH_alphaR_hist,[0:.001:2]);
            subplot(1,3,3), hold on, plot(xi,fx, 'LineStyle', '-', 'Color', colors(u), 'LineWidth', li_width, 'DisplayName', ['Spheroid ' num2str(u)]);
            xlabel('NADH \alpha_1 / \alpha_2 (-)'),
            xlim([0 2]);
            xticks([0 1 2]);
            ylabel('PDF'),
            set(gca,'FontSize',fntsiz/2,'LineWidth',ax_width),
            box off,

            NADH_alpha1_hist_line = [NADH_alpha1_hist_line; FLIM_parameters.NADH_alpha1_hist];
            NADH_alpha2_hist_line = [NADH_alpha2_hist_line; FLIM_parameters.NADH_alpha2_hist];
            NADH_alphaR_hist_line = [NADH_alphaR_hist_line; FLIM_parameters.NADH_alphaR_hist];


            
        end

        sgtitle(['PDFs of different samples. Day: ' days{i} ' Cell line: ' cell_lines{line}]);
        legend
        set(fig1,'PaperPositionMode','auto');
        print('-dtiff','-r300',[pdf_folder '/' figure_name '_params_spheroids_' days{i} '_' cell_lines{line}]);
        saveas(fig1,[pdf_folder '/' figure_name '_params_spheroids_' days{i} '_' cell_lines{line}]);        
        close(fig1)

        figure(fig2),
        [fx,xi]=ksdensity(NADH_alpha1_hist_line,[0:.001:1]);
        max_density = max(fx);
        subplot(1,3,1), hold on, plot(xi,fx, 'LineStyle', '-', 'Color', colors(line), 'LineWidth', li_width, 'DisplayName', cell_lines{line});
%         plot([xi(fx == max_density), xi(fx == max_density)], [max_density, 0], 'LineStyle', '--', 'Color', colors(line), 'LineWidth', li_width);
        xlabel('NADH \alpha_1 (-)'),
        xlim([0 1]);
        xticks([0 0.5 1]);
        ylabel('PDF'),
        set(gca,'FontSize',fntsiz/2,'LineWidth',ax_width),
        box off,
        
        figure(fig2),
        [fx,xi]=ksdensity(NADH_alpha2_hist_line,[0:.001:1]);
        subplot(1,3,2), hold on, plot(xi,fx, 'LineStyle', '-', 'Color', colors(line), 'LineWidth', li_width, 'DisplayName', cell_lines{line});
        xlabel('NADH \alpha_2 (-)'),
        xlim([0 1]);
        xticks([0 0.5 1]);
        ylabel('PDF'),
        set(gca,'FontSize',fntsiz/2,'LineWidth',ax_width),
        box off,
        
        figure(fig2),
        [fx,xi]=ksdensity(NADH_alphaR_hist_line,[0:.001:2]);
        subplot(1,3,3), hold on, plot(xi,fx, 'LineStyle', '-', 'Color', colors(line), 'LineWidth', li_width, 'DisplayName', cell_lines{line});
        xlabel('NADH \alpha_1 / \alpha_2 (-)'),
        xlim([0 2]);
        xticks([0 1 2]);
        ylabel('PDF'),
        set(gca,'FontSize',fntsiz/2,'LineWidth',ax_width),
        box off,

        if line == 1
            figure(fig4),
            [fx,xi]=ksdensity(NADH_alpha1_hist_line,[0:.001:1]);
            max_density = max(fx);
            subplot(1,3,1), hold on, plot(xi,fx, 'LineStyle', '-', 'Color', colors(i), 'LineWidth', li_width, 'DisplayName', days{i});
%             plot([xi(fx == max_density), xi(fx == max_density)], [max_density, 0], 'LineStyle', '--', 'Color', colors(line), 'LineWidth', li_width);
            xlabel('NADH \alpha_1 (-)'),
            xlim([0 1]);
            xticks([0 0.5 1]);
            ylabel('PDF'),
            set(gca,'FontSize',fntsiz/2,'LineWidth',ax_width),
            box off,
            
            figure(fig4),
            [fx,xi]=ksdensity(NADH_alpha2_hist_line,[0:.001:1]);
            subplot(1,3,2), hold on, plot(xi,fx, 'LineStyle', '-', 'Color', colors(i), 'LineWidth', li_width, 'DisplayName', days{i});
            xlabel('NADH \alpha_2 (-)'),
            xlim([0 1]);
            xticks([0 0.5 1]);
            ylabel('PDF'),
            set(gca,'FontSize',fntsiz/2,'LineWidth',ax_width),
            box off,
            
            figure(fig4),
            [fx,xi]=ksdensity(NADH_alphaR_hist_line,[0:.001:2]);
            subplot(1,3,3), hold on, plot(xi,fx, 'LineStyle', '-', 'Color', colors(i), 'LineWidth', li_width, 'DisplayName', days{i});
            xlabel('NADH \alpha_1 / \alpha_2 (-)'),
            xlim([0 2]);
            xticks([0 1 2]);
            ylabel('PDF'),
            set(gca,'FontSize',fntsiz/2,'LineWidth',ax_width),
            box off,

        else

            figure(fig5),
            [fx,xi]=ksdensity(NADH_alpha1_hist_line,[0:.001:1]);
            max_density = max(fx);
            subplot(1,3,1), hold on, plot(xi,fx, 'LineStyle', '-', 'Color', colors(i), 'LineWidth', li_width, 'DisplayName', days{i});
%             plot([xi(fx == max_density), xi(fx == max_density)], [max_density, 0], 'LineStyle', '--', 'Color', colors(line), 'LineWidth', li_width);
            xlabel('NADH \alpha_1 (-)'),
            xlim([0 1]);
            xticks([0 0.5 1]);
            ylabel('PDF'),
            set(gca,'FontSize',fntsiz/2,'LineWidth',ax_width),
            box off,
            
            figure(fig5),
            [fx,xi]=ksdensity(NADH_alpha2_hist_line,[0:.001:1]);
            subplot(1,3,2), hold on, plot(xi,fx, 'LineStyle', '-', 'Color', colors(i), 'LineWidth', li_width, 'DisplayName', days{i});
            xlabel('NADH \alpha_2 (-)'),
            xlim([0 1]);
            xticks([0 0.5 1]);
            ylabel('PDF'),
            set(gca,'FontSize',fntsiz/2,'LineWidth',ax_width),
            box off,
            
            figure(fig5),
            [fx,xi]=ksdensity(NADH_alphaR_hist_line,[0:.001:2]);
            subplot(1,3,3), hold on, plot(xi,fx, 'LineStyle', '-', 'Color', colors(i), 'LineWidth', li_width, 'DisplayName', days{i});
            xlabel('NADH \alpha_1 / \alpha_2 (-)'),
            xlim([0 2]);
            xticks([0 1 2]);
            ylabel('PDF'),
            set(gca,'FontSize',fntsiz/2,'LineWidth',ax_width),
            box off,
        end

%         figure(fig3),
%         [fx,xi]=ksdensity(NADH_alpha1_hist_line,[0:.001:1]);
%         subplot(1,3,1), hold on, plot(xi,fx, 'LineStyle', '-', 'Color', colors(2*(i-1)+line), 'LineWidth', li_width, 'DisplayName', [cell_lines{line} ' ' days{i}]);
%         xlabel('NADH \alpha_1 (-)'),
%         xlim([0 1]);
%         xticks([0 0.5 1]);
%         ylabel('PDF'),
%         set(gca,'FontSize',fntsiz/2,'LineWidth',ax_width),
%         box off,
%         
%         figure(fig3),
%         [fx,xi]=ksdensity(NADH_alpha2_hist_line,[0:.001:1]);
%         subplot(1,3,2), hold on, plot(xi,fx, 'LineStyle', '-', 'Color', colors(2*(i-1)+line), 'LineWidth', li_width, 'DisplayName', [cell_lines{line} ' ' days{i}]);
%         xlabel('NADH \alpha_2 (-)'),
%         xlim([0 1]);
%         xticks([0 0.5 1]);
%         ylabel('PDF'),
%         set(gca,'FontSize',fntsiz/2,'LineWidth',ax_width),
%         box off,
%         
%         figure(fig3),
%         [fx,xi]=ksdensity(NADH_alphaR_hist_line,[0:.001:2]);
%         subplot(1,3,3), hold on, plot(xi,fx, 'LineStyle', '-', 'Color', colors(2*(i-1)+line), 'LineWidth', li_width, 'DisplayName', [cell_lines{line} ' ' days{i}]);
%         xlabel('NADH \alpha_1 / \alpha_2 (-)'),
%         xlim([0 2]);
%         xticks([0 1 2]);
%         ylabel('PDF'),
%         set(gca,'FontSize',fntsiz/2,'LineWidth',ax_width),
%         box off,

    end

    figure(fig2),
    sgtitle(['PDFs of spheroids combined day: ' days{i}]);
    legend
    set(fig2,'PaperPositionMode','auto');
    print('-dtiff','-r300',[pdf_folder '/' figure_name '_params_line_' days{i}]);
    saveas(fig2,[pdf_folder '/' figure_name '_params_line_' days{i}]);
    close(fig2)



end

figure(fig4),
sgtitle(['PDFs of CNTRL over time']);
legend
set(fig4,'PaperPositionMode','auto');
print('-dtiff','-r300',[pdf_folder '/' figure_name '_params_cntrl_alldays']);
saveas(fig4,[pdf_folder '/' figure_name '_params_cntrl_alldays']);
close(fig4)

figure(fig5),
sgtitle(['PDFs of YAP over time']);
legend
set(fig5,'PaperPositionMode','auto');
print('-dtiff','-r300',[pdf_folder '/' figure_name '_params_yap_alldays']);
saveas(fig5,[pdf_folder '/' figure_name '_params_yap_alldays']);
close(fig5)

% figure(fig3),
% sgtitle('PDFs of spheroids combined over days and cell lines');
% legend
% set(fig3,'PaperPositionMode','auto');
% print('-dtiff','-r300',[pdf_folder '/' figure_name '_params_all']);
% saveas(fig3,[pdf_folder '/' figure_name '_params_all']);
% close all
