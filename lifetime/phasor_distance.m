function dist = phasor_distance(data_g,data_s,lifetimes_ns,met_factor,save_path, plot_phasor)

%% PLOT SETTINGS

% Define the size of plotting windows
size_cent = [0.25 0.1 0.5 0.75];

% Define the style of the plots
fntsiz = 25;
ax_width = 2;
mar_size = 10;
li_width = 1;

% Image size
imag_size = size(data_g);

%% PHASOR COORDINATES

% draw phasor unit circle
g_bound = 0:0.01:1;
s_bound = sqrt(0.25 - (g_bound - 0.5).^2);

% define phasor coordinates for pure substance
omega = 2*pi*80.2e6;                                                       % omega = 2*pi*160.4e6;
lifetimes_s = (1e-9) * lifetimes_ns';
omega_tau = omega*lifetimes_s;
denoms = (1+(omega_tau).^2);
g_coords = 1./denoms;
s_coords = omega_tau./denoms;
lt_coords = [g_coords s_coords];

% calculate the distance of each pixel from pure substance 
v_ref = [g_coords(1)-g_coords(2);s_coords(1)-s_coords(2)];
dist = NaN(imag_size(1),imag_size(2));
for i = 1:imag_size(1)
    for j = 1:imag_size(2)
        v_dat = [data_g(i,j)-g_coords(2);data_s(i,j)-s_coords(2)];
        dist(i,j) = dot(v_dat,v_ref) / norm(v_ref);
    end
end

% ANUP'S CODE
% fit=[[1;1] g_coords(:)]\s_coords(:); %y-.1478=m*(x-.9777);
% v=[1;fit(2)];
% p_ref=[g_coords(1)-g_coords(2);s_coords(1)-s_coords(2)];
% dist_ref=dot(p_ref,v)/dot(v,v);
% dist=(data_g-g_coords(2))+((data_s-s_coords(2)).*fit(2))/dot(v,v)/dist_ref;

%% PLOT RESULTS
if plot_phasor
    
    % plot unit circle
    fig = figure;
    set(fig,'Units','Normalized','OuterPosition',size_cent);
    hold on, plot(g_bound, s_bound, 'k', 'LineWidth', li_width)
    xlabel('G')
    ylabel('S')
    hold on, plot(g_coords,s_coords,'.','MarkerSize',2*mar_size)
    switch met_factor
        case 'NADH'
            text(g_coords-.05,s_coords-.05,{'Free','Bound'})
        case 'FAD'
            text(g_coords-.06,s_coords-.06,{'Bound','Free'})
    end
    hold on, scatter(reshape(data_g,imag_size(1)*imag_size(2),1),reshape(data_s,imag_size(1)*imag_size(2),1),mar_size,reshape(dist,imag_size(1)*imag_size(2),1),'filled')
    c = colorbar;
    c.Label.String = [met_factor ' \alpha_1 (-)'];
    c.Limits = [0 0.4];
    c.Ticks = [0 0.2 0.4];
    xlim([0 1]);
    xticks([0 0.5 1]);
    ylim([0 0.6]);
    yticks([0 0.3 0.6]);
    set(gca,'FontSize',fntsiz,'LineWidth',ax_width),
    box off,
    set(fig,'PaperPositionMode','auto');
    print('-dtiff','-r300',save_path);
    saveas(fig,save_path);
    close(fig)
end

end