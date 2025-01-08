% Script to compare GLEAM4 to alternative datasets
% Author: Diego G. Miralles, Ghent University, diego.miralles@ugent.be
% November 2024

% Define base paths for data files
basePathGLEAM4 = 'F:\GLEAM4\';
basePathGLEAMv3 = 'E:\Users\dmiralle\GLEAMv3\daily\';
basePathFLUXCOM = 'E:\Users\dmiralle\FLUXCOM\';
basePathERA5L = 'E:\Users\dmiralle\ERA5Land\';

% Flag to redo loading and aggregating data
redo = 1; % Set to 1 to redo loading and aggregating, 0 to load saved workspace

if redo == 1
    %% Define pixel coordinates for specific locations
    Spain = [510, 1780];
    Sahel = [780, 2100];
    Belgium = [380, 1890];
    Amazon = [895, 1240];
    Sweden = [220, 2000];
    Canada = [4265, 669];
    California = [513, 570];
    Australia = [1140, 3202];
    i = cat(1, Spain, Sahel, Belgium, Amazon, Sweden, Canada, California, Australia);
    
    %% Define resampling coordinates
    [lon, lat] = meshgrid(-180 + 0.25 / 2 : 0.25 : 180 - 0.25 / 2, -90 + 0.25 / 2 : 0.25 : 90 - 0.25 / 2);
    [long, lati] = meshgrid(-180 + 0.5 / 2 : 0.5 : 180 - 0.5 / 2, -90 + 0.5 / 2 : 0.5 : 90 - 0.5 / 2);
    [newLon, newLat] = meshgrid(-180 + 0.1 / 2 : 0.1 : 180 - 0.1 / 2, -90 + 0.1 / 2 : 0.1 : 90 - 0.1 / 2);
    
    %% Initialize data containers
    EE_GLEAM = []; EEE_GLEAM = [];
    EE_GLEAMv3 = []; EEE_GLEAMv3 = [];
    EE_FLUXCOM = []; EEE_FLUXCOM = [];
    EE_ERA5L = []; EEE_ERA5L = [];
    
    %% Load and aggregate data for each year
    for year = 1980:2020
        yr = num2str(year);
        disp(yr);
    
        % Load GLEAM4 data
        E_GLEAM = ncread([basePathGLEAM4, yr, '\E_', yr, '_GLEAM_v4.1a.nc'], 'E');
        E = nansum(E_GLEAM, 3);
        E = flipud(rot90(E));
        EEE_GLEAM = cat(3, EEE_GLEAM, E);
    
        % Load GLEAMv3 data
        E_GLEAMv3 = ncread([basePathGLEAMv3, yr, '\E_', yr, '_GLEAM_v3.8a.nc'], 'E');
        E = nansum(E_GLEAMv3, 3);
        E = flipud(rot90(E));
        E = interp2(lon, lat, E, newLon, newLat, 'nearest');
        EEE_GLEAMv3 = cat(3, EEE_GLEAMv3, E);
    
        % Load FLUXCOM data
        E_FLUXCOM = ncread([basePathFLUXCOM, 'LE.RS_METEO.EBC-ALL.MLM-ALL.METEO-era5.720_360.daily.', yr, '.nc'], 'LE');
        E_FLUXCOM = E_FLUXCOM * 0.408;
        E = nansum(E_FLUXCOM, 3);
        E = flipud(rot90(E));
        E = interp2(long, lati, E, newLon, newLat, 'nearest');
        EEE_FLUXCOM = cat(3, EEE_FLUXCOM, E);
    
        % Load ERA5L data
        E_ERA5L = ncread([basePathERA5L, 'TotalE_ERA5Land_010deg_', yr, '.nc'], 'Evaporation');
        E = nansum(E_ERA5L, 3);
        EEE_ERA5L = cat(3, EEE_ERA5L, E);
    
        % Extract data for specific locations
        egleam = []; egleam3 = []; efluxcom = []; eera5 = [];
        for w = 1:8
            % Ensure indices are within bounds
            if i(w, 2) <= size(E_GLEAM, 1) && i(w, 1) <= size(E_GLEAM, 2)
                egleam = cat(2, egleam, squeeze(E_GLEAM(i(w, 2), i(w, 1), :)));
            end
            if round(i(w, 2) / 2.5) <= size(E_GLEAMv3, 1) && round(i(w, 1) / 2.5) <= size(E_GLEAMv3, 2)
                egleam3 = cat(2, egleam3, squeeze(E_GLEAMv3(round(i(w, 2) / 2.5), round(i(w, 1) / 2.5), :)));
            end
            if round(i(w, 2) / 5) <= size(E_FLUXCOM, 1) && round(i(w, 1) / 5) <= size(E_FLUXCOM, 2)
                efluxcom = cat(2, efluxcom, squeeze(E_FLUXCOM(round(i(w, 2) / 5), round(i(w, 1) / 5), :)));
            end
            if i(w, 1) <= size(E_ERA5L, 1) && i(w, 2) <= size(E_ERA5L, 2)
                eera5 = cat(2, eera5, squeeze(E_ERA5L(i(w, 1), i(w, 2), :)));
            end
        end
        EE_GLEAM = cat(1, EE_GLEAM, egleam);
        EE_GLEAMv3 = cat(1, EE_GLEAMv3, egleam3);
        EE_FLUXCOM = cat(1, EE_FLUXCOM, efluxcom);
        EE_ERA5L = cat(1, EE_ERA5L, eera5);
        
        clear E_GLEAM E_GLEAMv3 E_ERA5L E_FLUXCOM;
    end
    
    %% Final processing of aggregated data
    EEE_GLEAM = nanmean(EEE_GLEAM, 3);
    EEE_GLEAMv3 = nanmean(EEE_GLEAMv3, 3);
    EEE_FLUXCOM = nanmean(EEE_FLUXCOM, 3);
    EEE_ERA5L = nanmean(EEE_ERA5L, 3);
    
    % Save the entire workspace
    save('workspace_compare.mat');
else
    % Load the saved workspace
    load('workspace_compare.mat');
end

%% Plotting and Analysis

% Define variables and locations for plotting
variables = {'E_GLEAM', 'E_GLEAMv3', 'E_ERA5L', 'E_FLUXCOM'};
locations = {'Spain', 'Sahel', 'Belgium', 'Amazon', 'Sweden', 'Canada', 'California', 'Australia'};
yLimits = [Inf, 6, 7, 12, Inf, 7, 7, 7]; 

num_locations = min(size(EE_GLEAM, 2), numel(locations));
figure('Position', [150, 100, 2200, 3500]); 

subplotHeight = 0.7 / num_locations; 
subplotWidth = 0.9; 
subplotYSpacing = 0.02; 
leftMargin = 0.05; 

for location_idx = 1:num_locations
    location_name = locations{location_idx};

    % Load the data for each variable at this location
    E_GLEAM_data = EE_GLEAM(:, location_idx);  
    E_GLEAMv3_data = EE_GLEAMv3(:, location_idx);  
    E_FLUXCOM_data = EE_FLUXCOM(:, location_idx); 
    E_ERA5L_data = EE_ERA5L(:, location_idx); 

    % Calculate the position of the current subplot
    subplotPosition = [leftMargin, 0.1 + (num_locations - location_idx) * (subplotHeight + subplotYSpacing), subplotWidth, subplotHeight];

    % Create the subplot
    subplot('Position', subplotPosition);

    % Plot the time series for each variable at this location    
    plot(E_GLEAMv3_data, 'LineWidth', 0.8, 'DisplayName', 'E\_GLEAMv3', 'Color', [0.4 0.4 0.8]);
    hold on;
    plot(E_FLUXCOM_data, 'LineWidth', 0.8, 'DisplayName', 'E\_FLUXCOM', 'Color', [0.4 0.8 0.4]);
    plot(E_ERA5L_data, 'LineWidth', 0.8, 'DisplayName', 'E\_ERA5L', 'Color', [0.8 0.4 0.4]);
    plot(E_GLEAM_data, 'LineWidth', 0.4, 'DisplayName', 'E\_GLEAM4', 'Color', [0 0 0]); 

    % Set labels, title, and y-axis limits
    ylabel('mm/day');
    title(location_name);
    grid on;
    if yLimits(location_idx) ~= Inf
        ylim([0 yLimits(location_idx)]);
    end

    % Remove X-axis label for all subplots except the bottom one
    if location_idx ~= num_locations
        xticklabels([]); 
    else
        xlabel('Time (days)'); 
    end

    % Store legend handle from the last subplot
    if location_idx == num_locations
        legend_handle = legend('Location', 'best', 'Orientation', 'horizontal', 'FontSize', 8); 
        legend_handle.Box = 'off'; 
        legend_handle.Color = 'none'; 
    end
end

% Position the legend centered at the bottom and make it smaller
subplot('Position', [leftMargin, 0.05, subplotWidth, 0.02]);
axis off;
legend_handle.Position = [0.5 - 0.125, 0.03, 0.25, 0.02]; 

% Save as PNG and vector PDF
printfigure_new([basePathGLEAM4, 'timeseries_comparison'])

%% Area covered by gridcell at model resolution [km^2]

Area=zeros(1800,3600);
circ_pol=40008; % polar circumference: 40,008 km.
lat_km=(circ_pol/2)/1800; % lat length [km].
circ_equ=40076; % equatorial circumference: 40,076 km.
lon_km=sin(((-90+0.1/2:0.1:90-0.1/2)+90)/180 *pi).*circ_equ/3600;
cell_size=lon_km*lat_km;
for i=1:3600, Area(:,i)=cell_size; end

load Msea.mat

%% Maps of evaporation

EEE_FLUXCOM(~Msea&EEE_FLUXCOM<-1000)=0;

A=EEE_GLEAM;
A(Msea)=NaN;
A(A>2000)=2000; A(A<-10)=-10;
A_area=(A/1000000).*Area; A_area=nansum(nansum(A_area)); disp(num2str(A_area));
A(isnan(A))=-200;
map_E(A,[-25 1400])
printfigure([basePathGLEAM4, 'E_GLEAM'])

A=EEE_GLEAMv3;
A(Msea)=NaN;
A(A>2000)=2000; A(A<-10)=-10;
A_area=(A/1000000).*Area; A_area=nansum(nansum(A_area)); disp(num2str(A_area));
A(isnan(A))=-200;
map_E(A,[-25 1400])
printfigure([basePathGLEAM4, 'E_GLEAMv3'])

A=EEE_ERA5L;
A(Msea)=NaN;
A(A>2000)=2000; A(A<-10)=-10;
A_area=(A/1000000).*Area; A_area=nansum(nansum(A_area)); disp(num2str(A_area));
A(isnan(A))=-200;
map_E(A,[-25 1400])
printfigure([basePathGLEAM4, 'E_ERA5L'])

A=EEE_FLUXCOM;
A(Msea)=NaN;
A(A>2000)=2000; A(A<-10)=-10;
A_area=(A/1000000).*Area; A_area=nansum(nansum(A_area)); disp(num2str(A_area));
A(isnan(A))=-200;
map_E(A,[-25 1400])
printfigure([basePathGLEAM4, 'E_FLUXCOM'])

%% Difference maps

A=EEE_GLEAM-EEE_GLEAMv3;
A(Msea)=NaN; A(A<-300)=-300;
A(isnan(A))=-305;
map_E_an(A,[-305 300])
printfigure([basePathGLEAM4, 'E_GLEAMv3_dif'])

A=EEE_GLEAM-EEE_ERA5L;
A(Msea)=NaN; A(A<-300)=-300;
A(isnan(A))=-305;
map_E_an(A,[-305 300])
printfigure([basePathGLEAM4, 'E_ERA5L_dif'])

A=EEE_GLEAM-EEE_FLUXCOM;
A(Msea)=NaN; A(A<-300)=-300;
A(isnan(A))=-305;
map_E_an(A,[-305 300])
printfigure([basePathGLEAM4, 'E_FLUXCOM_dif'])

%% Functions for Map Plotting

function map_E(A, cax)
    figure;
    A = double(A);
    lon_range = -179.95:0.1:179.95;
    lat_range = 89.95:-0.1:-89.95;
    [lon, lat] = meshgrid(lon_range, lat_range);

    load diecolor_EE.mat % load preferred colorscale

    % Set up the Plate Carrée (Equidistant Cylindrical) Projection
    m_proj('Equidistant', 'lon', [-180, 180], 'lat', [-75, 85]);

    % Plot the data using the custom colormap
    m_pcolor(lon, lat, A);
    shading flat;
    colormap(diecolor_EE);
    clim(cax);

    m_coast('line', 'Color', [.3, .3, .3]);
    m_grid('xaxis', 'middle', 'xticklabels', [], 'yticklabels', []);

    % Add a standard colorbar.
    colorbar('h');
end

function map_E_an(A, cax)
    figure;
    A = double(A);
    lon_range = -179.95:0.1:179.95;
    lat_range = 89.95:-0.1:-89.95;
    [lon, lat] = meshgrid(lon_range, lat_range);

    load diecolor_EE_an.mat % load preferred colorscale

    % Set up the Plate Carrée (Equidistant Cylindrical) Projection
    m_proj('Equidistant', 'lon', [-180, 180], 'lat', [-75, 85]);

    % Plot the data using the custom colormap
    m_pcolor(lon, lat, A);
    shading flat;
    colormap(diecolor_EE_an);
    clim(cax);

    m_coast('line', 'Color', [.3, .3, .3]);
    m_grid('xaxis', 'middle', 'xticklabels', [], 'yticklabels', []);

    % Add a standard colorbar.
    colorbar('h');
end
