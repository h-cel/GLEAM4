% Script to analyse global GLEAM4 patterns
% Author: Diego G. Miralles, Ghent University, diego.miralles@ugent.be
% January 2025

% Flag to redo loading and aggregating data
redo = 1;

% Paths
data_path = 'F:\GLEAM4\4.2\'; % Path to data files
workspace_file = 'workspace.mat'; % Path to save or load workspace

% Constants for calculating area (grid cell dimensions based on Earth circumference)
Area = zeros(1800, 3600);
circ_pol = 40008; % Polar circumference in km
lat_km = (circ_pol / 2) / 1800; % Latitude cell size in km
circ_equ = 40076; % Equatorial circumference in km
lon_km = sin(((-90 + 0.1 / 2 : 0.1 : 90 - 0.1 / 2) + 90) / 180 * pi) .* circ_equ / 3600;
cell_size = lon_km * lat_km;
for i = 1:3600, Area(:, i) = cell_size; end % Assign cell area across grid

if redo == 1
    %% Load and aggregate data
    variables = {'E', 'Ep', 'Ep_aero', 'Ep_rad', 'Et', 'Ec', 'Eb', 'Ew', 'Ei', 'Es', 'S', 'SMrz', 'H'};
    seasons = {'JJA', 'DJF'}; % Define seasons for seasonal averaging

    % Initialize empty arrays for each variable and season
    for var = variables
        varname = char(var);
        for season = seasons
            eval([varname '_' char(season) ' = [];']);
        end
        eval([varname ' = [];']);
    end

    % Load data for each year and calculate seasonal values
    for year = 1980:2023
        yr = num2str(year);
        disp(yr);

        for var = variables
            varname = char(var);
            filename = [data_path, yr, '\', varname, '_', yr, '_GLEAM_v4.2a.nc'];
            data = ncread(filename, varname);

            % Calculate seasonal means for each variable
            for season = seasons
                seasonname = char(season);
                if strcmp(seasonname, 'JJA')
                    data_season = mean(data(:, :, 154:154 + 91), 3, 'omitnan');
                else
                    data_season = (mean(data(:, :, 1:60), 3, 'omitnan') + mean(data(:, :, end - 30:end), 3, 'omitnan')) / 2;
                end
                data_season = flipud(rot90(data_season));
                eval([varname '_' seasonname ' = cat(3, ' varname '_' seasonname ', data_season);']);
            end

            % Calculate annual total and flip/rotate for each variable
            data = nansum(data, 3);
            data = flipud(rot90(data));
            eval([varname ' = cat(3, ' varname ', data);']);
        end
    end

    %% Calculate mean values over the time series
    for var = variables
        varname = char(var);
        eval([varname ' = mean(' varname ', 3, ''omitnan'');']);
        for season = seasons
            seasonname = char(season);
            eval([varname '_' seasonname ' = mean(' varname '_' seasonname ', 3, ''omitnan'');']);
        end
    end

    % Save the entire workspace
    save(workspace_file);
else
    % Load the saved workspace
    load(workspace_file);
end

%% Totals weighted per area

% Display totals of each component weighted by area
Et_area = calculate_area_weighted_total(Et, Area, S == 0);
disp(['Et_area: ', num2str(Et_area)]);

Ei_area = calculate_area_weighted_total(Ei, Area, S == 0);
disp(['Ei_area: ', num2str(Ei_area)]);

Eb_area = calculate_area_weighted_total(Eb, Area, S == 0);
disp(['Eb_area: ', num2str(Eb_area)]);

Eo = Ec + Ew + Es;
Eo_area = calculate_area_weighted_total(Eo, Area, S == 0);
disp(['Eo_area: ', num2str(Eo_area)]);

Ep_area = calculate_area_weighted_total(Ep, Area, S == 0);
disp(['Ep_area: ', num2str(Ep_area)]);

%% Component fractions

% Calculate total evaporation and fraction of each component
E_area = Et_area + Ei_area + Eb_area + Eo_area;
disp(['Et fraction: ', num2str(Et_area / E_area)]);
disp(['Ei fraction: ', num2str(Ei_area / E_area)]);
disp(['Eb fraction: ', num2str(Eb_area / E_area)]);
disp(['Eo fraction: ', num2str(Eo_area / E_area)]);

%% Latitudinal profiles of evaporation and components (Fig. 2)

f = S == 0;
et = mean_with_nan(Et, f);
ei = mean_with_nan(Ei, f);
eb = mean_with_nan(Eb, f);
ed = mean_with_nan(Ew + Es + Ec, f);
e = et + ei + eb + ed;

% Create figure for latitudinal profiles
figure1 = figure;

% Configure axes for plotting
axes1 = axes('Parent', figure1, ...
    'XTick', [], ...
    'YTick', [], ...
    'Position', [0.13 0.039 0.189 0.95]);
xlim(axes1, [50 1650]);
ylim(axes1, [0 1300]);
view(axes1, [90 90]);
box(axes1, 'on');
hold(axes1, 'on');

% Create stacked area plots for each component
area1 = area([et, eb, ei, ed], 'LineStyle', 'none', 'Parent', axes1);
set(area1(1), 'DisplayName', 'E_t', 'FaceColor', [0.4 0.9 0.6]);
set(area1(2), 'DisplayName', 'E_b', 'FaceColor', [0.8 0.6 0.5]);
set(area1(3), 'DisplayName', 'E_i', 'FaceColor', [0.8 0.9 1]);
set(area1(4), 'DisplayName', 'E_d', 'FaceColor', [0.6 0.6 0.6]);

% Plot total evaporation as a line
plot2 = plot(e, 'LineWidth', 1.2, 'Parent', axes1);
set(plot2, 'DisplayName', 'E', 'Color', [0 0 0]);

print(figure1, '-djpeg90', '-image', '-r800', '-cmyk', [data_path, 'latitudinal.jpg']);

%% Maps of components (Fig. 2)

% Map and save each component with specified color limits
map_and_save(E, [-25 1400], S == 0, [data_path, 'E.jpg']);
map_and_save(Et, [-25 1400], S == 0, [data_path, 'Et.jpg']);
map_and_save(Eb, [-25 1400], S == 0, [data_path, 'Eb.jpg']);
map_and_save(Ei, [-25 1400], S == 0, [data_path, 'Ei.jpg']);
map_and_save(Eo, [-25 1400], S == 0, [data_path, 'Eo.jpg']);

%% Seasonal maps (Fig. 3)

% Map and save seasonal evaporation and potential evaporation maps
map_and_save(E_JJA * 92, [-7 400], S == 0, [data_path, 'E_JJA.jpg']);
map_and_save(E_DJF * 90, [-7 400], S == 0, [data_path, 'E_DJF.jpg']);
map_and_save(Ep_JJA * 92, [-14 800], S == 0, [data_path, 'Ep_JJA.jpg']);
map_and_save(Ep_DJF * 90, [-14 800], S == 0, [data_path, 'Ep_DJF.jpg']);
map_and_save(S_JJA, [-0.017 1], S == 0, [data_path, 'S_JJA.jpg']);
map_and_save(S_DJF, [-0.017 1], S == 0, [data_path, 'S_DJF.jpg']);
map_and_save(SMrz_JJA, [-0.017 1], S == 0, [data_path, 'SMrz_JJA.jpg']);
map_and_save(SMrz_DJF, [-0.017 1], S == 0, [data_path, 'SMrz_DJF.jpg']);
H_JJA(H_JJA < 0) = 0; H_DJF(H_DJF < 0) = 0; % Set negative values to zero
map_and_save(H_JJA, [-2.2 120], S == 0, [data_path, 'H_JJA.jpg']);
map_and_save(H_DJF, [-2.2 120], S == 0, [data_path, 'H_DJF.jpg']);

%% Drivers of evaporation (Fig. 4)

% Calculate and display area-weighted totals for evaporation drivers
load Msea; 

DEF = Ep - E;
DEF(Msea) = NaN; % Mask coastline
DEF_area = calculate_area_weighted_total(DEF, Area, S);
disp(['DEF_area: ', num2str(DEF_area)]);

AERO = Ep_aero .* E ./ (Ep_aero + Ep_rad);
AERO(Msea) = NaN; % Mask coastline
AERO_area = calculate_area_weighted_total(AERO, Area, S);
disp(['AERO_area: ', num2str(AERO_area)]);

RAD = Ep_rad .* E ./ (Ep_aero + Ep_rad);
RAD(Msea) = NaN; % Mask coastline
RAD_area = calculate_area_weighted_total(RAD, Area, S);
disp(['RAD_area: ', num2str(RAD_area)]);

% Adjust and display the composite color map
DEF = clamp(DEF * 0.8, 0, 1800);
AERO = clamp(AERO * 1.2, 0, 1800);
RAD = clamp(RAD, 0, 1800);

triangle_colormap(DEF, AERO, RAD);
print(gcf, '-djpeg90', '-image', '-r800', '-cmyk', [data_path, 'drivers.jpg']);

%% Functions

% Function to map and save figures
function map_and_save(data, clim, S, filename)
data(S == 1) = -200;
map_E(data, clim);
print(gcf, '-djpeg90', '-image', '-r800', '-cmyk', filename);
end

% Function to create maps of evaporation and components
function map_E(A, cax)
figure;
A = double(A);
lon_range = -179.95:0.1:179.95;
lat_range = 89.95:-0.1:-89.95;
[lon, lat] = meshgrid(lon_range, lat_range);

load diecolor_EE.mat % Load preferred colormap

% Set up the Plate Carrée (Equidistant Cylindrical) Projection
m_proj('Equidistant', 'lon', [-180, 180], 'lat', [-75, 85]);

% Plot the data using the custom colormap
m_pcolor(lon, lat, A);
shading flat;
colormap(diecolor_EE);
clim(cax);

m_coast('line', 'Color', [.3, .3, .3]);
m_grid('xaxis', 'middle', 'xticklabels', [], 'yticklabels', []);

% Add a standard colorbar
colorbar('h');
end

% Function to create composite RGB image
function triangle_colormap(A, B, C)
RGB = cat(3, A, B, C);

% Normalize RGB values
RGB = RGB - min(RGB(:));
RGB = RGB / max(RGB(:));

% Create figure with custom colormap
figure;
axes;
imagesc(RGB);
colormap('jet');
colorbar;
end

% Function to clamp values
function result = clamp(data, lower_limit, upper_limit)
result = data;
result(result < lower_limit) = lower_limit;
result(result > upper_limit) = upper_limit;
end

% Function to calculate area-weighted total
function total = calculate_area_weighted_total(data, area, mask)
data_area = (data ./ 1000000) .* area;
data_area(isnan(mask)) = NaN;
total = nansum(nansum(data_area));
end

% Function to calculate mean with NaN handling
function mean_val = mean_with_nan(data, mask)
data(mask) = NaN;
mean_val = nanmean(data, 2);
end
