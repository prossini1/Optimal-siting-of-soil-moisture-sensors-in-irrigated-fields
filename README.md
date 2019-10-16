# Optimal-siting-of-soil-moisture-sensors-in-irrigated-fields

## Develop a quantitative protocol to decide the optimal location of a limited number of soil moisture sensors in irrigated fields









# Rank stability analysis
# Determine field locations (cells of the raster map of interpolated soil
# moisture) that have values similar to that of the average of the field.

%% Set script constants

% Colormap
cmap = flipud(parula);

% Measurement point
station_lat = 37.96637;
station_lon = -97.88042;
windturbin_lat = 37.967117;
windturbin_lon = -97.882092;

%% Load soil moisture data collected with the Hydrosense
data_1 = readtable('hutchinson_HyS_05_Nov_2018.csv');
data_2 = readtable('hutchinson_HyS_16_Nov_2018.csv');
data_3 = readtable('hutchinson_HyS_5_Jun_2019.csv');

%% Load field boundary
load Hutch_boundary.mat;

%% Interpolate irregular spatial observations into a regular grid
% I used data_1 to set the boudnary because is the day I went sampling with
% you and I know we covered the field pretty well.
lat_range = max(data_1.Latitude):-0.0001:min(data_1.Latitude);
lon_range = min(data_1.Longitude):0.0001:max(data_1.Longitude);
[lon_grid,lat_grid] = meshgrid(lon_range,lat_range);
IN = inpolygon(lon_grid,lat_grid,Hutch_boundary(:,1),Hutch_boundary(:,2));

%% Plot grid
f = figure('Color','w');
axesm('mercator')
nan_grid = ones(size(lon_grid))*0.9;
nan_grid(~IN) = nan;
smap = surfm(lat_grid,lon_grid,nan_grid); hold on
smap.EdgeColor = 'k';
smap.LineStyle = '-';
smap.EdgeColor = 'k';
colormap(gray)
hold on
plotm(Hutch_boundary(:,2),Hutch_boundary(:,1),'-k','Linewidth',2) 
framem('off')
f.CurrentAxes.XTick = [];
f.CurrentAxes.XColor = [1 1 1];
f.CurrentAxes.YTick = [];
f.CurrentAxes.YColor = [1 1 1];
axis tight
axis equal

%pcolorm(lat_grid,lon_grid,nan_grid,'MeshStyle','both');


%% Create Scattered Interpolant (point to raster)
% An interpolant is a surface model used to interpolate soil moisture
% values in the generated grid between observed points.
% The interpolant is based on the Natural Neighbor interpolation method, 
% which is between kriging and inverse distance weight.
F_1 = scatteredInterpolant(data_1.Longitude,data_1.Latitude,data_1.Moisture, 'natural');
z_grid_1 = F_1(lon_grid,lat_grid);
z_grid_1(~IN) = nan;

F_2 = scatteredInterpolant(data_2.Longitude,data_2.Latitude,data_2.Moisture, 'natural');
z_grid_2 = F_2(lon_grid,lat_grid);
IN = inpolygon(lon_grid,lat_grid,Hutch_boundary(:,1),Hutch_boundary(:,2));
z_grid_2(~IN) = nan;

F_3 = scatteredInterpolant(data_3.Longitude,data_3.Latitude,data_3.Moisture, 'natural');
z_grid_3 = F_3(lon_grid,lat_grid);
IN = inpolygon(lon_grid,lat_grid,Hutch_boundary(:,1),Hutch_boundary(:,2));
z_grid_3(~IN) = nan;

%% Plot gridded map of soil moisture
figure('Color','w');
axesm('mercator')
surfm(lat_grid,lon_grid,z_grid_1); hold on
plotm(Hutch_boundary(:,2),Hutch_boundary(:,1),'-k','Linewidth',1)
scatterm(data_1.Latitude,data_1.Longitude,'+r'); % Overlap measurements
colormap(cmap)
cb = colorbar;
cb.Limits = [0 50];
set(gca,'FontSize',18)
set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[])
axis tight
cb.Label.String = 'Percent water by soil volume';
title('5 Nov 2018', 'Fontsize', 15)

figure;
axesm('mercator')
surfm(lat_grid,lon_grid,z_grid_2); hold on
plotm(Hutch_boundary(:,2),Hutch_boundary(:,1),'-r','Linewidth',1) 
colormap(cmap)

figure;
axesm('mercator')
surfm(lat_grid,lon_grid,z_grid_3); hold on
plotm(Hutch_boundary(:,2),Hutch_boundary(:,1),'-r','Linewidth',1) 
colormap(cmap)

%% Leave-one-out cross-validation for a single date
% Main used to test the accuracy of the predictions of the interpolation
% model
 crossval = table(); % store results
for i=1:size(data_1,1)
    D = data_1; % generate copy of table to avoid erasing values in the original
    lat_true = D.Latitude(i);
    lon_true = D.Longitude(i);
    vwc_true = D.Moisture(i);
    D(i,:) = []; % remove row i
    F = scatteredInterpolant(D.Longitude, D.Latitude, D.Moisture, 'natural','linear');
    vwc_pred = F(lon_true,lat_true);
    
    crossval.vwc_true(i) = vwc_true;
    crossval.vwc_pred(i) = vwc_pred;
    crossval.error(i) = (vwc_true - vwc_pred);
end

% Compute crossvalidation error
r = corr(crossval.vwc_true,crossval.vwc_pred);
rmse = sqrt(mean((crossval.vwc_true - crossval.vwc_pred).^2));
mae = mean(abs(crossval.vwc_true - crossval.vwc_pred));

% Plot crossvalidation
scatter(crossval.vwc_true,crossval.vwc_pred);
hold on
plot([0 50],[0 50]);

%% Histograms
% Unmute to plot. Histograms are useful to show the range of values, the
% most predominant values, and the distribution of the values.
%figure;
%histogram(data_1.Moisture)
avg_1 = nanmean(z_grid_1(:));

%figure;
%histogram(data_2.Moisture)
avg_2 = nanmean(z_grid_2(:));

avg_3 = nanmean(z_grid_3(:));

%% Calculate the Relative Difference (RD) for each date

RD_1 = (z_grid_1 - avg_1)./avg_1; % [Eq 1]
RD_2 = (z_grid_2 - avg_2)./avg_2;
RD_3 = (z_grid_3 - avg_3)./avg_3;

%% Plot gridded map of relative differences
figure;
axesm('mercator')
surfm(lat_grid,lon_grid,RD_1); hold on
cb = colorbar;
cb.Limits = [-0.5 0.5];


%% Mean Relative Difference (MRD)
RD = cat(3,RD_1,RD_2,RD_3); % Concatenate maps of RD
MRD = nanmean(RD,3);
SRD = nanstd(RD,[],3); % standard deviation of the RD

%% Plot MRD
figure;
axesm('mercator')
surfm(lat_grid,lon_grid,MRD); hold on
cb = colorbar;
%cb.Limits = [-0.3 0.3];
%scatterm(station_lat, station_lon,'*r')
%scatterm(windturbin_lat,windturbin_lon,'*k');hold on
colormap(cmap)


%% Plot SRD
figure;
axesm('mercator')
surfm(lat_grid,lon_grid,SRD); hold on
cb = colorbar;
%scatterm(station_lat, station_lon,'*r')
colormap(cmap)

%% Soil texture
psa_data = readtable('PSA_data.xlsx');

F_sand = scatteredInterpolant(psa_data.longitud,psa_data.latitud,psa_data.sand, 'natural');
sand_grid = F_sand(lon_grid,lat_grid);
sand_grid(~IN) = nan;

F_clay = scatteredInterpolant(psa_data.longitud,psa_data.latitud,psa_data.clay, 'natural');
clay_grid = F_clay(lon_grid,lat_grid);
clay_grid(~IN) = nan;

[~,usdaclass] = soiltextureclass(sand_grid(:),clay_grid(:));
usdaclass = reshape(usdaclass,size(clay_grid));

%% Determine number of management zones
% rng('default');  % For reproducibility
% eval_CalinskiHarabasz = evalclusters(MRD(:),'kmeans','CalinskiHarabasz','KList',[1:10]);
% eval_gap = evalclusters(MRD(:),'kmeans','gap','KList',[1:10]);

%% Classification of management zones using K-Means

rng('default');
eval_silhouette = evalclusters(MRD(:),'kmeans','silhouette','klist',[1:6]);
Ngroups = eval_silhouette.OptimalK;
X = [MRD(:)];%,nanzscore(sand_grid(:)),nanzscore(clay_grid(:))];
[zones] = kmeans(X,Ngroups,'Replicates',10,'Distance','sqeuclidean');
zones = reshape(zones,size(MRD));

%% Distance transform of binary image
% For each pixel in BW, the distance transform assigns a number that is the 
% distance between that pixel and the nearest nonzero pixel.

zone_results = table();
f = figure('Color','w');
%cmap = colormap('gray');

for i = 1:Ngroups
   
    idx_zone = false(size(lat_grid));
    idx_zone(zones ~= i) = true;
    idx_zone_pad = padarray(idx_zone,[1 1],1,'both');
    
    % Add zero padding to surround the field with zeros.
    % Examine the matrix to see how it works.
    D = bwdist(idx_zone_pad,'euclidean');
    
    % Extract matrix without padding
    D = D(2:end-1,2:end-1); 
    
    % Affect the distance by the varaibility
    % For Hutch if you remove .* (1-SRD).^2 you can see the selection of
    % the wrong spot. Still centered, but on the unstable portion of the
    % zone.
    D = D .* (1-SRD).^2;
    [~,I] = nanmax(D(:));
    
    % Store lat and lon of pixel with highest value.
    zone_results.zone(i) = i;
    zone_results.lat(i) = lat_grid(I);
    zone_results.lon(i) = lon_grid(I);
    
    % Plot eroded zone with centroid
    subplot(1,3,i)
    axesm('mercator')
    surfm(lat_grid,lon_grid,D); hold on
    plotm(Hutch_boundary(:,2),Hutch_boundary(:,1),'-k','Linewidth',2) 
    hold on
    scatterm(zone_results.lat(i),zone_results.lon(i),'+r')
    hold on
    title(['Zone ',num2str(i)])
    f.CurrentAxes.XTick = [];
    f.CurrentAxes.XColor = [1 1 1];
    f.CurrentAxes.YTick = [];
    f.CurrentAxes.YColor = [1 1 1];
    axis tight
    axis equal
    %colormap(cmap)
end


%% Plot zones based on K-means
figure;
axesm('mercator')
surfm(lat_grid,lon_grid,zones); hold on
plotm(Hutch_boundary(:,2),Hutch_boundary(:,1),'-k','Linewidth',2) 
c = summer;
c = flipud(c);
colormap(c)
hold on
scatterm(zone_results.lat,zone_results.lon,'+r')

%% Show SRD inside management zone
ZOI = 1;
idx_current_zone = zones == ZOI;
Z = SRD .* double(idx_current_zone);
Z(~idx_current_zone) = nan;
figure;
axesm('mercator')
pcolorm(lat_grid, lon_grid, Z)
hold on
scatterm(zone_results.lat(ZOI),zone_results.lon(ZOI),50,'+r')
colormap(summer)
colorbar

%% Show K-Means cluster distance relative to cluster centroid inside management zone
figure;
zone_D_cluster = idx_zone .* D_cluster;
zone_D_cluster(~idx_zone) = nan;
axesm('mercator')
pcolorm(lat_grid,lon_grid,zone_D_cluster);
hold on
scatterm(zone_results.lat(i),zone_results.lon(i),100,'+r','LineWidth',3)
colormap(summer)
colorbar

% figure;scatter(zone_SRD(:),zone_D_cluster(:))


%% Plot soil texture
textures = {'Sand',...
          'Loamy sand',...
          'Sandy loam',...
          'Loam',...
          'Silt loam',...
          'Silt',...
          'Sandy clay loam',...
          'Clay loam',...
          'Silty clay loam',...
          'Sandy clay',...
          'Silty clay',...
          'Clay'};
figure;
axesm('mercator')
pcolorm(lat_grid,lon_grid,usdaclass);
hold on
%scatterm(zone_results.lat(i),zone_results.lon(i),100,'+r','LineWidth',3)
colormap(brewermap(12,'*Paired')) % Only few colormaps have 12 unique colors
cb = colorbar('Ticks',[1.5:12.5],...
         'TickLabels',textures,...
         'FontSize',12);
caxis([1 13])

%% Histogram soil texture
histogram(sand_grid(:))

%% Plot Soil Moisture Management Zones

figure1 = figure('color','white');
axesm('mercator')
[cc,h] = contourfm(lat_grid,lon_grid,MRD,4); % It will produce an extra zone
c = parula;
c = flipud(c);
colormap(c)

hold on
%scatterm(station_lat, station_lon,'*r')
axis equal
plotm(Hutch_boundary(:,2),Hutch_boundary(:,1),'-k','Linewidth',1)
box()
axis('off')
title(upper('Soil moisture management zones')); hold on

%load the data points from Hyprop sampling 

hyprop = readtable('Sampling_data.csv');
scatterm(hyprop.lat,hyprop.lon,'*k')

scatterm(zone_results.lat,zone_results.lon,'*r')

%%
kmlwritepoint('kml_sampling',hyprop.lat,hyprop.lon);

%% Plot SMRD
figure;
axesm('mercator')
surfm(lat_grid,lon_grid,SRD); hold on
cb = colorbar;
cb.Limits = [-0.5 0.5];

%% Ranks each grid cell
[pixel_rank,pixel_idx] = sort(MRD(:),'Ascend');

%% Plot rank
figure = figure('color','w');
scatter(1:numel(pixel_rank),pixel_rank)
box('on')
hold on
ylabel('MRD');
xlabel('Pixel');
title('Time stability analysis');
%%
figure;
errorbar(1:numel(pixel_rank),pixel_rank,SRD(pixel_idx))

%% Find field area exhibiting values similar to the field average
representative_area = double(MRD <= 0.05 & MRD >= -0.05);
representative_area(~IN) = nan;

%% Plot representative area for soil moisture
figure;
axesm('mercator')
surfm(lat_grid,lon_grid,representative_area); hold on
scatterm(station_lat,station_lon,'*r');hold on
%scatterm(windturbin_lat,windturbin_lon,'*w');hold on

%%
C = cc;
Cmod = mod(C(2,:),2);
Cidx = Cmod==0 | Cmod==1;
[~,col] = find(Cidx==1);

Clabels = C(1,col)';
ClabelsUnique = unique(C(1,col)');

zones = MRD >= ClabelsUnique(3) & MRD < ClabelsUnique(4);
B = bwboundaries(zones);

% axesm('mercator');
% plotm(Hutch_boundary(:,2),Hutch_boundary(:,1));hold on
% for k = 1:length(B)
%    boundary = B{k};
%    plotm(lat_grid(boundary(:,1),boundary(:,2)), lon_grid(boundary(:,1),boundary(:,2)), 'k', 'LineWidth', 2);hold on
% end

%% Select cordinates from each zone created in the contour map
 
zones = MRD >= 0.054875 & MRD <= 0.17707;
lat_idx = lat_grid(zones);
lon_idx = lon_grid(zones);

figure;
axesm('mercator');
plotm(Hutch_boundary(:,2),Hutch_boundary(:,1));hold on
scatterm(lat_idx,lon_idx)

%Create a kml file
%border = boundary(lon_idx,lat_idx,0.9);
kmlwritepoint('asampling_wet',lat_idx,lon_idx); 



