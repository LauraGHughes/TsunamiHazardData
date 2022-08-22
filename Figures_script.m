% Maps showing location of earthquakes:
% Load DEM data:
[bathymetry_Z,bathymetry_R] = readgeoraster('F:\PhD\DEM\nzsm_llbn_10sec_20140220.asc');
latlim = bathymetry_R.YWorldLimits; lat_points = [latlim(2):-0.0027778:latlim(1)];
lonlim = bathymetry_R.XWorldLimits; lon_points = [lonlim(1):0.0027778:lonlim(2)];
Z = 1*(bathymetry_Z>0); Z = Z./Z;
% New Zealand coastline
coastline = 'C:\Users\hughesla\OneDrive - Victoria University of Wellington - STAFF\Documents\Coastlines\MATLAB_coast.mat';
load(coastline);

event_data2 = csvread('C:\Users\hughesla\OneDrive - Victoria University of Wellington - STAFF\Documents\deformation_models\Bruce_NSHM_tsunami_energy.csv', 1,0);
nzsm_data = zeros(2595, 5);
for ii = 1:2595
    disp(ii)
    event_x = event_data2(ii,5); event_y = event_data2(ii,6); 
    event_z = event_data2(ii,7);
    [event_lon, event_lat] = UTMzoneX2Geo(event_x, event_y, 59);
    event_info = [event_data2(ii,1), event_lon, event_lat, event_z, event_data2(ii,4)];
    nzsm_data(ii, [1,2,3,4, 5]) = event_info;
end

pac = load('F:\PhD\figures_VUWcomputer\DataSets\pacific.mat'); lat = pac.pac(:,1); lon = pac.pac(:,2);

figure
tlo = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
nexttile(1, [2,1])
s = pcolor(lon_points,lat_points, Z); shading flat; s.FaceAlpha = 0.1; colormap parula; caxis([5 20]);
set(gca, 'FontSize',11);
hold on
plot(x,y, 'LineWidth', 1, 'Color', 'black');
ylim([-50 -30]); xlim([163 185]); text(159.5,-30,'a', 'FontSize', 20, 'FontWeight', 'bold');
plot(lon, lat, 'LineWidth', 2, 'Color', 'black')
t1 = text(178.2,-42,'Hikurangi'); set(t1,'Rotation',70);
t2 = text(163.5,-48,'Puysegur'); set(t2,'Rotation',70);
t3 = text(180.9, -34,'Kermadec');set(t3,'Rotation',70);
t4 = text(165,-32,'AUS');t5 = text(180,-48,'PAC');
t6 = text(175.5,-39,'NI');t7 = text(169,-45,'SI');
t8 = text(174.2,-41.7,'CS'); t9 = text(168.9, -43.4,'AF');set(t9,'Rotation',45);
hold off  

nexttile(5)
b = histogram(event_data2(:,4));
xlabel('Magnitude', 'FontSize',10); ylabel('Frequency', 'FontSize',10);
set(gca, 'YScale', 'log');  %yticklabels({'1','10','100', '1000'}); 
xlim([7 9.3]); ylim([1 1000]);
set(gca, 'FontSize',10)
b.FaceColor = '#000000'; b.EdgeColor = '#000000';hold on; text(6.65,1000,'b', 'FontSize', 20, 'FontWeight', 'bold')
text(8.8,600,'n = 2595', 'FontSize', 11, 'FontWeight', 'bold')


nexttile(2, [3,1]) % Dataset two
ss = pcolor(lon_points,lat_points, Z); shading flat; ss.FaceAlpha = 0.1; 
ylim([-50 -23]); xlim([163 187]); set(gca, 'FontSize',11);
hold on;
plot(x,y, 'LineWidth', 1, 'Color', 'black');
scatter(nzsm_data(:,2),nzsm_data(:,3),(10.^nzsm_data(:,5))/1000000,(nzsm_data(:,4)*-1)/1000,'filled','MarkerFaceAlpha',0.5)
c = colorbar('eastoutside', 'FontSize', 11);  colormap parula; caxis([5 20])
ylabel(c,'Depth (km)','FontSize',11,'Rotation',270);
set(gca, 'FontSize',11); text(161,-23,'c', 'FontSize', 20, 'FontWeight', 'bold')
%scatter(165.8,-29.8, 10^7/1000000, 'blue','filled','MarkerFaceAlpha',0.5); text(167.3,-29.9,'Mw7.0','FontSize',9);
scatter(165.8,-28, 10^7.5/1000000, 'blue','filled','MarkerFaceAlpha',0.5); text(167.7,-28.1,'Mw7.5','FontSize',10.5);
scatter(165.8,-27.5, 10^8/1000000, 'blue','filled','MarkerFaceAlpha',0.5); text(167.7,-27.5,'Mw8.0','FontSize',10.5);
scatter(165.8,-26.6, 10^8.5/1000000, 'blue','filled','MarkerFaceAlpha',0.5); text(167.7,-26.6,'Mw8.5','FontSize',10.5);
scatter(165.8,-25, 10^9/1000000, 'blue','filled','MarkerFaceAlpha',0.5); text(167.6,-25,'Mw9.0','FontSize',10.5);
hold off

%%-----------------------------------------------------------------------%%
% Maximum event - event 223533
mycolormap = customcolormap(linspace(0,1,11), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3','#d5e2f0','#93c5dc','#4295c1','#2265ad','#062e61'});
% NZ coast
coastline = 'C:\Users\hughesla\OneDrive - Victoria University of Wellington - STAFF\Documents\Coastlines\MATLAB_coast.mat';
load(coastline);

% Deformation information
event_path = 'F:\PhD\bruce2_NSHM\deformation_models\'; event = 'ev223533';
eventX_NZTM = ncread([event_path, event, '.grd'],'x'); eventY_NZTM = ncread([event_path, event, '.grd'],'y'); 
eventZ_NZTM = ncread([event_path, event, '.grd'],'z'); 
[NZTM_NN, NZTM_EE] = meshgrid(eventY_NZTM, eventX_NZTM); [Longitude, Latitude] = NZTM2Geo(NZTM_EE, NZTM_NN); 

%Grid 1 information
filename = ['F:\PhD\bruce2_NSHM\tsunamis\zmax_files\event223533\zmax_layer01.dat'];
layer01_max = importdata(filename); A1 = permute(layer01_max, [2 1]);
B1 = reshape(A1, 765, 826); B1(any(isnan(B1), 2), :) = []; maximum_grid1 = B1.';
filename = ['F:\PhD\bruce2_NSHM\tsunamis\zmax_files\event223533\layer01_x.dat']; layer01_x = importdata(filename);
filename = ['F:\PhD\bruce2_NSHM\tsunamis\zmax_files\event223533\layer01_y.dat']; layer01_y = importdata(filename);

%Grid 3 information
filename = ['F:\PhD\bruce2_NSHM\tsunamis\zmax_files\event223533\zmax_layer03.dat'];
layer03_max = importdata(filename); A1 = permute(layer03_max, [2 1]);
B1 = reshape(A1, 1560, 2238); B1(any(isnan(B1), 2), :) = []; maximum_grid3 = B1.';
filename = ['F:\PhD\bruce2_NSHM\tsunamis\zmax_files\event223533\layer03_x.dat']; layer03_x = importdata(filename);
filename = ['F:\PhD\bruce2_NSHM\tsunamis\zmax_files\event223533\layer03_y.dat']; layer03_y = importdata(filename);


figure
tlo = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
ax1 = nexttile; % Deformation model
pcolor(Longitude, Latitude,  eventZ_NZTM); colormap(ax1, mycolormap); caxis([max(max(eventZ_NZTM))*-1 max(max(eventZ_NZTM))]);set(gca, 'FontSize',10);
shading flat; a = colorbar('eastoutside', 'FontSize', 10); a.Label.String = 'Vertical displacement (m)';
a.Label.FontSize = 11.5;
hold on
plot(NIx,NIy, 'LineWidth', 1, 'Color', 'black'); plot(SIx,SIy, 'LineWidth', 1, 'Color', 'black');
xlim([170 186]); ylim([-42 -24]);text(167.5,-24,'a', 'FontSize', 20, 'FontWeight', 'bold');
hold off;

ax2 = nexttile; % Grid one
pcolor(layer01_x, layer01_y, maximum_grid1); shading flat; colormap(ax2, parula); 
set(gca,'ColorScale','log', 'CLim', [0.1 20]);xlim([160 190]); ylim([-50 -25]);set(gca, 'FontSize',10);
b = colorbar('FontSize', 10,'TickLabels', {'0.1','1','10'}); set(gca,'ColorScale','log', 'CLim', [0.1 20]);
b.Label.String = 'Wave height (m)'; b.Label.FontSize = 11.5; 
text(155,-25,'b', 'FontSize', 20, 'FontWeight', 'bold');
pos = get(b,'Position'); b.Label.Position = [2 pos(2)+1.2]; 

ax3 = nexttile; % Grid three
pcolor(layer03_x, layer03_y, maximum_grid3); shading flat; colormap(ax3, parula); 
set(gca,'ColorScale','log', 'CLim', [0.1 20]);xlim([166 179]); ylim([-48 -34]); set(gca, 'FontSize',10);
c = colorbar('FontSize', 10,'TickLabels', {'0.1','1','10'}); set(gca,'ColorScale','log', 'CLim', [0.1 20]);
c.Label.String = 'Wave height (m)'; c.Label.FontSize = 11.5; 
text(164,-34,'c', 'FontSize', 20, 'FontWeight', 'bold');
pos = get(c,'Position'); c.Label.Position = [2 pos(2)+1.2]; 

%%-----------------------------------------------------------------------%%
% Maximum wave height at the coast
bruce2_path = ['C:\Users\hughesla\OneDrive - Victoria University of Wellington '...
    '- STAFF\Documents\deformation_models\'...
    'Bruce_NSHM_tsunami_energy.csv'];
bruce2_tsunami_a = load('F:\PhD\bruce2_NSHM\tsunamis\max_height_locations.mat');
bruce2_tsunami_b = load('F:\PhD\bruce2_NSHM\tsunamis\max_height_locations2.mat');
bruce2_tsunami_c = load('F:\PhD\bruce2_NSHM\tsunamis\max_height_locations3.mat');

bruce2 = [bruce2_tsunami_a.max_height_locations; bruce2_tsunami_b.max_height_locations2; bruce2_tsunami_c.max_height_locations3];
bruce2 = sortrows(bruce2, 3);

[bathymetry_Z,bathymetry_R] = readgeoraster('F:\PhD\DEM\nzsm_llbn_10sec_20140220.asc');
latlim = bathymetry_R.YWorldLimits; lat_points = [latlim(2):-0.0027778:latlim(1)];
lonlim = bathymetry_R.XWorldLimits; lon_points = [lonlim(1):0.0027778:lonlim(2)];
Z = 1*(bathymetry_Z>0); Z = (Z./Z)-100;
% New Zealand coastline
coastline = 'C:\Users\hughesla\OneDrive - Victoria University of Wellington - STAFF\Documents\Coastlines\MATLAB_coast.mat';
load(coastline);

figure
tlo = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
nexttile(1, [2,1]) % Catalogue two
b = histogram(bruce2(:,3));
xlabel('Maximum wave height (m)', 'FontSize',10); ylabel('Number of events', 'FontSize',10);%title('DS2 maximum wave height distribution')
set(gca, 'YScale', 'log'); ylim([0.9 1000]); xlim([0 28]);
%yticklabels({'1','10','100', '1000'}); 
set(gca, 'FontSize',10)
b.FaceColor = '#000000'; b.EdgeColor = '#000000';hold on; text(-4,1000,'a', 'FontSize', 20, 'FontWeight', 'bold')
text(22,750,'n = 2595', 'FontSize', 11, 'FontWeight', 'bold')

nexttile(2,[3 1])
% Catalogue two
s = pcolor(lon_points,lat_points, Z); shading flat; s.FaceAlpha = 0.05;hold on
plot(x,y, 'LineWidth', 1, 'Color', 'black')
set(gca, 'FontSize',10);xlim([166 179]); ylim([-48 -34]);
scatter(bruce2(:,1), bruce2(:,2),50, bruce2(:,3),'filled', 'MarkerFaceAlpha', 0.5)
colormap parula;
a = colorbar('eastoutside','TickLabels', {'0.1','1','10'}, 'FontSize', 10); a.Label.String = 'Wave height (m)';
a.Label.FontSize = 11.5; set(gca,'ColorScale','log', 'CLim', [0.1 30]); %title('DS2 maximum wave height coastal location')
pos = get(a,'Position');
a.Label.Position = [2 pos(2)+2]; 
text(164,-34,'b', 'FontSize', 20, 'FontWeight', 'bold')
hold off


points_above5 = sum(bruce2(:,3)>5); percent_above5 = points_above5/2595*100;
points_above10 = sum(bruce2(:,3)>10); percent_above10 = points_above10/2595*100;
points_above15 = sum(bruce2(:,3)>15); percent_above15 = points_above15/2595*100;
points_above1 = sum(bruce2(:,3)>1); percent_above1 = points_above1/2595*100;

load('C:\Users\hughesla\OneDrive - Victoria University of Wellington - STAFF\Documents\Coastlines\MATLAB_coast.mat');

EC = polyshape([174.79 -41.61; 178.37 -37.44; 179.61 -37.85; 175.57 -41.96]);
US = polyshape([172 -43.5; 173.8 -41; 174.71 -41.78; 174 -43.5]);
LS = polyshape([165.47 -45.95; 168.41 -43; 169.43 -44; 167.09 -47]);


figure
plot(x,y); hold on;
plot(EC); plot(US); plot(LS); scatter(bruce2(:, 1), bruce2(:, 2))
hold off

bruce2_tsunami_a = load('F:\PhD\bruce2_NSHM\tsunamis\max_height_locations.mat');
bruce2_tsunami_b = load('F:\PhD\bruce2_NSHM\tsunamis\max_height_locations2.mat');
bruce2_tsunami_c = load('F:\PhD\bruce2_NSHM\tsunamis\max_height_locations3.mat');

bruce2 = [bruce2_tsunami_a.max_height_locations; bruce2_tsunami_b.max_height_locations2; bruce2_tsunami_c.max_height_locations3];
bruce2 = sortrows(bruce2, 3);

EC_e = zeros(length(bruce2), 3); US_e = zeros(length(bruce2), 3); LS_e = zeros(length(bruce2), 3); 


for e = 1:length(bruce2)
    if isinterior(EC, bruce2(e,1), bruce2(e,2)) == 1
        e_data = [bruce2(e,1), bruce2(e,2), bruce2(e,3)];
        EC_e(e, [1,2,3]) = e_data;

    elseif isinterior(US, bruce2(e,1), bruce2(e,2)) == 1
        e_data = [bruce2(e,1), bruce2(e,2), bruce2(e,3)];
        US_e(e, [1,2,3]) = e_data;

    elseif isinterior(LS, bruce2(e,1), bruce2(e,2)) == 1
        e_data = [bruce2(e,1), bruce2(e,2), bruce2(e,3)];
        LS_e(e, [1,2,3]) = e_data;        
    end
end

EC_events = EC_e(any(EC_e,2),:); US_events = US_e(any(US_e,2),:); LS_events = LS_e(any(LS_e,2),:);
figure
plot(x,y); hold on;
plot(EC); plot(US); plot(LS); 
scatter(EC_events(:, 1), EC_events(:, 2));scatter(US_events(:, 1), US_events(:, 2));scatter(LS_events(:, 1), LS_events(:, 2));
hold off


%%-----------------------------------------------------------------------%%
% Spatial distribution of maximum wave heights:
bruce2_tsunami_a = load('F:\PhD\bruce2_NSHM\tsunamis\wave_height_data.mat');
bruce2_tsunami_b = load('F:\PhD\bruce2_NSHM\tsunamis\wave_height_data2.mat');
bruce2_tsunami_c = load('F:\PhD\bruce2_NSHM\tsunamis\wave_height_data3.mat');
bruce2_waveheight = [bruce2_tsunami_a.wave_height_data bruce2_tsunami_b.wave_height_data2 bruce2_tsunami_c.wave_height_data];

coastal_length = zeros(2595, 4);
for ii = 1:2595
    coast = bruce2_waveheight(:, ii);
    
    points_above1 = sum(coast>1); length_above1 = points_above1*0.6;
    points_above5 = sum(coast>5); length_above5 = points_above5*0.6;
    points_above10 = sum(coast>10); length_above10 = points_above10*0.6;
    points_above15 = sum(coast>15); length_above15 = points_above15*0.6;
    
    coastal_length(ii, :) = [length_above1, length_above5, length_above10, length_above15];
end

figure
tlo = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;
a = histogram(coastal_length(:,1), 15);
xlabel('Length of coastline (km)', 'FontSize',10); ylabel('Number of events', 'FontSize',10);title('Wave height greater than 1m')
set(gca, 'YScale', 'log'); ylim([0.9 3000]); xlim([0 max(coastal_length(:,1))]);
set(gca, 'FontSize',10)
a.FaceColor = '#1BABDE'; a.EdgeColor = '#000000';
nexttile;
b = histogram(coastal_length(:,2), 15);
xlabel('Length of coastline (km)', 'FontSize',10); ylabel('Number of events', 'FontSize',10);title('Wave height greater than 5m')
set(gca, 'YScale', 'log'); ylim([0.9 3000]); xlim([0 max(coastal_length(:,2))]);
 set(gca, 'FontSize',10)
b.FaceColor = '#93CA4B'; b.EdgeColor = '#000000';
nexttile;
c = histogram(coastal_length(:,3), 15);
xlabel('Length of coastline (km)', 'FontSize',10); ylabel('Number of events', 'FontSize',10);title('Wave height greater than 10m')
set(gca, 'YScale', 'log'); ylim([0.9 3000]); xlim([0 max(coastal_length(:,3))]);
set(gca, 'FontSize',10)
c.FaceColor = '#EFBA35'; c.EdgeColor = '#000000';
nexttile;
d = histogram(coastal_length(:,4), 15);
xlabel('Length of coastline (km)', 'FontSize',10); ylabel('Number of events', 'FontSize',10);title('Wave height greater than 15m')
set(gca, 'YScale', 'log'); ylim([0.9 3000]); xlim([0 max(coastal_length(:,4))]);
set(gca, 'FontSize',10)
d.FaceColor = '#F7DD2A'; d.EdgeColor = '#000000';


points_above500 = sum(coastal_length(:,4)>50); 
med_above1 = median(coastal_length(:, 1));med_above5 = median(coastal_length(:, 2));
med_above10 = median(coastal_length(:, 3)); med_above15 = median(coastal_length(:, 4));

max_above1 = max(coastal_length(:, 1));max_above5 = max(coastal_length(:, 2));
max_above10 = max(coastal_length(:, 3)); max_above15 = max(coastal_length(:, 4));

q_above1 = quantile(coastal_length(:, 1), 0.75);q_above5 = quantile(coastal_length(:, 2), 0.75);
q_above10 = quantile(coastal_length(:, 3), 0.75); q_above15 = quantile(coastal_length(:, 4), 0.75);


event_data2 = csvread('C:\Users\hughesla\OneDrive - Victoria University of Wellington - STAFF\Documents\deformation_models\Bruce_NSHM_tsunami_energy.csv', 1,0);
event_magnitude = event_data2(1:2595, 4);

events_less75 = []; events_less8 = []; events_less85 = []; events_less9 = []; events_less95 = [];

for ii = 1:length(event_magnitude)
   mag = event_magnitude(ii);
   
   if mag < 7.5
       event_coastlength = coastal_length(ii, :); events_less75 = [events_less75,event_coastlength];
   elseif (7.5<=mag) && (mag<8.0)
       event_coastlength = coastal_length(ii, :); events_less8 = [events_less8, event_coastlength];     
   elseif (8.0<=mag) && (mag<8.5)
       event_coastlength = coastal_length(ii, :); events_less85 = [events_less85, event_coastlength];
   elseif (8.5<=mag) && (mag<9.0)
       event_coastlength = coastal_length(ii, :); events_less9 = [events_less9, event_coastlength];    
   elseif mag > 9.0
       event_coastlength = coastal_length(ii, :); events_less95 = [events_less95, event_coastlength];
   end
end


events_less75s = reshape(events_less75, [4,1510]); events_less8s = reshape(events_less8, [4, 767]);
events_less85s = reshape(events_less85, [4, 103]); events_less9s = reshape(events_less9, [4, 181]);
events_less95s = reshape(events_less95, [4, 34]);

med_above1_75 = median(events_less75s(1, :)); med_above5_75 = median(events_less75s(2, :));
med_above10_75 = median(events_less75s(3, :)); med_above15_75 = median(events_less75s(4, :));

med_above1_8 = median(events_less8s(1, :)); med_above5_8 = median(events_less8s(2, :));
med_above10_8 = median(events_less8s(3, :)); med_above15_8 = median(events_less8s(4, :));

med_above1_85 = median(events_less85s(1, :)); med_above5_85 = median(events_less85s(2, :));
med_above10_85 = median(events_less85s(3, :)); med_above15_85 = median(events_less85s(4, :));

med_above1_9 = median(events_less9s(1, :)); med_above5_9 = median(events_less9s(2, :));
med_above10_9 = median(events_less9s(3, :)); med_above15_9 = median(events_less9s(4, :));

med_above1_95 = median(events_less95s(1, :)); med_above5_95 = median(events_less95s(2, :));
med_above10_95 = median(events_less95s(3, :)); med_above15_95 = median(events_less95s(4, :));


max_above1_75 = max(events_less75s(1, :)); max_above5_75 = max(events_less75s(2, :));
max_above10_75 = max(events_less75s(3, :)); max_above15_75 = max(events_less75s(4, :));

max_above1_8 = max(events_less8s(1, :)); max_above5_8 = max(events_less8s(2, :));
max_above10_8 = max(events_less8s(3, :)); max_above15_8 = max(events_less8s(4, :));

max_above1_85 = max(events_less85s(1, :)); max_above5_85 = max(events_less85s(2, :));
max_above10_85 = max(events_less85s(3, :)); max_above15_85 = max(events_less85s(4, :));

max_above1_9 = max(events_less9s(1, :)); max_above5_9 = max(events_less9s(2, :));
max_above10_9 = max(events_less9s(3, :)); max_above15_9 = max(events_less9s(4, :));

max_above1_95 = max(events_less95s(1, :)); max_above5_95 = max(events_less95s(2, :));
max_above10_95 = max(events_less95s(3, :)); max_above15_95 = max(events_less95s(4, :));


q_above1_75 = quantile(events_less75s(1, :), 0.75); q_above5_75 = quantile(events_less75s(2, :), 0.75);
q_above10_75 = quantile(events_less75s(3, :), 0.75); q_above15_75 = quantile(events_less75s(4, :), 0.75);

q_above1_8 = quantile(events_less8s(1, :), 0.75); q_above5_8 = quantile(events_less8s(2, :), 0.75);
q_above10_8 = quantile(events_less8s(3, :), 0.75); q_above15_8 = quantile(events_less8s(4, :), 0.75);

q_above1_85 = quantile(events_less85s(1, :), 0.75); q_above5_85 = quantile(events_less85s(2, :), 0.75);
q_above10_85 = quantile(events_less85s(3, :), 0.75); q_above15_85 = quantile(events_less85s(4, :), 0.75);

q_above1_9 = quantile(events_less9s(1, :), 0.75); q_above5_9 = quantile(events_less9s(2, :), 0.75);
q_above10_9 = quantile(events_less9s(3, :), 0.75); q_above15_9 = quantile(events_less9s(4, :), 0.75);

q_above1_95 = quantile(events_less95s(1, :), 0.75); q_above5_95 = quantile(events_less95s(2, :), 0.75);
q_above10_95 = quantile(events_less95s(3, :), 0.75); q_above15_95 = quantile(events_less95s(4, :), 0.75);



event_data2 = csvread('F:\PhD\bruce2_NSHM\tsunamis\tsunami_events.csv', 1,0);
event_type = event_data2(:, 10);
events_k = []; events_h = []; events_p = []; events_c = []; 

for ii = 1:length(event_type)
   type = event_type(ii);
   
   if type == 1
       event_coastlength = coastal_length(ii, :); events_k = [events_k, event_coastlength];
   elseif type == 2
       event_coastlength = coastal_length(ii, :); events_h = [events_h, event_coastlength];     
   elseif type == 3
       event_coastlength = coastal_length(ii, :); events_p = [events_p, event_coastlength];
   elseif type == 4
       event_coastlength = coastal_length(ii, :); events_c = [events_c, event_coastlength];    
   end
end

events_ks = reshape(events_k, [4,370]); events_hs = reshape(events_h, [4, 125]);
events_ps = reshape(events_p, [4, 362]); events_cs = reshape(events_c, [4, 1738]);

med_above1_k = median(events_ks(1, :)); med_above5_k = median(events_ks(2, :));
med_above10_k = median(events_ks(3, :)); med_above15_k = median(events_ks(4, :));

med_above1_h = median(events_hs(1, :)); med_above5_h = median(events_hs(2, :));
med_above10_h = median(events_hs(3, :)); med_above15_h = median(events_hs(4, :));

med_above1_p = median(events_ps(1, :)); med_above5_p = median(events_ps(2, :));
med_above10_p = median(events_ps(3, :)); med_above15_p = median(events_ps(4, :));

med_above1_c = median(events_cs(1, :)); med_above5_c = median(events_cs(2, :));
med_above10_c = median(events_cs(3, :)); med_above15_c = median(events_cs(4, :));

max_above1_k = max(events_ks(1, :)); max_above5_k = max(events_ks(2, :));
max_above10_k = max(events_ks(3, :)); max_above15_k = max(events_ks(4, :));

max_above1_h = max(events_hs(1, :)); max_above5_h = max(events_hs(2, :));
max_above10_h = max(events_hs(3, :)); max_above15_h = max(events_hs(4, :));

max_above1_p = max(events_ps(1, :)); max_above5_p = max(events_ps(2, :));
max_above10_p = max(events_ps(3, :)); max_above15_p = max(events_ps(4, :));

max_above1_c = max(events_cs(1, :)); max_above5_c = max(events_cs(2, :));
max_above10_c = max(events_cs(3, :)); max_above15_c = max(events_cs(4, :));

q_above1_k = quantile(events_ks(1, :), 0.75); q_above5_k = quantile(events_ks(2, :), 0.75);
q_above10_k = quantile(events_ks(3, :), 0.75); q_above15_k = quantile(events_ks(4, :), 0.75);

q_above1_h = quantile(events_hs(1, :), 0.75); q_above5_h = quantile(events_hs(2, :), 0.75);
q_above10_h = quantile(events_hs(3, :), 0.75); q_above15_h = quantile(events_hs(4, :), 0.75);

q_above1_p = quantile(events_ps(1, :), 0.75); q_above5_p = quantile(events_ps(2, :), 0.75);
q_above10_p = quantile(events_ps(3, :), 0.75); q_above15_p = quantile(events_ps(4, :), 0.75);

q_above1_c = quantile(events_cs(1, :), 0.75); q_above5_c = quantile(events_cs(2, :), 0.75);
q_above10_c = quantile(events_cs(3, :), 0.75); q_above15_c = quantile(events_cs(4, :), 0.75);
%%-----------------------------------------------------------------------%%
% Return periods dataset one:
[bathymetry_Z,bathymetry_R] = readgeoraster('F:\PhD\DEM\nzsm_llbn_10sec_20140220.asc');
latlim = bathymetry_R.YWorldLimits; lat_points = [latlim(2):-0.0027778:latlim(1)];
lonlim = bathymetry_R.XWorldLimits; lon_points = [lonlim(1):0.0027778:lonlim(2)];
Z = 1*(bathymetry_Z>0); Z = (Z./Z)*-999;

% Return periods dataset two:
% Load the data:
bruce2_tsunami_a = load('F:\PhD\bruce2_NSHM\tsunamis\wave_height_data.mat');
bruce2_tsunami_b = load('F:\PhD\bruce2_NSHM\tsunamis\wave_height_data2.mat');
bruce2_tsunami_c = load('F:\PhD\bruce2_NSHM\tsunamis\wave_height_data3.mat');
bruce2_waveheight = [bruce2_tsunami_a.wave_height_data bruce2_tsunami_b.wave_height_data2 bruce2_tsunami_c.wave_height_data];

% 2500year return period 
t = maxk(bruce2_waveheight, 12, 2); min_exceed = mink(t, 1,2); rp2_2500 = round(min_exceed, 2); 
% 1000year return period 
t = maxk(bruce2_waveheight, 30, 2); min_exceed = mink(t, 1,2); rp2_1000 = round(min_exceed, 2); 
% 500year return period 
t = maxk(bruce2_waveheight, 60, 2); min_exceed = mink(t, 1,2); rp2_500 = round(min_exceed, 2); 
% 100year return period 
t = maxk(bruce2_waveheight, 300, 2); min_exceed = mink(t, 1,2); rp2_100 = round(min_exceed, 2); 


coastline = 'C:\Users\hughesla\OneDrive - Victoria University of Wellington - STAFF\Documents\Coastlines\MATLAB_coast.mat';
load(coastline);
coastal_points = load('C:\Users\hughesla\OneDrive - Victoria University of Wellington - STAFF\Documents\hazard_statistics\coastal_points.mat');
coastal_i = coastal_points.coastal_points(:,3); coastal_j = coastal_points.coastal_points(:,4);
coastal_x = coastal_points.coastal_points(:,1); coastal_y = coastal_points.coastal_points(:,2);

figure
tlo = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
% Dataset two
nexttile; % 2500year return period
s = pcolor(lon_points,lat_points, Z); shading flat; s.FaceAlpha = 0.05; hold on;
plot(x,y, 'LineWidth', 1, 'Color', 'black'); hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, rp2_2500)
colormap(gca,'parula'); set(gca,'ColorScale','log', 'CLim', [0.1 10]); 
set(gca, 'FontSize', 10); 
title('2500 year return period'); 

nexttile;% 1000year return period
s = pcolor(lon_points,lat_points, Z); shading flat; s.FaceAlpha = 0.05; hold on;
plot(x,y, 'LineWidth', 1, 'Color', 'black'); hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, rp2_1000)
colormap(gca,'parula'); set(gca,'ColorScale','log', 'CLim', [0.1 10]);
set(gca, 'FontSize', 10); title('1000 year return period')
c = colorbar('FontSize', 10,'TickLabels', {'0.1','1','10'}); set(gca,'ColorScale','log', 'CLim', [0.1 10]);
c.Label.String = 'Wave height (m)'; pos = get(c,'Position'); c.Label.Position = [2 pos(2)+0.6]; 

nexttile; % 500year return period
s = pcolor(lon_points,lat_points, Z); shading flat; s.FaceAlpha = 0.05; hold on;
plot(x,y, 'LineWidth', 1, 'Color', 'black'); hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, rp2_500)
colormap(gca,'parula'); set(gca,'ColorScale','log', 'CLim', [0.1 10]);
set(gca, 'FontSize', 10); title('500 year return period')

nexttile; % 100year return period
s = pcolor(lon_points,lat_points, Z); shading flat; s.FaceAlpha = 0.05; hold on;
plot(x,y, 'LineWidth', 1, 'Color', 'black'); hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, rp2_100)
colormap(gca,'parula'); c = colorbar('FontSize', 10,'TickLabels', {'0.1','1','10'}); set(gca,'ColorScale','log', 'CLim', [0.1 10]);
c.Label.String = 'Wave height (m)'; pos = get(c,'Position'); c.Label.Position = [2 pos(2)+1]; 
set(gca, 'FontSize', 10); title('100 year return period')


    
points_above_point5_rp100 = sum(rp2_100>=0.5); percent_above_point5_rp100 = points_above_point5_rp100/length(rp2_100);
points_above_1_rp500 = sum(rp2_500>=1); percent_points_above_1_rp500 = points_above_1_rp500/length(rp2_500);
points_above_2_rp1000 = sum(rp2_1000>=2); percent_points_above_2_rp1000 = points_above_2_rp1000/length(rp2_1000);
points_above_5_rp2500 = sum(rp2_2500>=5); percent_points_above_5_rp2500 = points_above_5_rp2500/length(rp2_2500);

NI_loc = polyshape([174.47 -41.33; 174.75 -40.27; 173.13 -39.46; 171.98 -33.19; 180.59 -36.35; 176.15 -42.55]);
%figure
%plot(x,y); hold on; plot(NI_loc); hold off

NI = zeros(length(rp2_2500), 3); 
for e = 1:length(rp2_2500)
    if isinterior(NI_loc, coastal_x(e),coastal_y(e)) == 1
        e_data = [coastal_x(e), coastal_y(e), rp2_2500(e)];
        NI(e, [1,2,3]) = e_data; 
    end
end
NI_events = NI(any(NI,2),:);
points_above5_rp2500 = sum(NI_events(:,3)>=5); percent_points_above5_rp2500NI = points_above5_rp2500/length(NI_events(:,3));

NI = zeros(length(rp2_1000), 3); 
for e = 1:length(rp2_1000)
    if isinterior(NI_loc, coastal_x(e),coastal_y(e)) == 1
        e_data = [coastal_x(e), coastal_y(e), rp2_1000(e)];
        NI(e, [1,2,3]) = e_data; 
    end
end
NI_events = NI(any(NI,2),:);
points_above2_rp1000 = sum(NI_events(:,3)>=2); percent_points_above2_rp1000NI = points_above2_rp1000/length(NI_events(:,3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
event_type = event_data2(:, 10);
events_k = []; events_h = []; events_p = []; events_c = []; 

for ii = 1:length(event_type)
   type = event_type(ii);
   
   if type == 1
       event_coastlength = bruce2_waveheight(:, ii); events_k = [events_k, event_coastlength];
   elseif type == 2
       event_coastlength = bruce2_waveheight(:, ii); events_h = [events_h, event_coastlength];     
   elseif type == 3
       event_coastlength = bruce2_waveheight(:, ii); events_p = [events_p, event_coastlength];
   elseif type == 4
       event_coastlength = bruce2_waveheight(:, ii); events_c = [events_c, event_coastlength];    
   end
   
end


t = maxk(bruce2_waveheight, 12, 2); min_exceed = mink(t, 1,2); all_2500 = round(min_exceed, 2); 
t = maxk(events_k, 12, 2); min_exceed = mink(t, 1,2); k_2500 = round(min_exceed, 2); 
t = maxk(events_h, 12, 2); min_exceed = mink(t, 1,2); h_2500 = round(min_exceed, 2); 
t = maxk(events_p, 12, 2); min_exceed = mink(t, 1,2); p_2500 = round(min_exceed, 2); 


coastline = 'C:\Users\hughesla\OneDrive - Victoria University of Wellington - STAFF\Documents\Coastlines\MATLAB_coast.mat';
load(coastline);
coastal_points = load('C:\Users\hughesla\OneDrive - Victoria University of Wellington - STAFF\Documents\hazard_statistics\coastal_points.mat');
coastal_i = coastal_points.coastal_points(:,3); coastal_j = coastal_points.coastal_points(:,4);
coastal_x = coastal_points.coastal_points(:,1); coastal_y = coastal_points.coastal_points(:,2);

[bathymetry_Z,bathymetry_R] = readgeoraster('F:\PhD\DEM\nzsm_llbn_10sec_20140220.asc');
latlim = bathymetry_R.YWorldLimits; lat_points = [latlim(2):-0.0027778:latlim(1)];
lonlim = bathymetry_R.XWorldLimits; lon_points = [lonlim(1):0.0027778:lonlim(2)];
Z = 1*(bathymetry_Z>0); Z = (Z./Z)*-999;

figure
tlo = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
% Dataset two
nexttile; % 2500year return period
s = pcolor(lon_points,lat_points, Z); shading flat; s.FaceAlpha = 0.05; hold on;
plot(x,y, 'LineWidth', 1, 'Color', 'black'); hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, all_2500)
colormap(gca,'parula'); set(gca,'ColorScale','log', 'CLim', [0.1 10]); 
set(gca, 'FontSize', 10); 
title('Entire catalogue'); 

nexttile;% 1000year return period
s = pcolor(lon_points,lat_points, Z); shading flat; s.FaceAlpha = 0.05; hold on;
plot(x,y, 'LineWidth', 1, 'Color', 'black'); hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, k_2500)
colormap(gca,'parula'); set(gca,'ColorScale','log', 'CLim', [0.1 10]);
set(gca, 'FontSize', 10); title('Tonga-Kermadec Subduction Zone')
c = colorbar('FontSize', 10,'TickLabels', {'0.1','1','10'}); set(gca,'ColorScale','log', 'CLim', [0.1 10]);
c.Label.String = 'Wave height (m)'; pos = get(c,'Position'); c.Label.Position = [2 pos(2)+0.6]; 

nexttile; % 500year return period
s = pcolor(lon_points,lat_points, Z); shading flat; s.FaceAlpha = 0.05; hold on;
plot(x,y, 'LineWidth', 1, 'Color', 'black'); hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, h_2500)
colormap(gca,'parula'); set(gca,'ColorScale','log', 'CLim', [0.1 10]);
set(gca, 'FontSize', 10); title('Hikurangi Subduction Margin')

nexttile; % 100year return period
s = pcolor(lon_points,lat_points, Z); shading flat; s.FaceAlpha = 0.05; hold on;
plot(x,y, 'LineWidth', 1, 'Color', 'black'); hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, p_2500)
colormap(gca,'parula'); c = colorbar('FontSize', 10,'TickLabels', {'0.1','1','10'}); set(gca,'ColorScale','log', 'CLim', [0.1 10]);
c.Label.String = 'Wave height (m)'; pos = get(c,'Position'); c.Label.Position = [2 pos(2)+1]; 
set(gca, 'FontSize', 10); title('Puysegur Subduction Zone')

%-------------------------%
event_number = event_data2(:, 1);
% Insert the event numbers into the overall matrix

time_events = sortrows(bruce2_waveheight.').';


t = maxk(time_events(2:end, :), 12, 2); min_exceed = mink(t, 1,2); all_2500 = round(min_exceed, 2); 
t = maxk(time_events(2:end,1:797), 12, 2); min_exceed = mink(t, 1,2); one_2500 = round(min_exceed, 2); 
t = maxk(time_events(2:end,798:1625), 12, 2); min_exceed = mink(t, 1,2); two_2500 = round(min_exceed, 2); 
t = maxk(time_events(2:end,1626:2595), 12, 2); min_exceed = mink(t, 1,2); three_2500 = round(min_exceed, 2); 

figure
tlo = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
% Dataset two
nexttile; % 2500year return period
s = pcolor(lon_points,lat_points, Z); shading flat; s.FaceAlpha = 0.05; hold on;
plot(x,y, 'LineWidth', 1, 'Color', 'black'); hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, all_2500)
colormap(gca,'parula'); set(gca,'ColorScale','log', 'CLim', [0.1 10]); 
set(gca, 'FontSize', 10); 
title('Entire catalogue'); 

nexttile;% 1000year return period
s = pcolor(lon_points,lat_points, Z); shading flat; s.FaceAlpha = 0.05; hold on;
plot(x,y, 'LineWidth', 1, 'Color', 'black'); hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, one_2500)
colormap(gca,'parula'); set(gca,'ColorScale','log', 'CLim', [0.1 10]);
set(gca, 'FontSize', 10); title('30,000-40,000years')
c = colorbar('FontSize', 10,'TickLabels', {'0.1','1','10'}); set(gca,'ColorScale','log', 'CLim', [0.1 10]);
c.Label.String = 'Wave height (m)'; pos = get(c,'Position'); c.Label.Position = [2 pos(2)+0.6]; 

nexttile; % 500year return period
s = pcolor(lon_points,lat_points, Z); shading flat; s.FaceAlpha = 0.05; hold on;
plot(x,y, 'LineWidth', 1, 'Color', 'black'); hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, two_2500)
colormap(gca,'parula'); set(gca,'ColorScale','log', 'CLim', [0.1 10]);
set(gca, 'FontSize', 10); title('40,000-50,000years')

nexttile; % 100year return period
s = pcolor(lon_points,lat_points, Z); shading flat; s.FaceAlpha = 0.05; hold on;
plot(x,y, 'LineWidth', 1, 'Color', 'black'); hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, three_2500)
colormap(gca,'parula'); c = colorbar('FontSize', 10,'TickLabels', {'0.1','1','10'}); set(gca,'ColorScale','log', 'CLim', [0.1 10]);
c.Label.String = 'Wave height (m)'; pos = get(c,'Position'); c.Label.Position = [2 pos(2)+1]; 
set(gca, 'FontSize', 10); title('50,000-60,000years')



coastline = 'C:\Users\hughesla\OneDrive - Victoria University of Wellington - STAFF\Documents\Coastlines\MATLAB_coast.mat';
load(coastline);

event2 = csvread('F:\PhD\bruce2_NSHM\tsunamis\tsunami_events.csv', 1,0);
events = sortrows(event2, 10); hikurangi_events = events(371:495);

bruce2_tsunami = ['F:\PhD\bruce2_NSHM\tsunamis\zmax_files\'];
% Load the coastal data points
coastal_points = load('C:\Users\hughesla\OneDrive - Victoria University of Wellington - STAFF\Documents\hazard_statistics\coastal_points.mat');
coastal_i = coastal_points.coastal_points(:,3); coastal_j = coastal_points.coastal_points(:,4);
coastal_x = coastal_points.coastal_points(:,1); coastal_y = coastal_points.coastal_points(:,2);

% Load layer 3 x and y reference
layer3_x = importdata('F:\PhD\bruce2_NSHM\tsunamis\layer03_x.dat');
layer3_y = importdata('F:\PhD\bruce2_NSHM\tsunamis\layer03_y.dat');


bruce2_maxheight = zeros(length(event2), 1);
max_locations = zeros(length(event2), 2);

for ii = 1:2%length(event2)
   disp(ii)
   filename = [bruce2_tsunami,'zmax_layer03_event', num2str(hikurangi_events(ii)),'.dat']; 
   if exist(filename, 'file')
      layer03_max = importdata(filename);
      A1 = permute(layer03_max, [2 1]); B1 = reshape(A1, 1560, 2238); 
      B1(any(isnan(B1), 2), :) = []; Z = B1.';
      for k=1:length(coastal_x)
          zm(k) = Z(coastal_j(k), coastal_i(k));
      end
      [max_height, I] = max(zm);
      bruce2_maxheight(ii, 1) = max_height;
      if max_height > 0
         x_point = coastal_x(I); y_point = coastal_y(I);
         max_locations(ii, [1,2]) = [x_point, y_point];
      end
   else
      disp('No max height file') 
   end  
end

bruce2_tsunami_info = [max_locations, bruce2_maxheight];