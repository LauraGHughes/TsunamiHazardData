% Maps showing location of earthquakes:
event_data2 = csvread('Earthquake_information.csv', 1,0);
nzsm_data = zeros(2585, 5);
for ii = 1:2585
    disp(ii)
    event_x = event_data2(ii,5); event_y = event_data2(ii,6); 
    event_z = event_data2(ii,7);
    [event_lon, event_lat] = UTMzoneX2Geo(event_x, event_y, 59);
    event_info = [event_data2(ii,1), event_lon, event_lat, event_z, event_data2(ii,4)];
    nzsm_data(ii, [1,2,3,4, 5]) = event_info;
end

pac = load('pacific.mat'); lat = pac.pac(:,1); lon = pac.pac(:,2);

f = figure('visible','off');
%figure
tlo = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
nexttile(1, [2,1])
s = pcolor(lon_points,lat_points, Z); shading flat; s.FaceAlpha = 0.1; colormap parula; caxis([5 20]);
set(gca, 'FontSize',11);
hold on
ylim([-50 -30]); xlim([163 185]); text(159.5,-30,'A', 'FontSize', 15, 'FontWeight', 'bold');
plot(lon, lat, 'LineWidth', 2, 'Color', 'black')
t1 = text(178.3,-42,'Hikurangi'); set(t1,'Rotation',67);
t2 = text(163.5,-48,'Puysegur'); set(t2,'Rotation',66);
t3 = text(180.9, -34,'Kermadec');set(t3,'Rotation',70);
t4 = text(165,-32,'AUS');t5 = text(180,-48,'PAC');
t6 = text(175.5,-39,'NI');t7 = text(169,-45,'SI');
t8 = text(174.2,-41.7,'CS'); t9 = text(168.9, -43.4,'AF');set(t9,'Rotation',45);

scatter(174.75, -36.9, 50, 'square', 'MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',1);
t10 = text(174.5, -36.9,'A','Color','red');
scatter(174.7, -41.3, 50, 'square', 'MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',1);
t11 = text(174.45, -41.3,'W','Color','red');
scatter(172.6, -43.5, 50, 'square', 'MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',1);
t12 = text(172.35, -43.5,'C','Color','red');
hold off  

nexttile(5)
b = histogram(event_data2(:,4));
xlabel('Magnitude', 'FontSize',10); ylabel('Frequency', 'FontSize',10);
xlim([7 9.3]); ylim([0 350]);
set(gca, 'FontSize',10)
b.FaceColor = '#000000'; b.EdgeColor = '#000000';hold on; text(6.65,350,'B', 'FontSize', 15, 'FontWeight', 'bold')
text(8.8,325,'n = 2585', 'FontSize', 11, 'FontWeight', 'bold')


nexttile(2, [3,1]) % Dataset two
ss = pcolor(lon_points,lat_points, Z); shading flat; ss.FaceAlpha = 0.1; 
ylim([-50 -23]); xlim([163 187]); set(gca, 'FontSize',11);
hold on;
scatter(nzsm_data(:,2),nzsm_data(:,3),(10.^nzsm_data(:,5))/1000000,(nzsm_data(:,4)*-1)/1000,'filled','MarkerFaceAlpha',0.5)
c = colorbar('eastoutside', 'FontSize', 11);  colormap parula; caxis([5 20])
ylabel(c,'Depth (km)','FontSize',11,'Rotation',270);
set(gca, 'FontSize',11); text(161,-23,'C', 'FontSize', 15, 'FontWeight', 'bold')
scatter(165.8,-28, 10^7.5/1000000, 'blue','filled','MarkerFaceAlpha',0.5); text(167.7,-28.1,'Mw7.5','FontSize',10.5);
scatter(165.8,-27.5, 10^8/1000000, 'blue','filled','MarkerFaceAlpha',0.5); text(167.7,-27.5,'Mw8.0','FontSize',10.5);
scatter(165.8,-26.6, 10^8.5/1000000, 'blue','filled','MarkerFaceAlpha',0.5); text(167.7,-26.6,'Mw8.5','FontSize',10.5);
scatter(165.8,-25, 10^9/1000000, 'blue','filled','MarkerFaceAlpha',0.5); text(167.6,-25,'Mw9.0','FontSize',10.5);
hold off

set(gcf,'position',[280,349,960,629]) 

fig = gcf;
exportgraphics(fig, 'earthquake_catalogue.png', 'BackgroundColor','white')


%%-----------------------------------------------------------------------%%
% Largest events
% Maximum event - event 223533
mycolormap = customcolormap(linspace(0,1,11), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3','#d5e2f0','#93c5dc','#4295c1','#2265ad','#062e61'});

% Deformation information
event_path = 'deformation_models\'; event = 'ev223533';
eventX_NZTM = ncread([event_path, event, '.grd'],'x'); eventY_NZTM = ncread([event_path, event, '.grd'],'y'); 
eventZ_NZTM223533 = ncread([event_path, event, '.grd'],'z'); 
[NZTM_NN, NZTM_EE] = meshgrid(eventY_NZTM, eventX_NZTM); [Longitude223533, Latitude223533] = NZTM2Geo(NZTM_EE, NZTM_NN); 

%Grid 3 information
filename = ['event223533\zmax_layer03.dat'];
layer03_max = importdata(filename); A1 = permute(layer03_max, [2 1]);
B1 = reshape(A1, 1560, 2238); B1(any(isnan(B1), 2), :) = []; maximum_grid223533 = B1.';
filename = ['event223533\layer03_x.dat']; layer03_x = importdata(filename);
filename = ['event223533\layer03_y.dat']; layer03_y = importdata(filename);

[event_lon223533, event_lat223533] = UTMzoneX2Geo(1903534.834, 7241987.405, 59);


event_path = 'deformation_models\'; event = 'ev199940';
eventX_NZTM = ncread([event_path, event, '.grd'],'x'); eventY_NZTM = ncread([event_path, event, '.grd'],'y'); 
eventZ_NZTM199940 = ncread([event_path, event, '.grd'],'z'); 
[NZTM_NN, NZTM_EE] = meshgrid(eventY_NZTM, eventX_NZTM); [Longitude199940, Latitude199940] = NZTM2Geo(NZTM_EE, NZTM_NN); 

%Grid 3 information
filename = ['zmax_layer03_event199940.dat'];
layer03_max = importdata(filename); A1 = permute(layer03_max, [2 1]);
B1 = reshape(A1, 1560, 2238); B1(any(isnan(B1), 2), :) = []; maximum_grid199940 = B1.';
filename = ['event223533\layer03_x.dat']; layer03_x = importdata(filename);
filename = ['event223533\layer03_y.dat']; layer03_y = importdata(filename);

[event_lon199940, event_lat199940] = UTMzoneX2Geo(1850844.943,7137532.401, 59);

f = figure('visible','off');
%figure
tlo = tiledlayout(6,3,'TileSpacing','compact','Padding','compact');
ax1 = nexttile(1,[3,1]); % Deformation model
pcolor(Longitude223533, Latitude223533,  eventZ_NZTM223533); colormap(ax1, mycolormap); caxis([-5 5]);set(gca, 'FontSize',10);
shading flat; a = colorbar('eastoutside', 'FontSize', 10); a.Label.String = 'Vertical displacement (m)';
a.Label.FontSize = 10;
hold on
scatter(event_lon223533, event_lat223533, 100, 'o', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1)
xlim([172 186]); ylim([-39 -23]);text(171,-23.5,'A', 'FontSize', 15, 'FontWeight', 'bold');
yticks([-45 -40 -35 -30 -25]); xticks([165 170 175 180 185]);
title('Mw 9.25 earthquake displacement', 'FontWeight', 'bold')
hold off;

ax2 = nexttile(2,[3,1]); % Grid one
pcolor(layer03_x, layer03_y, maximum_grid223533); shading flat; colormap(ax2, parula); 
hold on
set(gca,'ColorScale','log', 'CLim', [0.1 20]);set(gca, 'FontSize',10);
b = colorbar('FontSize', 10,'TickLabels', {'0.1','1','10'}); set(gca,'ColorScale','log', 'CLim', [0.1 10]);
b.Label.String = 'Wave height (m)'; b.Label.FontSize = 10;
text(165,-34,'B', 'FontSize', 15, 'FontWeight', 'bold');
pos = get(b,'Position'); b.Label.Position = [2.5 pos(2)+0.65]; 
yticks([-45 -40 -35 -30 -25]); xticks([165 170 175 180 185]);
title('Mw 9.25 tsunami', 'FontWeight', 'bold')
hold off

ax1 = nexttile(10,[3,1]); % Deformation model
pcolor(Longitude199940, Latitude199940,  eventZ_NZTM199940); colormap(ax1, mycolormap); caxis([-5 5]);set(gca, 'FontSize',10);
shading flat; a = colorbar('eastoutside', 'FontSize', 10); a.Label.String = 'Vertical displacement (m)';
a.Label.FontSize = 10;
hold on
scatter(event_lon199940, event_lat199940, 100, 'o', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1)
xlim([172 186]); ylim([-39 -23]);text(171,-23.5,'C', 'FontSize', 15, 'FontWeight', 'bold');
yticks([-45 -40 -35 -30 -25]); xticks([165 170 175 180 185]);
title('Mw 9.22 earthquake displacement', 'FontWeight', 'bold')
hold off;

ax2 = nexttile(11,[3,1]); % Grid one
pcolor(layer03_x, layer03_y, maximum_grid199940); shading flat; colormap(ax2, parula); 
hold on
set(gca,'ColorScale','log', 'CLim', [0.1 20]);set(gca, 'FontSize',10);
b = colorbar('FontSize', 10,'TickLabels', {'0.1','1','10'}); set(gca,'ColorScale','log', 'CLim', [0.1 10]);
b.Label.String = 'Wave height (m)'; b.Label.FontSize = 10;
text(165,-34,'D', 'FontSize', 15, 'FontWeight', 'bold');
pos = get(b,'Position'); b.Label.Position = [2.5 pos(2)+1.1]; 
yticks([-45 -40 -35 -30 -25]); xticks([165 170 175 180 185]);
title('Mw 9.22 tsunami', 'FontWeight', 'bold')
hold off


ax3 = nexttile(6,[4,1]); % Grid one
pcolor(layer03_x, layer03_y, maximum_grid223533 - maximum_grid199940); shading flat; colormap(ax3, parula); 
hold on
caxis([0 2.5]); set(gca, 'FontSize',10);
plot(NIx,NIy, 'LineWidth', 1, 'Color', 'white'); plot(SIx,SIy, 'LineWidth', 1, 'Color', 'white');
b = colorbar('FontSize', 10); 
b.Label.String = 'Wave height (m)'; b.Label.FontSize = 10;
text(165,-34,'E', 'FontSize', 15, 'FontWeight', 'bold');
pos = get(b,'Position'); b.Label.Position = [2.5 pos(2)+1]; 
yticks([-45 -40 -35 -30 -25]); xticks([165 170 175 180 185]);
title('Tsunami difference Mw9.25 - Mw9.22', 'FontWeight', 'bold')
hold off

%pos = get(gcf, 'Position')
set(gcf,'position',[196,347,1187,631]) 

fig = gcf;
exportgraphics(fig, 'largest_magnitude_comparison_alt_15_05_2023.png', 'BackgroundColor','white')

%%-----------------------------------------------------------------------%%
% Maximum wave height at the coast
load('coastal_maximum.mat');
bruce2 = sortrows(coastal_maximum, 4);

f = figure('visible','off');
%figure
tlo = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
nexttile(1, [2,1]) % Catalogue two
b = histogram(bruce2(:,4));
xlabel('Maximum wave height (m)', 'FontSize',10); ylabel('Number of events', 'FontSize',10);%title('DS2 maximum wave height distribution')
set(gca, 'YScale', 'log'); ylim([0.9 1000]); xlim([0 28]);
set(gca, 'FontSize',10)
b.FaceColor = '#000000'; b.EdgeColor = '#000000';hold on; text(-4,1000,'A', 'FontSize', 15, 'FontWeight', 'bold')
text(22,750,'n = 2585', 'FontSize', 11, 'FontWeight', 'bold')

nexttile(2,[3 1])
hold on
set(gca, 'FontSize',10);xlim([166 179]); ylim([-48 -34]);
scatter(bruce2(:,2), bruce2(:,3),50, bruce2(:,4),'filled', 'MarkerFaceAlpha', 0.5)
colormap parula;
a = colorbar('eastoutside', 'FontSize', 10); a.Label.String = 'Wave height (m)';
a.Label.FontSize = 11.5; set(gca,'ColorScale','log', 'CLim', [0.1 30]); %title('DS2 maximum wave height coastal location')
a.Ticks = unique([0.1, 1, 10, 20, 30]);  
a.TickLabels{1} = '0.1'; a.TickLabels{2} = '1'; a.TickLabels{3} = '10'; a.TickLabels{4} = '20'; a.TickLabels{5} = '30';
pos = get(a,'Position');
a.Label.Position = [2 pos(2)+2]; 
text(164,-34,'B', 'FontSize', 15, 'FontWeight', 'bold')
hold off

%pos = get(gcf, 'Position')
set(gcf,'position',[488, 587, 752, 391]) 

fig = gcf;
exportgraphics(fig, 'maximum_waveheight_to_coast_27_4_2023.png', 'BackgroundColor','white')


%%-----------------------------------------------------------------------%%
%%-----------------------------------------------------------------------%%
% Return periods dataset one:
load('F:\PhD\bruce2_NSHM\tsunamis\coastal_values.mat');
bruce2_waveheight = coastal_values(2:end, :);

% 2500year return period 
t = maxk(bruce2_waveheight, 12, 2); min_exceed = mink(t, 1,2); rp2_2500 = round(min_exceed, 2); 
% 1000year return period 
t = maxk(bruce2_waveheight, 30, 2); min_exceed = mink(t, 1,2); rp2_1000 = round(min_exceed, 2); 
% 500year return period 
t = maxk(bruce2_waveheight, 60, 2); min_exceed = mink(t, 1,2); rp2_500 = round(min_exceed, 2); 
% 100year return period 
t = maxk(bruce2_waveheight, 300, 2); min_exceed = mink(t, 1,2); rp2_100 = round(min_exceed, 2); 


coastal_points = load('coastal_points.mat');
coastal_i = coastal_points.coastal_points(:,3); coastal_j = coastal_points.coastal_points(:,4);
coastal_x = coastal_points.coastal_points(:,1); coastal_y = coastal_points.coastal_points(:,2);

f = figure('visible','off');
%figure
tlo = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
% Dataset two
nexttile; % 2500year return period
hold on;
xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, rp2_2500); 
scatter(174.75, -36.9, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1);
scatter(174.7, -41.3, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1);
scatter(172.6, -43.5, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1);
colormap(gca,'parula'); set(gca,'ColorScale','log', 'CLim', [0.1 10]); 
set(gca, 'FontSize', 10); yticks([-45 -40 -35]); xticks([165 170 175 180]);
text(163.5,-34,'A', 'FontSize', 15, 'FontWeight', 'bold')
title('2,500-year return period'); 

nexttile;% 1000year return period
hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, rp2_1000);
scatter(174.75, -36.9, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1);
scatter(174.7, -41.3, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1);
scatter(172.6, -43.5, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1);
colormap(gca,'parula'); set(gca,'ColorScale','log', 'CLim', [0.1 10]);
text(163.5,-34,'B', 'FontSize', 15, 'FontWeight', 'bold')
set(gca, 'FontSize', 10); title('1,000-year return period'); yticks([-45 -40 -35]); xticks([165 170 175 180]);
c = colorbar('FontSize', 10,'TickLabels', {'0.1','1','10'}); set(gca,'ColorScale','log', 'CLim', [0.1 10]);
c.Label.String = 'Wave height (m)'; pos = get(c,'Position'); c.Label.Position = [2 pos(2)+0.6]; 

nexttile; % 500year return period
hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, rp2_500); 
scatter(174.75, -36.9, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1);
scatter(174.7, -41.3, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1);
scatter(172.6, -43.5, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1);
yticks([-45 -40 -35]); xticks([165 170 175 180]);
colormap(gca,'parula'); set(gca,'ColorScale','log', 'CLim', [0.1 10]);
text(163.5,-34,'C', 'FontSize', 15, 'FontWeight', 'bold')
set(gca, 'FontSize', 10); title('500-year return period')

nexttile; % 100year return period
hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, rp2_100); 
scatter(174.75, -36.9, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1);
scatter(174.7, -41.3, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1);
scatter(172.6, -43.5, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1);
yticks([-45 -40 -35]); xticks([165 170 175 180]);
colormap(gca,'parula'); c = colorbar('FontSize', 10,'TickLabels', {'0.1','1','10'}); set(gca,'ColorScale','log', 'CLim', [0.1 10]);
c.Label.String = 'Wave height (m)'; pos = get(c,'Position'); c.Label.Position = [2 pos(2)+1]; 
text(163.5,-34,'D', 'FontSize', 15, 'FontWeight', 'bold')
set(gca, 'FontSize', 10); title('100-year return period')

%pos = get(gcf, 'Position')
set(gcf,'position',[604,409,636,569]) 

fig = gcf;
exportgraphics(fig, 'tsunami_hazard_return_periods_20_03_2023.png', 'BackgroundColor','white')
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Event nucleation:

event_data2 = csvread('Earthquake_information.csv', 1,0);
event_data2 = sortrows(event_data2, 1);
event_type = event_data2(:, 10);
events_k = []; events_h = []; events_p = []; events_c = []; 

for ii = 1:length(event_type)
   type = event_type(ii);
   disp(ii)
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


coastal_points = load('coastal_points.mat');
coastal_i = coastal_points.coastal_points(:,3); coastal_j = coastal_points.coastal_points(:,4);
coastal_x = coastal_points.coastal_points(:,1); coastal_y = coastal_points.coastal_points(:,2);

f = figure('visible','off');
%figure
tlo = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
% Dataset two
nexttile; % 2500year return period
hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, all_2500);
scatter(174.75, -36.9, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1);
scatter(174.7, -41.3, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1);
scatter(172.6, -43.5, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1);
colormap(gca,'parula'); set(gca,'ColorScale','log', 'CLim', [0.1 10]); 
set(gca, 'FontSize', 10); 
title('Full Dataset'); 
text(163.5,-34,'A', 'FontSize', 15, 'FontWeight', 'bold')
yticks([-45 -40 -35]); xticks([165 170 175 180]);

nexttile;% 1000year return period
hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, k_2500);
scatter(174.75, -36.9, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1);
scatter(174.7, -41.3, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1);
scatter(172.6, -43.5, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1);
colormap(gca,'parula'); set(gca,'ColorScale','log', 'CLim', [0.1 10]);
set(gca, 'FontSize', 10); title('Tonga-Kermadec Subduction Zone')
c = colorbar('FontSize', 10,'TickLabels', {'0.1','1','10'}); set(gca,'ColorScale','log', 'CLim', [0.1 10]);
c.Label.String = 'Wave height (m)'; pos = get(c,'Position'); c.Label.Position = [2 pos(2)+0.6]; 
text(163.5,-34,'B', 'FontSize', 15, 'FontWeight', 'bold')
yticks([-45 -40 -35]); xticks([165 170 175 180]);

nexttile; % 500year return period
hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, h_2500);
scatter(174.75, -36.9, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1);
scatter(174.7, -41.3, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1);
scatter(172.6, -43.5, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1);
colormap(gca,'parula'); set(gca,'ColorScale','log', 'CLim', [0.1 10]);
set(gca, 'FontSize', 10); title('Hikurangi Subduction Margin')
text(163.5,-34,'C', 'FontSize', 15, 'FontWeight', 'bold')
yticks([-45 -40 -35]); xticks([165 170 175 180]);

nexttile; % 100year return period
hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, p_2500);
scatter(174.75, -36.9, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1);
scatter(174.7, -41.3, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1);
scatter(172.6, -43.5, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','black','LineWidth',1);
colormap(gca,'parula'); c = colorbar('FontSize', 10,'TickLabels', {'0.1','1','10'}); set(gca,'ColorScale','log', 'CLim', [0.1 10]);
c.Label.String = 'Wave height (m)'; pos = get(c,'Position'); c.Label.Position = [2 pos(2)+1]; 
set(gca, 'FontSize', 10); title('Puysegur Subduction Zone')
text(163.5,-34,'D', 'FontSize', 15, 'FontWeight', 'bold')
yticks([-45 -40 -35]); xticks([165 170 175 180]);

set(gcf,'position',[604,409,636,569]) 

fig = gcf;
exportgraphics(fig, 'tsunami_hazard_return_periods_location_20_03_2023.png', 'BackgroundColor','white')

%%-----------------------------------------------------------------------%%
% Tsunami hazard curves:
load('coastal_values.mat');
BES2a_coastline = coastal_values;

% Coastal data
coastal_points = load('coastal_points.mat');
coastal_i = coastal_points.coastal_points(:,3); coastal_j = coastal_points.coastal_points(:,4);
coastal_x = coastal_points.coastal_points(:,1); coastal_y = coastal_points.coastal_points(:,2);

% New Zealand's coastline
Northland  = [172, 175, 175; -34.3, -36.5, -34.3]';
Eastcoast = [175.5, 175.5, 177, 177; -40.5, -42, -42, -40.5]';
HawkesBay = [176.5, 176.5, 178.1, 178.1; -39, -39.7, -39.7, -39]';
southland = [166.3, 167.9, 167.9, 166.3; -46.43, -46.43, -44.25, -44.25]';

BES2a_tsunami_data = zeros(2585, 5);
for event = 1:2585
    disp(event)
    event_id = BES2a_coastline(1,event);
    waveheight = BES2a_coastline(2:end, event);
    p100 = prctile(waveheight,99);
    p90 = prctile(waveheight,90); p10 = prctile(waveheight,10); p50 = prctile(waveheight,50);
    event_data = [event_id, p100, p90, p50, p10];
    BES2a_tsunami_data(event, [1,2,3,4,5]) = event_data; 
end


BES2a_tsunami_data = zeros(2585, 16);
for event = 1:2585
    disp(event)
    event_id = BES2a_coastline(1,event);
    waveheight = BES2a_coastline(2:end, event);
    
    %All of New Zealand
    waveheight = BES2a_coastline(2:end, event);
    p90 = prctile(waveheight,90); p10 = prctile(waveheight,10); p50 = prctile(waveheight,50);

    % For other regions
    temp_Northland = [];temp_Eastcoast = [];temp_southland = []; temp_HawkesBay = [];
    for point = 1:length(waveheight)
        x_point = coastal_x(point); y_point = coastal_y(point);
        if inpolygon(x_point, y_point, Northland(:, 1), Northland(:, 2)) == 1
           temp_Northland = [temp_Northland, waveheight(point)];
        elseif inpolygon(x_point, y_point, Eastcoast(:, 1), Eastcoast(:, 2)) == 1
           temp_Eastcoast = [temp_Eastcoast, waveheight(point)];
        elseif inpolygon(x_point, y_point, southland(:, 1), southland(:, 2)) == 1
           temp_southland = [temp_southland, waveheight(point)];
        elseif inpolygon(x_point, y_point, HawkesBay(:, 1), HawkesBay(:, 2)) == 1
           temp_HawkesBay = [temp_HawkesBay, waveheight(point)];           
        end    
    end
    p90_Northland = prctile(temp_Northland,90); p10_Northland = prctile(temp_Northland ,10); p50_Northland = prctile(temp_Northland,50);
    p90_Eastcoast = prctile(temp_Eastcoast,90); p10_Eastcoast = prctile(temp_Eastcoast,10); p50_Eastcoast = prctile(temp_Eastcoast,50);
    p90_southland = prctile(temp_southland,90); p10_southland = prctile(temp_southland,10); p50_southland = prctile(temp_southland,50);
    p90_HawkesBay = prctile(temp_HawkesBay,90); p10_HawkesBay = prctile(temp_HawkesBay,10); p50_HawkesBay = prctile(temp_HawkesBay,50);
    
    event_data = [event_id, p90, p50, p10, p90_Northland, p50_Northland, p10_Northland, p90_Eastcoast, p50_Eastcoast, p10_Eastcoast,...
        p90_southland, p50_southland, p10_southland, p90_HawkesBay, p50_HawkesBay, p10_HawkesBay];
    BES2a_tsunami_data(event, [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]) = event_data;  
end


% Calculate numbers for return periods
numbers = [1:1:2585]; return_period_numbers_BES2a = zeros(2585,1);
for ii = 1:length(numbers)
    disp(ii)
    num = 30000/ii; return_period_numbers_BES2a(ii, 1) = num;
end

x1 = [50, 100, 100, 50; 0, 0, 15, 15]';
x2 = [5000, 30000, 30000, 5000;0, 0, 15, 15]';

f = figure('visible','off');
%figure
tlo = tiledlayout(4,3,'TileSpacing','compact','Padding','compact');
%nexttile(2);
nexttile(1,[1,2]);
hold on;
plot(sort(return_period_numbers_BES2a)', sort(BES2a_tsunami_data(:,5)),'--','Color', 'k', 'LineWidth',2);
plot(sort(return_period_numbers_BES2a)', sort(BES2a_tsunami_data(:,6)),'Color', 'k', 'LineWidth',2);
plot(sort(return_period_numbers_BES2a)', sort(BES2a_tsunami_data(:,7)),'-.','Color', 'k', 'LineWidth',2);
xline(100, '-', {'100-year'}, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment','top');
xline(250, '-', {'250-year'}, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment','top'); xline(500, '-', {'500-year'}, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment','top'); 
xline(1000, '-', {'1000-year'}, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment','top'); xline(2500, '-', {'2500-year'}); 
xline(5000, '-', {'5000-year'}); 
plot(polyshape(x1),'FaceColor', 'black', 'FaceAlpha', 0.075, 'EdgeColor', 'none'); plot(polyshape(x2),'FaceColor','black', 'FaceAlpha', 0.075, 'EdgeColor', 'none'); box on;
set(gca, 'XScale', 'log');
ylabel('Wave height (m)'); xlabel('Return period');
xlim([50 30000]); ylim([0 15]); %ylim([0.001 30]);
xticks([100 500 1000 5000 10000 30000]); xticklabels({'100','500','1000','5000','10000','30000'});
title('Tsunami hazard curve - Northern North Island');box on;
text(30,16,'A', 'FontSize', 15, 'FontWeight', 'bold'); hold off; 

nexttile(3);
hold on; colormap(bone); caxis([0 1]);
plot(polyshape(Northland),'FaceColor', 'black', 'FaceAlpha', 0.6, 'EdgeColor', 'black');
scatter(174.75, -36.9, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','red','LineWidth',1);
scatter(174.7, -41.3, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','red','LineWidth',1);
scatter(172.6, -43.5, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','red','LineWidth',1);
ylim([-48 -34]); yticks([-45 -40 -35]); xlim([165 180]);  xticks([165 170 175 180]);
text(164,-33,'B', 'FontSize', 15, 'FontWeight', 'bold'); hold off; 

nexttile(4,[1,2]);
hold on;
plot(sort(return_period_numbers_BES2a)', sort(BES2a_tsunami_data(:,14)),'--','Color', 'k', 'LineWidth',2);
plot(sort(return_period_numbers_BES2a)', sort(BES2a_tsunami_data(:,15)),'Color', 'k', 'LineWidth',2);
plot(sort(return_period_numbers_BES2a)', sort(BES2a_tsunami_data(:,16)),'-.','Color', 'k', 'LineWidth',2);
xline(100, '-', {'100-year'}, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment','top');
xline(250, '-', {'250-year'}, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment','top'); xline(500, '-', {'500-year'}, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment','top'); 
xline(1000, '-', {'1000-year'}, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment','top'); xline(2500, '-', {'2500-year'}); 
xline(5000, '-', {'5000-year'}); 
plot(polyshape(x1),'FaceColor', 'black', 'FaceAlpha', 0.075, 'EdgeColor', 'none'); plot(polyshape(x2),'FaceColor','black', 'FaceAlpha', 0.075, 'EdgeColor', 'none'); box on;
set(gca, 'XScale', 'log');
ylabel('Wave height (m)'); xlabel('Return period');
xlim([50 30000]); ylim([0 15]); %ylim([0.001 30]);
xticks([100 500 1000 5000 10000 30000]); xticklabels({'100','500','1000','5000','10000','30000'});
title('Tsunami hazard curve - Hawkes Bay');box on;
text(30,16,'C', 'FontSize', 15, 'FontWeight', 'bold'); hold off; 

nexttile(6);
hold on; colormap(bone); caxis([0 1]);
plot(polyshape(HawkesBay),'FaceColor', 'black', 'FaceAlpha', 0.6, 'EdgeColor', 'black');
ylim([-48 -34]); yticks([-45 -40 -35]); xlim([165 180]);  xticks([165 170 175 180]);
scatter(174.75, -36.9, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','red','LineWidth',1);
scatter(174.7, -41.3, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','red','LineWidth',1);
scatter(172.6, -43.5, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','red','LineWidth',1);
text(164,-33,'D', 'FontSize', 15, 'FontWeight', 'bold'); hold off; 

nexttile(7,[1,2]);
hold on;
plot(sort(return_period_numbers_BES2a)', sort(BES2a_tsunami_data(:,8)),'--','Color', 'k', 'LineWidth',2);
plot(sort(return_period_numbers_BES2a)', sort(BES2a_tsunami_data(:,9)),'Color', 'k', 'LineWidth',2);
plot(sort(return_period_numbers_BES2a)', sort(BES2a_tsunami_data(:,10)),'-.','Color', 'k', 'LineWidth',2);
xline(100, '-', {'100-year'}, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment','top');
xline(250, '-', {'250-year'}, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment','top'); xline(500, '-', {'500-year'}, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment','top'); 
xline(1000, '-', {'1000-year'}, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment','top'); xline(2500, '-', {'2500-year'}); 
xline(5000, '-', {'5000-year'}); 
plot(polyshape(x1),'FaceColor', 'black', 'FaceAlpha', 0.075, 'EdgeColor', 'none'); plot(polyshape(x2),'FaceColor','black', 'FaceAlpha', 0.075, 'EdgeColor', 'none'); box on;
set(gca, 'XScale', 'log');
ylabel('Wave height (m)'); xlabel('Return period');
xlim([50 30000]); ylim([0 15]); %ylim([0.001 30]);
xticks([100 500 1000 5000 10000 30000]); xticklabels({'100','500','1000','5000','10000','30000'});
title('Tsunami hazard curve - Southeast North Island');box on;
text(30,16,'E', 'FontSize', 15, 'FontWeight', 'bold'); hold off; 

nexttile(9);
hold on; colormap(bone); caxis([0 1]);
plot(polyshape(Eastcoast),'FaceColor', 'black', 'FaceAlpha', 0.6, 'EdgeColor', 'black');
ylim([-48 -34]); yticks([-45 -40 -35]); xlim([165 180]);  xticks([165 170 175 180]);
scatter(174.75, -36.9, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','red','LineWidth',1);
scatter(174.7, -41.3, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','red','LineWidth',1);
scatter(172.6, -43.5, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','red','LineWidth',1);
text(164,-33,'F', 'FontSize', 15, 'FontWeight', 'bold'); hold off; 

nexttile(10,[1,2]);
hold on;
plot(sort(return_period_numbers_BES2a)', sort(BES2a_tsunami_data(:,11)),'--','Color', 'k', 'LineWidth',2);
plot(sort(return_period_numbers_BES2a)', sort(BES2a_tsunami_data(:,12)),'Color', 'k', 'LineWidth',2);
plot(sort(return_period_numbers_BES2a)', sort(BES2a_tsunami_data(:,13)),'-.','Color', 'k', 'LineWidth',2);
xline(100, '-', {'100-year'}, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment','top');
xline(250, '-', {'250-year'}, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment','top'); xline(500, '-', {'500-year'}, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment','top'); 
xline(1000, '-', {'1000-year'}, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment','top'); xline(2500, '-', {'2500-year'}, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment','top'); 
xline(5000, '-'); 
plot(polyshape(x1),'FaceColor', 'black', 'FaceAlpha', 0.075, 'EdgeColor', 'none'); plot(polyshape(x2),'FaceColor','black', 'FaceAlpha', 0.075, 'EdgeColor', 'none'); box on;
set(gca, 'XScale', 'log');
ylabel('Wave height (m)'); xlabel('Return period');
xlim([50 30000]); ylim([0 15]); %ylim([0.001 30]);
xticks([100 500 1000 5000 10000 30000]); xticklabels({'100','500','1000','5000','10000','30000'});
legend({'90th percentile', '50th percentile','10th percentile'}, 'Location','northeast');
title('Tsunami hazard curve - Southwest South Island'); box on;
text(30,16,'G', 'FontSize', 15, 'FontWeight', 'bold'); hold off; 

nexttile(12);
shold on; colormap(bone); caxis([0 1]);
plot(polyshape(southland),'FaceColor', 'black', 'FaceAlpha', 0.6, 'EdgeColor', 'black');
ylim([-48 -34]); yticks([-45 -40 -35]); xlim([165 180]);  xticks([165 170 175 180]);
scatter(174.75, -36.9, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','red','LineWidth',1);
scatter(174.7, -41.3, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','red','LineWidth',1);
scatter(172.6, -43.5, 50, 'square', 'MarkerEdgeColor','white','MarkerFaceColor','red','LineWidth',1);
text(164,-33,'H', 'FontSize', 15, 'FontWeight', 'bold'); hold off; 

set(gcf,'position',[522,42,718,954])
fig = gcf;
exportgraphics(fig, 'hazard_curves_26_04_2023.png', 'BackgroundColor','white')

%%-----------------------------------------------------------------------%%
% Length of earthquake catalogue
earthquake_catalogue = csvread(BES2_long.csv', 1,0);

bincentr = [7.05:0.1:9.25];

a_bincount = [780,866,771,562,419,430,179,62,180,44,42,27,13,15,6,8,28,79,37,29,12,16,6];
b_bincount = [767,803,703,520,385,405,176,84,150,49,27,22,14,16,7,8,32,69,36,28,6,16,6];
c_bincount = [523,534,477,327,253,259,126,49,104,19,21,19,11,8,3,4,19,39,29,16,2,12,5];


f = figure('visible','off');
%figure
tlo = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
nexttile([1,3]);
scatter(earthquake_catalogue(:,1)/(3.154*10^7), earthquake_catalogue(:,3),5, 'filled', ...
    'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black','MarkerFaceAlpha',0.5);
ylabel('Magnitude (M_W)'); xlabel('Time (years)');
xticks([20000 30000 40000 50000 60000 70000 80000 90000 100000 110000]); 
xticklabels({'20000','30000','40000','50000','60000','70000', ...
    '80000', '90000', '100000', '110000'});
yticks([7.0:0.2:9.3])
xlim([20000 112000]);ylim([6.95 9.3]);
hold on; 
scatter(earthquake_catalogue(1751:6361,1)/(3.154*10^7), earthquake_catalogue(1751:6361,3),5,'filled',...
    'MarkerEdgeColor', '#8dd3c7', 'MarkerFaceColor', '#8dd3c7','MarkerFaceAlpha',0.5);
text(14500,9.3,'A', 'FontSize', 15, 'FontWeight', 'bold'); hold off;

nexttile;
a = histogram(earthquake_catalogue(1751:6361,3));
xlabel('Magnitude (M_W)'); ylabel('Frequency');
a.FaceColor = '#8dd3c7'; a.EdgeColor = 'black'; 
xlim([7 9.3]);ylim([1 10000]);
set(gca, 'YScale', 'log');
hold on; 
plot(bincentr, cumsum(a_bincount, 'reverse'), '-o', 'LineWidth',2,'Color', 'black');
yticks([1 10 100 1000 10000]);yticklabels({'1','10','100','1000', '10000'});
title('32,000-63,000-years')
text(6.5,10000,'B', 'FontSize', 15, 'FontWeight', 'bold'); hold off;

nexttile;
b = histogram(earthquake_catalogue(6361:10689,3));
xlabel('Magnitude (M_W)'); ylabel('Frequency');
b.FaceColor = '#000000';
xlim([7 9.3]);ylim([1 10000]);
set(gca, 'YScale', 'log');
hold on; 
plot(bincentr, cumsum(b_bincount, 'reverse'), '-o', 'LineWidth',2,'Color', 'black');
yticks([1 10 100 1000 10000]);yticklabels({'1','10','100','1000', '10000'});
title('63,000-93,000-years')
text(6.5,10000,'C', 'FontSize', 15, 'FontWeight', 'bold'); hold off;

nexttile;
c = histogram(earthquake_catalogue(10689:end,3));
xlabel('Magnitude (M_W)'); ylabel('Frequency');
c.FaceColor = '#000000';
xlim([7 9.3]);ylim([1 10000]);
set(gca, 'YScale', 'log');
hold on; 
plot(bincentr, cumsum(c_bincount, 'reverse'), '-o', 'LineWidth',2,'Color', 'black');
yticks([1 10 100 1000 10000]);yticklabels({'1','10','100','1000', '10000'});
title('93,000-112,000-years')
text(6.5,10000,'D', 'FontSize', 15, 'FontWeight', 'bold'); hold off;

set(gcf,'position',[258,472,982,506]) 
fig = gcf;
exportgraphics(fig, 'E:\PhD\journal_article\PaperOne\figures\magnitude_frequency_alternative_15_05_2023.png', 'BackgroundColor','white') % this saves the figure

%%-----------------------------------------------------------------------%%
%K-S test for catalogue similarity:
earthquake_catalogue = csvread('BES2_long.csv', 1,0);

subset1 = earthquake_catalogue(1751:6361,3);
subset2 = earthquake_catalogue(6361:10689,3);
subset3 = earthquake_catalogue(10689:end,3); 

[h_subset1_2,p_subset1_2] = kstest2(subset1,subset2,'Alpha',0.05, 'Tail', 'unequal');
[h_subset1_3,p_subset1_3] = kstest2(subset1,subset3,'Alpha',0.05, 'Tail', 'unequal');
[h_subset2_3,p_subset2_3] = kstest2(subset2,subset3,'Alpha',0.05, 'Tail', 'unequal');
%%-----------------------------------------------------------------------%%
% Tsunami hazard to the coast for alternative time periods:
load('coastal_values.mat');

a_15000 = coastal_values(:, 1:1043); b_15000 = coastal_values(:, 1043:2585);
a_10000 = coastal_values(:, 1:642); b_10000 = coastal_values(:, 642:1467); c_10000 = coastal_values(:, 1467:2585);

a_5000 = coastal_values(:, 1:264); b_5000 = coastal_values(:, 264:642); c_5000 = coastal_values(:, 642:1043);
d_5000 = coastal_values(:, 1043:1467); e_5000 = coastal_values(:, 1467:1892); f_5000 = coastal_values(:, 1892:2585);

t = maxk(coastal_values(2:end, :), 6, 2); min_exceed = mink(t, 1,2); a_30000_2500 = round(min_exceed, 2); 

t = maxk(a_15000(2:end, :), 6, 2); min_exceed = mink(t, 1,2); a_15000_2500 = round(min_exceed, 2); 
t = maxk(b_15000(2:end,:), 6, 2); min_exceed = mink(t, 1,2); b_15000_2500 = round(min_exceed, 2); 

t = maxk(a_10000(2:end, :), 4, 2); min_exceed = mink(t, 1,2); a_10000_2500 = round(min_exceed, 2); 
t = maxk(b_10000(2:end, :), 4, 2); min_exceed = mink(t, 1,2); b_10000_2500 = round(min_exceed, 2); 
t = maxk(c_10000(2:end, :), 4, 2); min_exceed = mink(t, 1,2); c_10000_2500 = round(min_exceed, 2); 

t = maxk(a_5000(2:end, :), 2, 2); min_exceed = mink(t, 1,2); a_5000_2500 = round(min_exceed, 2); 
t = maxk(b_5000(2:end, :), 2, 2); min_exceed = mink(t, 1,2); b_5000_2500 = round(min_exceed, 2);
t = maxk(c_5000(2:end, :), 2, 2); min_exceed = mink(t, 1,2); c_5000_2500 = round(min_exceed, 2);
t = maxk(d_5000(2:end, :), 2, 2); min_exceed = mink(t, 1,2); d_5000_2500 = round(min_exceed, 2);
t = maxk(e_5000(2:end, :), 2, 2); min_exceed = mink(t, 1,2); e_5000_2500 = round(min_exceed, 2);
t = maxk(f_5000(2:end, :), 2, 2); min_exceed = mink(t, 1,2); f_5000_2500 = round(min_exceed, 2);

coastal_points = load('C:\Users\hughesla\OneDrive - Victoria University of Wellington - STAFF\Documents\hazard_statistics\coastal_points.mat');
coastal_i = coastal_points.coastal_points(:,3); coastal_j = coastal_points.coastal_points(:,4);
coastal_x = coastal_points.coastal_points(:,1); coastal_y = coastal_points.coastal_points(:,2);

mycolormap = customcolormap(linspace(0,1,11), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3','#d5e2f0','#93c5dc','#4295c1','#2265ad','#062e61'});

f = figure('visible','off');
%figure
tlo = tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
nexttile(1); % 2500year return period
hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, a_15000_2500); colormap(gca,'parula'); set(gca,'ColorScale','log', 'CLim', [0.1 10]); 
set(gca, 'FontSize', 10); title('30,000-45,000-years'); 
text(163.5,-34,'A', 'FontSize', 15, 'FontWeight', 'bold');yticks([-45 -40 -35]); xticks([165 170 175 180]);

nexttile(2); % 2500year return period
hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, b_15000_2500); colormap(gca,'parula'); set(gca,'ColorScale','log', 'CLim', [0.1 10]); 
set(gca, 'FontSize', 10); title('45,000-60,000-years'); 
text(163.5,-34,'B', 'FontSize', 15, 'FontWeight', 'bold');yticks([-45 -40 -35]); xticks([165 170 175 180]);

nexttile(3); % 2500year return period
hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, a_15000_2500 - b_15000_2500); colormap(gca,mycolormap); set(gca, 'CLim', [-2 2]);
set(gca, 'FontSize', 10); title('A-B difference'); 
text(163.5,-34,'C', 'FontSize', 15, 'FontWeight', 'bold');yticks([-45 -40 -35]); xticks([165 170 175 180]);

nexttile(4); % 2500year return period
hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, a_10000_2500); colormap(gca,'parula'); set(gca,'ColorScale','log', 'CLim', [0.1 10]); 
set(gca, 'FontSize', 10); title('30,000-40,000-years'); 
text(163.5,-34,'D', 'FontSize', 15, 'FontWeight', 'bold');yticks([-45 -40 -35]); xticks([165 170 175 180]);

nexttile(5); % 2500year return period
hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, b_10000_2500); colormap(gca,'parula'); set(gca,'ColorScale','log', 'CLim', [0.1 10]); 
set(gca, 'FontSize', 10); title('40,000-50,000-years'); 
text(163.5,-34,'E', 'FontSize', 15, 'FontWeight', 'bold');yticks([-45 -40 -35]); xticks([165 170 175 180]);

nexttile(6); % 2500year return period
hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, c_10000_2500); colormap(gca,'parula'); set(gca,'ColorScale','log', 'CLim', [0.1 10]); 
set(gca, 'FontSize', 10); title('50,000-60,000-years'); 
text(163.5,-34,'F', 'FontSize', 15, 'FontWeight', 'bold');yticks([-45 -40 -35]); xticks([165 170 175 180]);
c = colorbar('FontSize', 10,'TickLabels', {'0.1','1','10'}); set(gca,'ColorScale','log', 'CLim', [0.1 10]);
c.Label.String = 'Wave height (m)'; pos = get(c,'Position'); c.Label.Position = [2 pos(2)+1]; 

nexttile(7); % 2500year return period
hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, a_10000_2500 - b_10000_2500); colormap(gca,mycolormap); set(gca, 'CLim', [-2 2]);
set(gca, 'FontSize', 10); title('D-E difference'); 
text(163.5,-34,'G', 'FontSize', 15, 'FontWeight', 'bold');yticks([-45 -40 -35]); xticks([165 170 175 180]);

nexttile(8); % 2500year return period
hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, a_10000_2500 - c_10000_2500); colormap(gca,mycolormap); set(gca, 'CLim', [-2 2]);
set(gca, 'FontSize', 10); title('D-F difference'); 
text(163.5,-34,'H', 'FontSize', 15, 'FontWeight', 'bold');yticks([-45 -40 -35]); xticks([165 170 175 180]);

nexttile(9); % 2500year return period
hold on;xlim([165 180]); ylim([-48 -34]);
scatter(coastal_x,coastal_y, 5, b_10000_2500 - c_10000_2500); colormap(gca,mycolormap); set(gca, 'CLim', [-2 2]);
set(gca, 'FontSize', 10); title('E-F difference'); 
text(163.5,-34,'I', 'FontSize', 15, 'FontWeight', 'bold');yticks([-45 -40 -35]); xticks([165 170 175 180]);
d = colorbar('FontSize', 10); d.Label.String = 'Difference (m)'; pos = get(d,'Position'); d.Label.Position = [2 pos(2)+0.5]; 


set(gcf,'position',[510,214,744,746]) 

fig = gcf;
exportgraphics(fig, 'tsunami_hazard_data_subsets_27_04_2023.png', 'BackgroundColor','white')


%%-------------------------%%


