% Load the coastline data
coastal_points = load('coastal_points.mat');
coastal_x = coastal_points.coastal_points(:,1); coastal_y = coastal_points.coastal_points(:,2);
%%-----------------------------------------------------------------------%%
% Return periods for our study:
bruce2_tsunami_a = load('coastal_values1.mat'); bruce2_tsunami_b = load('coastal_values2.mat');
bruce2_tsunami_c = load('coastal_values3.mat');bruce2_tsunami_d = load('coastal_values4.mat');
bruce2_tsunami_e = load('coastal_values5.mat'); bruce2_tsunami_f = load('coastal_values6.mat');
bruce2_waveheight = [bruce2_tsunami_a.wave_height_data bruce2_tsunami_b.wave_height_data2 ...
    bruce2_tsunami_c.wave_height_data bruce2_tsunami_d.wave_height_data ...
    bruce2_tsunami_e.wave_height_data bruce2_tsunami_f.wave_height_data];

% Return periods
t = maxk(bruce2_waveheight, 12, 2); min_exceed = mink(t, 1,2); waveheight_2500 = round(min_exceed, 2);

%%-----------------------------------------------------------------------%%
% NTHM
NTHM_2021 = csvread('LocalOnly2500yr50pct.csv', 1,0);
%%-----------------------------------------------------------------------%%
% Load in the zone information and manipulate it to get a set of variables
zones = readtable('zones_small_name_coord.txt'); zones_new = removevars(zones, {'Var6', 'Var7'});
temp_zones_matrix = table2array(zones_new); temp_zones_matrix(37462, :) = [];
zones_matrix = temp_zones_matrix;


NaN_rows = find(all(isnan(zones_matrix),2)); zone_number = [];
% Isolate zone names
for ii = 1:length(NaN_rows)-1
    row_value = NaN_rows(ii)+1;
    zone_num = zones_matrix(row_value, 5);
    zone_number = [zone_number, zone_num];
end
zone_numbers = [42, zone_number];
% Remove unnecessary information
for ii = 1:length(NaN_rows)-1
   x = find(all(isnan(zones_matrix), 2)); zones_matrix(x(ii)+1, :) = [];
end
%Make a series of variables for the zones
for ii = 1:length(x)
   zone_num =  zone_numbers(ii);
   name = ['zone', num2str(zone_num), '_coord'];
   if ii == 1
      lat_lon = zones_matrix(1:x(ii)-1, 1:2); 
   else
      lat_lon = zones_matrix(x(ii-1)+1:x(ii)-1, 1:2); 
   end
   assignin('base',name,lat_lon)
end
%%-----------------------------------------------------------------------%%
% Find the datapoints that fit into each of the zones
waveheight = waveheight_100;

for zone = 1:252
    name = ['zone', num2str(zone), '_coord'];
    disp(['zone', num2str(zone)])
    eval(['d =', name, ';']);
    temp = [];
    for point = 1:length(waveheight)
        x_point = coastal_x(point); y_point = coastal_y(point); 
        if inpolygon(x_point, y_point, d(:, 1), d(:, 2)) == 1
           temp = [temp, waveheight(point)];
        end    
    end
    data_name = ['zone', num2str(zone), '_data'];
    assignin('base', data_name, temp)
end

waveheight = [];
for zone = 1:252
    disp(zone)
    name2 = ['zone', num2str(zone), '_data']; eval(['data =', name2, ';']);
    wave_height = quantile(data, 0.99);
    waveheight = [waveheight, wave_height];
end
%%-----------------------------------------------------------------------%%
% 2500 year return period hazard
figure
tlo = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; 
s = pcolor(lon_points,lat_points, Z); shading flat; s.FaceAlpha = 0.2; 
set(gca,'ColorScale','log', 'CLim', [1 10]); colormap gray
hold on
for zone = 1:252
    name = ['zone', num2str(zone), '_coord']; eval(['d =', name, ';']);
    wave_height = NTSH_2021(zone, 5);
    if wave_height < 2
        colour_poly = '#3AFF3D';
    elseif (wave_height > 2) && (wave_height < 4)
        colour_poly = '#4BFFFE';
    elseif (wave_height > 4) && (wave_height < 6)
        colour_poly = '#2F00F9';
    elseif (wave_height > 6) && (wave_height < 8)
        colour_poly = '#FCFF43';        
    elseif (wave_height > 8) && (wave_height < 10)
         colour_poly = '#F7001A';  
    elseif (wave_height > 10) && (wave_height < 12)
         colour_poly = '#FA00FA';  
    else
        colour_poly = '#000000';
    end
    plot(x,y, 'black', 'LineWidth',0.5)
    plot(polyshape(d(:,1),d(:,2)),'FaceColor',colour_poly, 'FaceAlpha',1)
end
xlim([165 180]); ylim([-48 -34]); title('Tsunami hazard: 2021 NTHM');

nexttile;
s = pcolor(lon_points,lat_points, Z); shading flat; s.FaceAlpha = 0.2; 
set(gca,'ColorScale','log', 'CLim', [1 10]); colormap gray
hold on
for zone = 1:252
    name = ['zone', num2str(zone), '_coord']; eval(['d =', name, ';']);
    wave_height = waveheight_2500(zone);
    if wave_height < 2
        colour_poly = '#3AFF3D';
    elseif (wave_height > 2) && (wave_height < 4)
        colour_poly = '#4BFFFE';
    elseif (wave_height > 4) && (wave_height < 6)
        colour_poly = '#2F00F9';
    elseif (wave_height > 6) && (wave_height < 8)
        colour_poly = '#FCFF43';        
    elseif (wave_height > 8) && (wave_height < 10)
         colour_poly = '#F7001A';  
    elseif (wave_height > 10) && (wave_height < 12)
         colour_poly = '#FA00FA';  
    else
        colour_poly = '#000000';
    end
    plot(x,y, 'black', 'LineWidth',0.5)
    plot(polyshape(d(:,1),d(:,2)),'FaceColor',colour_poly, 'FaceAlpha',1)
end
plot(polyshape([176.8 -42; 177.6 -42; 177.6 -42.5; 176.8 -42.5]),'FaceColor','#3AFF3D', 'FaceAlpha',1)
text(177.8,-42.2,'0-2m');
plot(polyshape([176.8 -42.7; 177.6 -42.7; 177.6 -43.2; 176.8 -43.2]),'FaceColor','#4BFFFE', 'FaceAlpha',1)
text(177.8,-42.9,'2-4m');
plot(polyshape([176.8 -43.4; 177.6 -43.4; 177.6 -43.9; 176.8 -43.9]),'FaceColor','#2F00F9', 'FaceAlpha',1)
text(177.8,-43.6,'4-6m');
plot(polyshape([176.8 -44.1; 177.6 -44.1; 177.6 -44.6; 176.8 -44.6]),'FaceColor','#FCFF43', 'FaceAlpha',1)
text(177.8,-44.3,'6-8m');
plot(polyshape([176.8 -44.8; 177.6 -44.8; 177.6 -45.3; 176.8 -45.3]),'FaceColor','#F7001A', 'FaceAlpha',1)
text(177.8,-45.0,'8-10m');
plot(polyshape([176.8 -45.5; 177.6 -45.5; 177.6 -46; 176.8 -46]),'FaceColor','#FA00FA', 'FaceAlpha',1)
text(177.8,-45.7,'10-12m');
plot(polyshape([176.8 -46.2; 177.6 -46.2; 177.6 -46.7; 176.8 -46.7]),'FaceColor','#000000', 'FaceAlpha',1)
text(177.8,-46.5,'12+m');
hold off
xlim([165 180]); ylim([-48 -34]); title('Tsunami hazard: this study');
sgtitle('2,500-year return period')

%%-----------------------------------------------------------------------%%
