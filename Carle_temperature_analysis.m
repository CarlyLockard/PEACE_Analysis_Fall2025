%% Code to analyze T2 and T1 maps acquired 8/21/2025 at 1.5 and 3T with
% 1 phantom per field strength, with phantom refrigerated 24 hours
% preceding scan and then imaged repeatedly without moving as returned to
% room temperature (based on temperature measurements of a water-filled 50
% mL conical centrifuge tubes that was placed adjaceent to the phantom). 

%% 1.5T
clear;
clc;
close all;

main_folder = 'S:\PHTx\Phantoms\Final phantoms\Carle scans and results\2023_FinalPhantoms_temperature\P22_1pt5TTemperatureTest';
cd(main_folder);

n_scans = 5;

for scan_number = 1:n_scans
    t1_long_series_folder = strcat(num2str(2+(scan_number-1)*9),'_',...
                                      't1map_long_t1_',num2str(scan_number),...
                                      '_Orig_MAPGRAY');
    cd(t1_long_series_folder);
    loop_file = [dir('*.DCM'),dir('*.IMA')];
    for i = 1:size(loop_file,1)
        filename_loop = string(loop_file(i).name);
        T1_long_series(:,:,i,scan_number) = dicomread(filename_loop);
    end
    cd(main_folder);

    t1_long_map_folder = strcat(num2str(3+(scan_number-1)*9),'_',...
                                      't1map_long_t1_',num2str(scan_number),...
                                      '_T1_Deriv_MAPGRAY');
    cd(t1_long_map_folder);
    loop_file = [dir('*.DCM'),dir('*.IMA')];
    for i = 1:size(loop_file,1)
        filename_loop = string(loop_file(i).name);
        T1_long_map(:,:,i,scan_number) = dicomread(filename_loop);
    end
    
    scan_info_from_DICOM = dicominfo(filename_loop);
    T1_info.resolution_in_plane = scan_info_from_DICOM.PerFrameFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing;
    cd(main_folder);

    t1_short_series_folder = strcat(num2str(5+(scan_number-1)*9),'_',...
                                      't1map_short_t1_',num2str(scan_number),...
                                      '_Orig_MAPGRAY');
    cd(t1_short_series_folder);
    loop_file = [dir('*.DCM'),dir('*.IMA')];
    for i = 1:size(loop_file,1)
        filename_loop = string(loop_file(i).name);
        T1_short_series(:,:,i,scan_number) = dicomread(filename_loop);
    end
    cd(main_folder);
    
    t1_short_map_folder = strcat(num2str(6+(scan_number-1)*9),'_',...
                                      't1map_short_t1_',num2str(scan_number),...
                                      '_T1_Deriv_MAPGRAY');
    cd(t1_short_map_folder);
    loop_file = [dir('*.DCM'),dir('*.IMA')];
    for i = 1:size(loop_file,1)
        filename_loop = string(loop_file(i).name);
        T1_short_map(:,:,i,scan_number) = dicomread(filename_loop);
    end
    scan_info_from_DICOM = dicominfo(filename_loop);
    T1short_info.resolution_in_plane = scan_info_from_DICOM.PerFrameFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing;
    cd(main_folder);
    
    t2_series_folder = strcat(num2str(8+(scan_number-1)*9),'_',...
                                      't2map_trufi_',num2str(scan_number),...
                                      '_Orig_MAPGRAY');
    cd(t2_series_folder);
    loop_file = [dir('*.DCM'),dir('*.IMA')];
    for i = 1:size(loop_file,1)
        filename_loop = string(loop_file(i).name);
        T2_series(:,:,i,scan_number) = dicomread(filename_loop);
    end
    cd(main_folder);
    
    t2_map_folder = strcat(num2str(9+(scan_number-1)*9),'_',...
                                      't2map_trufi_',num2str(scan_number),...
                                      '_T2_Deriv_MAPGRAY');
    cd(t2_map_folder);
    loop_file = [dir('*.DCM'),dir('*.IMA')];
    for i = 1:size(loop_file,1)
        filename_loop = string(loop_file(i).name);
        T2_map(:,:,i,scan_number) = dicomread(filename_loop);
    end
    scan_info_from_DICOM = dicominfo(filename_loop);
    T2_info.resolution_in_plane = scan_info_from_DICOM.PerFrameFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing;
    cd(main_folder);

end


% mask out loading solution bottles in selected series images 
for n = 1:n_scans
    %update below to pull out each individual scan for masking
    T1_matrix_r_c_s = squeeze(T1_long_series(:,:,:,n));
    T1short_matrix_r_c_s = squeeze(T1_short_series(:,:,:,n));
    T2_matrix_r_c_s = squeeze(T2_series(:,:,:,n));

    T1_map_matrix_r_c_s = squeeze(T1_long_map(:,:,:,n));
    T1short_map_matrix_r_c_s = squeeze(T1_short_map(:,:,:,n));
    T2_map_matrix_r_c_s = squeeze(T2_map(:,:,:,n));

    masked_images_T1long = mask_loading_phantom(T1_matrix_r_c_s,T1_info.resolution_in_plane); % mask out loading solution on T1-long 
    masked_images_T1short = mask_loading_phantom(T1short_matrix_r_c_s,T1short_info.resolution_in_plane); % mask out loading solution on T1-short
    masked_images_T2 = mask_loading_phantom(T2_matrix_r_c_s,T2_info.resolution_in_plane); % mask out loading solution on T2
    
    % create maximum intensity images to make finding the circular ROIs
    % most straightforward and to minimize variability between different
    % site protocols (e.g. different numbers of images or inversion times,
    % etc. leading to different optimum ROI-finding image numbers)
    masked_image_T1long = max(masked_images_T1long,[],3); 
    masked_image_T1short = max(masked_images_T1short,[],3);
    masked_image_T2 = max(masked_images_T2,[],3);
    
    % display to check visually
    figure; tiledlayout(1,3);
    nexttile;imshow(masked_image_T1long,[]);title('Loading soln masked - T1long');
    nexttile;imshow(masked_image_T1short,[]);title('Loading soln masked - T1short');
    nexttile;imshow(masked_image_T2,[]);title('Loading soln masked - T2');
    
    % now use the masked images from above for detection and localization
    % of the phantom tubes
    mm_radius_to_use = 6.5; % set radius to approximately match manual measurements performed at Vanderbilt in Medis software
    % extract values for each timepoint
    [timepoint_mean, timepoint_stdev, timepoint_median, timepoint_iqr] = find_rois(masked_image_T1long,T1_map_matrix_r_c_s,T1_info.resolution_in_plane,mm_radius_to_use);
    % organize timepoint values into an all-timepoints matrix
    T1long_means_all(:,n) =  timepoint_mean; 
    T1long_stdevs_all(:,n) =  timepoint_stdev; 
    T1long_medians_all(:,n) =  timepoint_median; 
    T1long_iqr_all(:,n) =  timepoint_iqr; 
    
    clearvars timepoint_mean timepoint_stdev timepoint_median timepoint_iqr 
    
    [timepoint_mean, timepoint_stdev, timepoint_median, timepoint_iqr] = find_rois(masked_image_T1short,T1short_map_matrix_r_c_s,T1short_info.resolution_in_plane,mm_radius_to_use);
    T1short_means_all(:,n) =  timepoint_mean; 
    T1short_stdevs_all(:,n) =  timepoint_stdev; 
    T1short_medians_all(:,n) =  timepoint_median; 
    T1short_iqr_all(:,n) =  timepoint_iqr; 
    
    clearvars timepoint_mean timepoint_stdev timepoint_median timepoint_iqr
    if ~isfield(scan_info_from_DICOM.PerFrameFunctionalGroupsSequence.Item_1.PixelValueTransformationSequence.Item_1,'RescaleIntercept')
        T2_map_matrix_r_c_s = T2_map_matrix_r_c_s/10; disp("rescale slope not applying;correcting in subsequent step")
    else
        T2_map_matrix_r_c_s = T2_map_matrix_r_c_s*scan_info_from_DICOM.PerFrameFunctionalGroupsSequence.Item_1.PixelValueTransformationSequence.Item_1.RescaleSlope;
    end
    [timepoint_mean, timepoint_stdev, timepoint_median, timepoint_iqr] = find_rois(masked_image_T2,T2_map_matrix_r_c_s,T2_info.resolution_in_plane,mm_radius_to_use);
    T2_means_all(:,n) =  timepoint_mean;
    T2_stdevs_all(:,n) =  timepoint_stdev; 
    T2_medians_all(:,n) =  timepoint_median; 
    T2_iqr_all(:,n) =  timepoint_iqr;
    
    clearvars timepoint_mean timepoint_stdev timepoint_median timepoint_iqr
end

%%
close all;

% Temperature array
temperature_array = [13.8, 16.1, 17.0, 18.0, 18.8];

% T1long, T1short, T2 for each concentration as f(temperature)
h1 = figure(1);
t = tiledlayout(1,6);
title('T_1 native')
for i = 1:6
    ax1 = nexttile;
    to_plot = squeeze(T1long_means_all(i,:));
    scatter(temperature_array,to_plot);hold on;
    plot([temperature_array(1) temperature_array(length(temperature_array))],[to_plot(1) to_plot(length(to_plot))],'k--');
    title(ax1,strcat('Tube ',num2str(i)));
    ylabel('Relaxation time [ms]');
    xlabel('Temperature [^{\circ} C]','Interpreter','tex');
    ylim([150 1600]);
    xlim([13 19]);
end

h2 = figure(2);
t2 = tiledlayout(1,6);
title('T_1 enhanced')
for i = 1:6
    ax1 = nexttile;
    to_plot = squeeze(T1short_means_all(i,:));
    scatter(temperature_array,to_plot);hold on;
    plot([temperature_array(1) temperature_array(length(temperature_array))],[to_plot(1) to_plot(length(to_plot))],'k--');
    title(ax1,strcat('Tube ',num2str(i)));
    ylabel('Relaxation time [ms]');
    xlabel('Temperature [^{\circ} C]','Interpreter','tex');
    ylim([150 1400]);
    xlim([13 19]);
end

h3 = figure(3);
t3 = tiledlayout(1,6);
title('T_2')
for i = 1:6
    ax1 = nexttile;
    to_plot = squeeze(T2_means_all(i,:));
    scatter(temperature_array,to_plot);hold on;
    plot([temperature_array(1) temperature_array(length(temperature_array))],[to_plot(1) to_plot(length(to_plot))],'k--');
    title(ax1,strcat('Tube ',num2str(i)));
    ylabel('Relaxation time [ms]');
    xlabel('Temperature [^{\circ} C]','Interpreter','tex');
    ylim([40 300]);
    xlim([13 19]);
end

cd('S:\PHTx\Phantoms\Final phantoms\Carle scans and results\2023_FinalPhantoms_temperature');
save('1pt5T_measurements.mat','temperature_array','T1long_means_all','T1short_means_all','T2_means_all');
%% 3T
clear;
clc;
close all;

main_folder = 'S:\PHTx\Phantoms\Final phantoms\Carle scans and results\2023_FinalPhantoms_temperature\P23_3TTemperatureTest';
cd(main_folder);

n_scans = 8; %scan 8 was incomplete so deleted and renumbered; removed that temperature from array in section below

for scan_number = 1:n_scans
    t1_long_series_folder = strcat(num2str(3+(scan_number-1)*12),'_',...
                                      '8_150cube_t1map_long_t1_50_',num2str(scan_number),...
                                      '_Orig_MAPGRAY');
    cd(t1_long_series_folder);
    loop_file = [dir('*.DCM'),dir('*.IMA')];
    for i = 1:size(loop_file,1)
        filename_loop = string(loop_file(i).name);
        T1_long_series(:,:,i,scan_number) = dicomread(filename_loop);
    end
    cd(main_folder);

    t1_long_map_folder = strcat(num2str(4+(scan_number-1)*12),'_',...
                                      '8_150cube_t1map_long_t1_50_',num2str(scan_number),...
                                      '_T1_Deriv_MAPGRAY');
    cd(t1_long_map_folder);
    loop_file = [dir('*.DCM'),dir('*.IMA')];
    for i = 1:size(loop_file,1)
        filename_loop = string(loop_file(i).name);
        T1_long_map(:,:,i,scan_number) = dicomread(filename_loop);
    end
    
    scan_info_from_DICOM = dicominfo(filename_loop);
    T1_info.resolution_in_plane = scan_info_from_DICOM.PerFrameFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing;
    cd(main_folder);

    t1_short_series_folder = strcat(num2str(6+(scan_number-1)*12),'_',...
                                      '8_150cube_t1map_short_t1_50_',num2str(scan_number),...
                                      '_Orig_MAPGRAY');
    cd(t1_short_series_folder);
    loop_file = [dir('*.DCM'),dir('*.IMA')];
    for i = 1:size(loop_file,1)
        filename_loop = string(loop_file(i).name);
        T1_short_series(:,:,i,scan_number) = dicomread(filename_loop);
    end
    cd(main_folder);
    
    t1_short_map_folder = strcat(num2str(7+(scan_number-1)*12),'_',...
                                      '8_150cube_t1map_short_t1_50_',num2str(scan_number),...
                                      '_T1_Deriv_MAPGRAY');
    cd(t1_short_map_folder);
    loop_file = [dir('*.DCM'),dir('*.IMA')];
    for i = 1:size(loop_file,1)
        filename_loop = string(loop_file(i).name);
        T1_short_map(:,:,i,scan_number) = dicomread(filename_loop);
    end
    scan_info_from_DICOM = dicominfo(filename_loop);
    T1short_info.resolution_in_plane = scan_info_from_DICOM.PerFrameFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing;
    cd(main_folder);
    
    t2_series_folder = strcat(num2str(9+(scan_number-1)*12),'_',...
                                      '8_150cube_t2map_flash_50_',num2str(scan_number),...
                                      '_Orig_MAPGRAY');
    cd(t2_series_folder);
    loop_file = [dir('*.DCM'),dir('*.IMA')];
    for i = 1:size(loop_file,1)
        filename_loop = string(loop_file(i).name);
        T2_series(:,:,i,scan_number) = dicomread(filename_loop);
    end
    cd(main_folder);
    
    t2_map_folder = strcat(num2str(10+(scan_number-1)*12),'_',...
                                      '8_150cube_t2map_flash_50_',num2str(scan_number),...
                                      '_T2_Deriv_MAPGRAY');
    cd(t2_map_folder);
    loop_file = [dir('*.DCM'),dir('*.IMA')];
    for i = 1:size(loop_file,1)
        filename_loop = string(loop_file(i).name);
        T2_map(:,:,i,scan_number) = dicomread(filename_loop);
    end
    scan_info_from_DICOM = dicominfo(filename_loop);
    T2_info.resolution_in_plane = scan_info_from_DICOM.PerFrameFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing;
    cd(main_folder);

end


% mask out loading solution bottles in selected series images 
for n = 1:n_scans
    %update below to pull out each individual scan for masking
    T1_matrix_r_c_s = squeeze(T1_long_series(:,:,:,n));
    T1short_matrix_r_c_s = squeeze(T1_short_series(:,:,:,n));
    T2_matrix_r_c_s = squeeze(T2_series(:,:,:,n));

    T1_map_matrix_r_c_s = squeeze(T1_long_map(:,:,:,n));
    T1short_map_matrix_r_c_s = squeeze(T1_short_map(:,:,:,n));
    T2_map_matrix_r_c_s = squeeze(T2_map(:,:,:,n));

    masked_images_T1long = mask_loading_phantom(T1_matrix_r_c_s,T1_info.resolution_in_plane); % mask out loading solution on T1-long 
    masked_images_T1short = mask_loading_phantom(T1short_matrix_r_c_s,T1short_info.resolution_in_plane); % mask out loading solution on T1-short
    masked_images_T2 = mask_loading_phantom(T2_matrix_r_c_s,T2_info.resolution_in_plane); % mask out loading solution on T2
    
    % create maximum intensity images to make finding the circular ROIs
    % most straightforward and to minimize variability between different
    % site protocols (e.g. different numbers of images or inversion times,
    % etc. leading to different optimum ROI-finding image numbers)
    masked_image_T1long = max(masked_images_T1long,[],3); 
    masked_image_T1short = max(masked_images_T1short,[],3);
    masked_image_T2 = max(masked_images_T2,[],3);
    
    % display to check visually
    figure; tiledlayout(1,3);
    nexttile;imshow(masked_image_T1long,[]);title('Loading soln masked - T1long');
    nexttile;imshow(masked_image_T1short,[]);title('Loading soln masked - T1short');
    nexttile;imshow(masked_image_T2,[]);title('Loading soln masked - T2');
    
    % now use the masked images from above for detection and localization
    % of the phantom tubes
    mm_radius_to_use = 6.5; % set radius to approximately match manual measurements performed at Vanderbilt in Medis software
    % extract values for each timepoint
    [timepoint_mean, timepoint_stdev, timepoint_median, timepoint_iqr] = find_rois(masked_image_T1long,T1_map_matrix_r_c_s,T1_info.resolution_in_plane,mm_radius_to_use);
    % organize timepoint values into an all-timepoints matrix
    T1long_means_all(:,n) =  timepoint_mean; 
    T1long_stdevs_all(:,n) =  timepoint_stdev; 
    T1long_medians_all(:,n) =  timepoint_median; 
    T1long_iqr_all(:,n) =  timepoint_iqr; 
    
    clearvars timepoint_mean timepoint_stdev timepoint_median timepoint_iqr 
    
    [timepoint_mean, timepoint_stdev, timepoint_median, timepoint_iqr] = find_rois(masked_image_T1short,T1short_map_matrix_r_c_s,T1short_info.resolution_in_plane,mm_radius_to_use);
    T1short_means_all(:,n) =  timepoint_mean; 
    T1short_stdevs_all(:,n) =  timepoint_stdev; 
    T1short_medians_all(:,n) =  timepoint_median; 
    T1short_iqr_all(:,n) =  timepoint_iqr; 
    
    clearvars timepoint_mean timepoint_stdev timepoint_median timepoint_iqr
    if ~isfield(scan_info_from_DICOM.PerFrameFunctionalGroupsSequence.Item_1.PixelValueTransformationSequence.Item_1,'RescaleIntercept')
        T2_map_matrix_r_c_s = T2_map_matrix_r_c_s/10; disp("rescale slope not applying;correcting in subsequent step")
    else
        T2_map_matrix_r_c_s = T2_map_matrix_r_c_s*scan_info_from_DICOM.PerFrameFunctionalGroupsSequence.Item_1.PixelValueTransformationSequence.Item_1.RescaleSlope;
    end
    [timepoint_mean, timepoint_stdev, timepoint_median, timepoint_iqr] = find_rois(masked_image_T2,T2_map_matrix_r_c_s,T2_info.resolution_in_plane,mm_radius_to_use);
    T2_means_all(:,n) =  timepoint_mean;
    T2_stdevs_all(:,n) =  timepoint_stdev; 
    T2_medians_all(:,n) =  timepoint_median; 
    T2_iqr_all(:,n) =  timepoint_iqr;
    
    clearvars timepoint_mean timepoint_stdev timepoint_median timepoint_iqr
end

%% Plot values vs. temperature
close all;

% Temperature array
temperature_array = [9.4 12.7 14.8 16.4 17.8 18.9 19.9 21.1]; %scan 8 was incomplete so deleted

% T1long, T1short, T2 for each concentration as f(temperature)
h1 = figure(1);
t = tiledlayout(1,6);
title('T_1 native')
for i = 1:6
    ax1 = nexttile;
    to_plot = squeeze(T1long_means_all(i,:));
    scatter(temperature_array,to_plot);hold on;
    plot([temperature_array(1) temperature_array(length(temperature_array))],[to_plot(1) to_plot(length(to_plot))],'k--');
    title(ax1,strcat('Tube ',num2str(i)));
    ylabel('Relaxation time [ms]');
    xlabel('Temperature [^{\circ} C]','Interpreter','tex');
    ylim([150 1700]);
    xlim([9 22]);
end

h2 = figure(2);
t2 = tiledlayout(1,6);
title('T_1 enhanced')
for i = 1:6
    ax1 = nexttile;
    to_plot = squeeze(T1short_means_all(i,:));
    scatter(temperature_array,to_plot);hold on;
    plot([temperature_array(1) temperature_array(length(temperature_array))],[to_plot(1) to_plot(length(to_plot))],'k--');
    title(ax1,strcat('Tube ',num2str(i)));
    ylabel('Relaxation time [ms]');
    xlabel('Temperature [^{\circ} C]','Interpreter','tex');
    ylim([150 1600]);
    xlim([9 22]);
end

h3 = figure(3);
t3 = tiledlayout(1,6);
title('T_2')
for i = 1:6
    ax1 = nexttile;
    to_plot = squeeze(T2_means_all(i,:));
    scatter(temperature_array,to_plot);hold on;
    plot([temperature_array(1) temperature_array(length(temperature_array))],[to_plot(1) to_plot(length(to_plot))],'k--');
    title(ax1,strcat('Tube ',num2str(i)));
    ylabel('Relaxation time [ms]');
    xlabel('Temperature [^{\circ} C]','Interpreter','tex');
    ylim([15 200]);
    xlim([9 22]);
end

cd('S:\PHTx\Phantoms\Final phantoms\Carle scans and results\2023_FinalPhantoms_temperature');
save('3T_measurements.mat','temperature_array','T1long_means_all','T1short_means_all','T2_means_all');

%%
% Quantify relationship for each of the 6 tubes at 1.5T and 3T for T1-long, T1-short,
% and T2 (linear relationship slope and intercept, and R^2). 
% For each measure and tube, plot the points and the fitted line to do a visual check.
% Using the line equation, predict for each field strength, tube, and measure, the expected min and max values for typical range of
% MRI scanner room temperatures (predict at min/max temperatures in range for 18 - 21 C, 20 - 22 C, and for wider range 18 - 25 C)
% [skyra_planning_guide.pdf, https://www.ismrm.org/smrt/files/con2033065.pdf]
%
% See MATLAB linear regression guide at https://www.mathworks.com/help/matlab/data_analysis/linear-regression.html

clear;
close all;
clc;

cd('S:\PHTx\Phantoms\Final phantoms\Carle scans and results\2023_FinalPhantoms_temperature');
%% For 1.5T
load('1pt5T_measurements.mat'); 

% add code to do the items above

%table for slope, y-intercept, and r^2 for T1 long means
quant_t1long_115t = quantify(T1long_means_all, temperature_array, "T1 long")
quant_t1short_115t = quantify(T1short_means_all, temperature_array, "T1 short")
quant_t2_115t = quantify(T2_means_all, temperature_array, "T2")

%% For 3T
load('3T_measurements.mat');

% add code to do the items above
quant_t1long_3t = quantify(T1long_means_all, temperature_array, "T1 long")
quant_t1short_3t = quantify(T1short_means_all, temperature_array, "T1 short")
quant_t2_3t = quantify(T2_means_all, temperature_array, "T2")
%%
function table_all = quantify(file_name, file_name2, ylab)
    tube = [];
    slope = [];
    y_intercept = [];
    Rsq2 = [];
    min_max1821 = [];
    min_max2022 = [];
    min_max1825 = [];
    for c = 1:6
        temp = file_name2';
        t1_long1 = file_name(c, :)';
     
        X = [ones(length(temp),1) temp];
        slope1 = X\t1_long1;
        slope2 = slope1(2,1);
        y_int = slope1(1,1);
        
        x_calc = (slope2*X(:, 2))+y_int;
        
        r2 = 1 - sum((t1_long1 - x_calc).^2)/sum((t1_long1 - mean(t1_long1)).^2);
        

        scatter(temp, t1_long1)
        plot(X(:, 2), x_calc)
        hold on;
        legend("Tube 1", "", "Tube 2", "", "Tube 3", "", "Tube 4", "", "Tube 5", "", "Tube 6", "", "Location", "northeastoutside")
        xlabel("temperature")
        ylabel(ylab)

    
        min_18 = (slope2*18) + y_int;
        max_21 = (slope2*21) + y_int;
        min_20 = (slope2*20) + y_int;
        max_22 = (slope2*22) + y_int;
        max_25 = (slope2*25) + y_int;
    
        tube = [tube; c];
        slope = [slope; slope2];
        y_intercept = [y_intercept; y_int];
        Rsq2 = [Rsq2; r2];
        min_max1821 = [min_max1821; min_18, max_21];
        min_max2022 = [min_max2022; min_20, max_22];
        min_max1825 = [min_max1825; min_18, max_25];
    end
    hold off;
    table_all = table(tube, slope, y_intercept, Rsq2, min_max1821, min_max2022, min_max1825)
end
