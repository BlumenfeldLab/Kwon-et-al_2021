%% ROI
clc;clear;
close all;

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

load r1_all.mat
fid=fopen('Subjects_list_all.txt','r');
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);

data_location='E:\RAM data set\RAM_Public_Data_all\FR1_final';
Montage_suffix='final_clean_dil';

mni_nifti_path = 'C:\yale\bioimagesuite30\images\MNI_T1_1mm_stripped.nii';
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
Tmni = mnit1_1mm_stripped.transform;
ROI_list=[84 121.3 109.6 ;
65.57 74.26 122.6;
27.93 114 56.79;
96 68.79 101.2;
148.9 105.6 57.95;
84.63 77.91 121.1;
96.04 97.36 124.8;
113.8 127.1 64.86;
61.08 115.4 59.8;
84.71 113.3 133.6;
44.03 136.3 118.6;
134.1 162.9 66.27;
100.4 166.1 79.38
146.3 134 84.63;
138.2 82.69 44.46;
140.2 80.06 82.72;
146 101.9 104.2;
150.2 114.3 66.33;
84.33 80.78 127.5;
82.02 53.45 98.45];

ROI_list=Tmni*[ROI_list ones(size(ROI_list,1),1)]';
ROI_list=ROI_list';
ROI_list=ROI_list(:,1:3);

radius=15;

for i=1:251
    try
        cd(data_location)
        cd([num2str(i),'_',r_sublist{i,1}]);
        
        for ROI_i=14:size(ROI_list,1)
            Buffer_MontageMap=[];
            L_MontageMap=[];
            R_MontageMap=[];
            if ROI_list(ROI_i,1) < 0
                side = 'L';
                load([side '_MontageMap_' Montage_suffix '.mat'])
                if ~isempty(L_MontageMap)
                    Buffer_MontageMap = L_MontageMap(:,2:4);
                    % calculate distance
                    Dist_MontageMap=[];
                    Dist_MontageMap=sqrt(sum((Buffer_MontageMap-ROI_list(ROI_i,:)).^2,2));
                    L_MontageMap(~(Dist_MontageMap<radius),:)=[];
                end
            else
                side = 'R';
                load([side '_MontageMap_' Montage_suffix '.mat'])
                if ~isempty(R_MontageMap)
                    Buffer_MontageMap = R_MontageMap(:,2:4);
                    % calculate distance
                    Dist_MontageMap=[];
                    Dist_MontageMap=sqrt(sum((Buffer_MontageMap-ROI_list(ROI_i,:)).^2,2));
                    R_MontageMap(~(Dist_MontageMap<radius),:)=[];
                end
            end
            save(['L_MontageMap_' Montage_suffix '_roi_' num2str(ROI_i) '.mat'], 'L_MontageMap');
            save(['R_MontageMap_' Montage_suffix '_roi_' num2str(ROI_i) '.mat'], 'R_MontageMap');
        end
        disp(sprintf('%s is completed! ',r_sublist{i,1}));
    catch
    end
end


patients=r_sublist;
overlay=4;
laterality=0;
views=4;
inflationstep=5;
electrodeDensity=0;
electrodeDisplay=1;
thr=2;
mean_flag=2;
color_range=[-6 8];
statistic_flag=0;

file_suffix='Wendy_merge';
for ROI_i=1:size(ROI_list,1)
    close all;
    Montage_suffix=['final_clean_dil_roi_' num2str(ROI_i) ];
    displayElectrodesInflated_new_ROI(data_location, patients, overlay, laterality, views, inflationstep,electrodeDensity, electrodeDisplay ,thr,mean_flag,color_range, file_suffix,Montage_suffix,statistic_flag) %, frames_set)
end



%% plot 
clc;clear;
close all;

mni_nifti_path = 'C:\yale\bioimagesuite30\images\MNI_T1_1mm_stripped.nii';
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
Tmni = mnit1_1mm_stripped.transform;
ROI_list=[84 121.3 109.6 ;
65.57 74.26 122.6;
27.93 114 56.79;
96 68.79 101.2;
148.9 105.6 57.95;
84.63 77.91 121.1;
96.04 97.36 124.8;
113.8 127.1 64.86;
61.08 115.4 59.8;
84.71 113.3 133.6;
44.03 136.3 118.6;
134.1 162.9 66.27;
100.4 166.1 79.38
146.3 134 84.63;
138.2 82.69 44.46;
140.2 80.06 82.72;
146 101.9 104.2;
150.2 114.3 66.33;
84.33 80.78 127.5;
82.02 53.45 98.45];

ROI_list=Tmni*[ROI_list ones(size(ROI_list,1),1)]';
ROI_list=ROI_list';
ROI_list=ROI_list(:,1:3);

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

load r1_all.mat
fid=fopen('Subjects_list_all.txt','r');
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);

file_suffix='Wendy_merge';
data_location='E:\RAM data set\RAM_Public_Data_all\FR1_final';
image_storage_folder = [data_location '/Storage_filtering_no_overlap/Image_storage/' file_suffix '/ROIs'];
mkdir(image_storage_folder)

for ROI_i=1:size(ROI_list,1)
    all_zscore_traces=[];
    Montage_suffix=['final_clean_dil_roi_' num2str(ROI_i) ];
    for i=1:251
        try
            cd(data_location)
            cd([num2str(i),'_',r_sublist{i,1}]);
            load(['stats_traces_' file_suffix '.mat']);
            load(['rejection_idx_' file_suffix '.mat']);
            
            trial_thr_ind=find(RAM_all_trials_total_number-noisy_RAM_all_trials_number<100);
            zscore_traces(trial_thr_ind,1,3,:)=NaN;
            if ROI_list(ROI_i,1) < 0
                load(['L_MontageMap_' Montage_suffix '.mat']);
                if ~isempty(L_MontageMap)
                    if size(L_MontageMap,1)==1
                        buffer_zscore_traces=squeeze(zscore_traces(L_MontageMap(:,1),1,3,:));
                        all_zscore_traces=[all_zscore_traces ; buffer_zscore_traces'];
                    else
                        buffer_zscore_traces=squeeze(zscore_traces(L_MontageMap(:,1),1,3,:));
                        all_zscore_traces=[all_zscore_traces ; buffer_zscore_traces];
                    end
                end
            else
                load(['R_MontageMap_' Montage_suffix '.mat']);
                if ~isempty(R_MontageMap)
                    if size(R_MontageMap,1)==1
                        buffer_zscore_traces=squeeze(zscore_traces(R_MontageMap(:,1),1,3,:));
                        all_zscore_traces=[all_zscore_traces ; buffer_zscore_traces'];
                    else
                        buffer_zscore_traces=squeeze(zscore_traces(R_MontageMap(:,1),1,3,:));
                        all_zscore_traces=[all_zscore_traces ; buffer_zscore_traces];
                    end
                end
            end
            disp(sprintf('%s is completed! ',r_sublist{i,1}));
        catch
        end
    end
    
    figure('position',[0 0 8000 4000]);
    set(gcf, 'color', [1 1 1]);
    set(gcf,'Visible','on');  
    
    subplot(3,2,1)
    plot(all_zscore_traces');
    
    hold on;
    sroot_all_zscore_traces=nansum(all_zscore_traces,1)/sqrt(size(all_zscore_traces,1));
%     sroot_all_zscore_traces=nansum(all_zscore_traces,1)/size(all_zscore_traces,1);
%     sroot_all_zscore_traces=nanmean(all_zscore_traces,1);
    plot(sroot_all_zscore_traces,'LineWidth',2,'color','k');
    
    set(gca,'xtick',[5 10 15 20 25 30 35]);
    set(gca,'xticklabel',[-250 0 250 500 750 1000 1250]);
    set(gca,'xlim',[1 39])
    set(gca,'ylim',[-12 16])
    line([1,39],[2,2],'color','b');
    line([1,39],[-2,-2],'color','b');
    line([10,10],[-12,16],'color','k');
    line([11,11],[-12,16],'color','k');
    line([12,12],[-12,16],'color','k');
    line([13,13],[-12,16],'color','k');
    set(gca,'Fontsize',15);
    title(['All trials,     ROI index : ' num2str(ROI_i) ', Total electrode # : ' num2str(size(all_zscore_traces,1)) ],'Fontsize',12)
    ylabel('Baseline z-score value');
    xlabel('Time (ms)');
    
    subplot(3,2,2)
    mean_all_zscore_traces=nanmean(all_zscore_traces,1);
    plot(mean_all_zscore_traces,'LineWidth',2,'color','k');
    
    set(gca,'xtick',[5 10 15 20 25 30 35]);
    set(gca,'xticklabel',[-250 0 250 500 750 1000 1250]);
    set(gca,'xlim',[1 39])
    line([10,10],[min(mean_all_zscore_traces),max(mean_all_zscore_traces)],'color','k');
    line([11,11],[min(mean_all_zscore_traces),max(mean_all_zscore_traces)],'color','k');
    line([12,12],[min(mean_all_zscore_traces),max(mean_all_zscore_traces)],'color','k');
    line([13,13],[min(mean_all_zscore_traces),max(mean_all_zscore_traces)],'color','k');
    set(gca,'Fontsize',15);
    title('Mean z-score value','Fontsize',12)
    ylabel('Baseline z-score value');
    xlabel('Time (ms)');
    
    subplot(3,2,3)
    histogram(all_zscore_traces(:,10),'BinLimits',[-50,50]);
    title('Histogram of zscore values, 0 ~ 50ms','Fontsize',12)
    
    subplot(3,2,4)
    histogram(all_zscore_traces(:,11),'BinLimits',[-50,50]);
    title('Histogram of zscore values, 50ms ~ 100ms','Fontsize',12)
    
    subplot(3,2,5)
    histogram(all_zscore_traces(:,12),'BinLimits',[-50,50]);
    title('Histogram of zscore values, 100ms ~ 150ms','Fontsize',12)
    
    subplot(3,2,6)
    histogram(all_zscore_traces(:,13),'BinLimits',[-50,50]);
    title('Histogram of zscore values, 150ms ~ 200ms','Fontsize',12)
    
    cd(image_storage_folder)
    saveas(gcf,['zscore_plot_roi_' num2str(ROI_i) '.png'])
    close all;
    
end


%% plot single zscores
clc;clear;
close all;

mni_nifti_path = 'C:\yale\bioimagesuite30\images\MNI_T1_1mm_stripped.nii';
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
Tmni = mnit1_1mm_stripped.transform;
ROI_list=[84 121.3 109.6 ;
65.57 74.26 122.6;
27.93 114 56.79;
96 68.79 101.2;
148.9 105.6 57.95;
84.63 77.91 121.1;
96.04 97.36 124.8;
113.8 127.1 64.86;
61.08 115.4 59.8;
84.71 113.3 133.6;
44.03 136.3 118.6;
134.1 162.9 66.27;
100.4 166.1 79.38
146.3 134 84.63;
138.2 82.69 44.46;
140.2 80.06 82.72;
146 101.9 104.2;
150.2 114.3 66.33;
84.33 80.78 127.5;
82.02 53.45 98.45];

ROI_list=Tmni*[ROI_list ones(size(ROI_list,1),1)]';
ROI_list=ROI_list';
ROI_list=ROI_list(:,1:3);

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

load r1_all.mat
fid=fopen('Subjects_list_all.txt','r');
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);

file_suffix='Wendy';
data_location='E:\RAM data set\RAM_Public_Data_all\FR1_final';
image_storage_folder_sigle = [data_location '/Storage_filtering_no_overlap/Image_storage/' file_suffix '_merge/Raw_power_plot_ROIs'];
mkdir(image_storage_folder_sigle)

% for ROI_i=5:size(ROI_list,1)
for ROI_i=[5 12]
    all_power_traces=[];
    Montage_suffix=['final_clean_dil_roi_' num2str(ROI_i) ];
    for i=1:251
        try
            cd(data_location)
            cd([num2str(i),'_',r_sublist{i,1}]);
            load(['stats_traces_' file_suffix '_raw.mat']);
            load(['rejection_idx_' file_suffix '.mat']);
            load(['rejection_power_idx_' file_suffix '.mat']);
            
            trial_thr_ind=find(RAM_all_trials_total_number-noisy_RAM_all_trials_number<100);
            mean_power_traces(trial_thr_ind,:,:)=NaN;
            
            for power_i=1:size(noisy_RAM_all_trials_power_idx,2)
                mean_power_traces(power_i,:,noisy_RAM_all_trials_power_idx{power_i})=NaN;
            end

            if ROI_list(ROI_i,1) < 0
                load(['L_MontageMap_' Montage_suffix '.mat']);
                if ~isempty(L_MontageMap)
                    if size(L_MontageMap,1)==1
%                         buffer_zscore_traces=squeeze(zscore_traces(L_MontageMap(:,1),1,3,:));
%                         all_zscore_traces=[all_zscore_traces ; buffer_zscore_traces'];
                        elec_mean_power_traces=[];
                        elec_mean_power_traces=squeeze(mean_power_traces(L_MontageMap(1,1),:,:));
                        elec_mean_power_traces(:,noisy_RAM_all_trials_idx{L_MontageMap(1,1)})=[];

                        figure('position',[0 0 8000 4000]);
                        set(gcf, 'color', [1 1 1]);
                        set(gcf,'Visible','off');  
                        subplot(2,1,1)
                        plot(elec_mean_power_traces)

                        set(gca,'xtick',[5 10 15 20 25 30 35]);
                        set(gca,'xticklabel',[-250 0 250 500 750 1000 1250]);
                        set(gca,'xlim',[1 39])
                        set(gca,'ylim',[0 100])
                        set(gca,'Fontsize',15);
                        title(['Subject index : ' num2str(i) ', ROI index : ' num2str(ROI_i) ', Electrode index : ' num2str(L_MontageMap(1,1)) ],'Fontsize',12)
                        ylabel('Gamma power');
                        xlabel('Time (ms)');

                        subplot(2,1,2)
                        title('Histogram of maximum powers (all trials)')
                        histogram(max(elec_mean_power_traces),'BinLimits',[1,200])

                        saveas(gcf,[image_storage_folder_sigle '/raw_power_plot_roi_' num2str(ROI_i) '_subj_' num2str(i) '_elec_' num2str(L_MontageMap(1,1)) '.png'])
                        close all;

                        buffer_elec_mean_power_traces=nanmean(elec_mean_power_traces,2);
                        all_power_traces=[all_power_traces ; buffer_elec_mean_power_traces'];
                    else
%                         buffer_zscore_traces=squeeze(zscore_traces(L_MontageMap(:,1),1,3,:));
%                         all_zscore_traces=[all_zscore_traces ; buffer_zscore_traces];
                        for elec_i=1:size(L_MontageMap,1)
                            elec_mean_power_traces=[];
                            elec_mean_power_traces=squeeze(mean_power_traces(L_MontageMap(elec_i,1),:,:));
                            elec_mean_power_traces(:,noisy_RAM_all_trials_idx{L_MontageMap(elec_i,1)})=[];
                            
                            figure('position',[0 0 8000 4000]);
                            set(gcf, 'color', [1 1 1]);
                            set(gcf,'Visible','off');  
                            subplot(2,1,1)
                            plot(elec_mean_power_traces)

                            set(gca,'xtick',[5 10 15 20 25 30 35]);
                            set(gca,'xticklabel',[-250 0 250 500 750 1000 1250]);
                            set(gca,'xlim',[1 39])
                            set(gca,'ylim',[0 100])
                            set(gca,'Fontsize',15);
                            title(['Subject index : ' num2str(i) ', ROI index : ' num2str(ROI_i) ', Electrode index : ' num2str(L_MontageMap(elec_i,1)) ],'Fontsize',12)
                            ylabel('Gamma power');
                            xlabel('Time (ms)');
                            
                            subplot(2,1,2)
                            title('Histogram of maximum powers (all trials)')
                            histogram(max(elec_mean_power_traces),'BinLimits',[1,200])
                            
                            saveas(gcf,[image_storage_folder_sigle '/raw_power_plot_roi_' num2str(ROI_i) '_subj_' num2str(i) '_elec_' num2str(L_MontageMap(elec_i,1)) '.png'])
                            close all;
                            
                            buffer_elec_mean_power_traces=nanmean(elec_mean_power_traces,2);
                            all_power_traces=[all_power_traces ; buffer_elec_mean_power_traces'];
                        end
                    end
                end
            else
                load(['R_MontageMap_' Montage_suffix '.mat']);
                if ~isempty(R_MontageMap)
                    if size(R_MontageMap,1)==1
%                         buffer_zscore_traces=squeeze(zscore_traces(R_MontageMap(:,1),1,3,:));
%                         all_zscore_traces=[all_zscore_traces ; buffer_zscore_traces'];
                        elec_mean_power_traces=[];
                        elec_mean_power_traces=squeeze(mean_power_traces(R_MontageMap(1,1),:,:));
                        elec_mean_power_traces(:,noisy_RAM_all_trials_idx{R_MontageMap(1,1)})=[];

                        figure('position',[0 0 8000 4000]);
                        set(gcf, 'color', [1 1 1]);
                        set(gcf,'Visible','off');  
                        subplot(2,1,1)
                        plot(elec_mean_power_traces)

                        set(gca,'xtick',[5 10 15 20 25 30 35]);
                        set(gca,'xticklabel',[-250 0 250 500 750 1000 1250]);
                        set(gca,'xlim',[1 39])
                        set(gca,'ylim',[0 100])
                        set(gca,'Fontsize',15);
                        title(['Subject index : ' num2str(i) ', ROI index : ' num2str(ROI_i) ', Electrode index : ' num2str(R_MontageMap(1,1)) ],'Fontsize',12)
                        ylabel('Gamma power');
                        xlabel('Time (ms)');

                        subplot(2,1,2)
                        title('Histogram of maximum powers (all trials)')
                        histogram(max(elec_mean_power_traces),'BinLimits',[1,200])

                        saveas(gcf,[image_storage_folder_sigle '/raw_power_plot_roi_' num2str(ROI_i) '_subj_' num2str(i) '_elec_' num2str(R_MontageMap(1,1)) '.png'])
                        close all;

                        buffer_elec_mean_power_traces=nanmean(elec_mean_power_traces,2);
                        all_power_traces=[all_power_traces ; buffer_elec_mean_power_traces'];
                    else
%                         buffer_zscore_traces=squeeze(zscore_traces(R_MontageMap(:,1),1,3,:));
%                         all_zscore_traces=[all_zscore_traces ; buffer_zscore_traces];
                        for elec_i=1:size(R_MontageMap,1)
                            elec_mean_power_traces=[];
                            elec_mean_power_traces=squeeze(mean_power_traces(R_MontageMap(elec_i,1),:,:));
                            elec_mean_power_traces(:,noisy_RAM_all_trials_idx{R_MontageMap(elec_i,1)})=[];
                            
                            figure('position',[0 0 8000 4000]);
                            set(gcf, 'color', [1 1 1]);
                            set(gcf,'Visible','off');  
                            subplot(2,1,1)
                            plot(elec_mean_power_traces)

                            set(gca,'xtick',[5 10 15 20 25 30 35]);
                            set(gca,'xticklabel',[-250 0 250 500 750 1000 1250]);
                            set(gca,'xlim',[1 39])
                            set(gca,'ylim',[0 100])
                            set(gca,'Fontsize',15);
                            title(['Subject index : ' num2str(i) ', ROI index : ' num2str(ROI_i) ', Electrode index : ' num2str(R_MontageMap(elec_i,1)) ],'Fontsize',12)
                            ylabel('Gamma power');
                            xlabel('Time (ms)');
                            
                            subplot(2,1,2)
                            title('Histogram of maximum powers (all trials)')
                            histogram(max(elec_mean_power_traces),'BinLimits',[1,200])
                            
                            saveas(gcf,[image_storage_folder_sigle '/raw_power_plot_roi_' num2str(ROI_i) '_subj_' num2str(i) '_elec_' num2str(R_MontageMap(elec_i,1)) '.png'])
                            close all;
                            
                            buffer_elec_mean_power_traces=nanmean(elec_mean_power_traces,2);
                            all_power_traces=[all_power_traces ; buffer_elec_mean_power_traces'];
                        end
                    end
                end
            end
            disp(sprintf('%s is completed! ',r_sublist{i,1}));
        catch
        end
    end
    
    figure('position',[0 0 8000 4000]);
    set(gcf, 'color', [1 1 1]);
    set(gcf,'Visible','on');  
    
    plot(all_power_traces');
    
    hold on;
    plot(nanmean(all_power_traces,1),'LineWidth',2,'color','k');
    
    set(gca,'xtick',[5 10 15 20 25 30 35]);
    set(gca,'xticklabel',[-250 0 250 500 750 1000 1250]);
%     set(gca,'xlim',[1 39])
%     set(gca,'ylim',[-12 16])
%     line([1,39],[2,2],'color','b');
%     line([1,39],[-2,-2],'color','b');
%     line([10,10],[-12,16],'color','k');
    set(gca,'Fontsize',15);
    title(['ROI index : ' num2str(ROI_i) ', Total electrode # : ' num2str(size(all_power_traces,1)) ],'Fontsize',12)
    ylabel('Gamma power');
    xlabel('Time (ms)');
    
    cd(image_storage_folder_sigle)
    saveas(gcf,['raw_power_plot_roi_' num2str(ROI_i) '.png'])
    close all;
    
end



%% plot (mean in each subject)
clc;clear;
close all;

mni_nifti_path = 'C:\yale\bioimagesuite30\images\MNI_T1_1mm_stripped.nii';
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
Tmni = mnit1_1mm_stripped.transform;
ROI_list=[84 121.3 109.6 ;
65.57 74.26 122.6;
27.93 114 56.79;
96 68.79 101.2;
148.9 105.6 57.95;
84.63 77.91 121.1;
96.04 97.36 124.8;
113.8 127.1 64.86;
61.08 115.4 59.8;
84.71 113.3 133.6;
44.03 136.3 118.6;
134.1 162.9 66.27;
100.4 166.1 79.38
146.3 134 84.63;
138.2 82.69 44.46;
140.2 80.06 82.72;
146 101.9 104.2;
150.2 114.3 66.33;
84.33 80.78 127.5;
82.02 53.45 98.45];


ROI_list=Tmni*[ROI_list ones(size(ROI_list,1),1)]';
ROI_list=ROI_list';
ROI_list=ROI_list(:,1:3);

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

load r1_all.mat
fid=fopen('Subjects_list_all.txt','r');
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);

file_suffix='Wendy_merge';
data_location='E:\RAM data set\RAM_Public_Data_all\FR1_final';
image_storage_folder = [data_location '/Storage_filtering_no_overlap/Image_storage/' file_suffix '/ROIs_mean_subject'];
mkdir(image_storage_folder)

for ROI_i=20:size(ROI_list,1)
    all_zscore_traces=[];
    Montage_suffix=['final_clean_dil_roi_' num2str(ROI_i) ];
    for i=1:251
        try
            cd(data_location)
            cd([num2str(i),'_',r_sublist{i,1}]);
            load(['stats_traces_' file_suffix '.mat']);
            if ROI_list(ROI_i,1) < 0
                load(['L_MontageMap_' Montage_suffix '.mat']);
                if ~isempty(L_MontageMap)
                    if size(L_MontageMap,1)==1
                        buffer_zscore_traces=squeeze(zscore_traces(L_MontageMap(:,1),1,3,:));
                        all_zscore_traces=[all_zscore_traces ; [i buffer_zscore_traces']];
                    else
                        buffer_zscore_traces=squeeze(zscore_traces(L_MontageMap(:,1),1,3,:));
                        mean_buffer_zscore_traces=nanmean(buffer_zscore_traces,1);
                        all_zscore_traces=[all_zscore_traces ; [i mean_buffer_zscore_traces]];
                    end
                end
            else
                load(['R_MontageMap_' Montage_suffix '.mat']);
                if ~isempty(R_MontageMap)
                    if size(R_MontageMap,1)==1
                        buffer_zscore_traces=squeeze(zscore_traces(R_MontageMap(:,1),1,3,:));
                        all_zscore_traces=[all_zscore_traces ; [i buffer_zscore_traces']];
                    else
                        buffer_zscore_traces=squeeze(zscore_traces(R_MontageMap(:,1),1,3,:));
                        mean_buffer_zscore_traces=nanmean(buffer_zscore_traces,1);
                        all_zscore_traces=[all_zscore_traces ; [i mean_buffer_zscore_traces]];
                    end
                end
            end
            disp(sprintf('%s is completed! ',r_sublist{i,1}));
        catch
        end
    end
    
    figure('position',[0 0 8000 4000]);
    set(gcf, 'color', [1 1 1]);
    set(gcf,'Visible','on');  
    
    plot(all_zscore_traces');
    
    hold on;
    sroot_all_zscore_traces=nansum(all_zscore_traces,1)/sqrt(size(all_zscore_traces,1));
    plot(sroot_all_zscore_traces,'LineWidth',2,'color','k');
    
    set(gca,'xtick',[5 10 15 20 25 30 35]);
    set(gca,'xticklabel',[-250 0 250 500 750 1000 1250]);
    set(gca,'xlim',[1 39])
    set(gca,'ylim',[-12 16])
    line([1,39],[2,2],'color','b');
    line([1,39],[-2,-2],'color','b');
    line([10,10],[-12,16],'color','k');
    set(gca,'Fontsize',15);
    title(['ROI index : ' num2str(ROI_i) ', Total subject # : ' num2str(size(all_zscore_traces,1)) ],'Fontsize',12)
    ylabel('Baseline z-score value');
    xlabel('Time (ms)');
    
    cd(image_storage_folder)
    saveas(gcf,['zscore_plot_roi_' num2str(ROI_i) '.png'])
    close all;
    
end




%% plot (aggregate in each subject)
clc;clear;
close all;

mni_nifti_path = 'C:\yale\bioimagesuite30\images\MNI_T1_1mm_stripped.nii';
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
Tmni = mnit1_1mm_stripped.transform;
ROI_list=[84 121.3 109.6 ;
65.57 74.26 122.6;
27.93 114 56.79;
96 68.79 101.2;
148.9 105.6 57.95;
84.63 77.91 121.1;
96.04 97.36 124.8;
113.8 127.1 64.86;
61.08 115.4 59.8;
84.71 113.3 133.6;
44.03 136.3 118.6;
134.1 162.9 66.27;
100.4 166.1 79.38];

ROI_list=Tmni*[ROI_list ones(size(ROI_list,1),1)]';
ROI_list=ROI_list';
ROI_list=ROI_list(:,1:3);

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

load r1_all.mat
fid=fopen('Subjects_list_all.txt','r');
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);

file_suffix='Wendy_merge';
data_location='E:\RAM data set\RAM_Public_Data_all\FR1_final';
image_storage_folder = [data_location '/Storage_filtering_no_overlap/Image_storage/' file_suffix '/ROIs_aggregate_subject'];
mkdir(image_storage_folder)

for ROI_i=1:size(ROI_list,1)
    all_zscore_traces=[];
    Montage_suffix=['final_clean_dil_roi_' num2str(ROI_i) ];
    for i=1:251
        try
            cd(data_location)
            cd([num2str(i),'_',r_sublist{i,1}]);
            load(['stats_traces_' file_suffix '.mat']);
            if ROI_list(ROI_i,1) < 0
                load(['L_MontageMap_' Montage_suffix '.mat']);
                if ~isempty(L_MontageMap)
                    if size(L_MontageMap,1)==1
                        buffer_zscore_traces=squeeze(zscore_traces(L_MontageMap(:,1),1,3,:));
                        all_zscore_traces=[all_zscore_traces ; buffer_zscore_traces'];
                    else
                        buffer_zscore_traces=squeeze(zscore_traces(L_MontageMap(:,1),1,3,:));
                        aggregate_buffer_zscore_traces=nansum(buffer_zscore_traces,1);
                        all_zscore_traces=[all_zscore_traces ; aggregate_buffer_zscore_traces];
                    end
                end
            else
                load(['R_MontageMap_' Montage_suffix '.mat']);
                if ~isempty(R_MontageMap)
                    if size(R_MontageMap,1)==1
                        buffer_zscore_traces=squeeze(zscore_traces(R_MontageMap(:,1),1,3,:));
                        all_zscore_traces=[all_zscore_traces ; buffer_zscore_traces'];
                    else
                        buffer_zscore_traces=squeeze(zscore_traces(R_MontageMap(:,1),1,3,:));
                        aggregate_buffer_zscore_traces=nansum(buffer_zscore_traces,1);
                        all_zscore_traces=[all_zscore_traces ; aggregate_buffer_zscore_traces];
                    end
                end
            end
            disp(sprintf('%s is completed! ',r_sublist{i,1}));
        catch
        end
    end
    
    figure('position',[0 0 8000 4000]);
    set(gcf, 'color', [1 1 1]);
    set(gcf,'Visible','on');  
    
    plot(all_zscore_traces');
    
    hold on;
    sroot_all_zscore_traces=nansum(all_zscore_traces,1)/sqrt(size(all_zscore_traces,1));
    plot(sroot_all_zscore_traces,'LineWidth',2,'color','k');
    
    set(gca,'xtick',[5 10 15 20 25 30 35]);
    set(gca,'xticklabel',[-250 0 250 500 750 1000 1250]);
    set(gca,'xlim',[1 39])
    set(gca,'ylim',[-12 16])
    line([1,39],[2,2],'color','b');
    line([1,39],[-2,-2],'color','b');
    line([10,10],[-12,16],'color','k');
    set(gca,'Fontsize',15);
    title(['ROI index : ' num2str(ROI_i) ', Total subject # : ' num2str(size(all_zscore_traces,1)) ],'Fontsize',12)
    ylabel('Baseline z-score value');
    xlabel('Time (ms)');
    
    cd(image_storage_folder)
    saveas(gcf,['zscore_plot_roi_' num2str(ROI_i) '.png'])
    close all;
    
end




