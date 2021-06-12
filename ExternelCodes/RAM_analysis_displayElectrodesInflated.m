%% Setting
clc;clear;
close all;

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

load r1_all.mat
load fs_DK_atlas.mat
fid=fopen('Subjects_list_all.txt','r');
cnt=1;
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);

% data_location='E:\RAM data set\RAM_Public_Data_all\FR1_final';
% data_location='E:\RAM data set\RAM_Public_Data_all\FR1_FARNAM';
data_location='E:\RAM data set\RAM_Public_Data_all\FR1_FS';

%% previous prefix
% file_suffix='Wendy_p005';
% file_suffix='Wendy_fdr';
% file_suffix='Wendy_ftest_fdr';
% file_suffix='Wendy';
% file_suffix='Wendy_fdr';
% file_suffix='final_clean';
% file_suffix='Wendy_thr_02';
% file_suffix='Wendy_thr_01';
% file_suffix='Wendy_merge';

% file_suffix='Wendy_hunki_std_5_20';
% file_suffix='Wendy_hunki_std_5_20_40';
% file_suffix='Wendy_hunki_std_5_20_50';
% file_suffix='Wendy_hunki_std_5_20_60';
% file_suffix='Wendy_hunki_std_5_20_base';
% file_suffix='Wendy_hunki_std_5_20_base_40';
% file_suffix='Wendy_hunki_std_5_20_base_50';
% file_suffix='Wendy_hunki_std_5_20_base_60';

% Montage_suffix='final_clean_dil_Wendy_fdr';
% Montage_suffix='final_clean_dil_Wendy';
% Montage_suffix='final_clean_dil';

%% new 
% file_suffix='Wendy_hunki_std_5_20_fs';
% file_suffix='Wendy_hunki_std_5_20_fs_xflip';
% file_suffix='Wendy_hunki_std_5_20_fs_Kate';
% file_suffix='Wendy_hunki_std_5_20_fs_Kate_xflip';

% file_suffix='Wendy_hunki_std_5_20_fs_Rec';
% file_suffix='Wendy_hunki_std_5_20_fs_nonRec';

% file_suffix='Wendy_hunki_std_5_20_fs_session_Rec';
% file_suffix='Wendy_hunki_std_5_20_fs_session_nonRec';
% file_suffix='Wendy_hunki_std_5_20_fs_session_all';

%% fs 500
% file_suffix='Final_std_5_20_fs_Rec';
% file_suffix='Final_std_5_20_fs_nonRec';
file_suffix='Final_std_5_20_fs_together';

Montage_suffix='final_roi_fs';
% Montage_suffix='final_roi_fs_xflip';

patients=r_sublist;
overlay=4;
laterality=0;
views=4;
inflationstep=5;
electrodeDensity=0;
electrodeDisplay=0;
thr=2;
mean_flag=2;
% color_range=[-8 10];
% color_range=[-15 15];
color_range=[-10 15];
% color_range=[-3 3];
% color_range=[-1 1];
statistic_flag=0;

all_sub_index=1:251;

displayElectrodesInflated_new_fs(data_location, patients, overlay, laterality, views, inflationstep,electrodeDensity, electrodeDisplay ,thr,mean_flag,color_range, file_suffix,Montage_suffix,statistic_flag,all_sub_index)


%% ROI mapping
% file_suffix='Final_std_5_20_fs_Rec';
% file_suffix='Final_std_5_20_fs_nonRec';
file_suffix='Final_std_5_20_fs_together';

Montage_suffix='final_roi_fs';
% Montage_suffix='final_roi_fs_xflip';

patients=r_sublist;
overlay=4;
laterality=0;
views=4;
inflationstep=5;
electrodeDensity=0;
electrodeDisplay=0;
thr=0.5;
color_range=[-2 2];
statistic_flag=0;

all_sub_index=1:251;
atlas_index_all=[3 4 6 8 9 10 11 12 13 15 16 17 22 24 26 28 29];
atlas_index_all_rest=1:36;
atlas_index_all_rest(atlas_index_all)=[];
atlas_index_all_rest(1)=[];
for atlas_index=atlas_index_all_rest
    displayElectrodesInflated_new_fs_Atlas(data_location, patients, overlay, laterality, views, inflationstep,electrodeDensity, electrodeDisplay ,thr,color_range, file_suffix,Montage_suffix,statistic_flag,all_sub_index,atlas_index)
    DK_atlas_names{atlas_index}
    close all;
end

%% ROI mapping (network level)
% file_suffix='Final_std_5_20_fs_Rec';
% file_suffix='Final_std_5_20_fs_nonRec';
file_suffix='Final_std_5_20_fs_together';

Montage_suffix='final_roi_fs';
% Montage_suffix='final_roi_fs_xflip';

patients=r_sublist;
overlay=4;
laterality=0;
views=4;
inflationstep=5;
electrodeDensity=0;
electrodeDisplay=0;
thr=0.5;
color_range=[-2 2];
statistic_flag=0;

all_sub_index=1:251;
network_index_all=[2 3 4 5 6 7 8];
for atlas_index=network_index_all
    displayElectrodesInflated_new_fs_Atlas(data_location, patients, overlay, laterality, views, inflationstep,electrodeDensity, electrodeDisplay ,thr,color_range, file_suffix,Montage_suffix,statistic_flag,all_sub_index,atlas_index)
    atlas_index
    close all;
end

displayElectrodesInflated_new_fs_Atlas_all(data_location, patients, overlay, laterality, views, inflationstep,electrodeDensity, electrodeDisplay ,thr,color_range, file_suffix,Montage_suffix,statistic_flag,all_sub_index)
   

%% ROI mapping (parcel level)
% file_suffix='Final_std_5_20_fs_Rec';
% file_suffix='Final_std_5_20_fs_nonRec';
file_suffix='Final_std_5_20_fs_together';

Montage_suffix='final_roi_fs';
% Montage_suffix='final_roi_fs_xflip';

patients=r_sublist;
overlay=4;
laterality=0;
views=4;
inflationstep=5;
electrodeDensity=0;
electrodeDisplay=0;
thr=0.5;
color_range=[-2 2];
statistic_flag=0;

all_sub_index=1:251;
network_index_all=[91 92 93];
for atlas_index=network_index_all
    displayElectrodesInflated_new_fs_Atlas(data_location, patients, overlay, laterality, views, inflationstep,electrodeDensity, electrodeDisplay ,thr,color_range, file_suffix,Montage_suffix,statistic_flag,all_sub_index,atlas_index)
    atlas_index
    close all;
end







%% stats
file_suffix='Final_std_5_20_fs_together_ttest';

Montage_suffix='final_roi_fs';
% Montage_suffix='final_roi_fs_xflip';

patients=r_sublist;
overlay=4;
laterality=0;
views=4;
inflationstep=5;
electrodeDensity=0;
electrodeDisplay=0;
thr=0.1;
stat_flag=1; % 1) uncorrected 0.05 / 2) FDR 0.05
color_range=[-6 6];

displayElectrodesInflated_new_fs_stat(data_location, patients, overlay, laterality, views, inflationstep,electrodeDensity, electrodeDisplay ,thr,stat_flag,color_range, file_suffix,Montage_suffix)




%% ROI plot zscore
% Montage_suffix='final_roi_fs';
% Montage_suffix='final_roi_fs_xflip';

atlas_index_all=[3 4 6 8 9 10 11 12 13 15 16 17 22 24 26 28 29];
for atlas_index=atlas_index_all_rest
    
    figure('position',[0 0 8000 8000]);
    set(gcf, 'color', [1 1 1]);
    set(gcf,'Visible','off');
    
    %% Recalled
    file_suffix='Final_std_5_20_fs_Rec';
    Atlas_storage_folder = [data_location '/Storage_filtering_overlap/Atlas/' file_suffix '/' ];
    cd(Atlas_storage_folder)
    atlas_all_vertex_val_both=[];
    
    % plot
    subplot(4,3,1)
    clearvars mean_atlas_all_vertex_val atlas_all_vertex_val sub_size
    load([num2str(atlas_index) '_' DK_atlas_names{atlas_index} '/R_vertex_values.mat' ])
    mean_atlas_all_vertex_val=squeeze(nanmean(atlas_all_vertex_val,1));
%     sub_size=size(find(~isnan(mean_atlas_all_vertex_val(:,1))>0),1);
    atlas_all_vertex_val_both=cat(1,atlas_all_vertex_val_both,atlas_all_vertex_val);
%     plot(nansum(mean_atlas_all_vertex_val)/sqrt(sub_size))
    plot_areaerrorbar(mean_atlas_all_vertex_val)
    test1=mean_atlas_all_vertex_val;
    title(['Right Hemisphere' ])
    set(gca,'xtick',[1 20 40 60 80]);
    set(gca,'xticklabel',[-500 0 500 1000 1500]);
    set(gca,'xlim',[1 98])
    set(gca,'ytick',[-6 -2 2 6]);
    set(gca,'yticklabel',[-6 -2 2 6 ]);
    set(gca,'ylim',[-8 10])
    set(gca,'Fontsize',15);
    line([1,98],[2,2],'Color','r','LineWidth',0.5)
    line([1,98],[-2,-2],'Color','r','LineWidth',0.5)
    ylabel({'Rec','Gamma power z-score'});
    xlabel('Time (ms)');

    subplot(4,3,2)
    clearvars mean_atlas_all_vertex_val atlas_all_vertex_val sub_size
    load([num2str(atlas_index) '_' DK_atlas_names{atlas_index} '/L_vertex_values.mat' ])
    mean_atlas_all_vertex_val=squeeze(nanmean(atlas_all_vertex_val,1));
%     sub_size=size(find(~isnan(mean_atlas_all_vertex_val(:,1))>0),1);
    atlas_all_vertex_val_both=cat(1,atlas_all_vertex_val_both,atlas_all_vertex_val);
%     plot(nansum(mean_atlas_all_vertex_val)/sqrt(sub_size))
    plot_areaerrorbar(mean_atlas_all_vertex_val)
    test2=mean_atlas_all_vertex_val;
    title(['Left Hemisphere'])
    set(gca,'xtick',[1 20 40 60 80]);
    set(gca,'xticklabel',[-500 0 500 1000 1500]);
    set(gca,'xlim',[1 98])
    set(gca,'ytick',[-6 -2 2 6 ]);
    set(gca,'yticklabel',[-6 -2 2 6 ]);
    set(gca,'ylim',[-8 10])
    set(gca,'Fontsize',15);
    line([1,98],[2,2],'Color','r','LineWidth',0.5)
    line([1,98],[-2,-2],'Color','r','LineWidth',0.5)
    ylabel('z-score');
    xlabel('Time (ms)');

    subplot(4,3,3)
    clearvars mean_atlas_all_vertex_val atlas_all_vertex_val sub_size
    mean_atlas_all_vertex_val_both=squeeze(nanmean(atlas_all_vertex_val_both,1));
%     sub_size=size(find(~isnan(mean_atlas_all_vertex_val_both(:,1))>0),1);
%     plot(nansum(mean_atlas_all_vertex_val_both)/sqrt(sub_size))
    plot_areaerrorbar(mean_atlas_all_vertex_val_both)
    test3=mean_atlas_all_vertex_val_both;
    title(['Both Hemispheres'])
    set(gca,'xtick',[1 20 40 60 80]);
    set(gca,'xticklabel',[-500 0 500 1000 1500]);
    set(gca,'xlim',[1 98])
    set(gca,'ytick',[-6 -2 2 6 ]);
    set(gca,'yticklabel',[-6 -2 2 6 ]);
    set(gca,'ylim',[-8 10])
    set(gca,'Fontsize',15);
    line([1,98],[2,2],'Color','r','LineWidth',0.5)
    line([1,98],[-2,-2],'Color','r','LineWidth',0.5)
    ylabel('z-score');
    xlabel('Time (ms)');

    
    
    %% Non recalled
    file_suffix='Final_std_5_20_fs_nonRec';
    Atlas_storage_folder = [data_location '/Storage_filtering_overlap/Atlas/' file_suffix '/' ];
    cd(Atlas_storage_folder)
    atlas_all_vertex_val_both=[];
    
    % plot
    subplot(4,3,4)
    clearvars mean_atlas_all_vertex_val atlas_all_vertex_val sub_size
    load([num2str(atlas_index) '_' DK_atlas_names{atlas_index} '/R_vertex_values.mat' ])
    mean_atlas_all_vertex_val=squeeze(nanmean(atlas_all_vertex_val,1));
%     sub_size=size(find(~isnan(mean_atlas_all_vertex_val(:,1))>0),1);
    atlas_all_vertex_val_both=cat(1,atlas_all_vertex_val_both,atlas_all_vertex_val);
%     plot(nansum(mean_atlas_all_vertex_val)/sqrt(sub_size))
    plot_areaerrorbar(mean_atlas_all_vertex_val)
    test4=mean_atlas_all_vertex_val;
%     title(['Right Hemisphere' ])
    set(gca,'xtick',[1 20 40 60 80]);
    set(gca,'xticklabel',[-500 0 500 1000 1500]);
    set(gca,'xlim',[1 98])
    set(gca,'ytick',[-6 -2 2 6]);
    set(gca,'yticklabel',[-6 -2 2 6 ]);
    set(gca,'ylim',[-8 10])
    set(gca,'Fontsize',15);
    line([1,98],[2,2],'Color','r','LineWidth',0.5)
    line([1,98],[-2,-2],'Color','r','LineWidth',0.5)
    ylabel({'NRec','z-score'});
    xlabel('Time (ms)');

    subplot(4,3,5)
    clearvars mean_atlas_all_vertex_val atlas_all_vertex_val sub_size
    load([num2str(atlas_index) '_' DK_atlas_names{atlas_index} '/L_vertex_values.mat' ])
    mean_atlas_all_vertex_val=squeeze(nanmean(atlas_all_vertex_val,1));
%     sub_size=size(find(~isnan(mean_atlas_all_vertex_val(:,1))>0),1);
    atlas_all_vertex_val_both=cat(1,atlas_all_vertex_val_both,atlas_all_vertex_val);
%     plot(nansum(mean_atlas_all_vertex_val)/sqrt(sub_size))
    plot_areaerrorbar(mean_atlas_all_vertex_val)
    test5=mean_atlas_all_vertex_val;
%     title(['Left Hemisphere'])
    set(gca,'xtick',[1 20 40 60 80]);
    set(gca,'xticklabel',[-500 0 500 1000 1500]);
    set(gca,'xlim',[1 98])
    set(gca,'ytick',[-6 -2 2 6 ]);
    set(gca,'yticklabel',[-6 -2 2 6 ]);
    set(gca,'ylim',[-8 10])
    set(gca,'Fontsize',15);
    line([1,98],[2,2],'Color','r','LineWidth',0.5)
    line([1,98],[-2,-2],'Color','r','LineWidth',0.5)
    ylabel('z-score');
    xlabel('Time (ms)');

    subplot(4,3,6)
    clearvars mean_atlas_all_vertex_val atlas_all_vertex_val sub_size
    mean_atlas_all_vertex_val_both=squeeze(nanmean(atlas_all_vertex_val_both,1));
%     sub_size=size(find(~isnan(mean_atlas_all_vertex_val_both(:,1))>0),1);
%     plot(nansum(mean_atlas_all_vertex_val_both)/sqrt(sub_size))
    plot_areaerrorbar(mean_atlas_all_vertex_val_both)
    test6=mean_atlas_all_vertex_val_both;
%     title(['Both Hemispheres'])
    set(gca,'xtick',[1 20 40 60 80]);
    set(gca,'xticklabel',[-500 0 500 1000 1500]);
    set(gca,'xlim',[1 98])
    set(gca,'ytick',[-6 -2 2 6 ]);
    set(gca,'yticklabel',[-6 -2 2 6 ]);
    set(gca,'ylim',[-8 10])
    set(gca,'Fontsize',15);
    line([1,98],[2,2],'Color','r','LineWidth',0.5)
    line([1,98],[-2,-2],'Color','r','LineWidth',0.5)
    ylabel('z-score');
    xlabel('Time (ms)');

    
    
    %% Recalled - Non recalled
    % plot
    subplot(4,3,7)
    clearvars mean_atlas_all_vertex_val atlas_all_vertex_val sub_size
    plot_areaerrorbar(test1-test4)
%     title(['Right Hemisphere' ])
    set(gca,'xtick',[1 20 40 60 80]);
    set(gca,'xticklabel',[-500 0 500 1000 1500]);
    set(gca,'xlim',[1 98])
    set(gca,'ytick',[-6 -2 2 6]);
    set(gca,'yticklabel',[-6 -2 2 6 ]);
    set(gca,'ylim',[-8 10])
    set(gca,'Fontsize',15);
    line([1,98],[2,2],'Color','r','LineWidth',0.5)
    line([1,98],[-2,-2],'Color','r','LineWidth',0.5)
    ylabel({'Rec - NRec','z-score'});
    xlabel('Time (ms)');

    subplot(4,3,8)
    clearvars mean_atlas_all_vertex_val atlas_all_vertex_val sub_size
    plot_areaerrorbar(test2-test5)
%     title(['Right Hemisphere' ])
    set(gca,'xtick',[1 20 40 60 80]);
    set(gca,'xticklabel',[-500 0 500 1000 1500]);
    set(gca,'xlim',[1 98])
    set(gca,'ytick',[-6 -2 2 6]);
    set(gca,'yticklabel',[-6 -2 2 6 ]);
    set(gca,'ylim',[-8 10])
    set(gca,'Fontsize',15);
    line([1,98],[2,2],'Color','r','LineWidth',0.5)
    line([1,98],[-2,-2],'Color','r','LineWidth',0.5)
    ylabel('z-score');
    xlabel('Time (ms)');
    
    subplot(4,3,9)
    clearvars mean_atlas_all_vertex_val atlas_all_vertex_val sub_size
    plot_areaerrorbar(test3-test6)
%     title(['Right Hemisphere' ])
    set(gca,'xtick',[1 20 40 60 80]);
    set(gca,'xticklabel',[-500 0 500 1000 1500]);
    set(gca,'xlim',[1 98])
    set(gca,'ytick',[-6 -2 2 6]);
    set(gca,'yticklabel',[-6 -2 2 6 ]);
    set(gca,'ylim',[-8 10])
    set(gca,'Fontsize',15);
    line([1,98],[2,2],'Color','r','LineWidth',0.5)
    line([1,98],[-2,-2],'Color','r','LineWidth',0.5)
    ylabel('z-score');
    xlabel('Time (ms)');

    
    
    %% Recalled + Non recalled
    % plot
    subplot(4,3,10)
    clearvars mean_atlas_all_vertex_val atlas_all_vertex_val sub_size
    plot_areaerrorbar(test1+test4)
%     title(['Right Hemisphere' ])
    set(gca,'xtick',[1 20 40 60 80]);
    set(gca,'xticklabel',[-500 0 500 1000 1500]);
    set(gca,'xlim',[1 98])
    set(gca,'ytick',[-6 -2 2 6]);
    set(gca,'yticklabel',[-6 -2 2 6 ]);
    set(gca,'ylim',[-8 10])
    set(gca,'Fontsize',15);
    line([1,98],[2,2],'Color','r','LineWidth',0.5)
    line([1,98],[-2,-2],'Color','r','LineWidth',0.5)
    ylabel({'Rec + NRec','z-score'});
    xlabel('Time (ms)');

    subplot(4,3,11)
    clearvars mean_atlas_all_vertex_val atlas_all_vertex_val sub_size
    plot_areaerrorbar(test2+test5)
%     title(['Right Hemisphere' ])
    set(gca,'xtick',[1 20 40 60 80]);
    set(gca,'xticklabel',[-500 0 500 1000 1500]);
    set(gca,'xlim',[1 98])
    set(gca,'ytick',[-6 -2 2 6]);
    set(gca,'yticklabel',[-6 -2 2 6 ]);
    set(gca,'ylim',[-8 10])
    set(gca,'Fontsize',15);
    line([1,98],[2,2],'Color','r','LineWidth',0.5)
    line([1,98],[-2,-2],'Color','r','LineWidth',0.5)
    ylabel('z-score');
    xlabel('Time (ms)');
    
    subplot(4,3,12)
    clearvars mean_atlas_all_vertex_val atlas_all_vertex_val sub_size
    plot_areaerrorbar(test3+test6)
%     title(['Right Hemisphere' ])
    set(gca,'xtick',[1 20 40 60 80]);
    set(gca,'xticklabel',[-500 0 500 1000 1500]);
    set(gca,'xlim',[1 98])
    set(gca,'ytick',[-6 -2 2 6]);
    set(gca,'yticklabel',[-6 -2 2 6 ]);
    set(gca,'ylim',[-8 10])
    set(gca,'Fontsize',15);
    line([1,98],[2,2],'Color','r','LineWidth',0.5)
    line([1,98],[-2,-2],'Color','r','LineWidth',0.5)
    ylabel('z-score');
    xlabel('Time (ms)');    

%     [ax,h1]=suplabel('super X label');
%     [ax,h2]=suplabel('super Y label','y');
    [ax,h1]=suplabel(['ROI : ' DK_atlas_names{atlas_index}] ,'t');
    set(h1,'FontSize',25)
    
    
    %% save figure
    cd([data_location '/Storage_filtering_overlap/Atlas/'])
    saveas(gcf,[num2str(atlas_index) '_' DK_atlas_names{atlas_index} '_zscore_plot.png'])
    close all;
    
    
end


%% Atlas plot
% file_suffix='atlas';
% 
% Montage_suffix='final_roi_fs';
% % Montage_suffix='final_roi_fs_xflip';
% 
% patients=r_sublist;
% overlay=4;
% laterality=0;
% views=4;
% inflationstep=5;
% electrodeDensity=0;
% electrodeDisplay=0;
% thr=0.5;
% color_range=[-2 2];
% statistic_flag=0;
% 
% all_sub_index=1:251;
% 
% displayElectrodesInflated_new_fs_Atlas(data_location, patients, overlay, laterality, views, inflationstep,electrodeDensity, electrodeDisplay ,thr,color_range, file_suffix,Montage_suffix,statistic_flag,all_sub_index)





% Montage_suffix='final_clean';
% displayElectrodesInflated_new(data_location, patients, overlay, laterality, views, inflationstep,electrodeDensity, electrodeDisplay ,thr,mean_flag,color_range, file_suffix,Montage_suffix,statistic_flag) %, frames_set)






%% fs electrodes plot
file_suffix='Electrodes_distribution_fs';
Montage_suffix='final_roi_fs_xflip';

patients=r_sublist;
overlay=4;
laterality=0;
views=4;
inflationstep=5;
electrodeDensity=0;
electrodeDisplay=1;
thr=0.01;
mean_flag=2;
color_range=[-1 3];
statistic_flag=0;

all_sub_index=1:251;

displayElectrodesInflated_new_fs(data_location, patients, overlay, laterality, views, inflationstep,electrodeDensity, electrodeDisplay ,thr,mean_flag,color_range, file_suffix,Montage_suffix,statistic_flag,all_sub_index)


%% single plot
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

data_location='E:\RAM data set\RAM_Public_Data_all\FR1_FS';
% file_suffix='final_clean';
% file_suffix='Wendy_thr_02';
file_suffix='Wendy_thr_01';
% file_suffix='Wendy_merge';
% file_suffix='Wendy';
patients=r_sublist;
overlay=4;
laterality=0;
views=4;
inflationstep=5;
electrodeDensity=0;
electrodeDisplay=1;
thr=2;
mean_flag=2;
color_range=[-8 10];
statistic_flag=0;

Montage_suffix='final_roi_fs';
for subject_i=1:251
    try
        displayElectrodesInflated_new_single(data_location,subject_i, patients{subject_i}, overlay, laterality, views, inflationstep,electrodeDensity, electrodeDisplay ,thr,mean_flag,color_range, file_suffix,Montage_suffix,statistic_flag) %, frames_set)
        close all;
    catch
    end
end

