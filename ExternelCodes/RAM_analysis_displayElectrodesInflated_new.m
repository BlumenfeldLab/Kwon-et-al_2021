
clc;clear;

addpath(genpath('E:\RAM data set\RAM_Public_Data_all\Parcellations'));
addpath(genpath('E:\RAM data set\RAM_Public_Data_all\fsaverage'));
addpath(genpath('E:\RAM data set\RAM_Public_Data_all\Codes'));

rootfolder = 'E:\RAM data set\RAM_Public_Data_all';
cd(rootfolder)

load r1_all.mat
load fs_DK_atlas.mat

fid=fopen('Subjects_list_all.txt','r');
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);

storage_suffix='base_agg';
% storage_suffix='base_avg';
% storage_suffix='base_sravg';

z_color_range=[-10 10];
t_color_range=[-6 6];

data_location='E:\RAM data set\RAM_Public_Data_all\FR1_FS';
% file_suffix_all={'rec'};
% file_suffix_all={'together'};
file_suffix_all={'together','rec','Nonrec','diff'};
Montage_suffix='final_roi_fs';
flip_flag=0;

% cluster_stats_all={'cluster_count','cluster_count_rand','cluster_sumt','cluster_sumt_rand'}; % 
% cluster_stats_all={'cluster_sumt_rand','cluster_sumt_rand_001','cluster_sumt'}; % 
% cluster_stats_all={'cluster_sumt_rand_005','cluster_sumt_rand_001','cluster_sumt_rand_0005'}; % 
cluster_stats_all={'cluster_sumt_rand_005'}; % 

% atlas_index_all={'400','200','DK'};
% atlas_index_all={'400','1000'}; % 'DK' % '200'
% atlas_index_all={'1000'}; % 'DK' % '200'
atlas_index_all={'400'}; % 'DK' % '200'

all_sub_index=1:251;
% all_sub_index=11;

%% mapping only zscore
zscore_mapping_flag=1;
Cluster_based_statistics_flag=0;
tscore_mapping_flag=0;
tscore_mapping_flag_raw=0;
zscore_mapping_stats_flag=0;
zscore_flag=0; % 0 or 1 
for file_suffix_i=1:size(file_suffix_all,2)
    file_suffix=['Final_std_20_fs_' file_suffix_all{file_suffix_i}];
    Final_pipeline_displayElectrodesInflated_HK(data_location,r_sublist,all_sub_index,Montage_suffix, ...
        file_suffix,storage_suffix,flip_flag,zscore_mapping_flag,Cluster_based_statistics_flag, ...
        tscore_mapping_flag,tscore_mapping_flag_raw,zscore_mapping_stats_flag,[],[],zscore_flag, ...
        z_color_range,t_color_range)
end

%% stats & mapping tscore & zscore
zscore_mapping_flag=0;
Cluster_based_statistics_flag=0;
tscore_mapping_flag=0;
tscore_mapping_flag_raw=0;
zscore_flag=0; % 0 or 1
zscore_mapping_stats_flag=1;
for file_suffix_i=1:size(file_suffix_all,2)
    for atlas_index_all_i=1:size(atlas_index_all,2)
        file_suffix=['Final_std_20_fs_' file_suffix_all{file_suffix_i}];
        Final_pipeline_displayElectrodesInflated_HK(data_location,r_sublist,all_sub_index,Montage_suffix, ...
            file_suffix,storage_suffix,flip_flag,zscore_mapping_flag,Cluster_based_statistics_flag, ...
            tscore_mapping_flag,tscore_mapping_flag_raw,zscore_mapping_stats_flag,atlas_index_all{atlas_index_all_i}, ...
            cluster_stats_all,zscore_flag,z_color_range,t_color_range)
    end
end

data_location='E:\RAM data set\RAM_Public_Data_all\FR1_FS';
overlay=4;
laterality=0; % 0 or 4
views=4;
electrodeDensity=0;
electrodeDisplay=0;
thr=1;
% atlas_index_all=[];
atlas_index_all=[1:36];
displayElectrodesInflated_new_fs_Atlas_all(data_location, overlay, laterality, views, ...
    electrodeDensity,electrodeDisplay, thr, atlas_index_all) 
close all




%% Plot electrodes only
data_location='E:\RAM data set\RAM_Public_Data_all\FR1_FS';
patients=r_sublist;
overlay=4;
laterality=0; % 0 or 4
views=4;
electrodeDensity=0;
electrodeDisplay=1;
mean_flag=2;
color_range=[];
thr=1;
%         colormap_mat='cmap_hunki_10to15.mat';
colormap_mat='cmap_hunki_new.mat';

for ii=1:251
    try
        patient_folder = [data_location '/' num2str(ii) '_' patients{ii}];
        cd(patient_folder)
        all_sub_index=ii;
        displayElectrodesInflated_new_fs_Electrodes(data_location, patients, overlay, laterality, views, electrodeDensity, electrodeDisplay ,thr,mean_flag,color_range, file_suffix,Montage_suffix,all_sub_index,colormap_mat)
        copyfile('*.tif',patient_folder)  
        close all;
    catch
    end
end



%% ROI mapping new
clc;clear;

rootfolder = 'E:\RAM data set\RAM_Public_Data_all';
cd(rootfolder)

load fs_DK_atlas.mat

data_location='E:\RAM data set\RAM_Public_Data_all\FR1_FS';

overlay=4;
laterality=0;
views=4;
electrodeDensity=0;
thr=0.5;

% atlas_all_index=[4 8 12 13 28 30];
atlas_all_index=1:36;
displayElectrodesInflated_new_fs_Atlas_all(data_location, overlay, laterality, views,electrodeDensity, thr, atlas_all_index) 
close all






























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

data_location='E:\RAM data set\RAM_Public_Data_all\FR1_FS';


%% fs new both hemi
close all

flip_flag=0;
Montage_suffix='final_roi_fs';

file_suffix='Final_std_20_fs_together';
% file_suffix='Final_std_20_fs_rec';
% file_suffix='Final_std_20_fs_Nonrec';
% file_suffix='Final_std_20_fs_diff';

patients=r_sublist;
overlay=4;
laterality=4;
views=4;
electrodeDensity=0; %
electrodeDisplay=0; % 
thr=2;
mean_flag=2;
% color_range=[-8 10];
% color_range=[-15 15];
% color_range=[-10 15];
% color_range=[-3 3];
% color_range=[-1 1];
color_range=[-20 30];
colormap_mat='cmap_hunki_20to30.mat';


all_sub_index=1:251;
% all_sub_index=find(session_n_all==1)'; % test

displayElectrodesInflated_new_fs(data_location, patients, overlay, laterality, views, electrodeDensity, electrodeDisplay ,thr,mean_flag,color_range, file_suffix,Montage_suffix,all_sub_index,colormap_mat)
close all;

colormap_label='Gamma power z-score';
display_dynamic_hunki(flip_flag,colormap_mat,color_range,colormap_label)







%% fs new fliped x
close all

flip_flag=0;
Montage_suffix='final_roi_fs_xflip';

% file_suffix='Final_std_20_fs_together_xflip';
% file_suffix='Final_std_20_fs_rec_xflip';
% file_suffix='Final_std_20_fs_Nonrec_xflip';
file_suffix='Final_std_20_fs_diff_xflip';

patients=r_sublist;
overlay=4;
laterality=4;
views=4;
electrodeDensity=0; %
electrodeDisplay=0; % 
thr=4;
mean_flag=2;
% color_range=[-8 10];
% color_range=[-15 15];
% color_range=[-10 15];
color_range=[-20 30];
% color_range=[-3 3];
% color_range=[-1 1];

all_sub_index=1:251;
% all_sub_index=find(session_n_all==1)'; % test

displayElectrodesInflated_new_fs(data_location, patients, overlay, laterality, views, electrodeDensity, electrodeDisplay ,thr,mean_flag,color_range, file_suffix,Montage_suffix,all_sub_index)
close all;


%% ROI mapping

clc;clear;

rootfolder = 'E:\RAM data set\RAM_Public_Data_all';
cd(rootfolder)

load r1_all.mat
load fs_DK_atlas.mat

fid=fopen('Subjects_list_all.txt','r');
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);

data_location='E:\RAM data set\RAM_Public_Data_all\FR1_FS';

patients=r_sublist;
overlay=4;
laterality=5;
views=4;
inflationstep=5;
electrodeDensity=0;
electrodeDisplay=0;
thr=0.5;
color_range=[-2 2];
statistic_flag=0;

file_suffix='Final_std_5_20_fs_together';
Montage_suffix='final_roi_fs';

all_sub_index=1:251;
atlas_index_all=[3 4 6 8 9 10 11 12 13 15 16 17 22 24 26 28 29];
atlas_index_all=2:36;
% atlas_index_all_rest=1:36;
% atlas_index_all_rest(atlas_index_all)=[];
% atlas_index_all_rest(1)=[];
for atlas_index=atlas_index_all
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
laterality=4;
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



%% ROI mapping (parcel level)
clc;clear;

rootfolder = 'E:\RAM data set\RAM_Public_Data_all';
cd(rootfolder)

load r1_all.mat
load fs_DK_atlas.mat

fid=fopen('Subjects_list_all.txt','r');
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);

data_location='E:\RAM data set\RAM_Public_Data_all\FR1_FS';

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
atlas_all_index=[4 8 12 13 28 30];
displayElectrodesInflated_new_fs_Atlas_all(data_location, patients, overlay, laterality, views,electrodeDensity, thr) 
close all



network_index_all=[91 92 93];
for atlas_index=network_index_all
    displayElectrodesInflated_new_fs_Atlas(data_location, patients, overlay, laterality, views, inflationstep,electrodeDensity, electrodeDisplay ,thr,color_range, file_suffix,Montage_suffix,statistic_flag,all_sub_index,atlas_index)
    atlas_index
    close all;
end



%% ROI mapping (neighbor ROIs)
overlay=4;
laterality=0;
views=4;
electrodeDensity=0;
thr=0.5;

network_index_all=randi([2,201],20,1)';
load Both_roi_neighbor_matrix_400_2.mat
atlas_suffix='Yeo_400_neighbor';

for atlas_index=network_index_all
    displayElectrodesInflated_new_fs_Atlas_n(data_location, overlay, laterality, views,electrodeDensity, thr,atlas_suffix,roi_neighbor_matrix,atlas_index)
    atlas_index
    close all;
end


%% stats
% file_suffix='Final_std_20_20_fs_together_ttest';
% Montage_suffix='final_roi_fs';

file_suffix='Final_std_20_20_fs_together_xflip_ttest';
Montage_suffix='final_roi_fs_xflip';

patients=r_sublist;
overlay=4;
laterality=0;
views=4;
inflationstep=5;
electrodeDensity=0;
electrodeDisplay=0;
thr=0.1;
stat_flag=1; % 1) uncorrected 0.05 / 2) FDR 0.05 / 3) FDR 0.2 / 4) 
color_range=[-5 5];

displayElectrodesInflated_new_fs_stat(data_location, patients, overlay, laterality, views, inflationstep,electrodeDensity, electrodeDisplay ,thr,stat_flag,color_range, file_suffix,Montage_suffix)



atlas_index_all=2:36;
atlas_index_all(find(atlas_index_all==5))=[];
% atlas_index_all=[3 4 6 8 9 10 11 12 13 15 16 17 22 24 26 28 29];
plot_test=nan(3,36,251,98);
[L_vertices, L_label, L_colortable] = read_annotation('lh.aparc.annot');
[R_vertices, R_label, R_colortable] = read_annotation('rh.aparc.annot');

for time_i=1:98
    load(['R_vertex_values_Final_std_20_20_fs_together_' num2str(time_i)]);
    for atlas_index=atlas_index_all
        plot_test(1,atlas_index,:,time_i)=nanmean(vertex_values_frame(find(R_label==R_colortable.table(atlas_index,5)),:),1);
    end
    time_i
end

for time_i=1:98
    load(['L_vertex_values_Final_std_20_20_fs_together_' num2str(time_i)]);
    for atlas_index=atlas_index_all
        plot_test(2,atlas_index,:,time_i)=nanmean(vertex_values_frame(find(L_label==L_colortable.table(atlas_index,5)),:),1);
    end
    time_i
end

plot_test=nan(36,251,98);
for time_i=1:98
    load(['R_vertex_values_Final_std_20_20_fs_together_xflip_' num2str(time_i)]);
    for atlas_index=atlas_index_all
        plot_test(3,atlas_index,:,time_i)=nanmean(vertex_values_frame(find(R_label==R_colortable.table(atlas_index,5)),:),1);
    end
    time_i
end


%% ROI plot zscore
% Montage_suffix='final_roi_fs';
% Montage_suffix='final_roi_fs_xflip';
load fs_DK_atlas.mat
atlas_index_all=[3 4 6 8 9 10 11 12 13 15 16 17 22 24 26 28 29];
for atlas_index=atlas_index_all
    
    
    figure('position',[0 0 8000 8000]);
    set(gcf, 'color', [1 1 1]);
    set(gcf,'Visible','off');

    % plot
    subplot(3,1,1)
    plot_buff=squeeze(plot_test(1,atlas_index,:,:));
    plot_buff(find(plot_buff(:,1)==0),:)=nan;
    plot_areaerrorbar(plot_buff)
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
    ylabel({'Gamma power z-score'});
    xlabel('Time (ms)');

    subplot(3,1,2)
    plot_buff=squeeze(plot_test(2,atlas_index,:,:));
    plot_buff(find(plot_buff(:,1)==0),:)=nan;
    plot_areaerrorbar(plot_buff)
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
    ylabel({'Gamma power z-score'});
    xlabel('Time (ms)');

    subplot(3,1,3)
    plot_buff=squeeze(plot_test(3,atlas_index,:,:));
    plot_buff(find(plot_buff(:,1)==0),:)=nan;
    plot_areaerrorbar(plot_buff)
    title(['Flip Hemispheres'])
    set(gca,'xtick',[1 20 40 60 80]);
    set(gca,'xticklabel',[-500 0 500 1000 1500]);
    set(gca,'xlim',[1 98])
    set(gca,'ytick',[-6 -2 2 6 ]);
    set(gca,'yticklabel',[-6 -2 2 6 ]);
    set(gca,'ylim',[-8 10])
    set(gca,'Fontsize',15);
    line([1,98],[2,2],'Color','r','LineWidth',0.5)
    line([1,98],[-2,-2],'Color','r','LineWidth',0.5)
    ylabel({'Gamma power z-score'});
    xlabel('Time (ms)');
    
    %% save figure
    [ax,h1]=suplabel(['ROI : ' DK_atlas_names{atlas_index}] ,'t');
    set(h1,'FontSize',25)
    cd('E:\RAM data set\RAM_Public_Data_all\FR1_FS\Storage_filtering_overlap\Atlas_DK')
    saveas(gcf,[num2str(atlas_index) '_' DK_atlas_names{atlas_index} '_zscore_plot.png'])
    close all;
    
    
end


%% zscore plot by ROIs
cd('E:\RAM data set\RAM_Public_Data_all\FR1_FS\Storage_filtering_overlap\Data_storage')
file_suffix='Final_std_20_fs';
Montage_suffix='final_roi_fs_xflip';
load fs_DK_atlas.mat

zscore_folder = ['E:\RAM data set\RAM_Public_Data_all\FR1_FS\Storage_filtering_overlap\zscore_plot\' Montage_suffix ];
mkdir(zscore_folder)

[L_vertices, L_label, L_colortable] = read_annotation('lh.aparc.annot');
[R_vertices, R_label, R_colortable] = read_annotation('rh.aparc.annot');

% plot_test_rec=nan(36,251,118);
% plot_test_Nonrec=nan(36,251,118);
% plot_test_together=nan(36,251,118);
% plot_test_diff=nan(36,251,118);
plot_test=nan(4,36,251,118);
task_name={'rec_xflip','Nonrec_xflip','together_xflip','diff_xflip'};

for time_i=1:118
    for task_i=1:4
        vertex_values_frame=[];
        load([file_suffix '_' task_name{task_i} '/' Montage_suffix '/L_vertex_values_' file_suffix '_' task_name{task_i} '_' num2str(time_i)]);
        for atlas_index=1:36
    %         plot_test(3,atlas_index,:,time_i)=nanmean(vertex_values_frame(find(R_label==R_colortable.table(atlas_index,5)),:),1);
            plot_test(task_i,atlas_index,:,time_i)=nanmean(vertex_values_frame(find(L_label==L_colortable.table(atlas_index,5)),:),1);

        end
    end
    time_i
end
save([file_suffix '_' Montage_suffix '_roi_val.mat'],'plot_test');
load([file_suffix '_' Montage_suffix '_roi_val.mat']);

task_name_fig={'Recall','NonRecall','Recall + NonRecall','Recall - NonRecall'};
for atlas_index=2:36
    max_val=max(max(max(squeeze(plot_test(:,atlas_index,:,:)))));
    min_val=min(min(min(squeeze(plot_test(:,atlas_index,:,:)))));
    if ~isnan(max_val)
    for task_i=1:4
        figure('position',[0 0 2000 2000]);
        set(gcf, 'color', [1 1 1]);
        set(gcf,'Visible','on');
        plot_test_b=squeeze(plot_test(task_i,atlas_index,:,:));
        plot_test_b(plot_test_b==0)=nan;
        plot_areaerrorbar(plot_test_b)
        
        set(gca,'xtick',[1 20 40 60 80 100 118]);
        set(gca,'xticklabel',[-500 0 500 1000 1500 2000 2500]);
        set(gca,'xlim',[1 118])
        set(gca,'box','off')

    %         set(gca,'ytick',[-6 -2 0 2 6 ]);
    %         set(gca,'yticklabel',[-6 -2 0 2 6]);
        set(gca,'ylim',[min_val max_val])
        set(gca,'Fontsize',30);

        x1=xline(20,'-.r',{'Word','ON'},'LabelVerticalAlignment','top','LabelHorizontalAlignment','center','LabelOrientation','horizontal','LineWidth',2,'Fontsize',30);
        x2=xline(84,'-.r',{'Word','OFF'},'LabelVerticalAlignment','top','LabelHorizontalAlignment','center','LabelOrientation','horizontal','LineWidth',2,'Fontsize',30);
        y1=yline(0,'k','LineWidth',2,'Fontsize',30);
        
        %         line([20,20],[-8,10],'Color','r','LineWidth',1)
    %         line([84,84],[-8,10],'Color','r','LineWidth',1)
        ylabel({'z-score'});
        xlabel('Time (ms)');

        %% save figure
%         [ax,h1]=suplabel(['ROI : ' DK_atlas_names{atlas_index} ', Condition : ' task_name_fig{task_i} ],'t');
        [ax,h1]=suplabel([task_name_fig{task_i} ],'t');
        set(h1,'FontSize',30)
        saveas(gcf,[zscore_folder '\' num2str(atlas_index) '_' DK_atlas_names{atlas_index} '_' task_name_fig{task_i} '_plot.png'])
        close all;
    end
    end
end

%% combine images
cd('E:\RAM data set\RAM_Public_Data_all\FR1_FS\Storage_filtering_overlap')
file_suffix='Final_std_20_fs';
Montage_suffix='final_roi_fs_xflip';
task_name_fig={'Recall','NonRecall','Recall + NonRecall','Recall - NonRecall'};
load fs_DK_atlas.mat

atlas_index_all=[2:4 6:36];
task_i_all=[1 3];
for atlas_index=atlas_index_all
for task_i=task_i_all
cd('E:\RAM data set\RAM_Public_Data_all\FR1_FS\Storage_filtering_overlap')
atlas_image=['Atlas\Final_std_5_20_fs_together\' num2str(atlas_index) '_' DK_atlas_names{atlas_index} '\combined_viewsL_1.tif' ];
imagea=imread(atlas_image);
zscore_image=['zscore_plot\' Montage_suffix '\' num2str(atlas_index) '_' DK_atlas_names{atlas_index} '_' task_name_fig{task_i} '_plot.png' ];
imageb=imread(zscore_image);
imageb_r=imresize(imageb,[651 size(imagea,2)]);

zscore_image=['zscore_plot\' Montage_suffix '\' num2str(atlas_index) '_' DK_atlas_names{atlas_index} '_' task_name_fig{task_i+1} '_plot.png' ];
imagec=imread(zscore_image);
imagec_r=imresize(imagec,[651 size(imagea,2)]);

imagebc_r=cat(1,imageb_r,255*ones(15,size(imagea,2),3),imagec_r);
multi=cat(2,imagea,imagebc_r);
imshow(multi);

text(550,1150,['ROI : ' DK_atlas_names{atlas_index}],'Fontsize',30)

set(gca,'position',[0 0 1 1],'units','normalized')
axis tight
print(gcf,[num2str(atlas_index) '_' DK_atlas_names{atlas_index} '_' num2str(task_i) '_plot.tif'],'-dtif','-r500')
close all;

end
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

