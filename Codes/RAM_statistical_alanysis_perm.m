
clc;clear;
% file_suffix='Final_std_20_20_fs_together';
% Montage_suffix='final_roi_fs';
file_suffix='Final_std_20_20_fs_together_xflip';
Montage_suffix='final_roi_fs_xflip';

data_folder=['E:\RAM data set\RAM_Public_Data_all\FR1_FS\Storage_filtering_overlap\Data_storage\' file_suffix '\' Montage_suffix];
cd(data_folder)

%% test 1 (elec)

start_side = 1;
sides_num = 2;
total_num_frames=118;

% file_suffix='Final_std_20_fs_together';
% Montage_suffix='final_roi_fs';

file_suffix='Final_std_20_fs_together_xflip';
% Montage_suffix='final_roi_fs_xflip';

mean_vertex_values=nan(2,201,98,251);
for side_index = start_side:sides_num    % Generate frames for L and R
    
    if side_index == 1
        side = 'L';
        [vertices, label, colortable]=read_annotation('lh.Schaefer2018_400Parcels_7Networks_order.annot');
%         [vertices, label, colortable] = read_annotation('lh.aparc.annot');
    else
        side = 'R';
        [vertices, label, colortable]=read_annotation('lh.Schaefer2018_400Parcels_7Networks_order.annot');
%         [vertices, label, colortable] = read_annotation('rh.aparc.annot');
    end

    for frame_index = 1:total_num_frames
        clearvars vertex_values_frame
        load([side '_vertex_values_' file_suffix '_' num2str(frame_index) '.mat']);
        vertex_values_frame(find(vertex_values_frame==0))=nan;
        
        for atlas_i=2:size(colortable.table,1)
            vertex_values_frame_buff=[];
            vertex_values_frame_buff=vertex_values_frame(label == colortable.table(atlas_i,5),:);
            mean_vertex_values(side_index,atlas_i,frame_index,:)=nanmean(vertex_values_frame_buff,1);
        end
        
        frame_index
    end
    
end

% load('cmap_hunki_new');
% figure('position',[0 0 8000 8000]);
% set(gcf, 'color', [1 1 1]);
% set(gcf,'Visible','on');
% 
% buffer_Left=squeeze(mean_vertex_values(:,:));
% % buffer_Left(isnan(buffer_Left(:,1)),:)=[];
% buffer_Left(isnan(buffer_Left(:,1)),:)=0;
% imagesc(buffer_Left)
% title(['Flip Hemisphere'])
% set(gca,'ydir','normal');
% set(gca,'xtick',[1 20 40 60 80]);
% set(gca,'xticklabel',[-500 0 500 1000 1500]);
% set(gca,'xlim',[1 98])
% set(gca,'Fontsize',15);
% line([20,20],ylim,'Color','r','LineWidth',0.5)
% ylabel('Subject');
% xlabel('Time (ms)');
% colormap(mycmap);
% colorbar
% caxis([-2 2])

% mean_vertex_values(:,11:end,:)=[];

%% version 1
mean_vertex_values(isnan(mean_vertex_values))=0;
test_values=zeros(size(mean_vertex_values,1),size(mean_vertex_values,2),size(mean_vertex_values,3));
[clusters, p_values, t_sums, permutation_distribution ] = permutest_ram(mean_vertex_values,test_values,true,0.05,10000,true);


rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)
save('cluster_based_400_vertex.mat','clusters','p_values','t_sums','permutation_distribution' ...
    ,'mean_vertex_values');
    

histogram(abs(permutation_distribution),100)
% H = histogram([log(abs(vertex_values_sumL(:))); log(abs(vertex_values_sumR(:)))], bins,...
%     'Normalization', 'probability');
xlabel('Permuted summed t value')
ylabel('Counts')
title('Histogram of permuted summed t value (10000 times)')


% find(p_values<0.05)
% 
% 
[a1 b1]=ind2sub([size(colortable.table,1) 98],clusters{1,1});
[a2 b2]=ind2sub([size(colortable.table,1) 98],clusters{1,2});
[a3 b3]=ind2sub([size(colortable.table,1) 98],clusters{1,3});
[a4 b4]=ind2sub([size(colortable.table,1) 98],clusters{1,4});




%% stats
clc;clear;
rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

load r1_all.mat
% load fs_DK_atlas.mat
fid=fopen('Subjects_list_all.txt','r');
cnt=1;
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);

data_location='E:\RAM data set\RAM_Public_Data_all\FR1_FS';


% file_suffix='Final_std_20_20_fs_together_ttest';
% Montage_suffix='final_roi_fs';

file_suffix='Final_std_20_20_fs_together_xflip_Cluster';
Montage_suffix='final_roi_fs_xflip';

patients=r_sublist;
overlay=4;
laterality=0;
views=4;
inflationstep=5;
electrodeDensity=0;
electrodeDisplay=0;
thr=0.1;
atlas_name=3; % 
color_range=[-5 5];

displayElectrodesInflated_new_fs_stat_cluster(data_location, patients, overlay, laterality, views, inflationstep,electrodeDensity, electrodeDisplay ,thr,color_range, file_suffix,Montage_suffix,atlas_name)




