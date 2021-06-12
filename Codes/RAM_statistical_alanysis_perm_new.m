
clc;clear;

atlas_folder='E:\RAM data set\RAM_Public_Data_all\FR1_FS\Atlas_info';
atlas_index='400'; % 'DK' % '200'
flip_flag=0; % 0:both 1:flip x

Montage_suffix='final_roi_fs';
file_suffix='Final_std_20_fs_together';
% file_suffix='Final_std_20_fs_rec';
% file_suffix='Final_std_20_fs_Nonrec';
% file_suffix='Final_std_20_fs_diff';

% Montage_suffix='final_roi_fs_xflip';
% file_suffix='Final_std_20_fs_together_xflip';
% file_suffix='Final_std_20_fs_rec_xflip';
% file_suffix='Final_std_20_fs_Nonrec_xflip';
% file_suffix='Final_std_20_fs_diff_xflip';

data_folder=['E:\RAM data set\RAM_Public_Data_all\FR1_FS\Storage_filtering_overlap\Data_storage\' file_suffix '\' Montage_suffix];
cd(data_folder)

start_side = 1;

if flip_flag==1
    sides_num = 1;
    flip_prefix='xflip';
else
    sides_num = 2;
    flip_prefix='both';
end

total_num_frames=118;

%% test 1 (vertex)
% switch atlas_index
%     case '400'
%         mean_vertex_values=nan(sides_num,201,total_num_frames,251);
%     case '200'
%         mean_vertex_values=nan(sides_num,101,total_num_frames,251);
%     case 'DK'
%         mean_vertex_values=nan(sides_num,36,total_num_frames,251);
%     otherwise
%         disp('Error, no such atlas')
% end
% 
% for side_index = start_side:sides_num    % Generate frames for L and R
%     
% % Determine what side we are working on
%     if side_index == 1
%         side = 'L';
%         switch atlas_index
%             case '400'
%                 [vertices, label, colortable]=read_annotation('lh.Schaefer2018_400Parcels_7Networks_order.annot');
%             case '200'
%                 [vertices, label, colortable]=read_annotation('lh.Schaefer2018_200Parcels_7Networks_order.annot');
%             case 'DK'
%                 [vertices, label, colortable] = read_annotation('lh.aparc.annot');
%             otherwise
%                 disp('Error, no such atlas.')
%         end
%     else
%         side = 'R';
%         switch atlas_index
%             case '400'
%                 [vertices, label, colortable]=read_annotation('rh.Schaefer2018_400Parcels_7Networks_order.annot');
%             case '200'
%                 [vertices, label, colortable]=read_annotation('rh.Schaefer2018_200Parcels_7Networks_order.annot');
%             case 'DK'
%                 [vertices, label, colortable] = read_annotation('rh.aparc.annot');
%             otherwise
%                 disp('Error, no such atlas.')
%         end
%     end
% 
%     for frame_index = 1:total_num_frames
%         clearvars vertex_values_frame
%         load([side '_vertex_values_' file_suffix '_' num2str(frame_index) '.mat']);
%         vertex_values_frame(find(vertex_values_frame==0))=nan;
%         
%         for atlas_i=2:size(colortable.table,1)
%             vertex_values_frame_buff=[];
%             vertex_values_frame_buff=vertex_values_frame(label == colortable.table(atlas_i,5),:);
%             mean_vertex_values(side_index,atlas_i,frame_index,:)=nanmean(vertex_values_frame_buff,1);
%         end
%         
%         frame_index
%     end
%     
% end

cd(atlas_folder)
% save(['./roi_time_subject_vertex_' atlas_index '_' file_suffix '_' flip_prefix '.mat'],'mean_vertex_values')
load(['./roi_time_subject_vertex_' atlas_index '_' file_suffix '_' flip_prefix '.mat'])

if flip_flag==1
    both_mean_vertex_values=[];
    both_mean_vertex_values=squeeze(mean_vertex_values(1,:,:,:));
    load(['./L_roi_neighbor_matrix_' atlas_index '_2.mat'])
else
    L_buff=[];R_buff=[];both_mean_vertex_values=[];
    L_buff=squeeze(mean_vertex_values(1,:,:,:));
    R_buff=squeeze(mean_vertex_values(2,:,:,:));
    both_mean_vertex_values=cat(1,L_buff,R_buff);
    load(['./Both_roi_neighbor_matrix_' atlas_index '_2.mat'])
end

% [pval, t_orig, clust_info, seed_state, est_alpha, mn_clust_mass]=clust_perm1_iceeg_sumt(both_mean_vertex_values,roi_neighbor_matrix);
% save(['cluster_sumt_' atlas_index '_' file_suffix '_' flip_prefix  '.mat'],'pval','t_orig','clust_info','seed_state','est_alpha','mn_clust_mass' ...
%     ,'mean_vertex_values');
% 
% [pval, t_orig, clust_info, seed_state, est_alpha, mn_clust_mass]=clust_perm1_iceeg_count(both_mean_vertex_values,roi_neighbor_matrix);
% save(['cluster_count_' atlas_index '_' file_suffix '_' flip_prefix  '.mat'],'pval','t_orig','clust_info','seed_state','est_alpha','mn_clust_mass' ...
%     ,'mean_vertex_values');

% [pval, t_orig, clust_info, seed_state, est_alpha, mn_clust_mass]=clust_perm1_iceeg_sumt_rand(both_mean_vertex_values,roi_neighbor_matrix);
% save(['cluster_sumt_rand_' atlas_index '_' file_suffix '_' flip_prefix  '.mat'],'pval','t_orig','clust_info','seed_state','est_alpha','mn_clust_mass' ...
%     ,'mean_vertex_values');
load(['cluster_sumt_rand_' atlas_index '_' file_suffix '_' flip_prefix  '.mat'])

% [pval, t_orig, clust_info, seed_state, est_alpha, mn_clust_mass]=clust_perm1_iceeg_count_rand(both_mean_vertex_values,roi_neighbor_matrix);
% save(['cluster_count_rand_' atlas_index '_' file_suffix '_' flip_prefix  '.mat'],'pval','t_orig','clust_info','seed_state','est_alpha','mn_clust_mass' ...
%     ,'mean_vertex_values');



%% stats
rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

load r1_all.mat

fid=fopen('Subjects_list_all.txt','r');
cnt=1;
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);
data_location='E:\RAM data set\RAM_Public_Data_all\FR1_FS';



patients=r_sublist;
overlay=4;
laterality=4;
views=4;
electrodeDensity=0;
electrodeDisplay=0;
thr=1;
% atlas_name_all={'cluster_count','cluster_count_rand','cluster_sumt','cluster_sumt_rand'}; % 
atlas_name_all={'cluster_sumt_rand'}; % 
color_range=[-6 6];

displayElectrodesInflated_new_fs_stat_cluster(data_location,patients, overlay, laterality, views,electrodeDensity, electrodeDisplay ,thr,color_range, ...
    file_suffix,Montage_suffix,atlas_name_all{1},atlas_folder,atlas_index,flip_prefix)
close all;





