function PipelineForPermutationTest(data_location,r_sublist,all_sub_index, ...
    Montage_suffix,file_suffix,storage_suffix,flip_flag,zscore_mapping_flag,Cluster_based_statistics_flag, ...
    zscore_mapping_stats_flag,atlas_index,cluster_stats_all, ...
    z_color_range)

%% common parameter setting
patients=r_sublist;
overlay=4;
figure_type=4; % 0 or 4
views=4;
electrodeDensity=0;
electrodeDisplay=0;
mean_flag=2;
% all_sub_index=1:251;
% all_sub_index=find(session_n_all==1)'; % test

if flip_flag==1
    sides_num = 1;
    flip_prefix='xflip';
else
    sides_num = 2;
    flip_prefix='both';
end

image_storage_folder = [data_location '\Storage_filtering_overlap_' storage_suffix '\Image_storage\' file_suffix '\' Montage_suffix];
data_storage_folder=[data_location '\Storage_filtering_overlap_' storage_suffix '\Data_storage\' file_suffix '\' Montage_suffix];
atlas_storage_folder=[data_location '\Storage_filtering_overlap_' storage_suffix '\Data_storage\' file_suffix '\' Montage_suffix '\Atlas_info'];

if flip_flag == 0
    %% zscore mapping on both hemispheres
    if zscore_mapping_flag == 1
        thr=0;
%         colormap_mat='cmap_hunki_10to15.mat';
%         colormap_mat='cmap_hunki_new.mat';
        colormap_mat='cmap_hunki_zmap_raw.mat';
        colormap_label='Gamma power z-score';
        
        % make figures
        displayElectrodesInflated_new_fs(data_location, patients, overlay, figure_type, ...
            views, electrodeDensity, electrodeDisplay ,thr,mean_flag,z_color_range, ...
            file_suffix,storage_suffix,Montage_suffix,all_sub_index,colormap_mat)
        if electrodeDensity~=1 | electrodeDensity~=1
            display_dynamic_HK(image_storage_folder,flip_flag,colormap_mat,z_color_range,colormap_label)
        end
        close all;
    end
    
    %% Cluster-based statistics on both hemispheres
    if Cluster_based_statistics_flag == 1
        RAM_average_zscore_by_ROI(data_storage_folder,atlas_storage_folder,file_suffix,atlas_index,flip_flag)
        for cluster_stats_i=1:size(cluster_stats_all,2)
            RAM_Cluster_based_stats(atlas_storage_folder,atlas_index,file_suffix, ...
                flip_prefix,cluster_stats_all{cluster_stats_i},flip_flag)
        end
    end
    
    %% zscore mapping with stats
    if zscore_mapping_stats_flag == 1
        thr=0;
        color_range=z_color_range;
        colormap_mat='cmap_hunki_zmap.mat';
        colormap_label='Gamma power z-score';
        % make figures
        for cluster_stats_i=1:size(cluster_stats_all,2)
            for time_i=5 %[3 5 7]
                displayElectrodesInflated_new_fs_stat_cluster_zscore(data_location,patients, overlay, figure_type, views,thr,color_range,colormap_mat, ...
                    file_suffix,storage_suffix,Montage_suffix,cluster_stats_all{cluster_stats_i},atlas_storage_folder,atlas_index,flip_prefix,time_i,all_sub_index,mean_flag)
                r_time=(time_i+1)*25;
                stats_storage_folder=[image_storage_folder '\' atlas_index '_' cluster_stats_all{cluster_stats_i} '_' num2str(r_time) 'ms_zscore'];
                display_dynamic_HK(stats_storage_folder,flip_flag,colormap_mat,color_range,colormap_label)
                close all;
            end
        end
    end

elseif flip_flag == 1
    % motage prefix should be changed
    % add 'flip' prefix

    %% zscore mapping on fliped hemisphere
    %% Cluster-based statistics on fliped hemisphere
    %% tscore mapping on fliped hemisphere
end
end


function RAM_average_zscore_by_ROI(data_location,atlas_storage_folder,file_suffix,atlas_index,flip_flag)

cd(data_location)

start_side = 1;
if flip_flag==1
    sides_num = 1;
    flip_prefix='xflip';
else
    sides_num = 2;
    flip_prefix='both';
end

total_num_frames=118;

switch atlas_index
    case '1000'
        mean_vertex_values=nan(sides_num,501,total_num_frames,251);
    case '400'
        mean_vertex_values=nan(sides_num,201,total_num_frames,251);
    case '200'
        mean_vertex_values=nan(sides_num,101,total_num_frames,251);
    case 'DK'
        mean_vertex_values=nan(sides_num,36,total_num_frames,251);
    otherwise
        disp('Error, no such atlas')
end

for side_index = start_side:sides_num    % Generate frames for L and R
    
% Determine what side we are working on
    if side_index == 1
        side = 'L';
        switch atlas_index
            case '1000'
                [vertices, label, colortable]=read_annotation('lh.Schaefer2018_1000Parcels_7Networks_order.annot');
            case '400'
                [vertices, label, colortable]=read_annotation('lh.Schaefer2018_400Parcels_7Networks_order.annot');
            case '200'
                [vertices, label, colortable]=read_annotation('lh.Schaefer2018_200Parcels_7Networks_order.annot');
            case 'DK'
                [vertices, label, colortable] = read_annotation('lh.aparc.annot');
            otherwise
                disp('Error, no such atlas.')
        end
    else
        side = 'R';
        switch atlas_index
            case '1000'
                [vertices, label, colortable]=read_annotation('rh.Schaefer2018_1000Parcels_7Networks_order.annot');
            case '400'
                [vertices, label, colortable]=read_annotation('rh.Schaefer2018_400Parcels_7Networks_order.annot');
            case '200'
                [vertices, label, colortable]=read_annotation('rh.Schaefer2018_200Parcels_7Networks_order.annot');
            case 'DK'
                [vertices, label, colortable] = read_annotation('rh.aparc.annot');
            otherwise
                disp('Error, no such atlas.')
        end
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
        
        disp(frame_index)
    end
    
end
mkdir(atlas_storage_folder)
save([atlas_storage_folder '/roi_time_subject_vertex_' atlas_index '_' file_suffix '_' flip_prefix '.mat'],'mean_vertex_values')

end

function RAM_Cluster_based_stats(atlas_storage_folder,atlas_index,file_suffix,flip_prefix,stats_method,flip_flag)

cd(atlas_storage_folder)
load(['roi_time_subject_vertex_' atlas_index '_' file_suffix '_' flip_prefix '.mat'])

if flip_flag==1
    both_mean_vertex_values=[];
    both_mean_vertex_values=squeeze(mean_vertex_values(1,:,:,:));
    load(['L_roi_neighbor_matrix_' atlas_index '_2.mat'])
else
    L_buff=[];R_buff=[];both_mean_vertex_values=[];
    L_buff=squeeze(mean_vertex_values(1,:,:,:));
    R_buff=squeeze(mean_vertex_values(2,:,:,:));
    both_mean_vertex_values=cat(1,L_buff,R_buff);
    load(['Both_roi_neighbor_matrix_' atlas_index '_2.mat'])
end
roi_n=size(both_mean_vertex_values,1);
time_n=size(both_mean_vertex_values,2);
sub_n=size(both_mean_vertex_values,3);
base_mean=squeeze(mean(both_mean_vertex_values(:,1:19,:),2));
for test_i=1:roi_n
    both_mean_vertex_values(test_i,:,:)=squeeze(both_mean_vertex_values(test_i,:,:))-base_mean(test_i,:);
end


switch stats_method
    
    case 'cluster_sumt_005'                                
        [pval, t_orig, clust_info, seed_state, est_alpha, mn_clust_mass]=clust_perm1_iceeg_sumt(both_mean_vertex_values,roi_neighbor_matrix,2000,0.05,0,0.05,2,[],0);
        save(['cluster_sumt_005_' atlas_index '_' file_suffix '_' flip_prefix  '.mat'],'pval','t_orig','clust_info','seed_state','est_alpha','mn_clust_mass' ...
        ,'both_mean_vertex_values');

    case 'cluster_sumt_001'
        [pval, t_orig, clust_info, seed_state, est_alpha, mn_clust_mass]=clust_perm1_iceeg_sumt(both_mean_vertex_values,roi_neighbor_matrix,2000,0.05,0,0.01,2,[],0);
        save(['cluster_sumt_001_' atlas_index '_' file_suffix '_' flip_prefix  '.mat'],'pval','t_orig','clust_info','seed_state','est_alpha','mn_clust_mass' ...
        ,'both_mean_vertex_values');
    
    case 'cluster_sumt_0005'
        [pval, t_orig, clust_info, seed_state, est_alpha, mn_clust_mass]=clust_perm1_iceeg_sumt(both_mean_vertex_values,roi_neighbor_matrix,2000,0.05,0,0.005,2,[],0);
        save(['cluster_sumt_0005_' atlas_index '_' file_suffix '_' flip_prefix  '.mat'],'pval','t_orig','clust_info','seed_state','est_alpha','mn_clust_mass' ...
        ,'both_mean_vertex_values');   
    
    otherwise
        disp('Error, no such method')
end

end




