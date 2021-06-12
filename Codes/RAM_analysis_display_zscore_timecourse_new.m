
%% zscore plot by ROIs
clc;clear;
rootfolder = 'E:\RAM data set\RAM_Public_Data_all';
cd([rootfolder '\FR1_FS\Storage_filtering_overlap_base_agg\Data_storage'])

file_suffix='Final_std_20_fs';
Montage_suffix='final_roi_fs';
load fs_DK_atlas.mat
flip_flag=0;

if flip_flag==1
    sides_num = 1;
    flip_prefix='xflip';
else
    sides_num = 2;
    flip_prefix='both';
end

% atlas_index_all=[4 30 8 12 13 15 28 29];
atlas_index_all_new={4,30,8,12,13,15,28,29,[19 20 21],[9,32],[10 16],17,[31 35]};

atlas_index_all_new_name={'Caudal Middle Frontal','Superior Parietal'...
    ,'Fusiform Gyrus','Lateral Occipital','Lateral Orbitofrontal','Medial Orbitofrontal' ...
    ,'Rostral Middle Frontal','Superior Frontal' ...
    ,'Inferior Frontal','Inferior Parietal','Medial Temporal' ...
    ,'Parahippocampal','Superior Temporal'};

zscore_folder = [rootfolder '\FR1_FS\Storage_filtering_overlap_base_agg\zscore_plot\' Montage_suffix ];
mkdir(zscore_folder)

atlas_color=distinguishable_colors(13);

% overlay=4;laterality=0;views=4;
% displayElectrodesInflated_new_fs_Atlas_all_final(zscore_folder, overlay, laterality, views, atlas_index_all_new,atlas_color) 
% close all


% % atlas_color=distinguishable_colors(200);
% overlay=4;laterality=4;views=4;
% % displayElectrodesInflated_new_fs_Atlas_all_final(zscore_folder, overlay, laterality, views, atlas_index_all_new,atlas_color) 
% displayElectrodesInflated_new_fs_Atlas_all(zscore_folder, overlay, laterality, views,0,0, 1, 21)
% 
% close all


if flip_flag ==0 
    task_name={'rec','Nonrec','together','diff'};
    plot_test=nan(2,4,size(atlas_index_all_new,2),251,98);
    plot_test_weighted=nan(2,4,size(atlas_index_all_new,2),251,98);
    plot_test_elec=nan(2,4,size(atlas_index_all_new,2),251,98);
    plot_test_elec_sig=nan(2,4,size(atlas_index_all_new,2),251,98);
    plot_test_sig=zeros(2,4,size(atlas_index_all_new,2),251,98);
    plot_test_stats_ratio=nan(2,4,size(atlas_index_all_new,2),98);
    plot_test_stats_avg_t=nan(2,4,size(atlas_index_all_new,2),98);
    plot_test_stats_sum_t=nan(2,4,size(atlas_index_all_new,2),98);
    plot_test_stats_avg_t_sig=nan(2,4,size(atlas_index_all_new,2),98);
    plot_test_stats_sum_t_sig=nan(2,4,size(atlas_index_all_new,2),98);
    
    plot_test_ROI_size=nan(2,4,size(atlas_index_all_new,2),98);
    
%     [L_vertices, L_label, L_colortable] = read_annotation('lh.aparc.annot');
%     [R_vertices, R_label, R_colortable] = read_annotation('rh.aparc.annot');
    [L_vertices, L_label, L_colortable] = read_annotation('lh.aparc.DKTatlas40.annot');
    [R_vertices, R_label, R_colortable] = read_annotation('rh.aparc.DKTatlas40.annot');
    [~, L_400_label, L_400_colortable]=read_annotation('lh.Schaefer2018_400Parcels_7Networks_order.annot');
    [~, R_400_label, R_400_colortable]=read_annotation('rh.Schaefer2018_400Parcels_7Networks_order.annot');
    [size_L_400_label, index_L_400_label] = hist(L_400_label,unique(L_400_label));
    [size_R_400_label, index_R_400_label] = hist(R_400_label,unique(R_400_label));
    
elseif flip_flag ==1
    task_name={'rec_xflip','Nonrec_xflip','together_xflip','diff_xflip'};
    plot_test=nan(1,4,size(atlas_index_all_new,2),251,98);
    plot_test_weighted=nan(1,4,size(atlas_index_all_new,2),251,98);
    plot_test_elec=nan(1,4,size(atlas_index_all_new,2),251,98);
    plot_test_elec_sig=nan(1,4,size(atlas_index_all_new,2),251,98);
    plot_test_sig=zeros(1,4,size(atlas_index_all_new,2),251,98);
    plot_test_stats_ratio=nan(1,4,size(atlas_index_all_new,2),98);
    plot_test_stats_avg_t=nan(1,4,size(atlas_index_all_new,2),98);
    plot_test_stats_sum_t=nan(1,4,size(atlas_index_all_new,2),98);
    plot_test_stats_avg_t_sig=nan(1,4,size(atlas_index_all_new,2),98);
    plot_test_stats_sum_t_sig=nan(1,4,size(atlas_index_all_new,2),98);
    
    plot_test_ROI_size=nan(1,4,size(atlas_index_all_new,2),98);
    
%     [L_vertices, L_label, L_colortable] = read_annotation('lh.aparc.annot');
    [L_vertices, L_label, L_colortable] = read_annotation('lh.aparc.DKTatlas40.annot');
    [~, L_400_label, L_400_colortable]=read_annotation('lh.Schaefer2018_400Parcels_7Networks_order.annot');
    [size_L_400_label, index_L_400_label] = hist(L_400_label,unique(L_400_label));
end


%% data extraction
% for task_i=1:4
%     vertex_values_sum_all=[];
%     load([file_suffix '_' task_name{task_i} '/' Montage_suffix '/Atlas_info/L_vertex_values_400_' file_suffix '_' task_name{task_i} '_both_150.mat']);
%     vertex_values_frame_L=vertex_values_sum_all;
%     vertex_values_sum_all=[];
%     load([file_suffix '_' task_name{task_i} '/' Montage_suffix '/Atlas_info/R_vertex_values_400_' file_suffix '_' task_name{task_i} '_both_150.mat']);
%     vertex_values_frame_R=vertex_values_sum_all;
% 
%     for time_i=1:98
%         if flip_flag == 0
%             %% Left
%             vertex_values_frame=[];
%             vertex_values_frame_e=[];
%             vertex_values_frame_L_buff=[];
%             vertex_values_frame_L_buff=vertex_values_frame_L(:,time_i);
%             load([file_suffix '_' task_name{task_i} '/' Montage_suffix '/L_vertex_values_' file_suffix '_' task_name{task_i} '_' num2str(time_i)]);
%             vertex_values_frame(find(vertex_values_frame==0))=nan;
%             
%             for atlas_set_index=1:size(atlas_index_all_new,2)
%                 
%                 atlas_index_sub_all=atlas_index_all_new{atlas_set_index};
%                 roi_index_L=[];
%                 for atlas_sub_index=atlas_index_sub_all
%                     roi_index_L_buff=find(L_label==L_colortable.table(atlas_sub_index,5));
%                     roi_index_L=[roi_index_L ; roi_index_L_buff];
%                 end
% 
%                 [size_L_parcel, index_L_parcel] = hist(L_400_label(roi_index_L),unique(L_400_label(roi_index_L)));
%                 
%                 for parcel_i=1:length(index_L_parcel)
%                     size_total=[]; size_parcel=[];
%                     size_total=size_L_400_label(find(index_L_400_label == index_L_parcel(parcel_i)));
%                     size_parcel=size_L_parcel(parcel_i);
%                     if size_parcel/size_total < 0.5
%                         idx=[];
%                         idx=ismember(roi_index_L,find(L_400_label == index_L_parcel(parcel_i)));
%                         roi_index_L(idx)=[];
%                     end
%                 end
%                 
%                 % total number of pacels in ROI
%                 plot_test_ROI_size(1,task_i,atlas_set_index,time_i)=length(unique(L_400_label(roi_index_L)));
%                 
%                 % number of significant parcels / total number
%                 all_values_L=[];
%                 all_values_L=unique(vertex_values_frame_L_buff(roi_index_L));
%                 plot_test_stats_ratio(1,task_i,atlas_set_index,time_i)=(length(find(all_values_L>0)) - length(find(all_values_L<0))) ...
%                     /length(unique(L_400_label(roi_index_L)));
%                 
%                 % average t values
%                 plot_test_stats_avg_t(1,task_i,atlas_set_index,time_i) = nanmean(vertex_values_frame_L_buff(roi_index_L));
%                 plot_test_stats_sum_t(1,task_i,atlas_set_index,time_i) = nansum(vertex_values_frame_L_buff(roi_index_L));
%                 
%                 % average sig t values
%                 if all_values_L == 0
%                     plot_test_stats_avg_t_sig(1,task_i,atlas_set_index,time_i) = 0;
%                 else
%                     plot_test_stats_avg_t_sig(1,task_i,atlas_set_index,time_i) = mean(all_values_L(all_values_L~=0));
%                     plot_test_stats_sum_t_sig(1,task_i,atlas_set_index,time_i) = nansum(all_values_L(all_values_L~=0));
%                 end
%                 
%                 % average z-score by all subjects
%                 plot_test(1,task_i,atlas_set_index,:,time_i)=nanmean(vertex_values_frame(roi_index_L,:),1);
%                 plot_test_weighted(1,task_i,atlas_set_index,:,time_i)=nansum(vertex_values_frame(roi_index_L,:),1) ...
%                     / sqrt(size(roi_index_L,1));
%                 plot_test_elec(1,task_i,atlas_set_index,:,time_i)=nanmean(vertex_values_frame_e(roi_index_L,:),1);
% 
%                 % average significant z-score by all subjects
%                 if plot_test_stats_ratio(1,task_i,atlas_set_index,time_i) ~= 0
%                     roi_index_L(vertex_values_frame_L_buff(roi_index_L)==0)=[];
%                     plot_test_sig(1,task_i,atlas_set_index,:,time_i)=nanmean(vertex_values_frame(roi_index_L,:),1);
%                     plot_test_elec_sig(1,task_i,atlas_set_index,:,time_i)=nanmean(vertex_values_frame_e(roi_index_L,:),1);
%                 end
%             end
%             
% 
%             %% Right
%             vertex_values_frame=[];
%             vertex_values_frame_e=[];
%             vertex_values_frame_R_buff=[];
%             vertex_values_frame_R_buff=vertex_values_frame_R(:,time_i);
%             load([file_suffix '_' task_name{task_i} '/' Montage_suffix '/R_vertex_values_' file_suffix '_' task_name{task_i} '_' num2str(time_i)]);
%             vertex_values_frame(find(vertex_values_frame==0))=nan;
%             
%             for atlas_set_index=1:size(atlas_index_all_new,2)
%                 
%                 atlas_index_sub_all=atlas_index_all_new{atlas_set_index};
%                 roi_index_R=[];
%                 for atlas_sub_index=atlas_index_sub_all
%                     roi_index_R_buff=find(R_label==R_colortable.table(atlas_sub_index,5));
%                     roi_index_R=[roi_index_R ; roi_index_R_buff];
%                 end
% 
%                 [size_R_parcel, index_R_parcel] = hist(R_400_label(roi_index_R),unique(R_400_label(roi_index_R)));
%                 
%                 for parcel_i=1:length(index_R_parcel)
%                     size_total=[]; size_parcel=[];
%                     size_total=size_R_400_label(find(index_R_400_label == index_R_parcel(parcel_i)));
%                     size_parcel=size_R_parcel(parcel_i);
%                     if size_parcel/size_total < 0.5
%                         idx=[];
%                         idx=ismember(roi_index_R,find(R_400_label == index_R_parcel(parcel_i)));
%                         roi_index_R(idx)=[];
%                     end
%                 end
%                 
%                 % total number of pacels in ROI
%                 plot_test_ROI_size(2,task_i,atlas_set_index,time_i)=length(unique(R_400_label(roi_index_R)));
%                 
%                 % number of significant parcels / total number
%                 all_values_R=[];
%                 all_values_R=unique(vertex_values_frame_R_buff(roi_index_R));
%                 plot_test_stats_ratio(2,task_i,atlas_set_index,time_i)=(length(find(all_values_R>0)) - length(find(all_values_R<0))) ...
%                     /length(unique(R_400_label(roi_index_R)));
%                 
%                 % average t values
%                 plot_test_stats_avg_t(2,task_i,atlas_set_index,time_i) = nanmean(vertex_values_frame_R_buff(roi_index_R));
%                 plot_test_stats_sum_t(2,task_i,atlas_set_index,time_i) = nansum(vertex_values_frame_R_buff(roi_index_R));
%                 
%                 % average sig t values
%                 if all_values_R == 0
%                     plot_test_stats_avg_t_sig(2,task_i,atlas_set_index,time_i) = 0;
%                 else
%                     plot_test_stats_avg_t_sig(2,task_i,atlas_set_index,time_i) = mean(all_values_R(all_values_R~=0));
%                     plot_test_stats_sum_t_sig(2,task_i,atlas_set_index,time_i) = nansum(all_values_R(all_values_R~=0));
%                 end
%                 
%                 % average z-score by all subjects
%                 plot_test(2,task_i,atlas_set_index,:,time_i)=nanmean(vertex_values_frame(roi_index_R,:),1);
%                 plot_test_weighted(2,task_i,atlas_set_index,:,time_i)=nansum(vertex_values_frame(roi_index_R,:),1) ...
%                     / sqrt(size(roi_index_R,1));
%                 plot_test_elec(2,task_i,atlas_set_index,:,time_i)=nanmean(vertex_values_frame_e(roi_index_R,:),1);
% 
%                 % average significant z-score by all subjects
%                 if plot_test_stats_ratio(2,task_i,atlas_set_index,time_i) ~= 0
%                     roi_index_R(vertex_values_frame_R_buff(roi_index_R)==0)=[];
%                     plot_test_sig(2,task_i,atlas_set_index,:,time_i)=nanmean(vertex_values_frame(roi_index_R,:),1);
%                     plot_test_elec_sig(2,task_i,atlas_set_index,:,time_i)=nanmean(vertex_values_frame_e(roi_index_R,:),1);
%                 end
%             end
%             
%         elseif flip_flag == 1
%            
%         end
%         time_i
%     end
% 
% end


cd(zscore_folder)

% save([file_suffix '_' Montage_suffix '_roi_val_' flip_prefix '_new.mat'],'plot_test','plot_test_sig' ...
%     ,'plot_test_stats_ratio','plot_test_stats_avg_t','plot_test_stats_sum_t','plot_test_stats_avg_t_sig' ...
%     ,'plot_test_stats_sum_t_sig','plot_test_elec','plot_test_elec_sig','plot_test_ROI_size','plot_test_weighted');
load([file_suffix '_' Montage_suffix '_roi_val_' flip_prefix '_new.mat']);

fig_order=[4 3 12 11 2 10 8 9 1 7 5 6 13];

%% plot all together
size_atlas_index_all=size(atlas_index_all_new,2);

% range_vals_max=nan(4,size_atlas_index_all);
% range_vals_min=nan(4,size_atlas_index_all);
% 
% load('xylim_vals_all.mat')
% range_vals_max=max(range_vals_max);
% range_vals_min=min(range_vals_min);

% range_vals_max=[2 2 12 12 2 2 2 2 2 2 4 4];
% range_vals_max=[1 1 10 10 1 1 1 1 1 1 3 3];
range_vals_max=[2 2 10 10 2 2 2 2 2 2 3 3 2];
% range_vals_max=[1 1 1 1 1 1 1 1 1 1 1 1];
range_vals_min=[-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1];
thr_p=0.05;
perm_p=0.05;
cluster_size_thr=6;

if flip_flag == 0
    task_name_fig={'Recalled','Not recalled','Recalled + Not recalled','Recalled - Not recalled'};
    
    for task_i=1:4
        figure('position',[0 0 4000 8000]);
%         figure('position',[0 0 1000 1000]);
        set(gcf, 'color', [1 1 1]);
        set(gcf,'Visible','off');
        set(gcf,'renderer','Painters')
        
        for atlas_index_i=1:size(atlas_index_all_new,2)

            subplot(4,4,atlas_index_i);
           
%             area_range=[20:32];
%             x_vector = [area_range, fliplr(area_range)];
%             patch = fill(x_vector, [ones(size(area_range,2),1)'*-range_vals(cnt), ...
%                 fliplr(ones(size(area_range,2),1)'*range_vals(cnt))],'g');
%             set(patch, 'edgecolor', 'none');
%             set(patch, 'FaceAlpha', 0.1);
%             hold on;
            
            plot_test_left=squeeze(plot_test(1,task_i,fig_order(atlas_index_i),:,:));
%             plot_test_left=squeeze(plot_test_weighted(1,task_i,atlas_index_i,:,:));
            plot_test_left(plot_test_left==0)=nan;
            [h1,output1]=plot_areaerrorbar(plot_test_left,1);
            hold on;

            plot_test_right=squeeze(plot_test(2,task_i,fig_order(atlas_index_i),:,:));
%             plot_test_right=squeeze(plot_test_weighted(2,task_i,atlas_index_i,:,:));
            plot_test_right(plot_test_right==0)=nan;
            [h2,output2]=plot_areaerrorbar(plot_test_right,2);
            hold on;
            
            %% simple ttest
%             [H_L,P_L,~,STATS_L]=ttest(plot_test_left-nanmean(plot_test_left(:,1:19),2),0,'Alpha',0.05,'tail','both');
%             [H_R,P_R,~,STATS_R]=ttest(plot_test_right-nanmean(plot_test_right(:,1:19),2),0,'Alpha',0.05,'tail','both');
%             [H_both,P_both,~,STATS_both]=ttest2(plot_test_left,plot_test_right,'tail','both');
            
            %% cluster based permutation analysis
            
            plot_test_left(isnan(plot_test_left(:,1)),:)=[]; % remove nan values
            plot_test_right(isnan(plot_test_right(:,1)),:)=[];
            
            %% left side
            plot_test_left_base=repmat(nanmean(plot_test_left(:,1:19),2),1,98);
            
            [clusters_left, p_values_left, ~ , ~ ] ...
                =permutest_ram(plot_test_left',plot_test_left_base',true,perm_p,10000,false);
            signicant_clusters=[];
            signicant_clusters=clusters_left(p_values_left<thr_p);

            if ~isempty(signicant_clusters)
                boundary_val_all=[];
                for c_i=1:size(signicant_clusters,2)
%                     if size(signicant_clusters{c_i},2)>cluster_size_thr
                        boundary_val_buf=[min(signicant_clusters{c_i}) max(signicant_clusters{c_i})];
                        boundary_val_all=[boundary_val_all boundary_val_buf];
%                     end
                end
                boundary_val_all=sort(boundary_val_all);
                
                b_flag=zeros(size(boundary_val_all,2),1);
                for b_i=2:2:size(boundary_val_all,2)-1
                    if (boundary_val_all(b_i+1)-boundary_val_all(b_i)<cluster_size_thr)
                        b_flag(b_i)=1;
                        b_flag(b_i+1)=1;
                    end
                end
                
                boundary_val_all(find(b_flag==1))=[];
                for b_i=1:2:size(boundary_val_all,2)
                    signicant_clusters_all=boundary_val_all(b_i):boundary_val_all(b_i+1);
                    if size(signicant_clusters_all,2) > cluster_size_thr
                        plot(signicant_clusters_all,ones(size(signicant_clusters_all,2),1) ...
                                *range_vals_max(fig_order(atlas_index_i))*0.8,'-b','LineWidth',4)
                    end
                end
            end
            
%             if(size(signicant_clusters,2)) > 1
%                 signicant_clusters=signicant_clusters(1);
%             end
%             if ~isempty(signicant_clusters)
%                 for c_i=1:size(signicant_clusters,2)
% %                     plot(signicant_clusters{c_i},ones(size(signicant_clusters{c_i})) ...
% %                         *range_vals_max(fig_order(atlas_index_i))*0.8,'*','Color',[52 148 255]./255)
%                     if size(signicant_clusters{c_i},2)>cluster_size_thr
%                         plot(signicant_clusters{c_i},ones(size(signicant_clusters{c_i})) ...
%                             *range_vals_max(fig_order(atlas_index_i))*0.8,'-b','LineWidth',4)
%                     end
%                 end
%             end

            [clusters_left, p_values_left, ~ , ~ ] ...
                =permutest_ram(plot_test_left_base',plot_test_left',true,perm_p,10000,false);
            signicant_clusters=[];
            signicant_clusters=clusters_left(p_values_left<thr_p);

            if ~isempty(signicant_clusters)
                boundary_val_all=[];
                for c_i=1:size(signicant_clusters,2)
%                     if size(signicant_clusters{c_i},2)>cluster_size_thr
                        boundary_val_buf=[min(signicant_clusters{c_i}) max(signicant_clusters{c_i})];
                        boundary_val_all=[boundary_val_all boundary_val_buf];
%                     end
                end
                boundary_val_all=sort(boundary_val_all);
                
                b_flag=zeros(size(boundary_val_all,2),1);
                for b_i=2:2:size(boundary_val_all,2)-1
                    if (boundary_val_all(b_i+1)-boundary_val_all(b_i)<cluster_size_thr)
                        b_flag(b_i)=1;
                        b_flag(b_i+1)=1;
                    end
                end
                
                boundary_val_all(find(b_flag==1))=[];
                for b_i=1:2:size(boundary_val_all,2)
                    signicant_clusters_all=boundary_val_all(b_i):boundary_val_all(b_i+1);
                    if size(signicant_clusters_all,2) > cluster_size_thr
                        plot(signicant_clusters_all,ones(size(signicant_clusters_all,2),1) ...
                                *range_vals_max(fig_order(atlas_index_i))*0.8,'-b','LineWidth',4)
                    end
                end
            end
            
%             if(size(signicant_clusters,2)) > 1
%                 signicant_clusters=signicant_clusters(1);
%             end            
%             if ~isempty(signicant_clusters)
%                 for c_i=1:size(signicant_clusters,2)
% %                     plot(signicant_clusters{c_i},ones(size(signicant_clusters{c_i})) ...
% %                         *range_vals_max(fig_order(atlas_index_i))*0.8,'*','Color',[52 148 255]./255)
%                     if size(signicant_clusters{c_i},2)>cluster_size_thr
%                         plot(signicant_clusters{c_i},ones(size(signicant_clusters{c_i})) ...
%                             *range_vals_max(fig_order(atlas_index_i))*0.8,'-b','LineWidth',4)
%                     end
%                 end
%             end
            
            %% right side
            plot_test_right_base=repmat(nanmean(plot_test_right(:,1:19),2),1,98);
            
            [clusters_right, p_values_right, ~ , ~ ] ...
                =permutest_ram(plot_test_right_base',plot_test_right',true,perm_p,10000,false);
            signicant_clusters=[];
            signicant_clusters=clusters_right(p_values_right<thr_p);
            
            if ~isempty(signicant_clusters)
                boundary_val_all=[];
                for c_i=1:size(signicant_clusters,2)
%                     if size(signicant_clusters{c_i},2)>cluster_size_thr
                        boundary_val_buf=[min(signicant_clusters{c_i}) max(signicant_clusters{c_i})];
                        boundary_val_all=[boundary_val_all boundary_val_buf];
%                     end
                end
                boundary_val_all=sort(boundary_val_all);
                
                b_flag=zeros(size(boundary_val_all,2),1);
                for b_i=2:2:size(boundary_val_all,2)-1
                    if (boundary_val_all(b_i+1)-boundary_val_all(b_i)<cluster_size_thr)
                        b_flag(b_i)=1;
                        b_flag(b_i+1)=1;
                    end
                end
                
                boundary_val_all(find(b_flag==1))=[];
                for b_i=1:2:size(boundary_val_all,2)
                    signicant_clusters_all=boundary_val_all(b_i):boundary_val_all(b_i+1);
                    if size(signicant_clusters_all,2) > cluster_size_thr
                        plot(signicant_clusters_all,ones(size(signicant_clusters_all,2),1) ...
                                *range_vals_max(fig_order(atlas_index_i))*0.7,'-r','LineWidth',4)
                    end
                end
            end
            
% %             if(size(signicant_clusters,2)) > 1
% %                 signicant_clusters=signicant_clusters(1);
% %             end
%             if ~isempty(signicant_clusters)
%                 for c_i=1:size(signicant_clusters,2)
% %                     plot(signicant_clusters{c_i},ones(size(signicant_clusters{c_i})) ...
% %                         *range_vals_max(fig_order(atlas_index_i))*0.7,'*','Color',[255 60 22]./255)
%                     if size(signicant_clusters{c_i},2)>cluster_size_thr                    
%                         plot(signicant_clusters{c_i},ones(size(signicant_clusters{c_i})) ...
%                             *range_vals_max(fig_order(atlas_index_i))*0.7,'-r','LineWidth',4)
%                     end
%                 end
%             end
            
            [clusters_right, p_values_right, ~ , ~ ] ...
                =permutest_ram(plot_test_right',plot_test_right_base',true,perm_p,10000,false);
            signicant_clusters=[];
            signicant_clusters=clusters_right(p_values_right<thr_p);
            
            if ~isempty(signicant_clusters)
                boundary_val_all=[];
                for c_i=1:size(signicant_clusters,2)
%                     if size(signicant_clusters{c_i},2)>cluster_size_thr
                        boundary_val_buf=[min(signicant_clusters{c_i}) max(signicant_clusters{c_i})];
                        boundary_val_all=[boundary_val_all boundary_val_buf];
%                     end
                end
                boundary_val_all=sort(boundary_val_all);
                
                b_flag=zeros(size(boundary_val_all,2),1);
                for b_i=2:2:size(boundary_val_all,2)-1
                    if (boundary_val_all(b_i+1)-boundary_val_all(b_i)<cluster_size_thr)
                        b_flag(b_i)=1;
                        b_flag(b_i+1)=1;
                    end
                end
                
                boundary_val_all(find(b_flag==1))=[];
                for b_i=1:2:size(boundary_val_all,2)
                    signicant_clusters_all=boundary_val_all(b_i):boundary_val_all(b_i+1);
                    if size(signicant_clusters_all,2) > cluster_size_thr
                        plot(signicant_clusters_all,ones(size(signicant_clusters_all,2),1) ...
                                *range_vals_max(fig_order(atlas_index_i))*0.7,'-r','LineWidth',4)
                    end
                end
            end            
            
% %             if(size(signicant_clusters,2)) > 1
% %                 signicant_clusters=signicant_clusters(1);
% %             end
%             if ~isempty(signicant_clusters)
%                 for c_i=1:size(signicant_clusters,2)
% %                     plot(signicant_clusters{c_i},ones(size(signicant_clusters{c_i})) ...
% %                         *range_vals_max(fig_order(atlas_index_i))*0.7,'*','Color',[255 60 22]./255)
%                     if size(signicant_clusters{c_i},2)>cluster_size_thr
%                         plot(signicant_clusters{c_i},ones(size(signicant_clusters{c_i})) ...
%                             *range_vals_max(fig_order(atlas_index_i))*0.7,'-r','LineWidth',4)
%                     end
%                 end
%             end
            
            %% left vs right
            [clusters_both, p_values_both, ~, ~ ] ...
                =permutest_ram(plot_test_right',plot_test_left',false,perm_p,10000,false);
            
            signicant_clusters=[];
            signicant_clusters=clusters_both(p_values_both<thr_p);
            
            if ~isempty(signicant_clusters)
                boundary_val_all=[];
                for c_i=1:size(signicant_clusters,2)
%                     if size(signicant_clusters{c_i},2)>cluster_size_thr
                        boundary_val_buf=[min(signicant_clusters{c_i}) max(signicant_clusters{c_i})];
                        boundary_val_all=[boundary_val_all boundary_val_buf];
%                     end
                end
                boundary_val_all=sort(boundary_val_all);
                
                b_flag=zeros(size(boundary_val_all,2),1);
                for b_i=2:2:size(boundary_val_all,2)-1
                    if (boundary_val_all(b_i+1)-boundary_val_all(b_i)<cluster_size_thr)
                        b_flag(b_i)=1;
                        b_flag(b_i+1)=1;
                    end
                end
                
                boundary_val_all(find(b_flag==1))=[];
                for b_i=1:2:size(boundary_val_all,2)
                    signicant_clusters_all=boundary_val_all(b_i):boundary_val_all(b_i+1);
                    if size(signicant_clusters_all,2) > cluster_size_thr
                        plot(signicant_clusters_all,ones(size(signicant_clusters_all,2),1) ...
                                *range_vals_max(fig_order(atlas_index_i))*0.9,'-g','LineWidth',4)
                    end
                end
            end    
            
% %             if(size(signicant_clusters,2)) > 1
% %                 signicant_clusters=signicant_clusters(1);
% %             end
%             if ~isempty(signicant_clusters)
%                 for c_i=1:size(signicant_clusters,2)
% %                     plot(signicant_clusters{c_i},ones(size(signicant_clusters{c_i})) ...
% %                         *range_vals_max(fig_order(atlas_index_i))*0.9,'g*')
%                     if size(signicant_clusters{c_i},2)>cluster_size_thr
%                         plot(signicant_clusters{c_i},ones(size(signicant_clusters{c_i})) ...
%                             *range_vals_max(fig_order(atlas_index_i))*0.9,'-g','LineWidth',4)
%                     end
%                 end
%             end
            
            [clusters_both, p_values_both, ~, ~ ] ...
                =permutest_ram(plot_test_left',plot_test_right',false,perm_p,10000,false);
            
            signicant_clusters=[];
            signicant_clusters=clusters_both(p_values_both<thr_p);
            
            if ~isempty(signicant_clusters)
                boundary_val_all=[];
                for c_i=1:size(signicant_clusters,2)
%                     if size(signicant_clusters{c_i},2)>cluster_size_thr
                        boundary_val_buf=[min(signicant_clusters{c_i}) max(signicant_clusters{c_i})];
                        boundary_val_all=[boundary_val_all boundary_val_buf];
%                     end
                end
                boundary_val_all=sort(boundary_val_all);
                
                b_flag=zeros(size(boundary_val_all,2),1);
                for b_i=2:2:size(boundary_val_all,2)-1
                    if (boundary_val_all(b_i+1)-boundary_val_all(b_i)<cluster_size_thr)
                        b_flag(b_i)=1;
                        b_flag(b_i+1)=1;
                    end
                end
                
                boundary_val_all(find(b_flag==1))=[];
                for b_i=1:2:size(boundary_val_all,2)
                    signicant_clusters_all=boundary_val_all(b_i):boundary_val_all(b_i+1);
                    if size(signicant_clusters_all,2) > cluster_size_thr
                        plot(signicant_clusters_all,ones(size(signicant_clusters_all,2),1) ...
                                *range_vals_max(fig_order(atlas_index_i))*0.9,'-g','LineWidth',4)
                    end
                end
            end                
            
% %             if(size(signicant_clusters,2)) > 1
% %                 signicant_clusters=signicant_clusters(1);
% %             end
%             if ~isempty(signicant_clusters)
%                 for c_i=1:size(signicant_clusters,2)
% %                     plot(signicant_clusters{c_i},ones(size(signicant_clusters{c_i})) ...
% %                         *range_vals_max(fig_order(atlas_index_i))*0.9,'g*')
%                     if size(signicant_clusters{c_i},2)>cluster_size_thr                   
%                         plot(signicant_clusters{c_i},ones(size(signicant_clusters{c_i})) ...
%                             *range_vals_max(fig_order(atlas_index_i))*0.9,'-g','LineWidth',4)
%                     end
%                 end
%             end

            %%
%             set(gca,'xtick',[1 20 40 60 80 98]);
%             set(gca,'xticklabel',[-500 0 500 1000 1500 2000]);
            if atlas_index_i==13            
                set(gca,'xtick',[8 20 32 44 56 68 80 92]);
                set(gca,'xticklabel',[-300 0 300 600 900 1200 1500 1800]);
                set(gca,'xlim',[1 98])
                xlabel('Time (ms)');
            else
                set(gca,'xtick',[]);
                set(gca,'xticklabel',[]);
            end

            set(gca,'box','off')

%             range_vals_min(task_i,atlas_index_i)=min([output1.min_val ; output1.max_val]);
%             range_vals_max(task_i,atlas_index_i)=max([output1.min_val ; output1.max_val]);
            set(gca,'ylim',[range_vals_min(fig_order(atlas_index_i)) range_vals_max(fig_order(atlas_index_i))])
            set(gca,'Fontsize',15);

%             x1=xline(20,'-.r',{'ON'},'LabelVerticalAlignment','top','LabelHorizontalAlignment','center','LabelOrientation','horizontal','LineWidth',3,'Fontsize',15);
%             x2=xline(84,'-.r',{'OFF'},'LabelVerticalAlignment','top','LabelHorizontalAlignment','center','LabelOrientation','horizontal','LineWidth',3,'Fontsize',15);
            x1=xline(20,':k','LineWidth',2,'Fontsize',15);
            x2=xline(84,':k','LineWidth',2,'Fontsize',15);
%             x3=xline(32,'--r','LineWidth',3,'Fontsize',15);
            y1=yline(0,'k','LineWidth',1,'Fontsize',15);

            %         line([20,20],[-8,10],'Color','r','LineWidth',1)
            %         line([84,84],[-8,10],'Color','r','LineWidth',1)


            if atlas_index_i==1
%                 ylabel({'Gamma power z-score'});
                legend([h1 h2],{'Left','Right'},'LineWidth',1.5,'FontSize',10)
            end

            title(atlas_index_all_new_name{fig_order(atlas_index_i)},'FontSize',20,'Color',atlas_color(fig_order(atlas_index_i),:))
        
        end
        
        [ax,h3]=suplabel(task_name_fig{task_i},'t');
        set(h3,'FontSize',30)
%         [ax2,h4]=suplabel({'Summation of significant t values'},'y');
        [ax2,h4]=suplabel({'Averaged gamma power z-score'},'y');
        set(h4,'FontSize',30)
        %% save figure
%         saveas(gcf,[zscore_folder '\' task_name_fig{task_i} '_all_plot.png'])

        export_fig([zscore_folder '\' task_name_fig{task_i} '_all_plot'],'-tiff')
        saveas(gcf,[zscore_folder '\' task_name_fig{task_i} '_all_plot'],'epsc')
        close all;
        
        
%         figure('position',[0 0 4000 8000]);
%         set(gcf, 'color', [1 1 1]);
%         set(gcf,'Visible','on');
%         set(gcf,'renderer','Painters')
%         
%         cnt=1;
%         for atlas_index_i=7:size(atlas_index_all_new,2)
% 
%             subplot(3,2,cnt);
%            
% %             area_range=[20:32];
% %             x_vector = [area_range, fliplr(area_range)];
% %             patch = fill(x_vector, [ones(size(area_range,2),1)'*-range_vals(cnt), ...
% %                 fliplr(ones(size(area_range,2),1)'*range_vals(cnt))],'g');
% %             set(patch, 'edgecolor', 'none');
% %             set(patch, 'FaceAlpha', 0.1);
% %             hold on;
%             
%             plot_test_left=squeeze(plot_test(1,task_i,fig_order(atlas_index_i),:,:));
% %             plot_test_left=squeeze(plot_test_weighted(1,task_i,atlas_index_i,:,:));
%             plot_test_left(plot_test_left==0)=nan;
%             [h1,output1]=plot_areaerrorbar(plot_test_left,1);
% 
%             hold on;
% 
%             plot_test_right=squeeze(plot_test(2,task_i,fig_order(atlas_index_i),:,:));
% %             plot_test_right=squeeze(plot_test_weighted(2,task_i,atlas_index_i,:,:));
%             plot_test_right(plot_test_right==0)=nan;
%             [h2,output2]=plot_areaerrorbar(plot_test_right,2);
% 
% %             set(gca,'xtick',[1 20 40 60 80 98]);
% %             set(gca,'xticklabel',[-500 0 500 1000 1500 2000]);
%             set(gca,'xtick',[8 20 32 44 56 68 80 92]);
%             set(gca,'xticklabel',[-300 0 300 600 900 1200 1500 1800]);
%             set(gca,'xlim',[1 98])
%             xlabel('Time (ms)');
% 
%             set(gca,'box','off')
% 
% %             range_vals_min(task_i,atlas_index_i)=min([output1.min_val ; output1.max_val]);
% %             range_vals_max(task_i,atlas_index_i)=max([output1.min_val ; output1.max_val]);
%             set(gca,'ylim',[range_vals_min(fig_order(atlas_index_i)) range_vals_max(fig_order(atlas_index_i))])
%             set(gca,'Fontsize',15);
% 
% %             x1=xline(20,'-.r',{'ON'},'LabelVerticalAlignment','top','LabelHorizontalAlignment','center','LabelOrientation','horizontal','LineWidth',3,'Fontsize',15);
% %             x2=xline(84,'-.r',{'OFF'},'LabelVerticalAlignment','top','LabelHorizontalAlignment','center','LabelOrientation','horizontal','LineWidth',3,'Fontsize',15);
%             x1=xline(20,':k','LineWidth',2,'Fontsize',15);
%             x2=xline(84,':k','LineWidth',2,'Fontsize',15);
% %             x3=xline(32,'--r','LineWidth',3,'Fontsize',15);
%             y1=yline(0,'k','LineWidth',1,'Fontsize',15);
% 
%             %         line([20,20],[-8,10],'Color','r','LineWidth',1)
%             %         line([84,84],[-8,10],'Color','r','LineWidth',1)
% 
% 
%             if cnt==1
% %                 ylabel({'Gamma power z-score'});
%                 legend([h1 h2],{'Left','Right'},'LineWidth',1.5,'FontSize',10)
%             end
% 
%             title(atlas_index_all_new_name{fig_order(atlas_index_i)},'FontSize',20,'Color',atlas_color(fig_order(atlas_index_i),:))
%             
%             cnt=cnt+1;
%         end
%         
%         [ax,h3]=suplabel(task_name_fig{task_i},'t');
%         set(h3,'FontSize',30)
% %         [ax2,h4]=suplabel({'Summation of significant t values'},'y');
%         [ax2,h4]=suplabel({'Averaged gamma power z-score'},'y');
%         set(h4,'FontSize',30)
%         %% save figure
% %         saveas(gcf,[zscore_folder '\' task_name_fig{task_i} '_all_plot.png'])
% 
%         export_fig([zscore_folder '\' task_name_fig{task_i} '_all_plot_2'],'-tiff')
%         saveas(gcf,[zscore_folder '\' task_name_fig{task_i} '_all_plot_1'],'epsc')
%         close all;

    end
elseif flip_flag == 1
    
    
    
    
end
% save('xylim_vals_all.mat','range_vals_min','range_vals_max')
