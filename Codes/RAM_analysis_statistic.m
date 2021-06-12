clc;clear;

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)
cd('E:\RAM data set\RAM_Public_Data_all\FR1_final\Storage_filtering_no_overlap\Data_storage\final_clean')


%% Ranksum test
% Left
load('L_vertex_values_sum_5_subjects.mat')
vertex_size=size(L_vertex_values_sum,1);
L_vertex_values_p=nan(vertex_size,1);
for vertex_i=1:vertex_size
%     [p,h,stats]=ranksum(abs(L_vertex_values_sum(vertex_i,1:9)),abs(L_vertex_values_sum(vertex_i,10:39)),'tail','both');
    [p,h,stats]=signrank(L_vertex_values_sum(vertex_i,10:39));
    L_vertex_values_p(vertex_i)=p;
    vertex_i
end

% Right
load('R_vertex_values_sum_5_subjects.mat')
vertex_size=size(R_vertex_values_sum,1);
R_vertex_values_p=nan(vertex_size,1);
for vertex_i=1:vertex_size
%     [p,h,stats]=ranksum(abs(R_vertex_values_sum(vertex_i,1:9)),abs(R_vertex_values_sum(vertex_i,10:39)),'tail','both');
    [p,h,stats]=signrank(R_vertex_values_sum(vertex_i,10:39));
    R_vertex_values_p(vertex_i)=p;
    vertex_i
end

%% FDR corrected
all_vertex_values_p=[L_vertex_values_p; R_vertex_values_p];
fdr_p=fdr(all_vertex_values_p,0.05);

L_vertex_values_p_fdr_mask=L_vertex_values_p<fdr_p;
R_vertex_values_p_fdr_mask=R_vertex_values_p<fdr_p;

vertex_values_p=[]; vertex_values_p_fdr_mask=[];
vertex_values_p=L_vertex_values_p; vertex_values_p_fdr_mask=L_vertex_values_p_fdr_mask;
save('L_vertex_values_sum_5_subjects_stats.mat', 'vertex_values_p', 'vertex_values_p_fdr_mask');

vertex_values_p=[]; vertex_values_p_fdr_mask=[];
vertex_values_p=R_vertex_values_p; vertex_values_p_fdr_mask=R_vertex_values_p_fdr_mask;
save('R_vertex_values_sum_5_subjects_stats.mat', 'vertex_values_p', 'vertex_values_p_fdr_mask');


% roi n = 235

%% Atlas based stats
clc;clear;close all;
load('E:\RAM data set\RAM_Public_Data_all\FR1_final\Atlas\Shen_atlas.mat')


rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)
cd('E:\RAM data set\RAM_Public_Data_all\FR1_final\Storage_filtering_no_overlap\Data_storage\final_clean')

vertex_values_p=nan(235,1);
vertex_values_t=nan(235,1);

% Left
load('L_vertex_values_sum_5_subjects.mat')
vertex_size=size(L_vertex_values_sum,1);
for ROI_i=1:235
    try
        base_mean=[];
        base_mean=nanmean(L_vertex_values_sum(L_vertex_values_frame==ROI_i,1:9),1); % baseline
        post_mean=[];
        post_mean=nanmean(L_vertex_values_sum(L_vertex_values_frame==ROI_i,10:end),1); % post        
        vertex_values_p(ROI_i,1)=signrank(post_mean);
    catch
    end
end

% Right
load('R_vertex_values_sum_5_subjects.mat')
vertex_size=size(R_vertex_values_sum,1);
for ROI_i=1:235
    try
        base_mean=[];
        base_mean=nanmean(R_vertex_values_sum(R_vertex_values_frame==ROI_i,1:9),1); % baseline
        post_mean=[];
        post_mean=nanmean(R_vertex_values_sum(R_vertex_values_frame==ROI_i,10:end),1); % post        
        vertex_values_p(ROI_i,1)=signrank(post_mean);
    catch
    end
end

%% FDR corrected (Atlas based)
fdr_p=fdr(vertex_values_p,0.05);
p_fdr_mask=vertex_values_p<fdr_p;
sig_rois=find(p_fdr_mask>0);

%left
vertex_values_p_fdr_mask=zeros(size(L_vertex_values_frame,1),1);
for i=1:size(sig_rois,1)
    vertex_values_p_fdr_mask(find(sig_rois(i)==L_vertex_values_frame))=1;
end
save('L_vertex_values_sum_5_subjects_stats_roi.mat', 'vertex_values_p_fdr_mask');

%right
vertex_values_p_fdr_mask=zeros(size(R_vertex_values_frame,1),1);
for i=1:size(sig_rois,1)
    vertex_values_p_fdr_mask(find(sig_rois(i)==R_vertex_values_frame))=1;
end
save('R_vertex_values_sum_5_subjects_stats_roi.mat', 'vertex_values_p_fdr_mask');







