
%% zscore plot by ROIs
clc;clear;

addpath(genpath('E:\RAM data set\RAM_Public_Data_all\Parcellations'));
addpath(genpath('E:\RAM data set\RAM_Public_Data_all\fsaverage'));
addpath(genpath('E:\RAM data set\RAM_Public_Data_all\Codes'));

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

atlas_index_all_new={4,30,8,12,13,[15 27],28,29,[19 20 21],[9,32],[10 16],17};

% Figure 4
% atlas_index_all_new={12,8,17,4,13,[15 27]};
% fig_order=[4 3 12 1 5 6];
% roi_result_folder='Figure4';
% roi_result_folder_map='Figure4_ROIs';

% % Figure 5
% atlas_index_all_new={30,[9,32],29,[19 20 21],28,[10 16]};
% fig_order=[2 10 8 9 7 11];
% roi_result_folder='Figure5';
% roi_result_folder_map='Figure5_ROIs';

atlas_index_all_new_name={'Caudal Middle Frontal','Superior Parietal'...
    ,'Fusiform Gyrus','Lateral Occipital','Lateral Orbitofrontal','Medial Orbitofrontal' ...
    ,'Rostral Middle Frontal','Superior Frontal' ...
    ,'Inferior Frontal','Inferior Parietal','Inferior Temporal' ...
    ,'Parahippocampal'};

zscore_folder = [rootfolder '\FR1_FS\Storage_filtering_overlap_base_agg\zscore_plot\' Montage_suffix ];
mkdir(zscore_folder)
mkdir([zscore_folder '\' roi_result_folder])

% atlas_color=distinguishable_colors(12);
% atlas_color=[237 28 36;247 148 29;39 170 225;33 64 154;127 63 152;0 148 68]/255; %OLD
% atlas_color=[237 28 36;247 148 29;39 170 225;254 222 0;127 63 152;0 148 68]/255; %NEW (yellow)
atlas_color=[237 28 36;247 148 29;39 170 225;41 42 116;127 63 152;0 148 68]/255; %NEW (blue)

overlay=4;figure_type=0;views=4;
displayElectrodesInflated_new_fs_Atlas_all_final(roi_result_folder_map,zscore_folder, overlay, figure_type, views, atlas_index_all_new,atlas_color) 
close all


% % atlas_color=distinguishable_colors(200);
% overlay=4;laterality=4;views=4;
% % displayElectrodesInflated_new_fs_Atlas_all_final(zscore_folder, overlay, laterality, views, atlas_index_all_new,atlas_color) 
% displayElectrodesInflated_new_fs_Atlas_all(zscore_folder, overlay, laterality, views,0,0, 1, 21)
% 
% close all

if flip_flag ==0 
    task_name={'rec','Nonrec','together','diff'};
    plot_test=nan(2,4,size(atlas_index_all_new,2),251,98);
   
elseif flip_flag ==1
    task_name={'rec_xflip','Nonrec_xflip','together_xflip','diff_xflip'};
    plot_test=nan(1,4,size(atlas_index_all_new,2),251,98);
end


%% data load/save
for task_i=1:4
    for time_i=1:98
        if flip_flag == 0
            %% Left
            vertex_values_frame=[];
            load([file_suffix '_' task_name{task_i} '/' Montage_suffix '/L_vertex_values_' file_suffix '_' task_name{task_i} '_' num2str(time_i)]);
            vertex_values_frame(find(vertex_values_frame==0))=nan;
            
            for atlas_set_index=1:size(atlas_index_all_new,2)
                
                atlas_index_sub_all=atlas_index_all_new{atlas_set_index};
                roi_index_L=[];
                for atlas_sub_index=atlas_index_sub_all
                    roi_index_L_buff=find(L_label==L_colortable.table(atlas_sub_index,5));
                    roi_index_L=[roi_index_L ; roi_index_L_buff];
                end
                plot_test(1,task_i,atlas_set_index,:,time_i)=nanmean(vertex_values_frame(roi_index_L,:),1);
            end
            

            %% Right
            vertex_values_frame=[];
            load([file_suffix '_' task_name{task_i} '/' Montage_suffix '/R_vertex_values_' file_suffix '_' task_name{task_i} '_' num2str(time_i)]);
            vertex_values_frame(find(vertex_values_frame==0))=nan;
            
            for atlas_set_index=1:size(atlas_index_all_new,2)
                
                atlas_index_sub_all=atlas_index_all_new{atlas_set_index};
                roi_index_R=[];
                for atlas_sub_index=atlas_index_sub_all
                    roi_index_R_buff=find(R_label==R_colortable.table(atlas_sub_index,5));
                    roi_index_R=[roi_index_R ; roi_index_R_buff];
                end
                plot_test(2,task_i,atlas_set_index,:,time_i)=nanmean(vertex_values_frame(roi_index_R,:),1); 
            end
            
        elseif flip_flag == 1
           
        end
        time_i
    end
end


cd(zscore_folder)

save([file_suffix '_' Montage_suffix '_roi_val_' flip_prefix '_new.mat'],'plot_test');
% load([file_suffix '_' Montage_suffix '_roi_val_' flip_prefix '_new.mat']);

%% plot all together
size_atlas_index_all=size(atlas_index_all_new,2);

% range_vals_max=nan(4,size_atlas_index_all);
% range_vals_min=nan(4,size_atlas_index_all);
% 
% load('xylim_vals_all.mat')
% range_vals_max=max(range_vals_max);
% range_vals_min=min(range_vals_min);

% range_vals_max=[2 2 12 12 2 2 2 2 2 2 4];
% range_vals_max=[1 1 10 10 1 1 1 1 1 1 2 2];
range_vals_max=[2 2 10 10 2 2 2 2 2 2 3 3];
% range_vals_max=[1 1 1 1 1 1 1 1 1 1 1];
range_vals_min=[-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1];
thr_p_1=0.05;
% thr_p_2=0.1;
perm_p=0.05;
cluster_size_thr=6; % 150ms
Total_patients_number=nan(4,12,2);
Total_patients_number_re=nan(4,12,2);

cd(zscore_folder)
load('cluster_stats_all.mat')
if flip_flag == 0
    task_name_fig={'Recalled','Not recalled','Recalled + Not recalled','Recalled - Not recalled'};
    
    for task_i=1:4
        figure('position',[0 0 1000 4000]);
        set(gcf,'color', [1 1 1]);
        set(gcf,'Visible','on');
        set(gcf,'renderer','Painters')
        
        for atlas_index_i=1:6%size(atlas_index_all_new,2)

            subplot(4,2,atlas_index_i);
%             subaxis(4,2,atlas_index_i,'SpacingVert',0.05,'MR',0)

            %% Left
            plot_test_left=squeeze(plot_test(1,task_i,fig_order(atlas_index_i),:,:));
%             plot_test_left=squeeze(plot_test_weighted(1,task_i,atlas_index_i,:,:));
            plot_test_left(plot_test_left==0)=nan;
            plot_test_left(isnan(plot_test_left(:,1)),:)=[]; % remove nan values
            
            Total_patients_number(task_i,atlas_index_i,1)=size(plot_test_left,1);
            
%             % additional step
%             meantrace = mean(plot_test_left)';
%             Sn=size(plot_test_left,1);
%             % calculate mean square error for each trial's trace
%             meansqerr = NaN(1,Sn);
%             samples_num=98;
%             for k = 1:Sn
%                 % first find the error from the mean of each of the 2048
%                 % samples in that trace
%                 sample_errors_squared = NaN(samples_num,1);
%                 for h = 1:samples_num
%                     sample_errors_squared(h,1) = (plot_test_left(k,h)-0)^2;
%                 end
%                 % the mean square error is the sum of all the errors for each
%                 % sample of that dtrace divided by the number of samples
%                 meansqerr(1,k) = nansum(sample_errors_squared)/samples_num;
%             end
% %             plot_test_left(meansqerr<1,:)=[];
%             plot_test_left(zscore(meansqerr)>1.96 | zscore(meansqerr)<-1.96,:)=[];
%             Total_patients_number_re(task_i,atlas_index_i,1)=size(plot_test_left,1);
%              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            [h1,output1]=plot_areaerrorbar_new(plot_test_left,1);
            hold on;

            %% Right
            plot_test_right=squeeze(plot_test(2,task_i,fig_order(atlas_index_i),:,:));
%             plot_test_right=squeeze(plot_test_weighted(2,task_i,atlas_index_i,:,:));
            plot_test_right(plot_test_right==0)=nan;
            plot_test_right(isnan(plot_test_right(:,1)),:)=[];
            
            Total_patients_number(task_i,atlas_index_i,2)=size(plot_test_right,1);
            
%             % additional step
%             meantrace = mean(plot_test_right)';
%             Sn=size(plot_test_right,1);
%             % calculate mean square error for each trial's trace
%             meansqerr = NaN(1,Sn);
%             samples_num=98;
%             for k = 1:Sn
%                 % first find the error from the mean of each of the 2048
%                 % samples in that trace
%                 sample_errors_squared = NaN(samples_num,1);
%                 for h = 1:samples_num
%                     sample_errors_squared(h,1) = (plot_test_right(k,h)-0)^2;
%                 end
%                 % the mean square error is the sum of all the errors for each
%                 % sample of that dtrace divided by the number of samples
%                 meansqerr(1,k) = nansum(sample_errors_squared)/samples_num;
%             end
%             plot_test_right(zscore(meansqerr)>1.96 | zscore(meansqerr)<-1.96,:)=[];
%             Total_patients_number_re(task_i,atlas_index_i,2)=size(plot_test_right,1);
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            [h2,output2]=plot_areaerrorbar_new(plot_test_right,2);
            hold on;

            %% simple ttest
%             [H_L,P_L,~,STATS_L]=ttest(plot_test_left-nanmean(plot_test_left(:,1:19),2),0,'Alpha',0.05,'tail','both');
%             [H_R,P_R,~,STATS_R]=ttest(plot_test_right-nanmean(plot_test_right(:,1:19),2),0,'Alpha',0.05,'tail','both');
%             [H_both,P_both,~,STATS_both]=ttest2(plot_test_left,plot_test_right,'tail','both');
            
            %% cluster based permutation analysis
            
%             plot_test_left(isnan(plot_test_left(:,1)),:)=[]; % remove nan values
%             plot_test_right(isnan(plot_test_right(:,1)),:)=[];
            
            %% left side
            plot_test_left_base=repmat(nanmean(plot_test_left(:,1:19),2),1,98);
            
%             [clusters_left, p_values_left, t_sums_left, permutation_distribution_left] ...
%                 =permutest_ram(plot_test_left',plot_test_left_base',true,perm_p,10000,true);
%             cluster_stats_left{task_i,fig_order(atlas_index_i)}.clusters=clusters_left;
%             cluster_stats_left{task_i,fig_order(atlas_index_i)}.p_values=p_values_left;
%             cluster_stats_left{task_i,fig_order(atlas_index_i)}.t_sums=t_sums_left;
%             cluster_stats_left{task_i,fig_order(atlas_index_i)}.permutation_distribution=permutation_distribution_left;
            
            clusters_left=cluster_stats_left{task_i,fig_order(atlas_index_i)}.clusters;
            p_values_left=cluster_stats_left{task_i,fig_order(atlas_index_i)}.p_values;

            signicant_clusters=[];
            signicant_clusters=clusters_left(p_values_left<thr_p_1);
%             signicant_clusters{:}
%             p_values_left
            
            if ~isempty(signicant_clusters)
                boundary_val_all=[];
                for c_i=1:size(signicant_clusters,2)
%                     if size(signicant_clusters{c_i},2)>cluster_size_thr
                        boundary_val_buf=[min(signicant_clusters{c_i}) max(signicant_clusters{c_i})];
                        boundary_val_all=[boundary_val_all boundary_val_buf];
%                     end
                end
                boundary_val_all=sort(boundary_val_all);
                
%                 b_flag=zeros(size(boundary_val_all,2),1);
%                 for b_i=2:2:size(boundary_val_all,2)-1
%                     if (boundary_val_all(b_i+1)-boundary_val_all(b_i)<cluster_size_thr)
%                         b_flag(b_i)=1;
%                         b_flag(b_i+1)=1;
%                     end
%                 end
%                 boundary_val_all(find(b_flag==1))=[];
                
%                 for b_i=1:2:size(boundary_val_all,2)
%                     signicant_cluster_size=[];
%                     signicant_clusters_all=boundary_val_all(b_i):1:boundary_val_all(b_i+1);
%                     signicant_cluster_size=size(signicant_clusters_all,2);
%                     if(mod(boundary_val_all(b_i),2)==1)
%                         boundary_val_all(b_i)=boundary_val_all(b_i)+1;
%                     end
%                     signicant_clusters_all=boundary_val_all(b_i):1:boundary_val_all(b_i+1);
%                     if signicant_cluster_size > cluster_size_thr
%                         plot(signicant_clusters_all,ones(size(signicant_clusters_all,2),1) ...
%                                 *range_vals_max(fig_order(atlas_index_i))*0.8,'*b','MarkerSize',8)
%                     end
%                 end
                for b_i=1:2:size(boundary_val_all,2)
                    signicant_clusters_all=boundary_val_all(b_i):1:boundary_val_all(b_i+1);
                    if size(signicant_clusters_all,2) > cluster_size_thr
                        plot(signicant_clusters_all,ones(size(signicant_clusters_all,2),1) ...
                                *range_vals_max(fig_order(atlas_index_i))*0.8,'-b','LineWidth',4)
                    end
                end
            end
            
            %% right side
            plot_test_right_base=repmat(nanmean(plot_test_right(:,1:19),2),1,98);

%             [clusters_right, p_values_right, t_sums_right, permutation_distribution_right] ...
%                 =permutest_ram(plot_test_right_base',plot_test_right',true,perm_p,10000,true);
%             cluster_stats_right{task_i,fig_order(atlas_index_i)}.clusters=clusters_right;
%             cluster_stats_right{task_i,fig_order(atlas_index_i)}.p_values=p_values_right;
%             cluster_stats_right{task_i,fig_order(atlas_index_i)}.t_sums=t_sums_right;
%             cluster_stats_right{task_i,fig_order(atlas_index_i)}.permutation_distribution=permutation_distribution_right;
            
            clusters_right=cluster_stats_right{task_i,fig_order(atlas_index_i)}.clusters;
            p_values_right=cluster_stats_right{task_i,fig_order(atlas_index_i)}.p_values;

            signicant_clusters=[];
            signicant_clusters=clusters_right(p_values_right<thr_p_1);
%             signicant_clusters{:}
%             p_values_right
            
            if ~isempty(signicant_clusters)
                boundary_val_all=[];
                for c_i=1:size(signicant_clusters,2)
%                     if size(signicant_clusters{c_i},2)>cluster_size_thr
                        boundary_val_buf=[min(signicant_clusters{c_i}) max(signicant_clusters{c_i})];
                        boundary_val_all=[boundary_val_all boundary_val_buf];
%                     end
                end
                boundary_val_all=sort(boundary_val_all);
                
%                 b_flag=zeros(size(boundary_val_all,2),1);
%                 for b_i=2:2:size(boundary_val_all,2)-1
%                     if (boundary_val_all(b_i+1)-boundary_val_all(b_i)<cluster_size_thr)
%                         b_flag(b_i)=1;
%                         b_flag(b_i+1)=1;
%                     end
%                 end
%                 
%                 boundary_val_all(find(b_flag==1))=[];
                
%                 for b_i=1:2:size(boundary_val_all,2)
%                     signicant_cluster_size=[];
%                     signicant_clusters_all=boundary_val_all(b_i):1:boundary_val_all(b_i+1);
%                     signicant_cluster_size=size(signicant_clusters_all,2);
%                     if(mod(boundary_val_all(b_i),2)==1)
%                         boundary_val_all(b_i)=boundary_val_all(b_i)+1;
%                     end
%                     signicant_clusters_all=boundary_val_all(b_i):1:boundary_val_all(b_i+1);
%                     if signicant_cluster_size > cluster_size_thr
%                         plot(signicant_clusters_all,ones(size(signicant_clusters_all,2),1) ...
%                                 *range_vals_max(fig_order(atlas_index_i))*0.7,'*r','MarkerSize',8)
%                     end
%                 end
                for b_i=1:2:size(boundary_val_all,2)
                    signicant_clusters_all=boundary_val_all(b_i):boundary_val_all(b_i+1);
                    if size(signicant_clusters_all,2) > cluster_size_thr
                        plot(signicant_clusters_all,ones(size(signicant_clusters_all,2),1) ...
                                *range_vals_max(fig_order(atlas_index_i))*0.7,'-r','LineWidth',4)
                    end
                end
            end                               
                           
            %% left vs right
%             [clusters_both, p_values_both, t_sums_both, permutation_distribution_both] ...
%                 =permutest_ram(plot_test_right',plot_test_left',false,perm_p,10000,true);
%             cluster_stats_both{task_i,fig_order(atlas_index_i)}.clusters=clusters_both;
%             cluster_stats_both{task_i,fig_order(atlas_index_i)}.p_values=p_values_both;
%             cluster_stats_both{task_i,fig_order(atlas_index_i)}.t_sums=t_sums_both;
%             cluster_stats_both{task_i,fig_order(atlas_index_i)}.permutation_distribution=permutation_distribution_both;
            
            clusters_both=cluster_stats_both{task_i,fig_order(atlas_index_i)}.clusters;
            p_values_both=cluster_stats_both{task_i,fig_order(atlas_index_i)}.p_values;
            
            signicant_clusters=[];
            signicant_clusters=clusters_both(p_values_both<thr_p_1);
%             signicant_clusters{:}
%             p_values_both
            
            
            if ~isempty(signicant_clusters)
                boundary_val_all=[];
                for c_i=1:size(signicant_clusters,2)
%                     if size(signicant_clusters{c_i},2)>cluster_size_thr
                        boundary_val_buf=[min(signicant_clusters{c_i}) max(signicant_clusters{c_i})];
                        boundary_val_all=[boundary_val_all boundary_val_buf];
%                     end
                end
                boundary_val_all=sort(boundary_val_all);
                
%                 b_flag=zeros(size(boundary_val_all,2),1);
%                 for b_i=2:2:size(boundary_val_all,2)-1
%                     if (boundary_val_all(b_i+1)-boundary_val_all(b_i)<cluster_size_thr)
%                         b_flag(b_i)=1;
%                         b_flag(b_i+1)=1;
%                     end
%                 end
%                 boundary_val_all(find(b_flag==1))=[];
                
%                 for b_i=1:2:size(boundary_val_all,2)
%                     signicant_cluster_size=[];
%                     signicant_clusters_all=boundary_val_all(b_i):1:boundary_val_all(b_i+1);
%                     signicant_cluster_size=size(signicant_clusters_all,2);
%                     if(mod(boundary_val_all(b_i),2)==1)
%                         boundary_val_all(b_i)=boundary_val_all(b_i)+1;
%                     end
%                     signicant_clusters_all=boundary_val_all(b_i):1:boundary_val_all(b_i+1);
%                     if signicant_cluster_size > cluster_size_thr
%                         plot(signicant_clusters_all,ones(size(signicant_clusters_all,2),1) ...
%                                 *range_vals_max(fig_order(atlas_index_i))*0.9,'*g','MarkerSize',8)
%                     end
%                 end
                for b_i=1:2:size(boundary_val_all,2)
                    signicant_clusters_all=boundary_val_all(b_i):boundary_val_all(b_i+1);
                    if size(signicant_clusters_all,2) > cluster_size_thr
                        plot(signicant_clusters_all,ones(size(signicant_clusters_all,2),1) ...
                                *range_vals_max(fig_order(atlas_index_i))*0.9,'-g','LineWidth',4)
                    end
                end
            end    

            %%%%%%%%%%%%%%%%%%%%%%% p<0.01
%             signicant_clusters=[];
%             signicant_clusters=clusters_both(p_values_both<thr_p_2 & p_values_both>thr_p_1);
% 
%             if ~isempty(signicant_clusters)
%                 boundary_val_all=[];
%                 for c_i=1:size(signicant_clusters,2)
% %                     if size(signicant_clusters{c_i},2)>cluster_size_thr
%                         boundary_val_buf=[min(signicant_clusters{c_i}) max(signicant_clusters{c_i})];
%                         boundary_val_all=[boundary_val_all boundary_val_buf];
% %                     end
%                 end
%                 boundary_val_all=sort(boundary_val_all);
%                 
%                 b_flag=zeros(size(boundary_val_all,2),1);
%                 for b_i=2:2:size(boundary_val_all,2)-1
%                     if (boundary_val_all(b_i+1)-boundary_val_all(b_i)<cluster_size_thr)
%                         b_flag(b_i)=1;
%                         b_flag(b_i+1)=1;
%                     end
%                 end
%                 
%                 boundary_val_all(find(b_flag==1))=[];
%                 for b_i=1:2:size(boundary_val_all,2)
%                     signicant_clusters_all=boundary_val_all(b_i):boundary_val_all(b_i+1);
%                     if size(signicant_clusters_all,2) > cluster_size_thr
%                         plot(signicant_clusters_all,ones(size(signicant_clusters_all,2),1) ...
%                                 *range_vals_max(fig_order(atlas_index_i))*0.9,':g','LineWidth',2)
%                     end
%                 end
%             end

            %%
%             if atlas_index_i==1
%                 legend([h1 h2],{'Left','Right'},'LineWidth',1.5,'FontSize',10)
%             end
            
            if atlas_index_i==12           
                set(gca,'xtick',[1 20 40 60 80 100]);
                set(gca,'xticklabel',[-500 0 500 1000 1500 2000]);
                set(gca,'xlim',[1 100])
                xlabel('Time (ms)');
                % short1
%                 set(gca,'xtick',[1 11 20 30 40]);
%                 set(gca,'xticklabel',[-500 -250 0 250 500]);
%                 set(gca,'xlim',[1 40])
%                 xlabel('Time (ms)');
                % short2
%                 set(gca,'xtick',[11 20 30]);
%                 set(gca,'xticklabel',[-250 0 250]);
%                 set(gca,'xlim',[11 30])
%                 xlabel('Time (ms)');
            else
                set(gca,'xtick',[1 20 40 60 80 100]);
                set(gca,'xticklabel',[]);
                set(gca,'xlim',[1 100])
                %short1
%                 set(gca,'xtick',[1 11 20 30 40]);
%                 set(gca,'xticklabel',[-500 -250 0 250 500]);
%                 set(gca,'xlim',[1 40])
                %short2
%                 set(gca,'xtick',[11 20 30]);
%                 set(gca,'xticklabel',[-250 0 250]);
%                 set(gca,'xlim',[11 29])
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
            y1=yline(0,'k','LineWidth',1,'Fontsize',15);

            %         line([20,20],[-8,10],'Color','r','LineWidth',1)
            %         line([84,84],[-8,10],'Color','r','LineWidth',1)
%             F_name_new{atlas_index_i,1}=atlas_index_all_new_name{fig_order(atlas_index_i)};
            title(atlas_index_all_new_name{fig_order(atlas_index_i)},'FontSize',20,'Color',atlas_color(atlas_index_i,:))
        
        end
        subplot(4,2,7)
%         subaxis(4,2,size(atlas_index_all_new,2)+1,'SpacingVert',0.05,'MR',0);
        ta=plot(1,7,'-b','LineWidth',1);hold on;
        tb=plot(1,8,'-r','LineWidth',1);hold on;
%         ha=plot(1,1,'-b','MarkerSize',8);hold on;
%         hc=plot(1,3,'-r','MarkerSize',8);hold on;
%         he=plot(1,5,'-g','MarkerSize',8);hold on;
        ha=plot(1,1,'-b','LineWidth',4);hold on;
        hc=plot(1,3,'-r','LineWidth',4);hold on;
        he=plot(1,5,'-g','LineWidth',4);hold on;

%         legend([ta tb he hf ha hb hc hd],{'Left','Right','Left vs Right, p<0.05','Left vs Right, p<0.1','Left vs Left baseline, p<0.05','Left vs Left baseline, p<0.1','Right vs Right baseline, p<0.05','Right vs Right baseline, p<0.1'}, ...
%             'LineWidth',1.5,'FontSize',15,'Location','best')
        legend([ta tb he ha hc],{'Left','Right','Left vs Right, p<0.05','Left vs Left baseline, p<0.05','Right vs Right baseline, p<0.05'}, ...
            'LineWidth',1.5,'FontSize',15,'Location','best')
        set(gca,'xtick',[]);
        set(gca,'xticklabel',[]);
        set(gca,'ytick',[]);
        set(gca,'yticklabel',[]);
        set(gca,'box','off')

        [ax,h3]=suplabel(task_name_fig{task_i},'t');
        set(h3,'FontSize',30)
%         [ax2,h4]=suplabel({'Summation of significant t values'},'y');
        [ax2,h4]=suplabel({'Gamma Power Z-score'},'y');
        set(h4,'FontSize',30)
        %% save figure
%         saveas(gcf,[zscore_folder '\' task_name_fig{task_i} '_all_plot.png'])

        export_fig([zscore_folder '\' roi_result_folder '\' task_name_fig{task_i} '_all_plot'],'-tiff')
        saveas(gcf,[zscore_folder '\' roi_result_folder '\' task_name_fig{task_i} '_all_plot'],'epsc')
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
%             [h1,output1]=plot_areaerrorbar_new(plot_test_left,1);
% 
%             hold on;
% 
%             plot_test_right=squeeze(plot_test(2,task_i,fig_order(atlas_index_i),:,:));
% %             plot_test_right=squeeze(plot_test_weighted(2,task_i,atlas_index_i,:,:));
%             plot_test_right(plot_test_right==0)=nan;
%             [h2,output2]=plot_areaerrorbar_new(plot_test_right,2);
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
%         [ax2,h4]=suplabel({'Gamma Power Z-score'},'y');
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
% save('cluster_stats_all.mat','cluster_stats_left','cluster_stats_right','cluster_stats_both')

% save('xylim_vals_all.mat','range_vals_min','range_vals_max')



