
%% zscore plot by ROIs
clc;clear;
rootfolder = 'E:\RAM data set\RAM_Public_Data_all';
cd([rootfolder '\FR1_FS\Storage_filtering_overlap_base_agg\Data_storage'])

file_suffix='Final_std_20_fs';
Montage_suffix='final_roi_fs';
load fs_DK_atlas.mat
flip_flag=0;
atlas_index_all=[4 30 8 12 13 15 28 29];

if flip_flag==1
    sides_num = 1;
    flip_prefix='xflip';
else
    sides_num = 2;
    flip_prefix='both';
end

zscore_folder = [rootfolder '\FR1_FS\Storage_filtering_overlap_base_agg\zscore_plot\' Montage_suffix ];
mkdir(zscore_folder)

if flip_flag ==0 
    task_name={'rec','Nonrec','together','diff'};
    plot_test=nan(2,4,36,251,118);
    plot_test_elec=nan(2,4,36,251,118);
    plot_test_elec_sig=nan(2,4,36,251,118);
    plot_test_sig=zeros(2,4,36,251,118);
    plot_test_stats_ratio=nan(2,4,36,118);
    plot_test_stats_avg_t=nan(2,4,36,118);
    plot_test_stats_avg_t_sig=nan(2,4,36,118);
    
    [L_vertices, L_label, L_colortable] = read_annotation('lh.aparc.annot');
    [R_vertices, R_label, R_colortable] = read_annotation('rh.aparc.annot');
    [~, L_400_label, L_400_colortable]=read_annotation('lh.Schaefer2018_400Parcels_7Networks_order.annot');
    [~, R_400_label, R_400_colortable]=read_annotation('rh.Schaefer2018_400Parcels_7Networks_order.annot');
    [size_L_400_label, index_L_400_label] = hist(L_400_label,unique(L_400_label));
    [size_R_400_label, index_R_400_label] = hist(R_400_label,unique(R_400_label));
    
elseif flip_flag ==1
    task_name={'rec_xflip','Nonrec_xflip','together_xflip','diff_xflip'};
    plot_test=nan(1,4,36,251,118);
    plot_test_elec=nan(1,4,36,251,118);
    plot_test_elec_sig=nan(1,4,36,251,118);
    plot_test_sig=zeros(1,4,36,251,118);
    plot_test_stats_ratio=nan(1,4,36,118);
    plot_test_stats_avg_t=nan(1,4,36,118);
    plot_test_stats_avg_t_sig=nan(1,4,36,118);

    [L_vertices, L_label, L_colortable] = read_annotation('lh.aparc.annot');
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
%     for time_i=1:118
%         if flip_flag == 0
%             %% Left
%             vertex_values_frame=[];
%             vertex_values_frame_e=[];
%             vertex_values_frame_L_buff=[];
%             vertex_values_frame_L_buff=vertex_values_frame_L(:,time_i);
%             load([file_suffix '_' task_name{task_i} '/' Montage_suffix '/L_vertex_values_' file_suffix '_' task_name{task_i} '_' num2str(time_i)]);
%             for atlas_index=2:36
%                 roi_index_L=[];
%                 roi_index_L=find(L_label==L_colortable.table(atlas_index,5));
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
%                 % number of significant parcels / total number
%                 all_values_L=[];
%                 all_values_L=unique(vertex_values_frame_L_buff(roi_index_L));
%                 plot_test_stats_ratio(1,task_i,atlas_index,time_i)=(length(find(all_values_L>0)) - length(find(all_values_L<0))) ...
%                     /length(unique(L_400_label(roi_index_L)));
%                 
%                 % average t values
%                 plot_test_stats_avg_t(1,task_i,atlas_index,time_i) = nanmean(vertex_values_frame_L_buff(roi_index_L));
%                 
%                 % average sig t values
%                 if all_values_L == 0
%                     plot_test_stats_avg_t_sig(1,task_i,atlas_index,time_i) = 0;
%                 else
%                     plot_test_stats_avg_t_sig(1,task_i,atlas_index,time_i) = mean(all_values_L(all_values_L~=0));
%                 end
%                 
%                 % average z-score by all subjects
%                 plot_test(1,task_i,atlas_index,:,time_i)=nanmean(vertex_values_frame(roi_index_L,:),1);
%                 plot_test_elec(1,task_i,atlas_index,:,time_i)=nanmean(vertex_values_frame_e(roi_index_L,:),1);
% 
%                 % average significant z-score by all subjects
%                 if plot_test_stats_ratio(1,task_i,atlas_index,time_i) ~= 0
%                     roi_index_L(vertex_values_frame_L_buff(roi_index_L)==0)=[];
%                     plot_test_sig(1,task_i,atlas_index,:,time_i)=nanmean(vertex_values_frame(roi_index_L,:),1);
%                     plot_test_elec_sig(1,task_i,atlas_index,:,time_i)=nanmean(vertex_values_frame_e(roi_index_L,:),1);
%                 end
%             end
%             
%             vertex_values_frame=[];
%             load([file_suffix '_' task_name{task_i} '/' Montage_suffix '/R_vertex_values_' file_suffix '_' task_name{task_i} '_' num2str(time_i)]);
%         
%             %% Right
%             vertex_values_frame=[];
%             vertex_values_frame_e=[];
%             vertex_values_frame_R_buff=[];
%             vertex_values_frame_R_buff=vertex_values_frame_R(:,time_i);
%             load([file_suffix '_' task_name{task_i} '/' Montage_suffix '/R_vertex_values_' file_suffix '_' task_name{task_i} '_' num2str(time_i)]);
%             for atlas_index=2:36
%                 roi_index_R=[];
%                 roi_index_R=find(R_label==R_colortable.table(atlas_index,5));
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
%                 % number of significant parcels / total number
%                 all_values_R=[];
%                 all_values_R=unique(vertex_values_frame_R_buff(roi_index_R));
%                 plot_test_stats_ratio(2,task_i,atlas_index,time_i)=(length(find(all_values_R>0)) - length(find(all_values_R<0))) ...
%                     /length(unique(R_400_label(roi_index_R)));
%                 
%                 % average t values
%                 plot_test_stats_avg_t(2,task_i,atlas_index,time_i) = nanmean(vertex_values_frame_R_buff(roi_index_R));
%                 
%                 % average sig t values
%                 if all_values_R == 0
%                     plot_test_stats_avg_t_sig(2,task_i,atlas_index,time_i) = 0;
%                 else
%                     plot_test_stats_avg_t_sig(2,task_i,atlas_index,time_i) = mean(all_values_R(all_values_R~=0));
%                 end
%                 
%                 % average z-score by all subjects
%                 plot_test(2,task_i,atlas_index,:,time_i)=nanmean(vertex_values_frame(roi_index_R,:),1);
%                 plot_test_elec(2,task_i,atlas_index,:,time_i)=nanmean(vertex_values_frame_e(roi_index_R,:),1);
% 
%                 % average significant z-score by all subjects
%                 if plot_test_stats_ratio(2,task_i,atlas_index,time_i) ~= 0
%                     roi_index_R(vertex_values_frame_R_buff(roi_index_R)==0)=[];
%                     plot_test_sig(2,task_i,atlas_index,:,time_i)=nanmean(vertex_values_frame(roi_index_R,:),1);
%                     plot_test_elec_sig(2,task_i,atlas_index,:,time_i)=nanmean(vertex_values_frame_e(roi_index_R,:),1);
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

% save([file_suffix '_' Montage_suffix '_roi_val_' flip_prefix '.mat'],'plot_test','plot_test_sig' ...
%     ,'plot_test_stats_ratio','plot_test_stats_avg_t','plot_test_stats_avg_t_sig','plot_test_elec','plot_test_elec_sig');
load([file_suffix '_' Montage_suffix '_roi_val_' flip_prefix '.mat']);

DK_atlas_names{4}='Caudal Middle Frontal';
DK_atlas_names{8}='Fusiform Gyrus';
DK_atlas_names{12}='Lateral Occipital';
DK_atlas_names{13}='Lateral Orbitofrontal';
DK_atlas_names{15}='Medial Orbitofrontal';
DK_atlas_names{28}='Rostral Middle Frontal';
DK_atlas_names{29}='Superior Frontal';
DK_atlas_names{30}='Superior Parietal';

font_color=L_colortable.table(:,1:3)/255;

%% plot signicant parcels
% 
% if flip_flag == 0
%     task_name_fig={'Recalled','Not recalled','Recalled + Not recalled','Recalled - Not recalled'};
%     
%     for task_i=1:4
%         cnt=1;
%         figure('position',[0 0 2000 8000]);
%         set(gcf, 'color', [1 1 1]);
%         set(gcf,'Visible','on');
%         for atlas_index=atlas_index_all
%             subplot(4,2,cnt);
%             plot_test_left=squeeze(plot_test_stats(1,task_i,atlas_index,:));
%             h1=plot(plot_test_left,'color',[52 148 186]./255,'LineWidth', 4);
%             hold on;
% 
%             plot_test_right=squeeze(plot_test_stats(2,task_i,atlas_index,:));
%             h2=plot(plot_test_right,'color',[236 112  22]./255,'LineWidth', 4);
% 
%             if cnt == 7 | cnt == 8
%                 set(gca,'xtick',[1 20 40 60 80 100 118]);
%                 set(gca,'xticklabel',[-500 0 500 1000 1500 2000 2500]);
%                 set(gca,'xlim',[1 118])
%                 xlabel('Time (ms)');
%             else
%                 set(gca,'xtick',[]);
%                 set(gca,'xticklabel',[]);
%             end
% 
%             set(gca,'box','off')
% 
%             set(gca,'ylim',[-1 1])
%             set(gca,'Fontsize',15);
% 
%             x1=xline(20,'-.r',{'ON'},'LabelVerticalAlignment','top','LabelHorizontalAlignment','center','LabelOrientation','horizontal','LineWidth',3,'Fontsize',15);
%             x2=xline(84,'-.r',{'OFF'},'LabelVerticalAlignment','top','LabelHorizontalAlignment','center','LabelOrientation','horizontal','LineWidth',3,'Fontsize',15);
%             y1=yline(0,'k','LineWidth',1,'Fontsize',15);
% 
%             %         line([20,20],[-8,10],'Color','r','LineWidth',1)
%             %         line([84,84],[-8,10],'Color','r','LineWidth',1)
% 
% 
%             if cnt==1
%                 ylabel({'number of significant parcels /','total number of parcels'});
%                 legend([h1 h2],{'Left','Right'},'LineWidth',1.5,'FontSize',10)
%             end
% 
%             title(DK_atlas_names{atlas_index},'FontSize',20,'Color',font_color(atlas_index,:))
%     %         [ax,h1]=suplabel(['ROI : ' DK_atlas_names{atlas_index} ', Condition : ' task_name_fig{task_i} ],'t');
% %                 [ax,h3]=suplabel([DK_atlas_names{atlas_index} ],'t');
% %                 set(h3,'FontSize',20,'Color',font_color(atlas_index,:))
%             cnt=cnt+1;
%         end  
%         [ax,h3]=suplabel(task_name_fig{task_i},'t');
%         set(h3,'FontSize',30)
%         %% save figure
% %         saveas(gcf,[zscore_folder '\' task_name_fig{task_i} '_all_plot.png'])
% 
%         export_fig([zscore_folder '\' task_name_fig{task_i} '_all_plot'],'-tiff')
%         saveas(gcf,[zscore_folder '\' task_name_fig{task_i} '_all_plot'],'epsc')
%         close all;
%     end
% elseif flip_flag == 1
%     
%     
%     
%     
% end

%% plot signicant parcels ratio
% 
% if flip_flag == 0
%     task_name_fig={'Recalled','Not recalled','Recalled + Not recalled','Recalled - Not recalled'};
%     
%     for task_i=1:4
%         cnt=1;
%         figure('position',[0 0 2000 8000]);
%         set(gcf, 'color', [1 1 1]);
%         set(gcf,'Visible','on');
%         for atlas_index=atlas_index_all
%             subplot(4,2,cnt);
% %             plot_test_left=squeeze(plot_test_stats_ratio(1,task_i,atlas_index,:));
%             plot_test_left=squeeze(plot_test_stats_avg_t(1,task_i,atlas_index,:));
% %             plot_test_left=squeeze(plot_test_stats_avg_t_sig(1,task_i,atlas_index,:));
%             
%             h1=plot(plot_test_left,'color',[52 148 186]./255,'LineWidth', 4);
%             hold on;
% 
% %             plot_test_right=squeeze(plot_test_stats_ratio(2,task_i,atlas_index,:));
%             plot_test_right=squeeze(plot_test_stats_avg_t(2,task_i,atlas_index,:));
% %             plot_test_right=squeeze(plot_test_stats_avg_t_sig(2,task_i,atlas_index,:));
%             
%             h2=plot(plot_test_right,'color',[236 112  22]./255,'LineWidth', 4);
% 
%             if cnt == 7 | cnt == 8
%                 set(gca,'xtick',[1 20 40 60 80 100 118]);
%                 set(gca,'xticklabel',[-500 0 500 1000 1500 2000 2500]);
%                 set(gca,'xlim',[1 118])
%                 xlabel('Time (ms)');
%             else
%                 set(gca,'xtick',[]);
%                 set(gca,'xticklabel',[]);
%             end
% 
%             set(gca,'box','off')
% 
% %             set(gca,'ylim',[-1 1])
%             set(gca,'ylim',[-6 6])
%             set(gca,'Fontsize',15);
% 
%             x1=xline(20,'-.r',{'ON'},'LabelVerticalAlignment','top','LabelHorizontalAlignment','center','LabelOrientation','horizontal','LineWidth',3,'Fontsize',15);
%             x2=xline(84,'-.r',{'OFF'},'LabelVerticalAlignment','top','LabelHorizontalAlignment','center','LabelOrientation','horizontal','LineWidth',3,'Fontsize',15);
%             y1=yline(0,'k','LineWidth',1,'Fontsize',15);
% 
%             %         line([20,20],[-8,10],'Color','r','LineWidth',1)
%             %         line([84,84],[-8,10],'Color','r','LineWidth',1)
% 
% 
%             if cnt==1
% %                 ylabel({'number of significant parcels /','total number of parcels'});
%                 ylabel({'Average t-value'});
%                 legend([h1 h2],{'Left','Right'},'LineWidth',1.5,'FontSize',10)
%             end
% 
%             title(DK_atlas_names{atlas_index},'FontSize',20,'Color',font_color(atlas_index,:))
%     %         [ax,h1]=suplabel(['ROI : ' DK_atlas_names{atlas_index} ', Condition : ' task_name_fig{task_i} ],'t');
% %                 [ax,h3]=suplabel([DK_atlas_names{atlas_index} ],'t');
% %                 set(h3,'FontSize',20,'Color',font_color(atlas_index,:))
%             cnt=cnt+1;
%         end  
%         [ax,h3]=suplabel(task_name_fig{task_i},'t');
%         set(h3,'FontSize',30)
%         %% save figure
% %         saveas(gcf,[zscore_folder '\' task_name_fig{task_i} '_all_plot.png'])
% 
%         export_fig([zscore_folder '\' task_name_fig{task_i} '_all_plot'],'-tiff')
%         saveas(gcf,[zscore_folder '\' task_name_fig{task_i} '_all_plot'],'epsc')
%         close all;
%     end
% elseif flip_flag == 1
%     
%     
%     
%     
% end

%% plot each
% if flip_flag == 0
%     task_name_fig={'Recall','NonRecall','Recall_P_NonRecall','Recall_N_NonRecall'};
%     for task_i=1:4
%         zscore_folder_n=[zscore_folder '\' task_name_fig{task_i}];
%         mkdir(zscore_folder_n)
%         cnt=1;
%         for atlas_index=atlas_index_all
%             max_val=max(max(max(max(squeeze(plot_test(:,:,atlas_index,:,:))))));
%             min_val=min(min(min(min(squeeze(plot_test(:,:,atlas_index,:,:))))));  
%             if ~isnan(max_val)
%                 figure('position',[0 0 2000 2000]);
%                 set(gcf, 'color', [1 1 1]);
%                 set(gcf,'Visible','on');
%                 plot_test_left=squeeze(plot_test(1,task_i,atlas_index,:,:));
%                 plot_test_left(plot_test_left==0)=nan;
%                 h1=plot_areaerrorbar(plot_test_left,1);
%                 hold on;
%                 
%                 plot_test_right=squeeze(plot_test(2,task_i,atlas_index,:,:));
%                 plot_test_right(plot_test_right==0)=nan;
%                 h2=plot_areaerrorbar(plot_test_right,2);
%                 
%                 set(gca,'xtick',[1 20 40 60 80 100 118]);
%                 set(gca,'xticklabel',[-500 0 500 1000 1500 2000 2500]);
%                 set(gca,'xlim',[1 118])
%                 set(gca,'box','off')
%                 
%                 
% 
%             %         set(gca,'ytick',[-6 -2 0 2 6 ]);
%             %         set(gca,'yticklabel',[-6 -2 0 2 6]);
% %                 set(gca,'ylim',[min_val max_val])
%                 set(gca,'Fontsize',40);
% 
%                 x1=xline(20,'-.r',{'Word','ON'},'LabelVerticalAlignment','top','LabelHorizontalAlignment','center','LabelOrientation','horizontal','LineWidth',2,'Fontsize',30);
%                 x2=xline(84,'-.r',{'Word','OFF'},'LabelVerticalAlignment','top','LabelHorizontalAlignment','center','LabelOrientation','horizontal','LineWidth',2,'Fontsize',30);
%                 y1=yline(0,'k','LineWidth',3,'Fontsize',30);
% 
%                 %         line([20,20],[-8,10],'Color','r','LineWidth',1)
%             %         line([84,84],[-8,10],'Color','r','LineWidth',1)
%                 ylabel({'Gamma power z-score'});
%                 xlabel('Time (ms)');
%                 
%                 legend([h1 h2],{'Left','Right'},'LineWidth',1,'FontSize',40)
%                 
%                 
%                 %% save figure
%         %         [ax,h1]=suplabel(['ROI : ' DK_atlas_names{atlas_index} ', Condition : ' task_name_fig{task_i} ],'t');
%                 [ax,h3]=suplabel([DK_atlas_names{atlas_index} ],'t');
%                 set(h3,'FontSize',50,'Color',font_color(atlas_index,:))
% %                 saveas(gcf,[zscore_folder_n '\' num2str(atlas_index) '_' DK_atlas_names{atlas_index} '_' task_name_fig{task_i} '_plot.png'])
%                 export_fig([zscore_folder_n '\' num2str(atlas_index) '_' DK_atlas_names{atlas_index} '_' task_name_fig{task_i} '_plot','-tiff'])
%                 
%                 close all;
%             end
%             cnt=cnt+1;
%         end
%     end
% elseif flip_flag == 1
%     
%     
%     
%     
% end

%% plot all together
size_atlas_index_all=size(atlas_index_all,2);

% max_vals=nan(4,size_atlas_index_all,2);
% min_vals=nan(4,size_atlas_index_all,2);
load('xylim_vals_all.mat')
for val_i=1:size_atlas_index_all
    ylim_max_all(val_i)=max(max(max_vals(:,val_i,:)));
    ylim_min_all(val_i)=min(min(min_vals(:,val_i,:)));
end


if flip_flag == 0
    task_name_fig={'Recalled','Not recalled','Recalled + Not recalled','Recalled - Not recalled'};
    
    for task_i=1:4
        cnt=1;
        figure('position',[0 0 2000 8000]);
        set(gcf, 'color', [1 1 1]);
        set(gcf,'Visible','on');
        for atlas_index=atlas_index_all
            max_val=max(max(max(max(squeeze(plot_test(:,:,atlas_index,:,:))))));
            min_val=min(min(min(min(squeeze(plot_test(:,:,atlas_index,:,:))))));  
            if ~isnan(max_val)
                subplot(4,2,cnt);
                plot_test_left=squeeze(plot_test(1,task_i,atlas_index,:,:));
                plot_test_left(plot_test_left==0)=nan;
                [h1,output1]=plot_areaerrorbar(plot_test_left,1);
                max_vals(task_i,cnt,1)=output1.max_val;
                min_vals(task_i,cnt,1)=output1.min_val;

                hold on;
                
                plot_test_right=squeeze(plot_test(2,task_i,atlas_index,:,:));
                plot_test_right(plot_test_right==0)=nan;
                [h2,output2]=plot_areaerrorbar(plot_test_right,2);
                max_vals(task_i,cnt,2)=output2.max_val;
                min_vals(task_i,cnt,2)=output2.min_val;
                

                if cnt == 7 | cnt == 8
%                     set(gca,'xtick',[1 20 40 60 80 100 118]);
%                     set(gca,'xticklabel',[-500 0 500 1000 1500 2000 2500]);
%                     set(gca,'xlim',[1 118])
                    set(gca,'xtick',[1 20 40 60 80 100]);
                    set(gca,'xticklabel',[-500 0 500 1000 1500 2000]);
                    set(gca,'xlim',[1 100])

                    xlabel('Time (ms)');
                else
                    set(gca,'xtick',[]);
                    set(gca,'xticklabel',[]);
                end
                
                set(gca,'box','off')

                set(gca,'ylim',[ylim_min_all(cnt) ylim_max_all(cnt)])
                set(gca,'Fontsize',15);

                x1=xline(20,'-.r',{'ON'},'LabelVerticalAlignment','top','LabelHorizontalAlignment','center','LabelOrientation','horizontal','LineWidth',3,'Fontsize',15);
                x2=xline(84,'-.r',{'OFF'},'LabelVerticalAlignment','top','LabelHorizontalAlignment','center','LabelOrientation','horizontal','LineWidth',3,'Fontsize',15);
                y1=yline(0,'k','LineWidth',1,'Fontsize',15);

                %         line([20,20],[-8,10],'Color','r','LineWidth',1)
                %         line([84,84],[-8,10],'Color','r','LineWidth',1)
                
                
                if cnt==1
                    ylabel({'Gamma power z-score'});
                    legend([h1 h2],{'Left','Right'},'LineWidth',1.5,'FontSize',10)
                end
                
                title(DK_atlas_names{atlas_index},'FontSize',20,'Color',font_color(atlas_index,:))
        %         [ax,h1]=suplabel(['ROI : ' DK_atlas_names{atlas_index} ', Condition : ' task_name_fig{task_i} ],'t');
%                 [ax,h3]=suplabel([DK_atlas_names{atlas_index} ],'t');
%                 set(h3,'FontSize',20,'Color',font_color(atlas_index,:))
            end
            cnt=cnt+1;
        end  
        [ax,h3]=suplabel(task_name_fig{task_i},'t');
        set(h3,'FontSize',30)
        %% save figure
%         saveas(gcf,[zscore_folder '\' task_name_fig{task_i} '_all_plot.png'])

        export_fig([zscore_folder '\' task_name_fig{task_i} '_all_plot'],'-tiff')
        saveas(gcf,[zscore_folder '\' task_name_fig{task_i} '_all_plot'],'epsc')
        close all;
    end
elseif flip_flag == 1
    
    
    
    
end

save('xylim_vals_all.mat','max_vals','min_vals')

%% plot all together (significant parcels)
% size_atlas_index_all=size(atlas_index_all,2);
% 
% % max_vals=nan(4,size_atlas_index_all,2);
% % min_vals=nan(4,size_atlas_index_all,2);
% load('xylim_vals_all.mat')
% for val_i=1:size_atlas_index_all
%     ylim_max_all(val_i)=max(max(max_vals(:,val_i,:)));
%     ylim_min_all(val_i)=min(min(min_vals(:,val_i,:)));
% end
% 
% 
% if flip_flag == 0
%     task_name_fig={'Recalled','Not recalled','Recalled + Not recalled','Recalled - Not recalled'};
%     
%     for task_i=1:4
%         cnt=1;
%         figure('position',[0 0 2000 8000]);
%         set(gcf, 'color', [1 1 1]);
%         set(gcf,'Visible','on');
%         for atlas_index=atlas_index_all
%             max_val=max(max(max(max(squeeze(plot_test_sig(:,:,atlas_index,:,:))))));
%             min_val=min(min(min(min(squeeze(plot_test_sig(:,:,atlas_index,:,:))))));  
%             if ~isnan(max_val)
%                 subplot(4,2,cnt);
%                 plot_test_left=squeeze(plot_test_sig(1,task_i,atlas_index,:,:));
%                 plot_test_left(plot_test_left==0)=nan;
%                 [h1,output1]=plot_areaerrorbar(plot_test_left,1);
%                 max_vals(task_i,cnt,1)=output1.max_val;
%                 min_vals(task_i,cnt,1)=output1.min_val;
% 
%                 hold on;
%                 
%                 plot_test_right=squeeze(plot_test_sig(2,task_i,atlas_index,:,:));
%                 plot_test_right(plot_test_right==0)=nan;
%                 [h2,output2]=plot_areaerrorbar(plot_test_right,2);
%                 max_vals(task_i,cnt,2)=output2.max_val;
%                 min_vals(task_i,cnt,2)=output2.min_val;
%                 
% 
%                 if cnt == 7 | cnt == 8
%                     set(gca,'xtick',[1 20 40 60 80 100 118]);
%                     set(gca,'xticklabel',[-500 0 500 1000 1500 2000 2500]);
%                     set(gca,'xlim',[1 118])
%                     xlabel('Time (ms)');
%                 else
%                     set(gca,'xtick',[]);
%                     set(gca,'xticklabel',[]);
%                 end
%                 
%                 set(gca,'box','off')
% 
%                 set(gca,'ylim',[ylim_min_all(cnt) ylim_max_all(cnt)])
%                 set(gca,'Fontsize',15);
% 
%                 x1=xline(20,'-.r',{'ON'},'LabelVerticalAlignment','top','LabelHorizontalAlignment','center','LabelOrientation','horizontal','LineWidth',3,'Fontsize',15);
%                 x2=xline(84,'-.r',{'OFF'},'LabelVerticalAlignment','top','LabelHorizontalAlignment','center','LabelOrientation','horizontal','LineWidth',3,'Fontsize',15);
%                 y1=yline(0,'k','LineWidth',1,'Fontsize',15);
% 
%                 %         line([20,20],[-8,10],'Color','r','LineWidth',1)
%                 %         line([84,84],[-8,10],'Color','r','LineWidth',1)
%                 
%                 
%                 if cnt==1
%                     ylabel({'Gamma power z-score'});
%                     legend([h1 h2],{'Left','Right'},'LineWidth',1.5,'FontSize',10)
%                 end
%                 
%                 title(DK_atlas_names{atlas_index},'FontSize',20,'Color',font_color(atlas_index,:))
%         %         [ax,h1]=suplabel(['ROI : ' DK_atlas_names{atlas_index} ', Condition : ' task_name_fig{task_i} ],'t');
% %                 [ax,h3]=suplabel([DK_atlas_names{atlas_index} ],'t');
% %                 set(h3,'FontSize',20,'Color',font_color(atlas_index,:))
%             end
%             cnt=cnt+1;
%         end  
%         [ax,h3]=suplabel(task_name_fig{task_i},'t');
%         set(h3,'FontSize',30)
%         %% save figure
% %         saveas(gcf,[zscore_folder '\' task_name_fig{task_i} '_all_plot.png'])
% 
%         export_fig([zscore_folder '\' task_name_fig{task_i} '_all_plot'],'-tiff')
%         saveas(gcf,[zscore_folder '\' task_name_fig{task_i} '_all_plot'],'epsc')
%         close all;
%     end
% elseif flip_flag == 1
%     
%     
%     
%     
% end
% 
% % save('xylim_vals_all.mat','max_vals','min_vals')

%% plot all together (t-test)
% abs_z_thr=0;
% if flip_flag == 0
%     task_name_fig={'Recalled','Not recalled','Recalled + Not recalled','Recalled - Not recalled'};
%     
%     for task_i=1:4
%         cnt=1;
%         figure('position',[0 0 2000 8000]);
%         set(gcf, 'color', [1 1 1]);
%         set(gcf,'Visible','on');
%         for atlas_index=atlas_index_all
%             max_val=max(max(max(max(squeeze(plot_test(:,:,atlas_index,:,:))))));
%             min_val=min(min(min(min(squeeze(plot_test(:,:,atlas_index,:,:))))));  
%             if ~isnan(max_val)
%                 subplot(4,2,cnt);
%                 plot_test_left=squeeze(plot_test(1,task_i,atlas_index,:,:));
%                 plot_test_left(plot_test_left==0)=nan;
%                 plot_test_left(max(abs(plot_test_left'))<abs_z_thr,:)=NaN;
%                 
%                 baseline_buf=[];
%                 baseline_buf=nanmean(plot_test_left(:,1:19),2);
%                 [H,P,CI,STATS]=ttest(plot_test_left-baseline_buf);
% %                 for time_i=1:118
% %                     stimul_buf=plot_test_left(:,time_i);
% %                     [H,P,CI,STATS]=ttest(baseline_buf,stimul_buf,'tail','both')
% %                 end
%                 h1=plot(STATS.tstat,'color',[52 148 186]./255,'LineWidth', 4);
%                 yline(tinv(0.975,STATS.df(1)),'color',[52 148 186]./255,'LineWidth',1);
%                 yline(tinv(0.025,STATS.df(1)),'color',[52 148 186]./255,'LineWidth',1);
%                 hold on;
% 
%                 plot_test_right=squeeze(plot_test(2,task_i,atlas_index,:,:));
%                 plot_test_right(plot_test_right==0)=nan;
%                 plot_test_right(max(abs(plot_test_right'))<abs_z_thr,:)=NaN;
%                 
%                 baseline_buf=[];
%                 baseline_buf=nanmean(plot_test_right(:,1:19),2);
%                 [H,P,CI,STATS]=ttest(plot_test_right-baseline_buf);
% %                 for time_i=1:118
% %                     stimul_buf=plot_test_left(:,time_i);
% %                     [H,P,CI,STATS]=ttest(baseline_buf,stimul_buf,'tail','both')
% %                 end
%                 h2=plot(STATS.tstat,'color',[236 112  22]./255,'LineWidth', 4);
%                 yline(tinv(0.975,STATS.df(1)),'color',[236 112  22]./255,'LineWidth',1);
%                 yline(tinv(0.025,STATS.df(1)),'color',[236 112  22]./255,'LineWidth',1);
%                 hold on;
%                 
% 
%                 if cnt == 7 | cnt == 8
%                     set(gca,'xtick',[1 20 40 60 80 100 118]);
%                     set(gca,'xticklabel',[-500 0 500 1000 1500 2000 2500]);
%                     set(gca,'xlim',[1 118])
%                     xlabel('Time (ms)');
%                 else
%                     set(gca,'xtick',[]);
%                     set(gca,'xticklabel',[]);
%                 end
%                 
%                 set(gca,'box','off')
% 
%                 set(gca,'ylim',[-6 6])
%                 set(gca,'Fontsize',15);
% 
%                 x1=xline(20,'-.r',{'ON'},'LabelVerticalAlignment','top','LabelHorizontalAlignment','center','LabelOrientation','horizontal','LineWidth',3,'Fontsize',15);
%                 x2=xline(84,'-.r',{'OFF'},'LabelVerticalAlignment','top','LabelHorizontalAlignment','center','LabelOrientation','horizontal','LineWidth',3,'Fontsize',15);
%                 y1=yline(0,'k','LineWidth',1,'Fontsize',15);
% 
%                 %         line([20,20],[-8,10],'Color','r','LineWidth',1)
%                 %         line([84,84],[-8,10],'Color','r','LineWidth',1)
%                 
%                 
%                 if cnt==1
%                     ylabel({'t value'});
%                     legend([h1 h2],{'Left','Right'},'LineWidth',1.5,'FontSize',10)
%                 end
%                 
%                 title([DK_atlas_names{atlas_index} ' dof:' num2str(STATS.df(1))],'FontSize',20,'Color',font_color(atlas_index,:))
%         %         [ax,h1]=suplabel(['ROI : ' DK_atlas_names{atlas_index} ', Condition : ' task_name_fig{task_i} ],'t');
% %                 [ax,h3]=suplabel([DK_atlas_names{atlas_index} ],'t');
% %                 set(h3,'FontSize',20,'Color',font_color(atlas_index,:))
%             end
%             cnt=cnt+1;
%         end  
%         [ax,h3]=suplabel([task_name_fig{task_i} ' ' num2str(abs_z_thr)],'t');
%         set(h3,'FontSize',30)
%         %% save figure
% %         saveas(gcf,[zscore_folder '\' task_name_fig{task_i} '_all_plot.png'])
% 
%         export_fig([zscore_folder '\' task_name_fig{task_i} '_all_plot'],'-tiff')
%         saveas(gcf,[zscore_folder '\' task_name_fig{task_i} '_all_plot'],'epsc')
%         close all;
%     end
% elseif flip_flag == 1
%     
%     
%     
%     
% end

% save('xylim_vals_all.mat','max_vals','min_vals')





%% electrode level
%% plot all together
% size_atlas_index_all=size(atlas_index_all,2);
% 
% % max_vals=nan(4,size_atlas_index_all,2);
% % min_vals=nan(4,size_atlas_index_all,2);
% load('xylim_vals_all.mat')
% for val_i=1:size_atlas_index_all
%     ylim_max_all(val_i)=max(max(max_vals(:,val_i,:)));
%     ylim_min_all(val_i)=min(min(min_vals(:,val_i,:)));
% end
% 
% 
% if flip_flag == 0
%     task_name_fig={'Recalled','Not recalled','Recalled + Not recalled','Recalled - Not recalled'};
%     
%     for task_i=1:4
%         cnt=1;
%         figure('position',[0 0 2000 8000]);
%         set(gcf, 'color', [1 1 1]);
%         set(gcf,'Visible','on');
%         for atlas_index=atlas_index_all
%             max_val=max(max(max(max(squeeze(plot_test_elec(:,:,atlas_index,:,:))))));
%             min_val=min(min(min(min(squeeze(plot_test_elec(:,:,atlas_index,:,:))))));  
%             if ~isnan(max_val)
%                 subplot(4,2,cnt);
%                 plot_test_left=squeeze(plot_test_elec(1,task_i,atlas_index,:,:));
%                 plot_test_left(plot_test_left==0)=nan;
%                 [h1,output1]=plot_areaerrorbar(plot_test_left,1);
%                 max_vals(task_i,cnt,1)=output1.max_val;
%                 min_vals(task_i,cnt,1)=output1.min_val;
% 
%                 hold on;
%                 
%                 plot_test_right=squeeze(plot_test_elec(2,task_i,atlas_index,:,:));
%                 plot_test_right(plot_test_right==0)=nan;
%                 [h2,output2]=plot_areaerrorbar(plot_test_right,2);
%                 max_vals(task_i,cnt,2)=output2.max_val;
%                 min_vals(task_i,cnt,2)=output2.min_val;
%                 
% 
%                 if cnt == 7 | cnt == 8
%                     set(gca,'xtick',[1 20 40 60 80 100 118]);
%                     set(gca,'xticklabel',[-500 0 500 1000 1500 2000 2500]);
%                     set(gca,'xlim',[1 118])
%                     xlabel('Time (ms)');
%                 else
%                     set(gca,'xtick',[]);
%                     set(gca,'xticklabel',[]);
%                 end
%                 
%                 set(gca,'box','off')
% 
%                 set(gca,'ylim',[ylim_min_all(cnt) ylim_max_all(cnt)])
%                 set(gca,'Fontsize',15);
% 
%                 x1=xline(20,'-.r',{'ON'},'LabelVerticalAlignment','top','LabelHorizontalAlignment','center','LabelOrientation','horizontal','LineWidth',3,'Fontsize',15);
%                 x2=xline(84,'-.r',{'OFF'},'LabelVerticalAlignment','top','LabelHorizontalAlignment','center','LabelOrientation','horizontal','LineWidth',3,'Fontsize',15);
%                 y1=yline(0,'k','LineWidth',1,'Fontsize',15);
% 
%                 %         line([20,20],[-8,10],'Color','r','LineWidth',1)
%                 %         line([84,84],[-8,10],'Color','r','LineWidth',1)
%                 
%                 
%                 if cnt==1
%                     ylabel({'Gamma power z-score'});
%                     legend([h1 h2],{'Left','Right'},'LineWidth',1.5,'FontSize',10)
%                 end
%                 
%                 title(DK_atlas_names{atlas_index},'FontSize',20,'Color',font_color(atlas_index,:))
%         %         [ax,h1]=suplabel(['ROI : ' DK_atlas_names{atlas_index} ', Condition : ' task_name_fig{task_i} ],'t');
% %                 [ax,h3]=suplabel([DK_atlas_names{atlas_index} ],'t');
% %                 set(h3,'FontSize',20,'Color',font_color(atlas_index,:))
%             end
%             cnt=cnt+1;
%         end  
%         [ax,h3]=suplabel(task_name_fig{task_i},'t');
%         set(h3,'FontSize',30)
%         %% save figure
% %         saveas(gcf,[zscore_folder '\' task_name_fig{task_i} '_all_plot.png'])
% 
%         export_fig([zscore_folder '\' task_name_fig{task_i} '_all_plot'],'-tiff')
%         saveas(gcf,[zscore_folder '\' task_name_fig{task_i} '_all_plot'],'epsc')
%         close all;
%     end
% elseif flip_flag == 1
%     
%     
%     
%     
% end
% 
% % save('xylim_vals_all.mat','max_vals','min_vals')


%% plot all together (t-test)

% if flip_flag == 0
%     task_name_fig={'Recalled','Not recalled','Recalled + Not recalled','Recalled - Not recalled'};
%     
%     for task_i=1:4
%         cnt=1;
%         figure('position',[0 0 2000 8000]);
%         set(gcf, 'color', [1 1 1]);
%         set(gcf,'Visible','on');
%         for atlas_index=atlas_index_all
%             max_val=max(max(max(max(squeeze(plot_test_elec(:,:,atlas_index,:,:))))));
%             min_val=min(min(min(min(squeeze(plot_test_elec(:,:,atlas_index,:,:))))));  
%             if ~isnan(max_val)
%                 subplot(4,2,cnt);
%                 plot_test_left=squeeze(plot_test_elec(1,task_i,atlas_index,:,:));
%                 plot_test_left(plot_test_left==0)=nan;
%                 
%                 baseline_buf=[];
%                 baseline_buf=nanmean(plot_test_left(:,1:19),2);
%                 [H,P,CI,STATS]=ttest(plot_test_left-baseline_buf);
% %                 for time_i=1:118
% %                     stimul_buf=plot_test_left(:,time_i);
% %                     [H,P,CI,STATS]=ttest(baseline_buf,stimul_buf,'tail','both')
% %                 end
%                 h1=plot(STATS.tstat,'color',[52 148 186]./255,'LineWidth', 4);
%                 yline(tinv(0.975,STATS.df(1)),'color',[52 148 186]./255,'LineWidth',1);
%                 yline(tinv(0.025,STATS.df(1)),'color',[52 148 186]./255,'LineWidth',1);
%                 hold on;
% 
%                 plot_test_right=squeeze(plot_test_elec(2,task_i,atlas_index,:,:));
%                 plot_test_right(plot_test_right==0)=nan;
% 
%                 
%                 baseline_buf=[];
%                 baseline_buf=nanmean(plot_test_right(:,1:19),2);
%                 [H,P,CI,STATS]=ttest(plot_test_right-baseline_buf);
% %                 for time_i=1:118
% %                     stimul_buf=plot_test_left(:,time_i);
% %                     [H,P,CI,STATS]=ttest(baseline_buf,stimul_buf,'tail','both')
% %                 end
%                 h2=plot(STATS.tstat,'color',[236 112  22]./255,'LineWidth', 4);
%                 yline(tinv(0.975,STATS.df(1)),'color',[236 112  22]./255,'LineWidth',1);
%                 yline(tinv(0.025,STATS.df(1)),'color',[236 112  22]./255,'LineWidth',1);
%                 hold on;
%                 
% 
%                 if cnt == 7 | cnt == 8
%                     set(gca,'xtick',[1 20 40 60 80 100 118]);
%                     set(gca,'xticklabel',[-500 0 500 1000 1500 2000 2500]);
%                     set(gca,'xlim',[1 118])
%                     xlabel('Time (ms)');
%                 else
%                     set(gca,'xtick',[]);
%                     set(gca,'xticklabel',[]);
%                 end
%                 
%                 set(gca,'box','off')
% 
%                 set(gca,'ylim',[-6 6])
%                 set(gca,'Fontsize',15);
% 
%                 x1=xline(20,'-.r',{'ON'},'LabelVerticalAlignment','top','LabelHorizontalAlignment','center','LabelOrientation','horizontal','LineWidth',3,'Fontsize',15);
%                 x2=xline(84,'-.r',{'OFF'},'LabelVerticalAlignment','top','LabelHorizontalAlignment','center','LabelOrientation','horizontal','LineWidth',3,'Fontsize',15);
%                 y1=yline(0,'k','LineWidth',1,'Fontsize',15);
% 
%                 %         line([20,20],[-8,10],'Color','r','LineWidth',1)
%                 %         line([84,84],[-8,10],'Color','r','LineWidth',1)
%                 
%                 
%                 if cnt==1
%                     ylabel({'t value'});
%                     legend([h1 h2],{'Left','Right'},'LineWidth',1.5,'FontSize',10)
%                 end
%                 
%                 title(DK_atlas_names{atlas_index},'FontSize',20,'Color',font_color(atlas_index,:))
%         %         [ax,h1]=suplabel(['ROI : ' DK_atlas_names{atlas_index} ', Condition : ' task_name_fig{task_i} ],'t');
% %                 [ax,h3]=suplabel([DK_atlas_names{atlas_index} ],'t');
% %                 set(h3,'FontSize',20,'Color',font_color(atlas_index,:))
%             end
%             cnt=cnt+1;
%         end  
%         [ax,h3]=suplabel(task_name_fig{task_i},'t');
%         set(h3,'FontSize',30)
%         %% save figure
% %         saveas(gcf,[zscore_folder '\' task_name_fig{task_i} '_all_plot.png'])
% 
%         export_fig([zscore_folder '\' task_name_fig{task_i} '_all_plot'],'-tiff')
%         saveas(gcf,[zscore_folder '\' task_name_fig{task_i} '_all_plot'],'epsc')
%         close all;
%     end
% elseif flip_flag == 1
%     
%     
%     
%     
% end
% 
% % save('xylim_vals_all.mat','max_vals','min_vals')

