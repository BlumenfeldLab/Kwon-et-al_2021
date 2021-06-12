
clc;clear;
% file_suffix='Final_std_20_20_fs_together';
% Montage_suffix='final_roi_fs';
file_suffix='Final_std_20_20_fs_together_xflip';
Montage_suffix='final_roi_fs_xflip';

data_folder=['E:\RAM data set\RAM_Public_Data_all\FR1_FS\Storage_filtering_overlap\Data_storage\' file_suffix '\' Montage_suffix];
cd(data_folder)

data_folder_new=['E:\RAM data set\RAM_Public_Data_all\FR1_FS\Storage_filtering_overlap\Data_storage\' file_suffix '_FDR\' Montage_suffix];
mkdir(data_folder_new)


%% baseline
start_side = 1;
sides_num = 2;
total_num_frames_base=19;

for side_index = start_side:sides_num    % Generate frames for L and R
    % Determine what side we are working on
    if side_index == 1
        side = 'L';
    else
        side = 'R';
    end
    vertex_values_frame_base=zeros(163842,251);  
    for frame_index = 1:total_num_frames_base
        clearvars vertex_values_frame
        load([side '_vertex_values_' file_suffix '_' num2str(frame_index) '.mat']);
        vertex_values_frame_base=vertex_values_frame_base+vertex_values_frame;
    end
    vertex_values_frame_base=vertex_values_frame_base/total_num_frames_base;
    save([data_folder_new '/' side '_vertex_values_' file_suffix '_base.mat'],'vertex_values_frame_base');
    
end

%% ttest
start_side = 1;
sides_num = 2;
total_num_frames=98;

for side_index = start_side:sides_num    % Generate frames for L and R
    % Determine what side we are working on
    if side_index == 1
        side = 'L';
    else
        side = 'R';
    end
    load([data_folder_new '/' side '_vertex_values_' file_suffix '_base.mat']);
    vertex_values_frame_base(find(vertex_values_frame_base==0))=nan;
    for frame_index = 1:total_num_frames
        clearvars vertex_values_frame_p vertex_values_frame_t
        vertex_values_frame_p=nan(163842,1);
        vertex_values_frame_t=nan(163842,1);
        vertex_values_frame_sub=nan(163842,1);
%         vertex_values_frame_rank_p=nan(163842,1);
%         vertex_values_frame_rank_t=nan(163842,1);
        
        load([side '_vertex_values_' file_suffix '_' num2str(frame_index) '.mat']);
        vertex_values_frame(find(vertex_values_frame==0))=nan;
        
        vertex_values_frame_new=vertex_values_frame-vertex_values_frame_base;
        
        for test_i=1:163842
            vertex_values_frame_sub(test_i,1)=length(find(~isnan(vertex_values_frame_new(test_i,:))>0));
        end
        
        % ttest
        [H,P,CI,STATS]=ttest(vertex_values_frame_new');
        vertex_values_frame_t=STATS.tstat';
        vertex_values_frame_p=P';
        
        vertex_values_frame_p(find(vertex_values_frame_sub<5),1)=nan;
        
        nonindex=isnan(vertex_values_frame_p);
        vertex_values_frame_t(nonindex)=0;
        vertex_values_frame_p(nonindex)=1;
        
        fdr(vertex_values_frame_p,0.05)
        
        save([data_folder_new '/' side '_vertex_values_' file_suffix '_ttest_' num2str(frame_index) '.mat'],'vertex_values_frame_t','vertex_values_frame_p');
        
        % signed ranksum test
%         signrank(vertex_values_frame_base',vertex_values_frame')

        frame_index
    end
    
end


%% FDR
start_side = 1;
sides_num = 2;
total_num_frames=98;

for side_index = start_side:sides_num    % Generate frames for L and R
    % Determine what side we are working on
    if side_index == 1
        side = 'L';
    else
        side = 'R';
    end
    vertex_values_frame_p_all=nan(163842,98);
    for frame_index = 1:total_num_frames
        clearvars vertex_values_frame_p vertex_values_frame_t
        load([data_folder_new '/' side '_vertex_values_' file_suffix '_ttest_' num2str(frame_index) '.mat']);
        
        vertex_values_frame_p_all(:,frame_index)=vertex_values_frame_p;
        
        frame_index
    end
    
    if side_index == 1
        vertex_values_frame_p_all_L=vertex_values_frame_p_all;
    else
        vertex_values_frame_p_all_R=vertex_values_frame_p_all;
    end
    
end

fdr([vertex_values_frame_p_all_L(:) ;vertex_values_frame_p_all_R(:)],0.05)




%% save data
start_side = 1;
sides_num = 2;
total_num_frames=98;
mean_vertex_values=nan(2,251,98);
for side_index = start_side:sides_num    % Generate frames for L and R
    % Determine what side we are working on
    if side_index == 1
        side = 'L';
    else
        side = 'R';
    end

    for frame_index = 1:total_num_frames
        clearvars vertex_values_frame

        load([side '_vertex_values_' file_suffix '_' num2str(frame_index) '.mat']);
        vertex_values_frame(find(vertex_values_frame==0))=nan;
        
        mean_vertex_values(side_index,:,frame_index)=nanmean(vertex_values_frame);
        
        frame_index
    end
    
end

load('cmap_hunki_new');
figure('position',[0 0 8000 8000]);
set(gcf, 'color', [1 1 1]);
set(gcf,'Visible','on');

subplot(2,1,1)
buffer_Left=squeeze(mean_vertex_values(1,:,:));
% buffer_Left(isnan(buffer_Left(:,1)),:)=[];
buffer_Left(isnan(buffer_Left(:,1)),:)=0;
imagesc(buffer_Left)
title(['Left Hemisphere'])
set(gca,'ydir','normal');
set(gca,'xtick',[1 20 40 60 80]);
set(gca,'xticklabel',[-500 0 500 1000 1500]);
set(gca,'xlim',[1 98])
set(gca,'Fontsize',15);
line([20,20],ylim,'Color','r','LineWidth',0.5)
ylabel('Subject');
xlabel('Time (ms)');
colormap(mycmap);
colorbar
caxis([-2 2])

subplot(2,1,2)
buffer_Right=squeeze(mean_vertex_values(2,:,:));
% buffer_Right(isnan(buffer_Right(:,1)),:)=[];
buffer_Right(isnan(buffer_Right(:,1)),:)=0;
imagesc(buffer_Right)
title(['Right Hemisphere'])
set(gca,'ydir','normal');
set(gca,'xtick',[1 20 40 60 80]);
set(gca,'xticklabel',[-500 0 500 1000 1500]);
set(gca,'xlim',[1 98])
set(gca,'Fontsize',15);
line([20,20],ylim,'Color','r','LineWidth',0.5)
ylabel('Subject');
xlabel('Time (ms)');
colormap(mycmap);
colorbar
caxis([-2 2])






%% test 1 (vertex)
clc;clear;

start_side = 1;
sides_num = 2;
total_num_frames=98;
mean_vertex_values=nan(2,251,98);

% file_suffix='Final_std_20_20_fs_together';
% Montage_suffix='final_roi_fs';
file_suffix='Final_std_20_20_fs_together_xflip';
% Montage_suffix='final_roi_fs_xflip';

for side_index = start_side:sides_num    % Generate frames for L and R
    % Determine what side we are working on
    if side_index == 1
        side = 'L';
        [vertices, label, colortable]=read_annotation('lh.Schaefer2018_400Parcels_7Networks_order.annot');
    else
        side = 'R';
        [vertices, label, colortable]=read_annotation('lh.Schaefer2018_400Parcels_7Networks_order.annot');
    end

    for frame_index = 1:total_num_frames
        clearvars vertex_values_frame

        load([side '_vertex_values_' file_suffix '_' num2str(frame_index) '.mat']);
        vertex_values_frame(find(vertex_values_frame==0))=nan;
        
        mean_vertex_values(side_index,:,frame_index)=nanmean(vertex_values_frame);
        
        frame_index
    end
    
end

load('cmap_hunki_new');
figure('position',[0 0 8000 8000]);
set(gcf, 'color', [1 1 1]);
set(gcf,'Visible','on');

% subplot(2,1,1)
buffer_Left=squeeze(mean_vertex_values(1,:,:));
% buffer_Left(isnan(buffer_Left(:,1)),:)=[];
buffer_Left(isnan(buffer_Left(:,1)),:)=0;
imagesc(buffer_Left)
title(['Flip Hemisphere'])
set(gca,'ydir','normal');
set(gca,'xtick',[1 20 40 60 80]);
set(gca,'xticklabel',[-500 0 500 1000 1500]);
set(gca,'xlim',[1 98])
set(gca,'Fontsize',15);
line([20,20],ylim,'Color','r','LineWidth',0.5)
ylabel('Subject');
xlabel('Time (ms)');
colormap(mycmap);
colorbar
caxis([-2 2])

% subplot(2,1,2)
buffer_Right=squeeze(mean_vertex_values(2,:,:));
% buffer_Right(isnan(buffer_Right(:,1)),:)=[];
buffer_Right(isnan(buffer_Right(:,1)),:)=0;
imagesc(buffer_Right)
title(['Right Hemisphere'])
set(gca,'ydir','normal');
set(gca,'xtick',[1 20 40 60 80]);
set(gca,'xticklabel',[-500 0 500 1000 1500]);
set(gca,'xlim',[1 98])
set(gca,'Fontsize',15);
line([20,20],ylim,'Color','r','LineWidth',0.5)
ylabel('Subject');
xlabel('Time (ms)');
colormap(mycmap);
colorbar
caxis([-2 2])

figure; plot(squeeze(mean_vertex_values(2,30,:)))



%% test 1 (elec)
clc;clear;

start_side = 1;
sides_num = 2;
total_num_frames=98;

% file_suffix='Final_std_20_20_fs_together';
% Montage_suffix='final_roi_fs';
file_suffix='Final_std_20_20_fs_together_xflip';
% Montage_suffix='final_roi_fs_xflip';

mean_vertex_values=nan(2,251,98);
for side_index = start_side:sides_num    % Generate frames for L and R
    
    if side_index == 1
        side = 'L';
%         [vertices, label, colortable]=read_annotation('lh.Schaefer2018_400Parcels_7Networks_order.annot');
    else
        side = 'R';
%         [vertices, label, colortable]=read_annotation('lh.Schaefer2018_400Parcels_7Networks_order.annot');
    end

    for frame_index = 1:total_num_frames
        clearvars vertex_values_frame_e

        load([side '_vertex_values_' file_suffix '_' num2str(frame_index) '.mat']);
        vertex_values_frame_e(find(vertex_values_frame_e==0))=nan;
        
        mean_vertex_values(side_index,:,frame_index)=nanmean(vertex_values_frame_e);
        
        frame_index
    end
    
end

load('cmap_hunki_new');
figure('position',[0 0 8000 8000]);
set(gcf, 'color', [1 1 1]);
set(gcf,'Visible','on');

buffer_Left=squeeze(mean_vertex_values(1,:,:));
% buffer_Left(isnan(buffer_Left(:,1)),:)=[];
buffer_Left(isnan(buffer_Left(:,1)),:)=0;
imagesc(buffer_Left)
title(['Flip Hemisphere'])
set(gca,'ydir','normal');
set(gca,'xtick',[1 20 40 60 80]);
set(gca,'xticklabel',[-500 0 500 1000 1500]);
set(gca,'xlim',[1 98])
set(gca,'Fontsize',15);
line([20,20],ylim,'Color','r','LineWidth',0.5)
ylabel('Subject');
xlabel('Time (ms)');
colormap(mycmap);
colorbar
caxis([-2 2])

figure; plot(squeeze(mean_vertex_values(1,60,:)))



%% test 2 (elec, stats)
clc;clear;

start_side = 1;
sides_num = 2;
total_num_frames=98;
total_num_frames_base=19;
load fs_DK_atlas.mat

% file_suffix='Final_std_20_20_fs_together';
% Montage_suffix='final_roi_fs';
file_suffix='Final_std_20_20_fs_together_xflip';
Montage_suffix='final_roi_fs_xflip';

mean_elec_values_t=zeros(2,36,total_num_frames);
mean_elec_values_p=zeros(2,36,total_num_frames);
mean_elec_values_df=zeros(2,36,total_num_frames);
mean_elec_values_z=zeros(2,36,total_num_frames);
mean_elec_values_base_all=nan(2,36,total_num_frames_base,251);
for side_index = start_side:1    % Generate frames for L and R
    
    if side_index == 1
        side = 'L';
%         [vertices, label, colortable]=read_annotation('lh.Schaefer2018_400Parcels_7Networks_order.annot');
        [vertices, label, colortable] = read_annotation('lh.aparc.annot');
    else
        side = 'R';
%         [vertices, label, colortable]=read_annotation('lh.Schaefer2018_400Parcels_7Networks_order.annot');
        [vertices, label, colortable] = read_annotation('rh.aparc.annot');
    end
    
    for frame_index = 1:total_num_frames_base
%         clearvars vertex_values_frame_e
        load([side '_vertex_values_' file_suffix '_' num2str(frame_index) '.mat']);
        vertex_values_frame_e(find(vertex_values_frame_e==0))=nan;
        for atlas_i=2:size(colortable.table,1)
            vertex_values_frame_e_buff=[];
            vertex_values_frame_e_buff=vertex_values_frame_e(label == colortable.table(atlas_i,5),:);
            
            mean_vertex_values_frame_e_buff=[];
            mean_vertex_values_frame_e_buff=nanmean(vertex_values_frame_e_buff);
            sub_n=length(find(~isnan(mean_vertex_values_frame_e_buff)>0));
            if sub_n > 5
                mean_elec_values_base_all(side_index,atlas_i,frame_index,:)=nanmean(mean_vertex_values_frame_e_buff,1);
            end
        end
    end
    mean_elec_values_base=squeeze(nanmean(mean_elec_values_base_all,3));
    
    for frame_index = 1:total_num_frames
        clearvars vertex_values_frame_e
        load([side '_vertex_values_' file_suffix '_' num2str(frame_index) '.mat']);
        vertex_values_frame_e(find(vertex_values_frame_e==0))=nan;
        for atlas_i=2:size(colortable.table,1)
            vertex_values_frame_e_buff=[];
            vertex_values_frame_e_buff=vertex_values_frame_e(label == colortable.table(atlas_i,5),:);
            
            mean_vertex_values_frame_e_buff=[];
            mean_vertex_values_frame_e_buff=nanmean(vertex_values_frame_e_buff);
            sub_n=length(find(~isnan(mean_vertex_values_frame_e_buff)>0));
            if sub_n > 5
                mean_elec_values_z(side_index,atlas_i,frame_index)=nansum(mean_vertex_values_frame_e_buff)/sqrt(sub_n);
                
                [H,P,CI,STATS]=ttest(mean_vertex_values_frame_e_buff-squeeze(mean_elec_values_base(side_index,atlas_i,:))');
                mean_elec_values_p(side_index,atlas_i,frame_index)=P;
                mean_elec_values_t(side_index,atlas_i,frame_index)=STATS.tstat;
                mean_elec_values_df(side_index,atlas_i,frame_index)=STATS.df;
            end
        end
        frame_index
    end
    
end


% plot 
figure('position',[0 0 8000 8000]);
set(gcf, 'color', [1 1 1]);
set(gcf,'Visible','on');

subplot(1,2,1)
imagesc(squeeze(mean_elec_values_z(1,:,:)))

title(['Square root mean z-score'])
set(gca,'ydir','normal');
set(gca,'ytick',[1:36]);
set(gca,'yticklabel',DK_atlas_names');
set(gca,'xtick',[1 20 40 60 80]);
set(gca,'xticklabel',[-500 0 500 1000 1500]);
set(gca,'xlim',[1 98])
set(gca,'Fontsize',15);
line([20,20],ylim,'Color','r','LineWidth',0.5)
ylabel('ROI index');
xlabel('Time (ms)');
load('cmap_hunki_10to15.mat');
colormap(mycmap);
colorbar
caxis([-10 10])

subplot(1,2,2)
imagesc(squeeze(mean_elec_values_t(1,:,:)))

title(['T-value'])
set(gca,'ydir','normal');
set(gca,'ytick',[1:36]);
set(gca,'yticklabel',DK_atlas_names');
set(gca,'xtick',[1 20 40 60 80]);
set(gca,'xticklabel',[-500 0 500 1000 1500]);
set(gca,'xlim',[1 98])
set(gca,'Fontsize',15);
line([20,20],ylim,'Color','r','LineWidth',0.5)
ylabel('ROI index');
xlabel('Time (ms)');
load('cmap_hunki_new.mat');
colormap(mycmap);
colorbar
caxis([-5 5])



%% test 4 (elec, stats, each ROI)
clc;clear;
close all;

rootfolder='E:/RAM data set/RAM_Public_Data_all';
cd(rootfolder)

load r1_all.mat
load fs_DK_atlas.mat
fid=fopen('Subjects_list_all.txt','r');
cnt=1;
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);

total_num_frames=98;
total_num_frames_base=19;
load fs_DK_atlas.mat

% file_suffix='Final_std_20_20_fs_together';
% Montage_suffix='final_roi_fs';
file_suffix='Final_std_20_20_fs_together_xflip';
Montage_suffix='final_roi_fs_xflip';

for subject_i=1:251
    try
        data_location='E:/RAM data set/RAM_Public_Data_all/FR1_FS';
        cd(data_location)

        cd([num2str(subject_i) '_' r_sublist{subject_i,1}]);

        start_side = 1;
        sides_num = 1;

        mean_elec_values_z=zeros(2,36,total_num_frames);
        for side_index = start_side:sides_num    % Generate frames for L and R

            if side_index == 1
                side = 'L';
        %         [vertices, label, colortable]=read_annotation('lh.Schaefer2018_400Parcels_7Networks_order.annot');
                [vertices, label, colortable] = read_annotation('lh.aparc.annot');
            else
                side = 'R';
        %         [vertices, label, colortable]=read_annotation('lh.Schaefer2018_400Parcels_7Networks_order.annot');
                [vertices, label, colortable] = read_annotation('rh.aparc.annot');
            end

            clearvars vertex_values_frame_e
            load([side '_vertex_values_' file_suffix '.mat']);
            vertex_values_e(find(vertex_values_e==0))=nan;
            for atlas_i=2:size(colortable.table,1)
                vertex_values_e_buff=[];
                vertex_values_e_buff=vertex_values_e(label == colortable.table(atlas_i,5),:);
                mean_vertex_values_frame_e_buff=[];
                elec_n=length(find(~isnan(vertex_values_e_buff(:,1))>0));
                if elec_n > 0
%                     mean_vertex_values_frame_e_buff=nansum(vertex_values_e_buff); % mean zscore across electrodes
%                     mean_elec_values_z(side_index,atlas_i,:)=mean_vertex_values_frame_e_buff./sqrt(elec_n);
                    mean_elec_values_z(side_index,atlas_i,:)=nanmean(vertex_values_e_buff);
                end
            end
            

        end
        
        figure('position',[0 0 8000 8000]);
        set(gcf, 'color', [1 1 1]);
        set(gcf,'Visible','off');

        imagesc(squeeze(mean_elec_values_z(1,:,:)))

        title(['Square root mean z-score, ' num2str(subject_i) '_' r_sublist{subject_i,1}]) % subject name
        set(gca,'ydir','normal');
        set(gca,'ytick',[1:36]);
        set(gca,'yticklabel',DK_atlas_names');
        set(gca,'xtick',[1 20 40 60 80]);
        set(gca,'xticklabel',[-500 0 500 1000 1500]);
        set(gca,'xlim',[1 98])
        set(gca,'Fontsize',15);
        line([20,20],ylim,'Color','r','LineWidth',0.5)
        ylabel('ROI index');
        xlabel('Time (ms)');
        load('cmap_hunki_10to15.mat');
        colormap(mycmap);
        colorbar
        caxis([-10 15])
%         caxis([-5 7.5])
        saveas(gcf,[data_location '/all_zscore_ROI_plot/subject/' num2str(subject_i) '_' r_sublist{subject_i,1} '.png'])
        close all;
        
        disp(['Completed : ' num2str(subject_i) '_' r_sublist{subject_i,1}]);
    catch
        disp(['Not completed : ' num2str(subject_i) '_' r_sublist{subject_i,1}]);
    end
end



%% test 5 (elec, stats, each roi)
clc;clear;
close all;

rootfolder='E:/RAM data set/RAM_Public_Data_all';
cd(rootfolder)

load r1_all.mat
load fs_DK_atlas.mat
fid=fopen('Subjects_list_all.txt','r');
cnt=1;
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);

total_num_frames=98;
total_num_frames_base=19;
load fs_DK_atlas.mat

% file_suffix='Final_std_20_20_fs_together';
% Montage_suffix='final_roi_fs';
file_suffix='Final_std_20_20_fs_together_xflip';
Montage_suffix='final_roi_fs_xflip';
for atlas_i=2:size(DK_atlas_names,1)
    mean_elec_values_z=[];
    sub_index_all=[];
    for subject_i=1:251
        try
            data_location='E:/RAM data set/RAM_Public_Data_all/FR1_FS';
            cd(data_location)

            cd([num2str(subject_i) '_' r_sublist{subject_i,1}]);

            start_side = 1;
            sides_num = 1;

            
            for side_index = start_side:sides_num    % Generate frames for L and R

                if side_index == 1
                    side = 'L';
            %         [vertices, label, colortable]=read_annotation('lh.Schaefer2018_400Parcels_7Networks_order.annot');
                    [vertices, label, colortable] = read_annotation('lh.aparc.annot');
                else
                    side = 'R';
            %         [vertices, label, colortable]=read_annotation('lh.Schaefer2018_400Parcels_7Networks_order.annot');
                    [vertices, label, colortable] = read_annotation('rh.aparc.annot');
                end

                clearvars vertex_values_frame_e
                load([side '_vertex_values_' file_suffix '.mat']);
                vertex_values_e(find(vertex_values_e==0))=nan;


                vertex_values_e_buff=[];
                vertex_values_e_buff=vertex_values_e(label == colortable.table(atlas_i,5),:);
                mean_vertex_values_frame_e_buff=[];
                elec_n=length(find(~isnan(vertex_values_e_buff(:,1))>0));
                if elec_n > 0
%                     mean_vertex_values_frame_e_buff=nansum(vertex_values_e_buff); % mean zscore across electrodes
%                     mean_elec_values_z = [mean_elec_values_z ; mean_vertex_values_frame_e_buff./sqrt(elec_n)];
                    mean_vertex_values_frame_e_buff=nanmean(vertex_values_e_buff); % mean zscore across electrodes
                    mean_elec_values_z = [mean_elec_values_z ; mean_vertex_values_frame_e_buff];
                    sub_index_all=[sub_index_all subject_i];
                end


            end

            disp(['Completed : ' num2str(subject_i) '_' r_sublist{subject_i,1}]);
        catch
            disp(['Not completed : ' num2str(subject_i) '_' r_sublist{subject_i,1}]);
        end
    end
    
    figure('position',[0 0 8000 8000]);
    set(gcf, 'color', [1 1 1]);
    set(gcf,'Visible','on');

    imagesc(mean_elec_values_z)

    title(['ROI : ' DK_atlas_names{atlas_i} ', Subjects : ' num2str(size(mean_elec_values_z,1))]) % subject name
    set(gca,'ydir','normal');
    set(gca,'ytick',1:size(sub_index_all,2));
    set(gca,'yticklabel',sub_index_all);
    set(gca,'xtick',[1 20 40 60 80]);
    set(gca,'xticklabel',[-500 0 500 1000 1500]);
    set(gca,'xlim',[1 98])
    set(gca,'Fontsize',15);
    line([20,20],ylim,'Color','r','LineWidth',0.5)
    ylabel('Subject index');
    xlabel('Time (ms)');
    load('cmap_hunki_10to15.mat');
    colormap(mycmap);
    colorbar
    caxis([-10 15])
%     caxis([-5 7.5])
    saveas(gcf,[data_location '/all_zscore_ROI_plot/roi/' num2str(atlas_i) '_' DK_atlas_names{atlas_i} '.png'])
    close all;
    
end

