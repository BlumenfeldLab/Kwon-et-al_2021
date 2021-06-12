% Arrival Times Script

% Modified by Hunki 12/11/2018 for RAM data set


function displayElectrodesInflated_new_fs_stat_cluster_zscore(data_location,patients,overlay, figure_type, views, thr,color_range,colormap_mat,file_suffix,storage_suffix,Montage_suffix,atlas_name,atlas_folder,atlas_index,flip_prefix,time_thr,all_sub_index,mean_flag)
r_time=(time_thr+1)*25;
image_storage_folder = [data_location '/Storage_filtering_overlap_' storage_suffix '/Image_storage/' file_suffix '/' Montage_suffix '/' atlas_index '_' atlas_name '_' num2str(r_time) 'ms_zscore'];
data_storage_folder = [data_location '/Storage_filtering_overlap_' storage_suffix '/Data_storage/' file_suffix '/' Montage_suffix ];

xtick_spacing_filename = 'time_50ms_overlap.mat';
mni_nifti_path = 'MNI_T1_1mm_stripped.nii';

% fsaverage_dir = 'E:/RAM data set/RAM_Public_Data_all/fsaverage';

% mni_nifti_path = 'MNI_T1_1mm_stripped.nii';
% mni2fs_dir = 'Functions for Inflated Brain Display/mni2fs-master';

load(xtick_spacing_filename);

mkdir(image_storage_folder)
% mkdir(data_storage_folder)

mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
Tmni = mnit1_1mm_stripped.transform;
alpha = 0.7;

%% Side figure_type - which side of the brain to draw
start_side = 1;
sides_num = 2;

figHandle = gobjects(1, sides_num);

load('FS_ori_avg_surface.mat')
load('FS_avg_surface.mat')
%% Generating all frames for each side of the brain
for side_index = start_side:sides_num    % Generate frames for L and R

    % Determine what side we are working on
    if side_index == 1
        side = 'L';
        hem = 'lh';
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
        hem = 'rh';
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
             
    s.tri=inflated_surf(side_index).faces;
    s.coord=inflated_surf(side_index).vertices';
    s.coord=double(s.coord);
%     [nbr, deg] = findnbr(s.tri);
    
    %% hard coded number of frames for each overlay setting, need to change later
    
    if overlay == 1 || overlay == 0 % arrival times
        total_frames = 1;
    elseif overlay == 2 || overlay == 3 || overlay == 4%
%         total_frames=[20 22 24 26 28 30 9 20 32 44 56 68 80 92];
        total_frames = [1:length(T)];
    end
    
    %------------------------------------------------------------------
    %   Set-up some local variables
    %------------------------------------------------------------------
    
%     s = SurfStatReadSurf( {...
%     ['E:/RAM data set/RAM_Public_Data_all/fsaverage/fsaverage/surf/' hem '.orig_avg']} );
%     s_buff.faces=s.tri;
%     s_buff.vertices=s.coord';
% %     surf_fn = fullfile(fsaverage_dir,['/inflated_' hem '.gii']);
% %     fs_surf(side_index) = export(gifti(surf_fn));
%     fs_surf(side_index)=s_buff;
    
    v_num = size(fs_surf(side_index).vertices, 1);
    
    fs_surf(side_index).vertices = [fs_surf(side_index).vertices, ones(v_num, 1)] / Tmni';
    fs_surf(side_index).vertices = fs_surf(side_index).vertices(:, 1:3);
    
%     ver_test=fs_surf(side_index).vertices;
%     test=fmris_read_nifti('MNI_T1_1mm_graywhite.nii');
%     fmris_write_image(test)
%     test=fmris_read_nifti('fs_surf_ori.nii');
% %     new_mat=zeros(181,217,181);
% %     for test_i=1:size(ver_test,1)
% %         ver_x=round(ver_test(test_i,1));
% %         ver_y=round(ver_test(test_i,2));
% %         ver_z=round(ver_test(test_i,3));
% %         new_mat(ver_x,ver_y,ver_z)=1;
% %     end
% %     test.data=new_mat;
%     test.file_name='fs_surf_ori.nii';
%     test.byte=4;
%     fmris_write_nifti(test);

    inflated_surf(side_index).vertices = [inflated_surf(side_index).vertices, ones(v_num, 1)] / Tmni';
    inflated_surf(side_index).vertices = inflated_surf(side_index).vertices(:, 1:3);
    

    %% create a figure of the left or right side of the brain
    % you will be working with this for different electrodes, overlays,
    % views, and light
    
    figHandle(side_index) = figure('Position',[70,70,700,700]); % Original = [50,50,600,600]Creates a figure graphic object
    %     g=trisurf(cortex.faces,cortex.vertices(:,1),cortex.vertices(:,2),cortex.vertices(:,3));
    
    temp_surf = [];
    temp_surf.hem = hem; % choose the hemesphere 'lh' or 'rh'
    temp_surf.inflationstep = 1; % 1 no inflation, 6 fully inflated
    temp_surf.decimation = 0;
    temp_surf = mni2fs_brain(temp_surf);
    set(temp_surf.p, 'Faces', inflated_surf(side_index).faces, 'Vertices', inflated_surf(side_index).vertices)
    
    set(gcf, 'color', [1 1 1]);
    hold on

    axis off;
    axis equal;

    if ~ismember(overlay,1)
        %------------------------------------------------------------------
        %   Frame Generation with overlaying color onto the brain surface
        %------------------------------------------------------------------
%         load('cmap6to8WRX.mat');
%         load('cmap_hunki_new');
%         load('cmap_hunki_5to15.mat');       
%         load('cmap_hunki_10to15.mat'); 
        load(colormap_mat);
        
        cd(atlas_folder)
        load([atlas_name '_' atlas_index '_' file_suffix '_' flip_prefix '.mat']);
        
        cluster_based_sig_roi=zeros(size(colortable.table,1),size(T,2));
        pos_clust_ids_buff=[];
        neg_clust_ids_buff=[];
        t_orig_buff=[];
        
        if side_index==1
            pos_clust_ids_buff = clust_info.pos_clust_ids(1:size(colortable.table,1),:);
            neg_clust_ids_buff = clust_info.neg_clust_ids(1:size(colortable.table,1),:);
            t_orig_buff = t_orig(1:size(colortable.table,1),:);
        else
            pos_clust_ids_buff = clust_info.pos_clust_ids(size(colortable.table,1)+1:end,:);
            neg_clust_ids_buff = clust_info.neg_clust_ids(size(colortable.table,1)+1:end,:);
            t_orig_buff = t_orig(size(colortable.table,1)+1:end,:);
        end
        
        pos_cluster_index=find(clust_info.pos_clust_pval<est_alpha);
        if ~isempty(pos_cluster_index)
            for c_i=1:size(pos_cluster_index,2)
                buff_index=[];
                buff_index=pos_clust_ids_buff==pos_cluster_index(c_i);
                [a_roi b_time]=find(buff_index);
                if length(unique(b_time)) > time_thr
                    cluster_based_sig_roi(buff_index)=t_orig_buff(buff_index);
                end
            end
        end
        
        neg_cluster_index=find(clust_info.neg_clust_pval<est_alpha);
        if ~isempty(neg_cluster_index)
            for c_i=1:size(neg_cluster_index,2)
                buff_index=[];
                buff_index=neg_clust_ids_buff==neg_cluster_index(c_i);
                [a_roi b_time]=find(buff_index);
                if length(unique(b_time)) > time_thr
                    cluster_based_sig_roi(buff_index)=t_orig_buff(buff_index);
                end
            end
        end
        
%         test=nansum(logical(cluster_based_sig_roi),2);
%         
%         if side_index==1
% 
%         else
% 
%         end
        
        cd(image_storage_folder)
        for frame_index = total_frames
            
            vertex_values_sum=zeros(v_num,1);
            disp(['>> Rendering ' side ' Cortex, Frame ' num2str(frame_index)]);  % Output Status
            
            vertex_values_frame=zeros(v_num,length(all_sub_index));
            vertex_values_frame_masked=zeros(v_num,length(all_sub_index));
            cd(data_storage_folder)
            load([side '_vertex_values_' file_suffix '_' num2str(frame_index) '.mat']); 
            vertex_values_frame=vertex_values_frame(:,all_sub_index);
            
            
            %% sum the count of vertices across patients
            if (overlay == 2 || overlay == 3 || overlay == 4) % the surface is colored according to zscores of power
                cluster_all_index=find(cluster_based_sig_roi(:,frame_index)~=0);
                if ~isempty(cluster_all_index)
                    empty_flag=0;
                    for roi_i=1:size(cluster_all_index,1)
                        vertex_values_frame_masked(label == colortable.table(cluster_all_index(roi_i),5),:) ...
                            = vertex_values_frame(label == colortable.table(cluster_all_index(roi_i),5),:);
                    end  
                    
                    if mean_flag == 1
                        vertex_values_sum = squeeze(nansum(vertex_values_frame_masked, 2));
                        divisor = squeeze(sum(vertex_values_frame_masked ~= 0, 2));
%                         divisor(divisor == 0) = 1;

%                         % subject threshold
%                         divisor(divisor <5 ) = 1;
%                         vertex_values_sum(divisor<5)=0;

                        vertex_values_sum = vertex_values_sum ./ divisor;
%                         vertex_values_sum(~subject_mask)=0;


                    elseif mean_flag ==2
                        vertex_values_sum = squeeze(nansum(vertex_values_frame_masked, 2));
                        divisor = squeeze(sum(vertex_values_frame_masked ~= 0, 2));
%                         divisor(divisor == 0) = 1;

%                         % subject threshold
%                         divisor(divisor <5 ) = 1;
%                         vertex_values_sum(divisor<5)=0;

                        vertex_values_sum = vertex_values_sum ./ sqrt(divisor);
%                         vertex_values_sum(~subject_mask)=0;

                    end
                else
                    empty_flag=1;
                end
                
                % load the colormap and set axis of colormap
                colormap(mycmap)
                caxis(color_range);
                
%                 vertex_values_sum = surface_smooth_hk_ini(vertex_values_sum,s,2,20,nbr,deg);
                
                electrodes_present(1,1)=1;
                electrodes_present(2,1)=1;
                
%                 if side_index ==2
%                     vertex_values_sum=zeros(v_num,1);
%                 end
                
                if any(electrodes_present(side_index,:) == 1)
                    if empty_flag == 1
                        h = patch('Faces',inflated_surf(side_index).faces, 'Vertices', inflated_surf(side_index).vertices,...
                            'FaceVertexCData',vertex_values_sum(:,1),'FaceColor','interp',...
                            'LineStyle', 'none', 'FaceAlpha', 'flat');
                        alpha_data = zeros(numel(h.FaceVertexCData), 1);
                        lightset = [0.6 0.5 0.1];
                        material(lightset);
                        set(h, 'FaceVertexAlphaData', alpha_data)
                    else
                        h = patch('Faces',inflated_surf(side_index).faces, 'Vertices', inflated_surf(side_index).vertices,...
                            'FaceVertexCData',vertex_values_sum(:,1),'FaceColor','interp',...
                            'LineStyle', 'none', 'FaceAlpha', alpha);
                        lightset = [0.6 0.5 0.1];
                        material(lightset);
                    end
                end

                %             elseif overlay == 0 % the surface is colored according to overlap of space around the electrodes
                %                 vertex_values_sum = NaN(cortex_v_num(side_index),1);
                %                 for i = 1:size(vertex_values,1)
                %                     vertex_values_sum(i,1) = sum(vertex_values(i,:));
                %                 end
                %                 colormap jet
                %                 caxis([-1 7])
            else
                lightset = [0.6 0.5 0.1];
                material(lightset);
            end

            %% Save each side of the brain in 3 views
            frames_set(side_index,frame_index,:) = savebrainimages(side_index,frame_index,views);

            if any(electrodes_present(side_index,:) == 1) && exist('h','var')
                % get rid of the color patch so you can lay on the color patch for the
                % next frame
                delete(h)
            end
%             vertex_values_sum_all(:,frame_index)=vertex_values_sum;
        end
    end
    
%     cd(data_storage_folder)
%     if side_index == 1
%         save(['L_vertex_values_' file_suffix '_all.mat'],'vertex_values_sum_all','-v7.3')
%     else
%         save(['R_vertex_values_' file_suffix '_all.mat'],'vertex_values_sum_all','-v7.3')
%     end
end


cd(image_storage_folder)


%% Put the 3 views into one image

if overlay == 1
    all_frames_num = electrode_mapping_times;
else
    all_frames_num = total_frames;
end

% for frame_index = 1:all_frames_num
for frame_index = total_frames    
    ViewAngles={'Lateral','Medial','Ventral','Posterior'};
    for side_index = start_side:sides_num
        if side_index == 1
            side = 'L';
        else
            side = 'R';
        end
        for view_index = 1:numel(ViewAngles)
            eval([side '_' ViewAngles{view_index} ' = [];'])
        end
    end
    
    %% name the frames for each side
    for side_index = start_side:sides_num
        % Determine what side we are working on
        if(side_index == 1)
            side = 'L';
            ViewAngles = {'Medial'};
        else
            side = 'R';
            ViewAngles = {'Lateral'};
        end
        
        if views == 4
            ViewAngles={'Lateral','Medial','Ventral','Posterior'};
        elseif views == 2
            ViewAngles={'Lateral','Medial'};
        elseif strcmp(views, 'Lateral')
            ViewAngles={'Lateral'};
        elseif strcmp(views, 'Medial')
            ViewAngles={'Medial'};
        end
        
        for view_index=1:length(ViewAngles)
            Frame = frames_set(side_index,frame_index,view_index).cdata;
            % Do cropping here
            Frame(all(all(Frame == 255, 3), 2), :, :) = [];
            Frame(:, all(all(Frame == 255, 3), 1), :) = [];
            frames_set(side_index,frame_index,view_index).cdata=Frame;
            eval([side '_' ViewAngles{view_index} ' = Frame;'])
        end
    end
    
    %% create the frames for each side
    if views == 4
        for side_index = start_side:sides_num
            if side_index == 1
                
                %             figure
                L_Lateral_Padded = [L_Lateral, ones(size(L_Lateral, 1),...
                    max([0, size(L_Medial, 2) - size(L_Lateral, 2)]), 3) * 255];
                %             imshow(L_Lateral_Padded)
                
                L_Medial_Padded = [L_Medial, ones(size(L_Medial, 1),...
                    max([0, size(L_Lateral, 2) - size(L_Medial, 2)]), 3) * 255];
                %             figure
                %             imshow(L_Medial_Padded)
                
                L_Ventral_Padded = [L_Ventral, ones(size(L_Ventral, 1),...
                    max([0, size(R_Posterior, 2) - size(L_Ventral, 2)]), 3) * 255];

                L_Lateral_Medial = [L_Medial_Padded; L_Lateral_Padded];
                
                L_Medial_Padded_n=[L_Medial_Padded ; ones(100,size(L_Lateral_Padded,2),3)*255];
                L_Lateral_Padded_n=[L_Lateral_Padded ; ones(100,size(L_Lateral_Padded,2),3)*255];

                resize_ratio=size(L_Medial_Padded_n,1)/size(L_Ventral,1);
                L_Ventral_resize=imresize(L_Ventral,[size(L_Ventral,1)*resize_ratio round(size(L_Ventral,2)*resize_ratio)]);
                L_Lateral_Medial_new = [L_Ventral_resize ones(size(L_Lateral_Padded_n,1),40,3)*255 L_Medial_Padded_n ones(size(L_Lateral_Padded_n,1),40,3)*255 L_Lateral_Padded_n];

%                 figure
%                 imshow(L_Lateral_Medial_new)
                
                if figure_type == 0
                    R_Posterior_Padded = [R_Posterior, ones(size(R_Posterior, 1),...
                        max([0, size(L_Ventral, 2) - size(R_Posterior, 2)]), 3) * 255];
                    %               figure
                    %               imshow(R_Posterior_Crop)
                    
                    L_Ventral_R_Posterior = [L_Ventral_Padded; R_Posterior_Padded];
                    
                    L_Lateral_Medial_Padded = [L_Lateral_Medial; ones(max([0, size(L_Ventral_R_Posterior, 1) - ...
                        size(L_Lateral_Medial, 1)]), size(L_Lateral_Medial, 2), 3) * 255];
                    
                    L_Ventral_R_Posterior_Padded = [L_Ventral_R_Posterior; ones(max([0, size(L_Lateral_Medial, 1) - ...
                        size(L_Ventral_R_Posterior, 1)]), size(L_Ventral_R_Posterior, 2), 3) * 255];
                    %                     L_Ventral_R_Posterior_Padded = [L_Ventral_R_Posterior_Padded; ones(155, 545, 3)*255];
                    %             figure
                    %             imshow(L_Ventral_R_Posterior)
                    combined_views_left{frame_index} = [L_Ventral_R_Posterior_Padded, L_Lateral_Medial_Padded];
                elseif figure_type == 3
                    L_Lateral_Medial_Padded = L_Lateral_Medial;
                    L_Posterior_Padded = L_Posterior; %(80:495,155:390,:);
                    %                     figure
                    %                     imshow(L_Posterior_Crop)
                    
                    L_Ventral_L_Posterior = [L_Ventral_Padded;L_Posterior_Padded];
                    %                     figure
                    %                     imshow(L_Ventral_L_Posterior)
                    combined_views_left_f{frame_index} = [L_Ventral_L_Posterior, L_Lateral_Medial_Padded];
                    combined_views_left{frame_index} = [L_Ventral_L_Posterior, L_Lateral_Medial_Padded,...
                        255*ones(957, 1, 3)];
                elseif figure_type == 4
                    combined_views_left_f{frame_index} = L_Lateral_Medial_new;
                    combined_views_left{frame_index} = [L_Lateral_Medial_new 255*ones(size(L_Lateral_Medial_new,1), 400, 3)];
                end
                
                figure
                imshow(combined_views_left{frame_index});
                set(gca,'position',[0 0 1 1],'units','normalized')
                
                % adding time stamp text to the figure
                if figure_type == 4 && ismember(overlay,[2 3 4])
                    load(xtick_spacing_filename)
                    time_txt = {[num2str(round((T(frame_index)-0.5)*1000)) 'ms ~ ' ],
                        ['  ' num2str(round((T(frame_index)+0.05-0.5)*1000)) 'ms']};
                    text(size(combined_views_left{frame_index}, 2) -350,...
                        size(combined_views_left{frame_index}, 1)/2 , time_txt, 'Fontsize', 35)

                end
                
                fn = ['combined_viewsL' '_' num2str(frame_index)];
                export_fig(fn, '-tiff');
%                 print(fn, '-dtiff');
                close

            elseif side_index == 2
                
                R_Lateral_Padded = [R_Lateral, ones(size(R_Lateral, 1),...
                    max([0, size(R_Medial, 2) - size(R_Lateral, 2)]), 3) * 255];
                
                R_Medial_Padded = [R_Medial, ones(size(R_Medial, 1),...
                    max([0, size(R_Lateral, 2) - size(R_Medial, 2)]), 3) * 255];

%                 figure
%                 imshow(R_Lateral_Medial_new)
                
                R_Ventral_Padded = [R_Ventral, ones(size(R_Ventral, 1),...
                    max([0, size(L_Posterior, 2) - size(R_Ventral, 2)]), 3) * 255];
%                 figure
%                 imshow(R_Ventral_Padded)
                
                R_Lateral_Medial = [R_Medial_Padded ; R_Lateral_Padded];

%                 figure
%                 imshow(R_Lateral_Medial_new)  
                
                R_Medial_Padded_n=[R_Medial_Padded ; ones(100,size(R_Lateral_Padded,2),3)*255];
                R_Lateral_Padded_n=[R_Lateral_Padded ; ones(100,size(R_Lateral_Padded,2),3)*255];

                resize_ratio=size(R_Medial_Padded_n,1)/size(R_Ventral,1);
                R_Ventral_resize=imresize(R_Ventral,[size(R_Ventral,1)*resize_ratio round(size(R_Ventral,2)*resize_ratio)]);
                R_Lateral_Medial_new = [R_Lateral_Padded_n ones(size(R_Lateral_Padded_n,1),40,3)*255 R_Medial_Padded_n ones(size(R_Lateral_Padded_n,1),40,3)*255 R_Ventral_resize];  

%                 figure
%                 imshow(R_Lateral_Medial_new)  

                
                if figure_type == 0
                    L_Posterior_Padded = [L_Posterior, ones(size(L_Posterior, 1),...
                        max([0, size(R_Ventral, 2) - size(L_Posterior, 2)]), 3) * 255];
                    %               figure
                    %               imshow(L_Posterior_Crop)
                    
                    R_Ventral_L_Posterior = [R_Ventral_Padded; L_Posterior_Padded];
                    %               figure
                    %               imshow(R_Ventral_R_Posterior)
                    
                    R_Lateral_Medial_Padded = [R_Lateral_Medial; ones(max([0, size(R_Ventral_L_Posterior, 1) - ...
                        size(R_Lateral_Medial, 1)]), size(R_Lateral_Medial, 2), 3) * 255];
                    
                    R_Ventral_L_Posterior_Padded = [R_Ventral_L_Posterior; ones(max([0, size(R_Lateral_Medial, 1) - ...
                        size(R_Ventral_L_Posterior, 1)]), size(R_Ventral_L_Posterior, 2), 3) * 255];
                    
                    combined_views_right{frame_index} = [R_Lateral_Medial_Padded, R_Ventral_L_Posterior_Padded];
                elseif figure_type == 2
                    R_Posterior_Padded = R_Posterior; %(80:495,155:395,:);
                    %               figure
                    %               imshow(R_Posterior_Crop)
                    
                    R_Lateral_Medial_Padded = [R_Lateral_Medial; ones(116,545,3)*255];
                    
                    R_Ventral_R_Posterior = [R_Ventral_Padded;R_Posterior_Padded];
                    %               figure
                    %               imshow(R_Ventral_R_Posterior)
                    
                    combined_views_right{frame_index} = [R_Lateral_Medial_Padded, R_Ventral_R_Posterior];
                    combined_views_right{frame_index} = [combined_views_right{frame_index}; 255*ones(15, 776, 3)];
                elseif figure_type == 4
                    combined_views_right{frame_index} = R_Lateral_Medial_new;
                end
                
                figure
                imshow(combined_views_right{frame_index});
%                 set(gca,'position',[0 0 1 1],'units','normalized')
                fn = ['combined_viewsR' '_' num2str(frame_index)];
                export_fig(fn, '-tiff');
%                 saveas(gcf,[fn '.tif'])
                close
            end
        end
    elseif views == 2
        L_Medial_Padded = L_Medial; %(75:485,:,:);
        R_Lateral_Padded = R_Lateral; %(95:505,:,:);
    end
    
    %------------------------------------------------------------------
    %   Movie Assembly Procedure
    %------------------------------------------------------------------
    % Add all 6 components of each movie frame together onto one plot after
    % cropping. Then after assembling the frame, add to an avi movie buffer
    % output
    
    if figure_type == 0
        figure
        if views == 4
        combined_views = [[combined_views_right{frame_index}; ones(max([0, size(combined_views_left{frame_index}, 1) - ...
            size(combined_views_right{frame_index}, 1)]), size(combined_views_right{frame_index}, 2), 3) * 255]...
            [combined_views_left{frame_index}; ones(max([0, size(combined_views_right{frame_index}, 1) - ...
            size(combined_views_left{frame_index}, 1)]), size(combined_views_left{frame_index}, 2), 3) * 255]];
        elseif views == 2
            combined_views = [L_Medial_Padded R_Lateral_Padded];
        end
        imshow(combined_views);
        set(gca,'position',[0 0 1 1],'units','normalized')
        axis tight
%         frame(frame_index,:,:,:) = combined_views;
        if exist('j','var')
            fn = ['combined_views_full_' num2str(frame_index) '_' num2str(j)];
        else
            fn = ['combined_views_full_' num2str(frame_index)];
        end
        
        if ismember(overlay,[2 3 4])
            load(xtick_spacing_filename)
            time_txt = ['Time = ' num2str(round((T(frame_index)-0.5)*1000)) 'ms, ' num2str(round((T(frame_index)+0.05-0.5)*1000)) 'ms '];
            if views == 4
                uc = uicontrol('Style', 'text', 'Visible', 'off', 'FontName', 'Helvetica',...
                    'FontSize', 24, 'String', time_txt);
                text(size(combined_views, 2) - 1.5 * uc.Extent(3),...
                    size(combined_views, 1) - uc.Extent(4) - 20, time_txt, 'Fontsize', 24)
            else
%                 uc = uicontrol('Style', 'text', 'Visible', 'off', 'FontName', 'Helvetica',...
%                     'FontSize', 20, 'String', time_txt);
%                 text(size(combined_views, 2) - 1.5 * uc.Extent(3),...
%                     size(combined_views, 1) - uc.Extent(4) - 20, time_txt, 'Fontsize', 20)
            end
        elseif overlay == 1
            T = [100 200 300 400 500 600 700 800 900 1000];
            time_txt = ['Time = ' num2str(round((T(frame_index)-0.5)*1000)) 'ms, ' num2str(round((T(frame_index)+0.05-0.5)*1000)) 'ms '];
            text(txt_pos_x,txt_pos_y,time_txt,'Fontsize',20)
        end
        
        export_fig(fn, '-tiff');
        close;
        
    elseif figure_type == 4
        figure
        if views == 4
            combined_views = [[combined_views_right{frame_index}; ones(max([0, size(combined_views_left{frame_index}, 1) - ...
            size(combined_views_right{frame_index}, 1)]), size(combined_views_right{frame_index}, 2), 3) * 255]...
            ones(size(combined_views_right{frame_index}, 1),20,3)*255 ...
            [combined_views_left{frame_index}; ones(max([0, size(combined_views_right{frame_index}, 1) - ...
            size(combined_views_left{frame_index}, 1)]), size(combined_views_left{frame_index}, 2), 3) * 255]];
        end
        imshow(combined_views);
        set(gca,'position',[0 0 1 1],'units','normalized')
        axis tight
        fig=gcf;
        fig.PaperSize=[fig.PaperPosition(3) fig.PaperPosition(4)];
%         frame(frame_index,:,:,:) = combined_views;
        if exist('j','var')
            fn = ['combined_views_full_RL_' num2str(frame_index) '_' num2str(j)];
        else
            fn = ['combined_views_full_RL_' num2str(frame_index)];
        end
        
        if ismember(overlay,[2 3 4])
            load(xtick_spacing_filename)
            time_txt = {[num2str(round((T(frame_index)-0.5)*1000)) 'ms ~ ' ],
                ['  ' num2str(round((T(frame_index)+0.05-0.5)*1000)) 'ms']};
            if views == 4
                uc = uicontrol('Style', 'text', 'Visible', 'off', 'FontName', 'Helvetica',...
                    'FontSize', 24, 'String', time_txt);
                text(size(combined_views, 2) -350,...
                    size(combined_views, 1)/2 , time_txt, 'Fontsize', 22)
            end
        elseif overlay == 1
            T = [100 200 300 400 500 600 700 800 900 1000];
            time_txt = ['Time = ' num2str(round((T(frame_index)-0.5)*1000)) 'ms, ' num2str(round((T(frame_index)+0.05-0.5)*1000)) 'ms '];
            text(txt_pos_x,txt_pos_y,time_txt,'Fontsize',20)
        end

%         saveas(gca,[fn '.tiff'])
        export_fig(fn, '-tiff');
        %         print(fn, '-dtiff');
        % clear combined_views_with_colorbar; clear combined_views_left; clear combined_colorbars;
        % clear Beta_L_Lateral_Crop; clear Beta_L_Medial_Crop; clear Beta_L_Ventral_Crop;
        % clear Beta_R_Lateral_Crop; clear Beta_R_Medial_Crop; clear Beta__Ventral_Crop;
        % clear combined_views;
        
        close
    end
end

%% test figure

if figure_type == 4
%     T2=T*1000;
%     for frame_index=total_frames
%         combined_views_f=[];
%         combined_views_f = [[combined_views_right{frame_index}; ones(max([0, size(combined_views_left_f{frame_index}, 1) - ...
%         size(combined_views_right{frame_index}, 1)]), size(combined_views_right{frame_index}, 2), 3) * 255]...
%         ones(size(combined_views_right{frame_index}, 1),20,3)*255 ...
%         [combined_views_left_f{frame_index}; ones(max([0, size(combined_views_right{frame_index}, 1) - ...
%         size(combined_views_left_f{frame_index}, 1)]), size(combined_views_left_f{frame_index}, 2), 3) * 255]];
%         imshow(combined_views_f)
%         set(gca,'position',[0 0 1 1],'units','normalized')
%         axis tight
%         saveas(gcf,['figure_' file_suffix '_time_' num2str(T2(frame_index)-500)],'epsc')
%         close all   
%     end
    
    % short
    image_set_index=[20 22 24 26 28 30];
    figure('position',[0 0 4000 4000]);
    set(gcf, 'color', [1 1 1]);
    set(gcf,'Visible','on');
    set(gcf,'renderer','Painters') 
    cnt=1;
    for frame_index=image_set_index
        %right
        view_index_all=[1 2 3 4];
        for view_index=view_index_all
            subaxis(size(image_set_index,2),8,cnt,'SpacingVert',0,'MR',0)
            imshow(frames_set(2,frame_index,view_index).cdata)
            cnt=cnt+1;
        end
        
        %left
        view_index_all=[4 3 2 1];
        for view_index=view_index_all
            subaxis(size(image_set_index,2),8,cnt,'SpacingVert',0,'MR',0)
            imshow(frames_set(1,frame_index,view_index).cdata)
            cnt=cnt+1;
        end
    end
    
    saveas(gcf,['figure_' file_suffix '_short'],'epsc')
    close all
    
%     % all
    image_set_index=[32 44 56 68 80 92];
    figure('position',[0 0 4000 4000]);
    set(gcf, 'color', [1 1 1]);
    set(gcf,'Visible','on');
    set(gcf,'renderer','Painters') 
    cnt=1;
    for frame_index=image_set_index
        %right
        view_index_all=[1 2 3 4];
        for view_index=view_index_all
            subaxis(size(image_set_index,2),8,cnt,'SpacingVert',0,'MR',0)
            imshow(frames_set(2,frame_index,view_index).cdata)
            cnt=cnt+1;
        end
        
        %left
        view_index_all=[4 3 2 1];
        for view_index=view_index_all
            subaxis(size(image_set_index,2),8,cnt,'SpacingVert',0,'MR',0)
            imshow(frames_set(1,frame_index,view_index).cdata)
            cnt=cnt+1;
        end
    end
    
    saveas(gcf,['figure_' file_suffix '_all'],'epsc')
    close all




end
%%


% if you want to create a montage of evenly space post-stimulus frames,
% then uncomment this code

if views == 2
    %indices of frames to capture in montage
    frame_idx = [33 37 41 45 49 53 57];
    
    large_image = [];
    for i = 1:length(frame_idx)
        large_image = cat(1,squeeze(frame(frame_idx(i),:,:,:)),large_image);
    end
    figure
    imshow(large_image)
    set(gca,'position',[0 0 1 1],'units','normalized')
    saveas(gcf,'frames_montage.tiff')
    saveas(gcf,'frames_montage.eps')
end


end
%% Below are functions called upon within this function

%==================================================================
%
%Title:         savebrainimages
%Description:   This function takes a brain you have prepared and shines a
%               light on it in 3 views and saves them
%
%Inputs:        "side"          [str] left = 'L', right = 'R'
%               "frame_index"   [int] the "frame" of the "movie" you are
%                               on, which can be #clusters, time period, or
%                               bin number for power analysis
%==================================================================

function frames_set = savebrainimages(side_index,frame_index,views)

if views == 4
    % For each frame, take a snapshot for each of these views
    ViewAngles={'Lateral','Medial','Ventral','Posterior'};
    % These arrays define the viewing and light perspective of the frame
    zoom_factors = [1.25, 1.25, 1.25, 1.25];
    if side_index == 1
        side = 'L';
        Light_Pos = [1200 -200 500; -700 -800 -100;-200, -100, -800; 200, 800, 0];
        View_Pos = [90 0;-90, 0;0 -90;180 0];
        %         View_Pos_xyz = [200 0 0;-200 0 0;0 0 -200;0 1000 0];
    elseif side_index == 2
        side = 'R';
        Light_Pos = [-1200 -300 700; 1000 -1100 -200; -200, -100, -800;-200, 800, 0];
        View_Pos = [-90 0;90, 0;0 -90;180 0];
    else
        Light_Pos=0;
        View_Pos=0;
    end
elseif views == 2
    zoom_factors = [1.25, 1.25];
    if side_index == 1
        side = 'L';
        ViewAngles = {'Lateral', 'Medial'};
        Light_Pos = [1200 -200 500; -700 -800 -100];
        View_Pos = [90 0;-90, 0];
%         Light_Pos = [-700 -800 -300];
%         View_Pos = [-90 -30];
    elseif side_index == 2
        side = 'R';
        ViewAngles = {'Lateral', 'Medial'};
        Light_Pos = [-1200 -300 700; 1000 -1100 -200];
        View_Pos = [-90 0;90, 0];
%         Light_Pos = [-1200 -300 -200];
%         View_Pos = [-90 -15];
    else
        Light_Pos=0;
        View_Pos=0;
    end
elseif strcmp(views, 'Lateral')
    zoom_factors = 1.25;
    if side_index == 1
        side = 'L';
        ViewAngles = {'Lateral'};
        Light_Pos = [1200 -200 500];
        View_Pos = [90 0];
%         Light_Pos = [-700 -800 -300];
%         View_Pos = [-90 -30];
    elseif side_index == 2
        side = 'R';
        ViewAngles = {'Lateral'};
        Light_Pos = [-1200 -300 700];
        View_Pos = [-90 0];
%         Light_Pos = [-1200 -300 -200];
%         View_Pos = [-90 -15];
    else
        Light_Pos=0;
        View_Pos=0;
    end
elseif strcmp(views, 'Medial')
    zoom_factors = 1.25;
    if side_index == 1
        side = 'L';
        ViewAngles = {'Medial'};
        Light_Pos = [-700 -800 -100];
        View_Pos = [-90, 0];
%         Light_Pos = [-700 -800 -300];
%         View_Pos = [-90 -30];
    elseif side_index == 2
        side = 'R';
        ViewAngles = {'Medial'};
        Light_Pos = [1000 -1100 -200];
        View_Pos = [90, 0];
%         Light_Pos = [-1200 -300 -200];
%         View_Pos = [-90 -15];
    else
        Light_Pos=0;
        View_Pos=0;
    end
end

for i = 1:length(ViewAngles)
    lightHandle = light('Position', Light_Pos(i,:));
    
    view(View_Pos(i,:));
    if i == 1
        % with this command, matlab will freeze the aspect ratio of
        % your image instead of changing the object size to fill the
        % window each time
        axis vis3d
    end
    camzoom(zoom_factors(i))
    drawnow
    fn = [side '_' num2str(frame_index) '_' ViewAngles{i}];
    frames_set(i) = getframe;
    %         frames_set(side_index,frame_index,i) = getframe;
    
    % before you save this picture, add a colorbar
    if frame_index == 1
        colorbar
        print(gcf, fn, '-dtiff');       
        %             saveas(gcf, fn);
    end
    
    camzoom(1/zoom_factors(i))
    
    % you have to turn the light off before applying the new light for
    % the next position. The light will be the most recent axes child
    axes_children = get(gca,'children');
    delete(axes_children(1))
    
    % also take off the colorbar so that the next image will not
    % have a colorbar until right before you save it
    if frame_index == 1
        delete(colorbar)
    end
end
end

%==================================================================
%
%Title:         FormatMeshSurface
%Description:   This function loads a surface txt file output from Bioimage
%               Suite into Matlab, and arranges it into a readable format.
%               Repeat twice (once for cortex and once for smooth). The
%               output is a reorganized matrix of vertices.
%
%Inputs:        "vertices_num"  [int] The number of vertices in the file
%               "faces_num"     [int] The number of faces in the file
%               "data_location" [string] should point to the root folder
%                               for the seizure
%                               i.e. 'C:/IcEEGProject_ProcessedData/25/25_1'
%               "filename"      [string] Location of the vertice data
%==================================================================

function [mesh_vertices, mesh_faces]=FormatMeshSurface(vertices_num, faces_num, data_location, filename)
% cd(patient_folder);
% cd ..
cd(data_location)
data_location_previous=pwd;
vertices_range_str=strcat('A1..I',num2str(ceil(vertices_num/3)));
raw_mesh_vertices=dlmread(fullfile( data_location_previous, filename),'',vertices_range_str);
faces_range_str=strcat('B',num2str(ceil(vertices_num/3)+1),'..D',num2str(ceil(vertices_num/3)+faces_num));
raw_mesh_faces=dlmread(fullfile( data_location_previous, filename),'',faces_range_str );
mesh_faces=raw_mesh_faces+1;

vertices_range=1:ceil(vertices_num/3);
mesh_vertices=zeros(ceil(vertices_num/3)*3,3);
mesh_vertices(3*vertices_range-2,1:3)=raw_mesh_vertices(vertices_range, 1:3);
mesh_vertices(3*vertices_range-1,1:3)=raw_mesh_vertices(vertices_range, 4:6);
mesh_vertices(3*vertices_range,1:3)=raw_mesh_vertices(vertices_range, 7:9);
mesh_vertices(ismember(mesh_vertices,zeros(1,3),'rows'),:)=[];
end



%==================================================================
%
%Title:         ProjectElectrode2TransSurf
%Description:
%
%Inputs:
%==================================================================

% Project electrodes to the surface
%Input: number of vertices of the projected surface, matrix of vertices, all
%electrodes, electrodes that is going to be projected, Powerdata, frame index
%Output:electtrode vertices on the new surface, the color/power value of that electrode
function [elecVertices elecVertices_min_dis]=ProjectElectrode2TransSurf(vertices_num, vertices, all_electrodes, electrode)
elecVertices = zeros(size(electrode,1),1);
elecVertices_min_dis = zeros(size(electrode,1),1);
for elec_index = 1:size(electrode,1)
    xCoord = ones(vertices_num,1)*electrode(elec_index,2);
    yCoord = ones(vertices_num,1)*electrode(elec_index,3);
    zCoord = ones(vertices_num,1)*electrode(elec_index,4);
    [minVal, minInd] = min(sqrt((vertices(:,1) - xCoord).^2 + (vertices(:,2) - yCoord).^2 + (vertices(:,3) - zCoord).^2));
    elecVertices(elec_index, 1) = minInd;
    elecVertices_min_dis(elec_index, 1) = minVal;
end
% load('labels');
% overlay_frame_data = overlay_frame_data';
% electrode_vertex_values = overlay_frame_data(frame_index, all_electrodes);
end


% %==================================================================
% %
% %Title:         electrode_data_overlay_ram
% %Description:   This function loads a surface txt file output from Bioimage
% %               Suite into Matlab, and arranges it into a readable format.
% %               Repeat twice (once for cortex and once for smooth). The
% %               output is a reorganized matrix of vertices.
% %
% %Inputs:        "vertices_num"  [int] The number of vertices in the file
% %               "faces_num"     [int] The number of faces in the file
% %               "data_location" [string] should point to the root folder
% %                               for the seizure
% %                               i.e. 'C:/IcEEGProject_ProcessedData/25/25_1'
% %               "filename"      [string] Location of the vertice data
% %==================================================================

function [electrode_vertex_values] = electrode_data_overlay_ram(overlay,frame_index,all_electrodes, file_suffix)
overlay_frame_data = [];
if overlay == 0 || overlay == 1

elseif overlay == 2 || overlay == 3 || overlay == 4
    %% getting power zscores from mean power of each electrode
    load(['stats_traces_', file_suffix, '.mat'])
    overlay_frame_data=squeeze(zscore_traces(:,1,3,:));
end

% once you have the overlay data, you need to match it to an electrode from
% labels

overlay_frame_data = overlay_frame_data';
electrode_vertex_values = overlay_frame_data(frame_index, all_electrodes);

end




%% Determine the transparency for every face/vertex
%==================================================================
%
%Title:         CalculateTransparency
%Description:   This function takes the input values from overlay frame
%               data for each electrode and gives it to all the vertices
%               within radius2 distance from each electrode's vertex
%
%Inputs:        "vertices_num"  The number of vertices in the file
%               "elecVertices"  vertex indices of each electrode to be
%                               plotted on this side of the brain
%               "electrode_vertex_values"
%                               the overlay value assigned to each of those
%                               vertices in elecVertices
%               "overlay"      	[int] this number is an input to the parent
%                               function. If overlay = 0 then it does not
%                               do any weighing based on other electrodes
%                               within a radius2 proximity. All vertices
%                               within radius 2 receive a 1
%
%Outputs:       "vertexCdata"   matrix the length of the number of vertices
%                               with the electrode_vertex_values assigned
%                               to all vertices within radius2 of each
%                               electrode
%==================================================================

%This function should be combined with OverlayColor to speed up processing
%time...

function [ vertexCdata]= CalculateTransparency(vertices_num,vertices, elecVertices, electrode_vertex_values, overlay)
maxAlpha    = 1;
radius1     = 1;
radius2     = 15;
electrodes_num      = length(elecVertices);
electrode_vector    = vertices(elecVertices,:);
vertexCdata = zeros(vertices_num, 1);
j=1;
k=1;

for vertIndex = 1 : vertices_num
    sharedTransparency = [];
    sharedElectrode = [];
    
    xCoord = ones(electrodes_num,1)*vertices(vertIndex, 1);
    yCoord = ones(electrodes_num,1)*vertices(vertIndex, 2);
    zCoord = ones(electrodes_num,1)*vertices(vertIndex, 3);
    
    % Calculate "Transparency" in other words the "weight" based on
    % distance from the vertex
    [distanceToElectrode, electrodeIndice] = sort(sqrt((electrode_vector(:,1) - xCoord).^2 + (electrode_vector(:,2) - yCoord).^2 + (electrode_vector(:,3) - zCoord).^2));
    
    if overlay == 0
        if any(distanceToElectrode <= radius2)
            vertexCdata(vertIndex,1) = 1;
        end
    else
        
        for n = 1:length(electrodeIndice)
            if distanceToElectrode(n) <= radius1
                % Calculates the transparencies for each electrode within a
                % distance less than 15mm from this vertex
                sharedTransparency = [sharedTransparency; maxAlpha];
                % Saves the indices/index of the electrodes/electrode that contribute transparency to this vertex
                sharedElectrode = [sharedElectrode; electrodeIndice(n)];
            elseif distanceToElectrode(n) <= radius2
                % Calculates the transparencies for each electrode within a
                % distance less than 15mm from this vertex
                sharedTransparency = [sharedTransparency; (maxAlpha - ((distanceToElectrode(n) - radius1)/(radius2-radius1))*maxAlpha)];
                % Saves the indices/index of the electrodes/electrode that contribute transparency to this vertex
                sharedElectrode = [sharedElectrode; electrodeIndice(n)];
            end
        end
        % Grabs the zScore values associated with the electrodes that
        % contribute power to this vertex
        zIn = electrode_vertex_values(sharedElectrode)';
        weightedZ = []; %transparencySum = [];
        % Calculates the weight of the z score according to the transparencies and zScores associated with this vertex
        for h = 1:length(zIn)
            %                 weightedZ = [weightedZ zIn(h)*(sharedTransparency(h).^2)]
            weightedZ = [weightedZ zIn(h)*(sharedTransparency(h))];
            %                 transparencySum = [transparencySum sharedTransparency(h)];
        end
        
        
        %             indDelete = [];
        %             for y = 1:length(zIn)
        %                 if isnan(weightedZ(y))
        %                     indDelete = [indDelete y];
        %                 end
        %             end
        %             weightedZ(indDelete) = [];
        %             transparencySum(indDelete) = [];
        
        
        weightedZ = nansum(weightedZ);
        %             transparencySum = sum(transparencySum);
        %             weightedZ = weightedZ/transparencySum;
        % If there is no weighted zScore that means the vertex was greater than
        % 15mm from eevery electrode and has no power or zScore to display in
        % this frame
        if isempty(zIn)
            vertexCdata(vertIndex,1) = 0;
        elseif isnan(zIn)
            vertexCdata(vertIndex,1) = 0;
        else
            vertexCdata(vertIndex,1) = weightedZ;
        end
    end
end
end


