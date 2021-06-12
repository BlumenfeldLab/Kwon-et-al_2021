
% Modified by Hunki 12/11/2018 for RAM data set


function displayElectrodesInflated_new_fs_Atlas_all(atlas_result_location,atlas_index,overlay, figure_type, views, thr)

mni_nifti_path = 'MNI_T1_1mm_stripped.nii';
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
Tmni = mnit1_1mm_stripped.transform;
alpha = 0.7;
% colors = distinguishable_colors(length(patients)); % generate different colors

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
    
    %% hard coded number of frames for each overlay setting, need to change later
    total_num_frames = 1;

    v_num = size(fs_surf(side_index).vertices, 1);
    
    fs_surf(side_index).vertices = [fs_surf(side_index).vertices, ones(v_num, 1)] / Tmni';
    fs_surf(side_index).vertices = fs_surf(side_index).vertices(:, 1:3);
    
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

    electrode_mapping_times = 1;
    
    mkdir(atlas_result_location)
    cd(atlas_result_location)
    if ~ismember(overlay,1)
        
        atlas_vertex_index_all=colortable.table(:,5);

        %% make Atlas data
        atlas_color = distinguishable_colors(size(colortable.table,1)); % generate different colors
        for frame_index = 1
            disp(['>> Rendering ' side ' Cortex' ]);  % Output Status
            
            vertex_values_frame=zeros(v_num,1);
            for atlas_i=2:size(atlas_vertex_index_all,1)           
                roi_set_vertex=find(label == atlas_vertex_index_all(atlas_i));
                vertex_values_frame(roi_set_vertex,1)=atlas_i;
            end

            if ((overlay == 2 || overlay == 3 || overlay == 4) && ndims(vertex_values_frame) == 2)% the surface is colored according to zscores of power
                vertex_values_sum = vertex_values_frame;

                % load the colormap and set axis of colormap
                colormap(atlas_color)
                caxis([1 length(atlas_color)]); % can modify

                electrodes_present(1,1)=1;
                electrodes_present(2,1)=1;
                
                if any(electrodes_present(side_index,:) == 1)
                    h = patch('Faces',inflated_surf(side_index).faces, 'Vertices', inflated_surf(side_index).vertices,...
                        'FaceVertexCData',vertex_values_sum(:, 1),'FaceColor','flat',...
                        'LineStyle', 'none', 'FaceAlpha', 'flat');
                    alpha_data = ones(numel(h.FaceVertexCData), 1) * alpha;
                    alpha_data(abs(h.FaceVertexCData) <= 0.01) = 0;
                    lightset = [0.6 0.5 0.1];
                    material(lightset);
                    set(h, 'FaceVertexAlphaData', alpha_data)
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
        end
    end
    
  
    
end


cd(atlas_result_location)
total_num_frames=1;


%% Put the 3 views into one image

if overlay == 1
    all_frames_num = electrode_mapping_times;
else
    all_frames_num = total_num_frames;
end

for frame_index = 1:all_frames_num
    
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
                
                L_Medial_Padded = [L_Medial, ones(size(L_Medial, 1),...
                    max([0, size(L_Lateral, 2) - size(L_Medial, 2)]), 3) * 255];
%                             figure
%                             imshow(L_Lateral_Medial)
                
                L_Ventral_Padded = [L_Ventral, ones(size(L_Ventral, 1),...
                    max([0, size(R_Posterior, 2) - size(L_Ventral, 2)]), 3) * 255];

                L_Lateral_Medial = [L_Medial_Padded; L_Lateral_Padded];
                
                L_Medial_Padded_n=[L_Medial_Padded ; ones(100,size(L_Lateral_Padded,2),3)*255];
                L_Lateral_Padded_n=[L_Lateral_Padded ; ones(100,size(L_Lateral_Padded,2),3)*255];
                if figure_type == 4
                    resize_ratio=size(L_Medial_Padded_n,1)/size(L_Ventral,1);
                    L_Ventral_resize=imresize(L_Ventral,[size(L_Ventral,1)*resize_ratio round(size(L_Ventral,2)*resize_ratio)]);
                    L_Lateral_Medial_new = [L_Ventral_resize ones(size(L_Lateral_Padded_n,1),40,3)*255 L_Medial_Padded_n ones(size(L_Lateral_Padded_n,1),40,3)*255 L_Lateral_Padded_n];
                end
%                 figure
%                 imshow(L_Ventral_Padded)
                
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
%                                 figure
%                                 imshow(combined_views_left{frame_index})
                    combined_views_left{frame_index} = [L_Ventral_R_Posterior_Padded, L_Lateral_Medial_Padded];
                elseif figure_type == 3
                    L_Lateral_Medial_Padded = L_Lateral_Medial;
                    L_Posterior_Padded = L_Posterior; %(80:495,155:390,:);
                    %                     figure
                    %                     imshow(L_Posterior_Crop)
                    
                    L_Ventral_L_Posterior = [L_Ventral_Padded;L_Posterior_Padded];
                    %                     figure
                    %                     imshow(L_Ventral_L_Posterior)
                    combined_views_left{frame_index} = [L_Ventral_L_Posterior, L_Lateral_Medial_Padded,...
                        255*ones(957, 1, 3)];
                elseif figure_type == 4
                    combined_views_left_f{frame_index} = L_Lateral_Medial_new;
                    combined_views_left{frame_index} = [L_Lateral_Medial_new 255*ones(size(L_Lateral_Medial_new,1), 400, 3)];
                end
                
                figure
                imshow(combined_views_left{frame_index});
                set(gca,'position',[0 0 1 1],'units','normalized')
                
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
                if figure_type == 4
                    resize_ratio=size(R_Medial_Padded_n,1)/size(R_Ventral,1);
                    R_Ventral_resize=imresize(R_Ventral,[size(R_Ventral,1)*resize_ratio round(size(R_Ventral,2)*resize_ratio)]);
                    R_Lateral_Medial_new = [R_Lateral_Padded_n ones(size(R_Lateral_Padded_n,1),40,3)*255 R_Medial_Padded_n ones(size(R_Lateral_Padded_n,1),40,3)*255 R_Ventral_resize];  
                end
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

        export_fig(fn, '-tiff');
%         saveas(gcf,fn,'pdf');
        close;
        
    elseif figure_type == 4
        figure
        if views == 4
            combined_views_f = [[combined_views_right{frame_index}; ones(max([0, size(combined_views_left_f{frame_index}, 1) - ...
            size(combined_views_right{frame_index}, 1)]), size(combined_views_right{frame_index}, 2), 3) * 255]...
            ones(size(combined_views_right{frame_index}, 1),20,3)*255 ...
            [combined_views_left_f{frame_index}; ones(max([0, size(combined_views_right{frame_index}, 1) - ...
            size(combined_views_left_f{frame_index}, 1)]), size(combined_views_left_f{frame_index}, 2), 3) * 255]];
        end
        imshow(combined_views_f);
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
%         colorbar
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
function [elecVertices]=ProjectElectrode2TransSurf(vertices_num, vertices, all_electrodes, electrode)
elecVertices = zeros(size(electrode,1),1);
for elec_index = 1:size(electrode,1)
    xCoord = ones(vertices_num,1)*electrode(elec_index,2);
    yCoord = ones(vertices_num,1)*electrode(elec_index,3);
    zCoord = ones(vertices_num,1)*electrode(elec_index,4);
    [minVal, minInd] = min(sqrt((vertices(:,1) - xCoord).^2 + (vertices(:,2) - yCoord).^2 + (vertices(:,3) - zCoord).^2));
    elecVertices(elec_index, 1) = minInd;
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


