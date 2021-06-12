% Arrival Times Script

% Inputs: overlay_frame_data: 0 = nothing, just draw electrodes on blank brain
%                             1 = arrival_times
%                             2 = zScores confirmed perceived
%                             3 = zscores confirmed not perceived
%         createMontage: 1 = create each patient's montage and map for EEG channels, 0 = nothing
%         laterality: collect data and draw on: 0 = bilateral, 1 = left only, 2 = right only, 3 = project both left and right electrodes on the left brain
%         views: 4 = 4 views (lateral, medial, ventral, and posterior)
%                2 = views (Left medial and right lateral)
%
% Modified by Hunki 12/11/2018 for RAM data set


function displayElectrodesInflated_new(data_location, patients, overlay, laterality, views, inflationstep ,electrodeDensity, electrodeDisplay, thr,mean_flag,color_range,file_suffix,Montage_suffix,statistic_flag,all_sub_index)


image_storage_folder = [data_location '/Storage_filtering_no_overlap/Image_storage/' file_suffix '/' Montage_suffix];
data_storage_folder = [data_location '/Storage_filtering_no_overlap/Data_storage/' file_suffix '/' Montage_suffix ];
xtick_spacing_filename = 'time_50ms_no_overlap.mat';
mni_nifti_path = 'C:\yale\bioimagesuite30\images\MNI_T1_1mm_stripped.nii';
mni2fs_dir = 'Functions for Inflated Brain Display\mni2fs-master';
load([mni2fs_dir, filesep, 'surf', filesep, 'transmats.mat']);
load(xtick_spacing_filename);

mkdir(image_storage_folder)
mkdir(data_storage_folder)

mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
Tmni = mnit1_1mm_stripped.transform;
alpha = 0.7;
colors = distinguishable_colors(length(patients)); % generate different colors

%% Side laterality - which side of the brain to draw
start_side = 1;
sides_num = 2;

figHandle = gobjects(1, sides_num);

%% Generating all frames for each side of the brain
for side_index = start_side:sides_num    % Generate frames for L and R
    
    % Determine what side we are working on
    if side_index == 1
        side = 'L';
        hem = 'lh';
    else
        side = 'R';
        hem = 'rh';
    end
    
    %% hard coded number of frames for each overlay setting, need to change later
    
    if overlay == 1 || overlay == 0 % arrival times
        total_num_frames = 1;
    elseif overlay == 2 || overlay == 3 || overlay == 4%
%         total_num_frames = 121;
        total_num_frames = length(T);
    end
    
    %------------------------------------------------------------------
    %   Set-up some local variables
    %------------------------------------------------------------------
    
    surf_fn = fullfile(mni2fs_dir,['/surf/' hem '.surf.gii']);
    fs_surf(side_index) = export(gifti(surf_fn));
    
    v_num = size(fs_surf(side_index).vertices, 1);
    
    fs_surf(side_index).vertices = [fs_surf(side_index).vertices, ones(v_num, 1)] *...
        Tfstovox_rcor' * Trsvoxtomni_rcor' / Tmni';
    fs_surf(side_index).vertices = fs_surf(side_index).vertices(:, 1:3);
    
    surfrender_fn = fullfile(mni2fs_dir,['/surf/' hem '.inflated' num2str(inflationstep) '.surf.gii']);
    inflated_surf(side_index) = export(gifti(surfrender_fn));
    inflated_surf(side_index).vertices = [inflated_surf(side_index).vertices, ones(v_num, 1)] *...
        Tfstovox_rcor' * Trsvoxtomni_rcor' / Tmni';
    inflated_surf(side_index).vertices = inflated_surf(side_index).vertices(:, 1:3);
    
   
    %% create a figure of the left or right side of the brain
    % you will be working with this for different electrodes, overlays,
    % views, and light
    
    figHandle(side_index) = figure('Position',[70,70,700,700]); % Original = [50,50,600,600]Creates a figure graphic object
    %     g=trisurf(cortex.faces,cortex.vertices(:,1),cortex.vertices(:,2),cortex.vertices(:,3));
    
    temp_surf = [];
    temp_surf.hem = hem; % choose the hemesphere 'lh' or 'rh'
    temp_surf.inflationstep = inflationstep; % 1 no inflation, 6 fully inflated
    temp_surf.decimation = 0;
    temp_surf = mni2fs_brain(temp_surf);
    set(temp_surf.p, 'Faces', inflated_surf(side_index).faces, 'Vertices', inflated_surf(side_index).vertices)
    
    set(gcf, 'color', [1 1 1]);
    hold on

    axis off;
    axis equal;
    % if you are plotting arrival times, then each frame is a different set
    % of electrodes
    if overlay == 1
        electrode_mapping_times = 10;
    else
        electrode_mapping_times = 1;
    end
    
    for e = 1:electrode_mapping_times
        scatter_count = 0;
        for p = 1:length(patients)
            try
                %% Get electrode info for the patient

                patient_folder = [data_location '/' num2str(p) '_' patients{p}];

                cd(patient_folder);
                display(patients{p})
                load([side '_MontageMap_' Montage_suffix '.mat']);
                electrode = eval([side '_MontageMap' ]);

                if ~isempty(electrode)
                    buf_electrode = [electrode(:,2:4), ones(size(electrode,1), 1)] / Tmni';
                    electrode(:,2:4) = buf_electrode(:,1:3);
                end



                %% Plot patient's electrodes
                % assign this patient's electrodes one of the colors
                colorElectrode = colors(p,:);

                if isempty(electrode) % Some files only have electrodes on the one side
                    disp('-- No electrode data!');
                    electrodes_present(side_index,p) = 0;
                else
                    all_electrodes = electrode(:,1);
                    electrodes_present(side_index,p) = 1;
                    [elecVertices] = ProjectElectrode2TransSurf(v_num, fs_surf(side_index).vertices, all_electrodes, electrode);

                    % Plot the electrodes in a specified color and count them
                    % in case you need to remove them later
                    if electrodeDisplay == 1
                        scatter3(inflated_surf(side_index).vertices(elecVertices,1),...
                            inflated_surf(side_index).vertices(elecVertices,2),...
                            inflated_surf(side_index).vertices(elecVertices,3), 40, colorElectrode,'filled');
                        scatter_count = scatter_count + 1;
                    end
                end

%                 %% cycle through all the frames to aggregate zscores for each vertex
%                 vertex_values=zeros(v_num,total_num_frames);
%                 if overlay == 2 || overlay == 3 || overlay == 4
%                     for frame_index = 1:total_num_frames
% 
%                         if isempty(electrode) % Some files only have electrodes on one hemisphere, hence the if statement
%                             disp('-- No electrode data!');
% %                             FaceVertexAlphaData    = zeros(v_num, 1);
% %                             vertexCdata            = zeros(v_num, 1);
%                         else
%                             display([patients{p} ' Aggregating Frame ' num2str(frame_index)])
% 
%                             [electrode_vertex_values] = electrode_data_overlay_ram(overlay,frame_index,all_electrodes, file_suffix);
%                             [vertexCdata]        = CalculateTransparency(v_num, inflated_surf(side_index).vertices, elecVertices, electrode_vertex_values, overlay);
% 
%                             vertex_values(:,frame_index) = vertexCdata;
%                         end
%                     end
%                 end
%                 save([side '_vertex_values_' file_suffix '.mat'],'vertex_values');
                
                
            catch
            end
        end
    end
    
    
    
    if ~ismember(overlay,1)
        %------------------------------------------------------------------
        %   Frame Generation with overlaying color onto the brain surface
        %------------------------------------------------------------------
%         load('cmap6to8WRX.mat');
%         load('cmap_hunki_new');
        load('cmap_hunki_5to15.mat');       
%         load('cmap_hunki_10to15.mat');   

        if statistic_flag==1
            cd(data_storage_folder)
            load([side '_vertex_values_sum_5_subjects_stats.mat']);
%             load([side '_vertex_values_sum_5_subjects_stats_roi.mat']);
        end
                    
        for frame_index = 1:total_num_frames

            disp(['>> Rendering ' side ' Cortex, Frame ' num2str(frame_index)]);  % Output Status

            %% read RAM data

            vertex_values_frame=zeros(v_num,length(patients));

            for p = 1:length(patients)
               try
                   patient_folder = [data_location '/' num2str(p) '_' patients{p}];
                   vertex_values=[];
                   load([patient_folder '/' side '_vertex_values_' file_suffix '.mat']);
                   vertex_values_frame(:,p)=vertex_values(:,frame_index);
                   
%                    patient_folder = [data_location '/' num2str(p) '_' patients{p}];
%                    vertex_values_time=[];
%                    load([patient_folder '/' file_suffix '/' side '_vertex_values_' num2str(frame_index) '.mat']);
%                    vertex_values_frame(:,p)=vertex_values_time(:,1);

               catch
               end
            end
            cd(data_storage_folder)
            save([side '_vertex_values_' file_suffix '_' num2str(frame_index) '.mat'],'vertex_values_frame');  
%             load([side '_vertex_values_' file_suffix '_' num2str(frame_index) '.mat']);  
            
            vertex_values_frame=vertex_values_frame(:,all_sub_index);

            cd(image_storage_folder)
            %% sum the count of vertices across patients
            if ((overlay == 2 || overlay == 3 || overlay == 4) && ndims(vertex_values_frame) == 2) || ~electrodeDensity% the surface is colored according to zscores of power
                if ~electrodeDensity
                    if mean_flag == 1
                        vertex_values_sum = squeeze(nansum(vertex_values_frame, 2));
                        divisor = squeeze(sum(vertex_values_frame ~= 0, 2));
%                         divisor(divisor == 0) = 1;

%                         % subject threshold
%                         divisor(divisor <5 ) = 1;
%                         vertex_values_sum(divisor<5)=0;

                        vertex_values_sum = vertex_values_sum ./ divisor;
%                         vertex_values_sum(~subject_mask)=0;

                        % load the colormap and set axis of colormap
                        colormap(mycmap)
    %                     caxis([-6 8]); % can modify
                        caxis(color_range); % can modify

                    elseif mean_flag ==2
                        vertex_values_sum = squeeze(nansum(vertex_values_frame, 2));
                        divisor = squeeze(sum(vertex_values_frame ~= 0, 2));
%                         divisor(divisor == 0) = 1;

%                         % subject threshold
%                         divisor(divisor <5 ) = 1;
%                         vertex_values_sum(divisor<5)=0;

                        vertex_values_sum = vertex_values_sum ./ sqrt(divisor);
%                         vertex_values_sum(~subject_mask)=0;

                        % load the colormap and set axis of colormap
                        colormap(mycmap)
    %                     caxis([-6 8]); % can modify
                        caxis(color_range); % can modify

                        %% save final data
%                         if side_index == 1
%                             L_vertex_values_sum(:,frame_index)=vertex_values_sum;
%                         else
%                             R_vertex_values_sum(:,frame_index)=vertex_values_sum;
%                         end

                    end

                else
                    divisor = squeeze(sum(vertex_values_frame ~= 0, 2));
                    colormap jet
                    caxis([min(divisor(:)) max(divisor(:))])

                end
                
                if statistic_flag==1
                    vertex_values_sum=vertex_values_sum.*vertex_values_p_fdr_mask;
                end
                
                electrodes_present(1,1)=1;
                electrodes_present(2,1)=1;
                
                if any(electrodes_present(side_index,:) == 1)
                    if ~electrodeDensity
                        % for now just overlay the color as number of patients overlapping in
                        % this area without any transparency
                        h = patch('Faces',inflated_surf(side_index).faces, 'Vertices', inflated_surf(side_index).vertices,...
                            'FaceVertexCData',vertex_values_sum(:, 1),'FaceColor','interp',...
                            'LineStyle', 'none', 'FaceAlpha', 'flat');
                        alpha_data = ones(numel(h.FaceVertexCData), 1) * alpha;
                        alpha_data(abs(h.FaceVertexCData) <= thr) = 0;
                        lightset = [0.6 0.5 0.1];
                        material(lightset);
                        set(h, 'FaceVertexAlphaData', alpha_data)
                    else
                        h = patch('Faces',inflated_surf(side_index).faces, 'Vertices', inflated_surf(side_index).vertices,...
                            'FaceVertexCData',divisor(:, 1),'FaceColor','interp',...
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
        end
    end
    
end


cd(image_storage_folder)
%% Put the 3 views into one image

if overlay == 1
    all_frames_num = electrode_mapping_times;
else
    all_frames_num = total_num_frames;
end

for frame_index = 1:all_frames_num-1
    
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
                %             imshow(L_Lateral_Crop)
                
                L_Medial_Padded = [L_Medial, ones(size(L_Medial, 1),...
                    max([0, size(L_Lateral, 2) - size(L_Medial, 2)]), 3) * 255];
                %             figure
                %             imshow(L_Medial_Crop)
                
                L_Ventral_Padded = [L_Ventral, ones(size(L_Ventral, 1),...
                    max([0, size(R_Posterior, 2) - size(L_Ventral, 2)]), 3) * 255];
                %             figure
                %             imshow(L_Ventral_Crop)
                
                %             L_Ventral_padded = [ones(15,1090,3)*255; L_Ventral_Crop,ones(572,839,3)*255];
                %             L_Lateral_Medial = [L_Medial_Crop, L_Lateral_Crop];
                %             combined_views_left{frame_index} = [L_Lateral_Medial; L_Ventral_padded];
                L_Lateral_Medial = [L_Medial_Padded; L_Lateral_Padded];
                
                %             figure
                %             imshow(L_Lateral_Medial)
                
                if laterality == 0
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
                elseif laterality == 3
                    L_Lateral_Medial_Padded = L_Lateral_Medial;
                    L_Posterior_Padded = L_Posterior; %(80:495,155:390,:);
                    %                     figure
                    %                     imshow(L_Posterior_Crop)
                    
                    L_Ventral_L_Posterior = [L_Ventral_Padded;L_Posterior_Padded];
                    %                     figure
                    %                     imshow(L_Ventral_L_Posterior)
                    combined_views_left{frame_index} = [L_Ventral_L_Posterior, L_Lateral_Medial_Padded,...
                        255*ones(957, 1, 3)];
                end
                
                figure
                imshow(combined_views_left{frame_index});
                set(gca,'position',[0 0 1 1],'units','normalized')
                
                % adding time stamp text to the figure
                if laterality == 3 && ismember(overlay,[2 3 4])
                    load(xtick_spacing_filename)
                    time_txt = ['Time = ' num2str(round((T(frame_index)-0.5)*1000)) 'ms, ' num2str(round((T(frame_index+1)-0.5)*1000)) 'ms '];
                    text(400,900,time_txt,'Fontsize',24)
                end
                
                fn = ['combined_viewsL' '_' num2str(frame_index)];
                %             export_fig(fn, '-tiff');
                print(fn, '-dtiff');
                close
            elseif side_index == 2
                
                R_Lateral_Padded = [R_Lateral, ones(size(R_Lateral, 1),...
                    max([0, size(R_Medial, 2) - size(R_Lateral, 2)]), 3) * 255];
                %             figure
                %             imshow(R_Lateral_Crop)
                
                R_Medial_Padded = [R_Medial, ones(size(R_Medial, 1),...
                    max([0, size(R_Lateral, 2) - size(R_Medial, 2)]), 3) * 255];
                %             figure
                %             imshow(R_Medial_Crop)
                R_Lateral_Medial = [R_Medial_Padded; R_Lateral_Padded];
                
                %             figure
                %             imshow(L_Lateral_Medial_padded)
                
                R_Ventral_Padded = [R_Ventral, ones(size(R_Ventral, 1),...
                    max([0, size(L_Posterior, 2) - size(R_Ventral, 2)]), 3) * 255];
                %             figure
                %             imshow(R_Ventral_Crop)
                
                %             R_Ventral_padded = [ones(572,849,3)*255, R_Ventral_Crop];
                %             R_Lateral_Medial = [R_Lateral_Crop, R_Medial_Crop];
                %             combined_views_right{frame_index} = [R_Lateral_Medial; R_Ventral_padded];
                
                if laterality == 0
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
                elseif laterality == 2
                    R_Posterior_Padded = R_Posterior; %(80:495,155:395,:);
                    %               figure
                    %               imshow(R_Posterior_Crop)
                    
                    R_Lateral_Medial_Padded = [R_Lateral_Medial; ones(116,545,3)*255];
                    
                    R_Ventral_R_Posterior = [R_Ventral_Padded;R_Posterior_Padded];
                    %               figure
                    %               imshow(R_Ventral_R_Posterior)
                    
                    combined_views_right{frame_index} = [R_Lateral_Medial_Padded, R_Ventral_R_Posterior];
                    combined_views_right{frame_index} = [combined_views_right{frame_index}; 255*ones(15, 776, 3)];
                end
                
                figure
                imshow(combined_views_right{frame_index});
                set(gca,'position',[0 0 1 1],'units','normalized')
                fn = ['combined_viewsR' '_' num2str(frame_index)];
                %             export_fig(fn, '-tiff');
                saveas(gcf,[fn '.tif'])
                %             print(fn, '-dtiff');
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
    
    if laterality == 0
        figure
%         if views == 4
        combined_views = [[combined_views_right{frame_index}; ones(max([0, size(combined_views_left{frame_index}, 1) - ...
            size(combined_views_right{frame_index}, 1)]), size(combined_views_right{frame_index}, 2), 3) * 255]...
            [combined_views_left{frame_index}; ones(max([0, size(combined_views_right{frame_index}, 1) - ...
            size(combined_views_left{frame_index}, 1)]), size(combined_views_left{frame_index}, 2), 3) * 255]];
%         elseif views == 2
%             combined_views = [L_Medial_Padded R_Lateral_Padded];
%         end
        imshow(combined_views);
        set(gca,'position',[0 0 1 1],'units','normalized')
        axis tight
        frame(frame_index,:,:,:) = combined_views;
        if exist('j','var')
            fn = ['combined_views_full_' num2str(frame_index) '_' num2str(j)];
        else
            fn = ['combined_views_full_' num2str(frame_index)];
        end
        
        if ismember(overlay,[2 3 4])
            load(xtick_spacing_filename)
            time_txt = ['Time = ' num2str(round((T(frame_index)-0.5)*1000)) 'ms, ' num2str(round((T(frame_index+1)-0.5)*1000)) 'ms '];
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
            time_txt = ['Time = ' num2str(round((T(frame_index)-0.5)*1000)) 'ms, ' num2str(round((T(frame_index+1)-0.5)*1000)) 'ms '];
            text(txt_pos_x,txt_pos_y,time_txt,'Fontsize',20)
        end
        
        saveas(gcf,[fn '.tiff'])
        %         export_fig(fn, '-tiff');
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
%                               i.e. 'C:\IcEEGProject_ProcessedData\25\25_1'
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
% %                               i.e. 'C:\IcEEGProject_ProcessedData\25\25_1'
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


