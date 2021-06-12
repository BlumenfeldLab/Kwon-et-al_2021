
function RAM_FARNAM_AR_Power_all(i,subject_name)

% addpath(genpath('/home/hk582/project/Codes'));
cd([num2str(i),'_',subject_name]);
load('Raw_data_extended_4s_all.mat')

actual_time=257:768;
node_n=size(location_data_pair,2);
session_n=size(RAM_all_trials_sessions,2);
prefix='Wendy_hunki_std_5_20';

folder_name=['Power_plot_std_5_20_test'];
mkdir(folder_name)

fs=256;
load('time_50ms_no_overlap.mat')
window_size=round(fs*0.05); % 50ms
samples=round(T*fs);
base_sample_idx=9;

% fs=256;
% load('time_31ms_no_overlap.mat') % 31.125ms
% base_sample_idx=16;
% window_size=0.03125;

% load('time_63ms_no_overlap.mat') % 63.25ms
% base_sample_idx=8;
% window_size=0.0626;


frequency_bands=[45 95; 3 8;40 115;13 30];
% design filter
band_i=3;
d_filter= designfilt('bandpassiir','FilterOrder',40, ...
'HalfPowerFrequency1',frequency_bands(band_i,1),'HalfPowerFrequency2',frequency_bands(band_i,2), ...
'SampleRate',fs);

RAM_all_trials=[];
for session_i=1:session_n
    try
        if ~isempty(RAM_all_trials_sessions{session_i})
            RAM_all_trials=cat(3,RAM_all_trials,RAM_all_trials_sessions{session_i});
        end
    catch
    end
end

% artifact rejection
% noisy_RAM_all_trials_number=nan(node_n,1);
% RAM_all_trials_total_number=size(RAM_all_trials,3);
% thr_i=0.2;
% 
% for node_i=1:node_n
%     try
%         all_trials_node=[];
%         all_trials_node=squeeze(RAM_all_trials(node_i,actual_time,:));
%         
%         [meansqerr,noisy_idx_buffer]=artifact_rejection_all_together(all_trials_node,fs,thr_i);                        
%         noisy_idx_buffer=unique(noisy_idx_buffer);
% %         histogram(meansqerr,200,'BinLimits',[0,20000],'Normalization', 'probability')
%         
%         noisy_RAM_all_trials_idx{node_i} = noisy_idx_buffer;
%         noisy_RAM_all_trials_mse{node_i} = meansqerr;
%         noisy_RAM_all_trials_number(node_i) = length(unique(noisy_idx_buffer));
%         disp(sprintf('electrode #%d : artifact rejection is completed! ',node_i));
%     catch
%         disp(sprintf('electrode #%d : artifact rejection is failed! ',node_i));
%     end
% end
% save(['rejection_idx_' prefix '.mat'],'noisy_RAM_all_trials_idx','noisy_RAM_all_trials_mse',...
%     'noisy_RAM_all_trials_number','RAM_all_trials_total_number','-v7.3')
% load(['rejection_idx_' prefix '.mat'])

% save(['rejection_idx_Wendy_hunki_2.mat'],'noisy_RAM_all_trials_idx','noisy_RAM_all_trials_mse',...
%     'noisy_RAM_all_trials_number','RAM_all_trials_total_number','-v7.3')
load(['rejection_idx_' prefix '.mat'])


% power
zscore_traces=NaN(size(location_data_pair,2),1,4,size(T,2));
noisy_RAM_all_trials=nan(node_n,2);

for node_i=1:node_n
    try
        channel_epochs=squeeze(RAM_all_trials(node_i,:,:));
        TIME = linspace(-0.5,1.5,size(T,2));
        TIME_raw = linspace(-0.5,1.5,size(actual_time,2));  
        
        figure('position',[0 0 4000 8000]);
        set(gcf, 'color', [1 1 1]);
        set(gcf,'Visible','off');  
        
        subplot(3,2,1)
        plot(TIME_raw,channel_epochs(actual_time,:)); hold on;
        title(['<Raw voltages>   total number : ' num2str(size(channel_epochs,2))])
        xlabel('Time (s)','Fontsize',12);
        ylabel('Voltage','Fontsize',12);
        ylim([-500 500])
        set(gca,'fontsize',12);
        
        
        % remove bad trials 
        rest_num_voltage=size(channel_epochs,2)-size(noisy_RAM_all_trials_idx{node_i},2);
        channel_epochs(:,noisy_RAM_all_trials_idx{node_i})=[];

        subplot(3,2,2)
        plot(TIME_raw,channel_epochs(actual_time,:)); hold on;
        
        title(['<Remaining voltages>   remaining total number : ' num2str(rest_num_voltage)])
        xlabel('Time (s)','Fontsize',12);
        ylabel('Voltage','Fontsize',12);
        ylim([-500 500])
        set(gca,'fontsize',12);

        fpower = [];
%         trialpowers_window = [];
        for k = 1:size(channel_epochs,2)
            data = squeeze(channel_epochs(:,k));
            data_buf = filtfilt(d_filter,data);
            % power
            data_buf=data_buf.^2;
            fpower(:,k)= data_buf(actual_time,:);
%                 trialpowers_window(:,k)=nanmean(reshape(data_buf(actual_time,:),[fs*window_size,2*1/window_size]),1);
        end

        trialpowers_window=[];
        %% method 1
        for k = 1:size(fpower,2)
            fdata_new_trial=fpower(:,k);
            for window_i=1:39
                trialpowers_window(window_i,k)=nanmean(fdata_new_trial(samples(window_i):samples(window_i)+window_size-1,1));
            end
        end

        subplot(3,2,3)
        plot(TIME,trialpowers_window); hold on;
        title(['<Raw power>  total number : ' num2str(size(trialpowers_window,2))])
        xlabel('Time (s)','Fontsize',12);
        ylabel('Gamma power (40-115Hz)','Fontsize',12);
        ylim([0 100])
        set(gca,'fontsize',12);
        
        trialpowers_window_buffer=trialpowers_window;
        
        %% method 1
        std_num =5;
        mean_mean_all_trials = nanmean(trialpowers_window_buffer(:));
        std_mean_all_trials = nanstd(trialpowers_window_buffer(:));
        
        trialpowers_window_buffer(trialpowers_window_buffer > mean_mean_all_trials + std_num*std_mean_all_trials | trialpowers_window_buffer < mean_mean_all_trials - std_num*std_mean_all_trials) = NaN;
        
        noisy_trials_power_idx=[];
        for h = 1:size(trialpowers_window_buffer,2)
            if ~isempty(find(isnan(trialpowers_window_buffer(:,h)),1))       
                noisy_trials_power_idx = [noisy_trials_power_idx h];       
            end
        end
        trialpowers_window_buffer(:,noisy_trials_power_idx) = [];

        %% method 2 (baseline)
        std_num =20;
        mean_all_trials_base = nanmean(trialpowers_window_buffer(1:base_sample_idx,:),1);
        mean_mean_all_trials_base = nanmean(mean_all_trials_base);
        std_mean_all_trials_base = nanstd(mean_all_trials_base);
        
        trialpowers_window_buffer(trialpowers_window_buffer > mean_mean_all_trials_base + std_num*std_mean_all_trials_base | trialpowers_window_buffer < mean_mean_all_trials_base - std_num*std_mean_all_trials_base) = NaN;
        
        noisy_trials_power_idx=[];
        for h = 1:size(trialpowers_window_buffer,2)
            if ~isempty(find(isnan(trialpowers_window_buffer(:,h)),1))       
                noisy_trials_power_idx = [noisy_trials_power_idx h];       
            end
        end
        trialpowers_window_buffer(:,noisy_trials_power_idx) = [];
        
%         figure;
%         plot(TIME,trialpowers_window(:,noisy_trials_power_idx)); hold on;
        
        trialpowers_window=trialpowers_window_buffer;
%         trialpowers_window(:,noisy_trials_power_idx) = [];
        remaining_ratio=[];
        remaining_ratio=size(trialpowers_window,2)/RAM_all_trials_total_number;
        noisy_RAM_all_trials(node_i,1)=size(trialpowers_window,2);
        noisy_RAM_all_trials(node_i,2)=remaining_ratio;
        
%         subplot(3,2,4)
%         plot(TIME,trialpowers_window); hold on;
%         title(['<Raw powers>  remaining total number : ' num2str(size(trialpowers_window,2))])
%         xlabel('Time (s)','Fontsize',12);
%         ylabel('Gamma power (40-115Hz)','Fontsize',12);
%         ylim([0 100])
%         set(gca,'fontsize',12);

        subplot(3,2,4)
        plot(TIME,trialpowers_window); hold on;
        title(['<Raw powers>  remaining total number : ' num2str(size(trialpowers_window,2))])
        xlabel('Time (s)','Fontsize',12);
        ylabel('Gamma power (40-115Hz)','Fontsize',12);
        ylim([0 100])
        set(gca,'fontsize',12);

%         %z-scoring 1
%         mean_trialpowers_window=[];
%         mean_trialpowers_window = squeeze(nanmean(trialpowers_window,2));
%         base_std=std(mean_trialpowers_window(1:base_sample_idx));
%         base_mean=mean(mean_trialpowers_window(1:base_sample_idx));
%         z_mean_trialpowers_window=[];
%         z_mean_trialpowers_window=(mean_trialpowers_window-base_mean)/base_std;

        %z-scoring 2
        base_std=nanstd(nanmean(trialpowers_window(1:base_sample_idx,:)));
        base_mean=nanmean(nanmean(trialpowers_window(1:base_sample_idx,:)));

        z_trialpowers_window=[];
        z_trialpowers_window=(trialpowers_window-base_mean)/base_std;

        z_mean_trialpowers_window=[];
        z_mean_trialpowers_window=nanmean(z_trialpowers_window,2);

        subplot(3,2,5)
        plot(TIME,z_trialpowers_window); hold on;
        title(['<zscore powers>  remaining total number : ' num2str(size(trialpowers_window,2))])
        xlabel('Time (s)','Fontsize',12);
        ylabel('Gamma power (40-115Hz)','Fontsize',12);
        ylim([-20 20])
        set(gca,'fontsize',12);
        
%         subplot(3,2,5)
%         plot(TIME,mean_trialpowers_window); hold on;
%         title('<Mean raw gamma power>')
%         xlabel('Time (s)','Fontsize',12);
%         ylabel('Gamma power (40-115Hz)','Fontsize',12);
%         set(gca,'fontsize',12);
        

        
%         subplot(3,2,6)
%         plot(TIME,z_mean_trialpowers_window);
%         ylabel('','Fontsize',12);
%         title(['<Gamma Power z-score>   remaining trial ratio : ' num2str(remaining_ratio*100) '%'])
%         ylim([-10,30]);
%         set(gca,'fontsize',12);
%         line([-0.5,1.5],[2,2],'color','r');
%         line([-0.5,1.5],[-2,-2],'color','r');
        
        subplot(3,2,6)
        plot(TIME,z_mean_trialpowers_window); hold on;
        title('<Mean zscore gamma power>')
        xlabel('Time (s)','Fontsize',12);
        ylabel('Gamma power (40-115Hz)','Fontsize',12);
        ylim([-2 2])
        set(gca,'fontsize',12);

        
        zscore_traces(node_i,1,band_i,:)=z_mean_trialpowers_window';

        saveas(gcf,[folder_name '/plot_power_elec_'  num2str(node_i) '_test.png'])
        close all;

        disp(sprintf('electrode #%d : power extraction is completed! ',node_i));
    catch
        disp(sprintf('electrode #%d : power extraction is not completed! ',node_i));
        
    end
end

% patient_filename = ['./stats_traces_' prefix '.mat'];
% save(patient_filename,'zscore_traces','noisy_RAM_all_trials','location_data_pair','T','-v7.3')



%% mapping zscore to surface
% new_fdr_val=fdr(ranksum_p_traces(Atlas_subject_elec_val),0.05);
% Atlas_subject_elec_val_fdr=[];
% if ~isempty(new_fdr_val)
%     Atlas_subject_elec_val_fdr=find(ranksum_p_traces<=new_fdr_val)';
% end
% patient_filename = ['./stats_traces_' prefix '_fdr.mat'];
% save(patient_filename,'zscore_traces','ranksum_p_traces','Atlas_subject_elec_val_fdr','location_data_pair','T','-v7.3')
% 
% % Make the Montage
% M=[0.9975 -0.0073 0.0176 -0.0429 ; 
%    0.0146 1.0009 -0.0024 1.5496 ; 
%   -0.0130 -0.0093 0.9971 1.1840]; % https://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems
% 
% % Montage new
% left_index=1;
% right_index=1;
% L_MontageMap=[];
% R_MontageMap=[];
% Montage_prefix='final_clean_dil_Wendy';
% 
% masked_elec_n=size(Atlas_subject_elec_val,2);
% for elec_i=1:masked_elec_n
%     masked_elec=Atlas_subject_elec_val(elec_i);
%     channel_inf=location_data_pair{1,masked_elec};
%     fsavg_coord=[channel_inf.atlases.avg.x channel_inf.atlases.avg.y channel_inf.atlases.avg.z]';
%     mni_coord = M*[fsavg_coord ; 1];
% 
%     if(mni_coord(1)<0)
%         L_MontageMap(left_index,:)=[masked_elec mni_coord(1) mni_coord(2) mni_coord(3)];
%         left_index=left_index+1;
%     else
%         R_MontageMap(right_index,:)=[masked_elec mni_coord(1) mni_coord(2) mni_coord(3)]; 
%         right_index=right_index+1;
%     end
% end
% save(['R_MontageMap_' Montage_prefix '.mat'], 'R_MontageMap');
% save(['L_MontageMap_' Montage_prefix '.mat'], 'L_MontageMap');
% 
% % Montage new fdr
% left_index=1;
% right_index=1;
% L_MontageMap=[];
% R_MontageMap=[];
% Montage_prefix_fdr='final_clean_dil_Wendy_fdr';
% masked_elec_n=size(Atlas_subject_elec_val_fdr,2);
% 
% if masked_elec_n>0
%     for elec_i=1:masked_elec_n
%         masked_elec=Atlas_subject_elec_val(elec_i);
%         channel_inf=location_data_pair{1,masked_elec};
%         fsavg_coord=[channel_inf.atlases.avg.x channel_inf.atlases.avg.y channel_inf.atlases.avg.z]';
%         mni_coord = M*[fsavg_coord ; 1];
% 
%         if(mni_coord(1)<0)
%             L_MontageMap(left_index,:)=[masked_elec mni_coord(1) mni_coord(2) mni_coord(3)];
%             left_index=left_index+1;
%         else
%             R_MontageMap(right_index,:)=[masked_elec mni_coord(1) mni_coord(2) mni_coord(3)]; 
%             right_index=right_index+1;
%         end
%     end
% end
% save(['R_MontageMap_' Montage_prefix_fdr '.mat'], 'R_MontageMap');
% save(['L_MontageMap_' Montage_prefix_fdr '.mat'], 'L_MontageMap');
% 
% %% elec to vertex
% overlay=4;
% inflationstep=5;
% file_suffix=prefix;
% 
% displayElectrodesInflated_new_FARNAM(overlay, inflationstep ,file_suffix,Montage_prefix) %, frames_set)
% 
% file_suffix=[prefix '_fdr'];
% displayElectrodesInflated_new_FARNAM(overlay, inflationstep ,file_suffix,Montage_prefix_fdr) %, frames_set)

end





function displayElectrodesInflated_new_FARNAM(overlay, inflationstep ,file_suffix,Montage_prefix)

xtick_spacing_filename = 'time_31ms_no_overlap.mat';

mni_nifti_path = 'MNI_T1_1mm_stripped.nii';
mni2fs_dir = 'Functions_Inflated/mni2fs-master';
% mni_nifti_path = 'C:\yale\bioimagesuite30\images\MNI_T1_1mm_stripped.nii';
% mni2fs_dir = 'E:\RAM data set\RAM_Public_Data_all\Functions for Inflated Brain Display\Functions for Inflated Brain Display\mni2fs-master';

load([mni2fs_dir, filesep, 'surf', filesep, 'transmats.mat']);
load(xtick_spacing_filename);

mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
Tmni = mnit1_1mm_stripped.transform;

%% Side laterality - which side of the brain to draw
start_side = 1;
sides_num = 2;

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
    total_num_frames = length(T);
    
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

    
    % if you are plotting arrival times, then each frame is a different set
    % of electrodes
    if overlay == 1
        electrode_mapping_times = 10;
    else
        electrode_mapping_times = 1;
    end
    
    for e = 1:electrode_mapping_times
            
        load([side '_MontageMap_' Montage_prefix '.mat']);
        electrode = eval([side '_MontageMap' ]);

        if ~isempty(electrode)
            buf_electrode = [electrode(:,2:4), ones(size(electrode,1), 1)] / Tmni';
            electrode(:,2:4) = buf_electrode(:,1:3);
        end

        %% Plot patient's electrodes

        if isempty(electrode) % Some files only have electrodes on the one side
            disp('-- No electrode data!');
        else
            all_electrodes = electrode(:,1);
            [elecVertices] = ProjectElectrode2TransSurf(v_num, fs_surf(side_index).vertices, all_electrodes, electrode);

        end

        %% cycle through all the frames to aggregate zscores for each vertex
        vertex_values=zeros(v_num,total_num_frames);

        if overlay == 2 || overlay == 3 || overlay == 4
            for frame_index = 1:total_num_frames

                if isempty(electrode) % Some files only have electrodes on one hemisphere, hence the if statement
                    disp('-- No electrode data!');
%                             FaceVertexAlphaData    = zeros(v_num, 1);
%                             vertexCdata            = zeros(v_num, 1);
                else
                    display([' Aggregating Frame ' num2str(frame_index)])

                    [electrode_vertex_values] = electrode_data_overlay_ram(overlay,frame_index,all_electrodes, file_suffix);
                    [vertexCdata]        = CalculateTransparency(v_num, inflated_surf(side_index).vertices, elecVertices, electrode_vertex_values, overlay);

                    vertex_values(:,frame_index) = vertexCdata;
                end
            end
        end
        save([side '_vertex_values_' file_suffix '.mat'],'vertex_values');
    end
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




