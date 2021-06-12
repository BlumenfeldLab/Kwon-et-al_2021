function [trials_rate]=RAM_surface_mapping_filtering_FARNAM_new(i,subject_name,prefix2)

% addpath(genpath('/home/hk582/project/Codes'));

fs=256;
load('time_63ms_no_overlap.mat')
prefix='63ms_no_overlap_together';
window_size=fs*0.0625; % 62.5ms
base_sample_idx=7;

% fs=256;
% load('time_31ms_no_overlap.mat')
% prefix='31ms_no_overlap_together';
% window_size=fs*0.03125; % 31.25ms
% base_sample_idx=15;

% fs=256;
% load('spectrogram_times_121bins.mat')
% base_sample_idx=28;
% % prefix='121bin_both_baseline';
% % prefix='121bin';
% prefix='121bin_together';
% window_size=fs*0.125; % 125ms

% fs=256;
% window_time=0.0625;
% T=[window_time:window_time/4:0.5-window_time 0.5:window_time:2-window_time];
% prefix='63ms_no_overlap_with_overlaped_base';
% window_size=fs*window_time; % 63ms
% base_sample_idx=25;

% fs=256;
% window_time=0.0625;
% T=[window_time:window_time/4:0.5-window_time 0.5:window_time/4:2-window_time];
% prefix='63ms_overlap_together';
% window_size=fs*window_time; % 63ms
% base_sample_idx=25;

% fs=256;
% window_time=0.0625;
% T=[window_time:window_time/4:0.5-window_time 0.5:window_time/4:2-window_time];
% prefix='63ms_overlap_together';
% window_size=fs*window_time; % 63ms
% base_sample_idx=25;

samples=T*fs;

cd([num2str(i),'_',subject_name]);

frequency_bands=[45 95; 3 8;40 115;13 30];

load Raw_data_extended_5s_recentered_all.mat
load(['rejection_idx_extended_recentered_together_' prefix2 '.mat'])

% prefix2=[prefix2 '_40_95'];
% mkdir(prefix2);

thr_trials_rate=0.2;
meanpower_traces=NaN(size(location_data_pair,2),1,4,size(T,2));
zscore_meanpower_traces=NaN(size(location_data_pair,2),1,4,size(T,2));

% basestim=1:127;
for band_i=3
    % design filter
    d_filter= designfilt('bandpassiir','FilterOrder',40, ...
    'HalfPowerFrequency1',frequency_bands(band_i,1),'HalfPowerFrequency2',frequency_bands(band_i,2), ...
    'SampleRate',fs);
%     fvtool(d_filter)
    
    for node_i=1:size(location_data_pair,2)
        try
            % set the trials you will be working with
            channel_epochs = squeeze(RAM_together_trials_final(node_i,:,:));

            % figure out the number of trials of that type
            trial_num = size(RAM_together_trials_final,3);

            % trial rejection
            % set the indicies of rejected trials
            noisy_trials_idx = noisy_RAM_together_trials_final_idx{node_i};

            % change rejected trials NaNs
            channel_epochs(:,noisy_trials_idx) = [];
            
            fERP = [];
            for k = 1:size(channel_epochs,2)
                data = squeeze(channel_epochs(:,k));
                data_buf = filtfilt(d_filter,data);
                fERP(:,k)= data_buf(257:1024,:); % -1500 to 3500(1 to  1280 samples) --> -500 to 2500 (257 to 1024)  
            end
            fERP_recenter = fERP-repmat(nanmean(fERP(1:128,:),1),768,1);
            fpower = fERP_recenter.^2;

            trials_rate(node_i)=size(fpower,2)/trial_num;
            trialpowers_window=[];
            
            if(trials_rate(node_i)>thr_trials_rate)
                for k = 1:size(fpower,2)
                    fdata_new_trial=fpower(:,k);
                    for window_i=1:size(T,2)
%                              trialpowers_window(window_i,k)=nanmean(fdata_new_trial(samples(window_i)-window_size/2+1:samples(window_i)+window_size/2,1));
                         trialpowers_window(window_i,k)=nanmean(fdata_new_trial(samples(window_i):samples(window_i)+window_size-1,1));
                    end
                end
                mean_trialpowers_window=nanmean(trialpowers_window,2);
                meanpower_traces(node_i,1,band_i,:)=mean_trialpowers_window;                    

                all_mean=nanmean(mean_trialpowers_window(1:base_sample_idx,1)); %% check edge effect
                all_std=nanstd(mean_trialpowers_window(1:base_sample_idx,1));   %% check
                zpower=(mean_trialpowers_window-all_mean)/all_std;
                zscore_meanpower_traces(node_i,1,band_i,:)=zpower;
                
            end

            disp(sprintf(' %s/%s electrode completed!! ',num2str(node_i),num2str(size(location_data_pair,2))));

        catch
            disp(sprintf(' %s/%s electrode not completed!! ',num2str(node_i),num2str(size(location_data_pair,2))));
        end
    end
end

patient_filename_new = ['./meanpower_traces_zscore_filtering_' prefix '_' prefix2 '.mat'];
save(patient_filename_new,'zscore_meanpower_traces','meanpower_traces','T')
% load(patient_filename_new)

% % plot zscore for all electdoes
% for plot_i=1:size(location_data_pair,2)
%     plot(T-0.5,squeeze(meanpower_traces_zscore(plot_i,1,1,:)));
%     hold on;
% end
% ylabel('Gamma power z-score 63ms no overlap');
% xlabel('Time (ms)');
% set(gca,'ylim',[-10 10])
% export_fig([prefix '_' prefix2 '_' subject_name '_zscore_plot' ], '-tiff');
% close all;


%% make Montage
% left_index=1;
% right_index=1;
% MontageMap=[];
% L_MontageMap=[];
% R_MontageMap=[];
%             
% for elec_i=1:size(location_data_pair,2)
%     channer_inf=location_data_pair{1,elec_i};
%     mni_coord=tal2mni([channer_inf.atlases.tal.x channer_inf.atlases.tal.y channer_inf.atlases.tal.z]');
% %     mni_coord=[channer_inf.atlases.mni.x channer_inf.atlases.mni.y channer_inf.atlases.mni.z];
%     MontageMap(elec_i,:)=[elec_i mni_coord(1) mni_coord(2) mni_coord(3)];
% 
%     if(mni_coord(1)<0)
% %                     L_MontageMap(left_index,:)=[elec_i channer_inf.atlases.mni.x channer_inf.atlases.mni.y channer_inf.atlases.mni.z];
%         L_MontageMap(left_index,:)=[elec_i mni_coord(1) mni_coord(2) mni_coord(3)];
% 
%         left_index=left_index+1;
%     else
% %                     R_MontageMap(right_index,:)=[elec_i channer_inf.atlases.mni.x channer_inf.atlases.mni.y channer_inf.atlases.mni.z];
%         R_MontageMap(right_index,:)=[elec_i mni_coord(1) mni_coord(2) mni_coord(3)];                                      
%         right_index=right_index+1;
%     end
% end
% save('R_MontageMap.mat', 'R_MontageMap' );
% save('L_MontageMap.mat', 'L_MontageMap' ); 



%% mapping zscore
patient_name=[num2str(i),'_',subject_name];
prefix_all=[prefix '_' prefix2];

displayElectrodesInflated_FARNAM_elec(patient_name,4, 0, 0, 4, 0, prefix_all,T)

% displayElectrodesInflated_FARNAM(patient_name,4, 0, 0, 4, 0, prefix,T)

end
% 
% overlay=4;
% createMontage=0;
% laterality=0;
% inflationstep=4;
% subtraction=0;
% file_suffix='121bin_together';



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

function displayElectrodesInflated_FARNAM_elec(patient_name, overlay, createMontage, laterality, inflationstep, subtraction, file_suffix,T) %, frames_set)

% image_storage_folder = ['/home/hk582/project/RAM/Storage_vertex_filtering'];
% xtick_spacing_filename = '/home/hk582/project/Codes/spectrogram_times_121bins.mat';
% xtick_spacing_filename = '/home/hk582/project/Codes/spectrogram_times_57bins.mat';

% mni_nifti_path = '/home/hk582/project/Codes/MNI_T1_1mm_stripped.nii';
% mni2fs_dir = '/home/hk582/project/Codes/Functions_Inflated/mni2fs-master';

mni_nifti_path = 'C:\yale\bioimagesuite30\images\MNI_T1_1mm_stripped.nii';
mni2fs_dir = 'E:\RAM data set\RAM_Public_Data_all\Functions for Inflated Brain Display\Functions for Inflated Brain Display/mni2fs-master';

% mni2fs_dir='Functions for Inflated Brain Display\mni2fs-master'
load([mni2fs_dir, filesep, 'surf', filesep, 'transmats.mat']);

% file_suffix = '';
% file_suffix = '_sounds_restricted';

mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
% mnit1_1mm_stripped = mni2fs_load_nii('MNI_T1_1mm_stripped.nii');

Tmni = mnit1_1mm_stripped.transform;

alpha = 0.7;

% list of patients with only right sided contacts

% VTK_folder_name = 'VTK Leah';

% the colors for the different patients' electrodes
%
%
%

%% Side laterality - which side of the brain to draw
start_side = 1;
if laterality == 0 % bilateral
    sides_num = 2;
elseif laterality == 1 % left only
    sides_num = 1;
elseif laterality == 2 % right only
    start_side = 2;
    sides_num = 2;
elseif laterality == 3 % do both left and right on the left brain
    sides_num = 1;
end

% figHandle = gobjects(1, sides_num);

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
    elseif overlay == 2 || overlay == 3 || overlay == 4
%         total_num_frames = 57;
%         total_num_frames = 121;
        total_num_frames = size(T,2);
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

    % if you are plotting arrival times, then each frame is a different set
    % of electrodes
    if overlay == 1
        electrode_mapping_times = 10;
    else
        electrode_mapping_times = 1;
    end
    
    for e = 1:electrode_mapping_times
%         scatter_count = 0;
           
        %% Make the Montage
        % if the left/right electrode montages have not been created, do so
        % now, indicated by createMontage = 0 or = 1
        if createMontage == 1
            % in this case, just do the normal thing: create the normal
            % (non-X-flipped) montage for left and right
            flipX = 0;
            create_electrode_montage(patient_folder,flipX)
            if  laterality == 3 % for x-flipped in addition to regular
                % in this case, create the x-flipped montage for the left and
                % right electrodes but you will later only use the RIGHT
                % electrodes with the coordinates flipped to the LEFT
                flipX = 1;
                create_electrode_montage(patient_folder,flipX)
            end
        end

        %         % Combine EEG Channel with Electrode Position
        %         cd(data_location)
        %
        %         [overlay_frame_data] = OrganizeElectrodePosition(data_location, patient_folder, p, overlay, createMontage, overlay_frame_data);

        % Now load some data files needed for electrode placement and surface shading
        
        load([side '_MontageMap.mat']);
        electrode = eval([side '_MontageMap' ]);

       if ~isempty(electrode)
%                 buf_electrode = [electrode(:,2:4), ones(size(electrode,1), 1)] *...
%                     Tfstovox_rcor' * Trsvoxtomni_rcor' / Tmni';
                buf_electrode = [electrode(:,2:4), ones(size(electrode,1), 1)] / Tmni';
                electrode(:,2:4) = buf_electrode(:,1:3);
       end
            %         electrode_num = size(electrode,1);
            % load the X-flipped montage and add to list of electrodes
        if laterality == 3
            load('R_MontageMap_flipX.mat');
            electrode = cat(1,electrode,R_MontageMap);
            if~isempty(electrode)
    %                     buf_electrode = [electrode(:,2:4), ones(size(electrode,1), 1)] *...
    %                         Tfstovox_rcor' * Trsvoxtomni_rcor' / Tmni';
                buf_electrode = [electrode(:,2:4), ones(size(electrode,1), 1)] / Tmni';
                electrode(:,2:4) = buf_electrode(:,1:3);
            end
        end

        %% Plot patient's electrodes

        %             electrode = [electrode(:, 1), [electrode(:, 2:end), ones(size(electrode, 1), 1)] * Tmni'];
        %             electrode = electrode(:, 1:4);

        % assign this patient's electrodes one of the colors
%         colorElectrode = colors{1};

        if isempty(electrode) % Some files only have electrodes on the one side
            disp('-- No electrode data!');
            electrodes_present(side_index,1) = 0;
            elecVertices = [];
        else
            all_electrodes = electrode(:,1);
            electrodes_present(side_index,1) = 1;

            %                 [elecVertices2]= ProjectElectrode2TransSurf(smooth_v_num, smooth(side_index).vertices, all_electrodes, electrode, side_index);
            %                 electrode(:,2:4)                       = smooth(side_index).vertices(elecVertices2,1:3);
            [elecVertices] = ProjectElectrode2TransSurf(v_num, fs_surf(side_index).vertices, all_electrodes, electrode);

            %Plot the electrodes in a specified color and count them
            %in case you need to remove them later
%             scatter3(inflated_surf(side_index).vertices(elecVertices,1),...
%                 inflated_surf(side_index).vertices(elecVertices,2),...
%                 inflated_surf(side_index).vertices(elecVertices,3), 40, colorElectrode,'filled');
%             scatter_count = scatter_count + 1;
        end

        %% cycle through all the frames to aggregate zscores for each vertex
        if overlay == 2 || overlay == 3 || overlay == 4

            if isempty(electrode) % Some files only have electrodes on one hemisphere, hence the if statement
                vertexCdata            = zeros(v_num, total_num_frames);
                disp('-- No electrode data!');

            else
                vertexCdata            = zeros(v_num, total_num_frames);
                [electrode_vertex_values_zscore] = electrode_data_overlay(all_electrodes,file_suffix);
                vertexCdata(elecVertices,:)=electrode_vertex_values_zscore';
            end
        end
    end

    if side_index == 1
        vertex_infoL = vertexCdata;
    elseif side_index == 2
        vertex_infoR = vertexCdata;
    end
end

if laterality == 0 & overlay ==4
%     delete([patient_name '_vertex_values_elec' file_suffix '_filtering.mat'])
    save([patient_name '_vertex_values_elec_' file_suffix '_filtering.mat'],'vertex_infoL','vertex_infoR','-v7.3')
end


end

%% Below are functions called upon within this function

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


%==================================================================
%
%Title:         create_electrode_montage
%Description:   This function takes a comma delimited sheet with all the
%               names of the electrodes and coordinates and reads it into
%               Matlab. It then saves the electrodes according to left and
%               right side with the (x,y,z) coordinates and the index of
%               that electrode in 'labels'. Essentially the electrode names
%               in 'labels' comes from the EEG .edf file and the coordinates
%               come from the .mgrid file so this reconciles the naming of
%               electrodes between the two
%
%Inputs:        "patient_folder"    [str] location of this patient's files
%               "flipX"             [int] 0 means it's the regular montage
%                                         1 means use the X-flipped montage
%==================================================================

function create_electrode_montage(patient_folder,flipX)
cd (patient_folder);

% if you are creating a montage for the x-flipped electrodes, you want to
% open and then save everything with a "_flipX" attached to it
if flipX == 0
    flipX_str = '';
elseif flipX == 1
    flipX_str = '_flipX';
end

% Read in Electrode locations (Map.xls or Map_flipX.xls) and Montage file (Montage.xls).
% load('labels.mat');
try
    load('labels_all_gray_matter.mat');
catch
    load('labels_first_surgery_all_gray_matter.mat');
end
if size(labels, 1) == 1
    labels = labels';
end
b1 = linspace(1, length(labels), length(labels))';
b1 = num2cell(b1);
b2 = labels(:,1);
montage = [b1 b2];
% [a1, a2, map]=xlsread('Map.xls');
[a1, a2, map]=xlsread('Map_all_gray_matter.xlsx');
% Convert the long channel names (i.e "A_L_Most_Mid_Frontal_Polar_12" into the abbreviated names found in the montage (i.e. "A12")
try
    for i=1:size(map,1)
        str             =   map{i,1};
        str_delim       =   regexp(str,'_','split');
        if length(str_delim)==1 || length(str_delim) == 2
            if startsWith(str_delim{1,1}, {'L', 'R'})
                map{i,5}=strtrim(str_delim{1,1}(1));
                map{i,1}=strtrim(str_delim{1,2});
                map{i,1}(map{i,1} == ';') = [];
            else
                map{i,5}=strtrim(str_delim{1,1}(end));
            end
        else
            firstLetter     =   str_delim{1};
            positionLetter  =   strtrim(str_delim{2});
            digitStr        =   str(isstrprop(str,'digit'));
            replaceStr      =   strcat(firstLetter, digitStr);
            map{i,1}        =   replaceStr; % i.e. "A12"
            map{i,5}        =   positionLetter; % Either L or R
        end
    end
catch
    disp(['Error reading Map file! There was a problem on line ' num2str(i) '. Verify that file is in the correct format (i.e each row is Letter_Side*Number. Unacceptable entries: J_Sup_Parietal8,PEG1)']);
    error('MATLAB:myCode:dimensions', '');
end

Montage=montage;
%%
%     try
% Combine montage with 3-d electrode locations
for i=1:size(montage,1)
    MontageMap(i,1)=montage{i,1}; % Abbreviated Channel ("A12")
    for j=1:size(map,1)
        if strcmp(strtrim(montage{i,2}), map{j,1})
            if i == 113
                leah = [];
            end
            MontageMap(i,2)=map{j,2}; % X Coordinates
            MontageMap(i,3)=map{j,3}; % Y Coordinates
            MontageMap(i,4)=map{j,4}; % Z Coordinates
            %                     MontageMap(i,5) = newTimes{i,1}; % Arrival Time
            Position(i,1)=map{j,5}; % Side (R or L)
        end
    end
end
%         % catch
%         %     disp(['Error on MontageMap write! There was a problem on line ' num2str(i) '. Verify that file is in the correct format (i.e each row is Letter_Side*Number. Unacceptable entries: J_Sup_Parietal8)']);
%         %     error('MATLAB:myCode:dimensions', '');
%     end
j=1;k=1;

% check MontageMap to see if there are any errors that result in
% (0,0,0) coordinates for (x,y,z)
for i = 1:size(MontageMap,2)
    if MontageMap(i,1) == 0 && MontageMap(i,2) == 0 && MontageMap(i,3) == 0
        display(MontageMap{i,1});
        error('This electrode name in Map or Map_flipX is not matching, double check the original Map sheet for misnamed electrodes compared to labels');
    end
end


% Break the montage into L and R sides
for i=1: length(Position)
    if (Position(i)=='L')
        L_MontageMap(j,1:4)=MontageMap(i,1:4);
        j=j+1;
    else
        R_MontageMap(k,1:4)=MontageMap(i,1:4);
        k=k+1;
    end
end

if (j==1)
    L_MontageMap=[];
elseif (k==1)
    R_MontageMap=[];
end
% Save data into mat files
MontageMap=[L_MontageMap; R_MontageMap];
MontageMap=sortrows(MontageMap,1);
save(['Montage' flipX_str '.mat'], 'Montage');
save(['MontageMap' flipX_str '.mat'], 'MontageMap');
save(['L_MontageMap' flipX_str '.mat'], 'L_MontageMap');
save(['R_MontageMap' flipX_str '.mat'], 'R_MontageMap');

end

%==================================================================
%
%Title:         electrode_data_overlay
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

function [electrode_vertex_values_zscore] = electrode_data_overlay(all_electrodes,file_suffix)
    %% getting power zscores from mean power of each electrode
load(['./meanpower_traces_zscore_filtering_' file_suffix '.mat']);
overlay_frame_data = squeeze(zscore_meanpower_traces(:,1,3,:));
    
% once you have the overlay data, you need to match it to an electrode from
% labels

overlay_frame_data = overlay_frame_data';
electrode_vertex_values_zscore = overlay_frame_data(:, all_electrodes);

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



