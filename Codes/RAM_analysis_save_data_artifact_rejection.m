%% artifact rejection
clc;clear;

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

load r1_all.mat
fid=fopen('Subjects_list_all.txt','r');
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);

fs=256;
trialtypes = {'RAM_Rec_trials','RAM_nonRec_trials'};
task={'FR1' 'FR2' 'FR3'};
session={'x0x30_' 'x0x31_' 'x0x32_' 'x0x33_' 'x0x34_' 'x0x35_' 'x0x36_' 'x0x37_' 'x0x38_' 'x0x39_'};
session_folder={'0' '1' '2' '3' '4' '5' '6' '7' '8' '9'};

load('time_50ms_no_overlap.mat')
prefix='50ms_no_overlap_together';
window_size=round(fs*0.05); % 50ms
samples=round(T*fs);

% load('time_63ms_no_overlap.mat')
% window_size=round(fs*0.063); % 50ms
% samples=round(T*fs);

base_sample_idx=9;

% atlas

mni_nifti_path = 'C:\yale\bioimagesuite30\images\MNI_T1_1mm_stripped.nii';
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
Tmni = mnit1_1mm_stripped.transform;

M=[0.9975 -0.0073 0.0176 -0.0429 ; 
   0.0146 1.0009 -0.0024 1.5496 ; 
  -0.0130 -0.0093 0.9971 1.1840]; % https://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems

% mask_nifti_path = 'E:\RAM data set\RAM_Public_Data_all\Atlases\yale_broadmann_cortical_dil.nii';
network_nifti_path = 'E:\RAM data set\RAM_Public_Data_all\Atlases\shen_1mm_268_parcellation_network_cortical.nii';
network_1mm_image = mni2fs_load_nii(network_nifti_path);
network_1mm_image_data=network_1mm_image.img;

atlas_nifti_path = 'E:\RAM data set\RAM_Public_Data_all\Atlases\shen_1mm_268_parcellation_cortical.nii';
atlas_1mm_image = mni2fs_load_nii(atlas_nifti_path);
atlas_1mm_image_data=atlas_1mm_image.img;
Atlas_index=unique(atlas_1mm_image_data);
Atlas_index(1,:)=[];
test=find(Atlas_index>0);
Atlas_all_val=cell(max(Atlas_index),251);

% unique(atlas_1mm_image_data(find(network_1mm_image_data==1)))

% Network_name_all={'MedialFrontal','FrontoParietal','DefaultMode' ... 
%     ,'SubcorticalCerebellum','Motor','Visual1','Visual2','VisualAssociation'};

Network_name_all={'FrontoParietal(FEF)','FrontoParietal(IPL)','FrontoParietal(MTG)' ...
    ,'DefaultMode', 'VisualPrimary','VisualAssociation'};


actual_time=257:768;
% actual_time=257:512;
frequency_bands=[45 95; 3 8;40 115;13 30];

% design filter
band_i=3;
d_filter= designfilt('bandpassiir','FilterOrder',40, ...
'HalfPowerFrequency1',frequency_bands(band_i,1),'HalfPowerFrequency2',frequency_bands(band_i,2), ...
'SampleRate',fs);

prefix_s='all_baseline_zscore';
% absolute_pow_val=300;
% folder_name=['plot_' prefix_s '_thr_' num2str(absolute_pow_val)];

folder_name='Voltage_plot_session';
for j=1:1 % FR1 or FR2 or FR3

    for i=1:251 % 251 subject
        try
            clearvars RAM_Rec_trials_sessions RAM_nonRec_trials_sessions location_data_pair

            cd(rootfolder)
            cd('FR1_final')
            cd([num2str(i),'_',r_sublist{i,1}]);
            
            mkdir(folder_name)
            
            load('Raw_data_extended_4s.mat')
            node_n=size(location_data_pair,2);
            
            % grey matter masking
            Atlas_subject_elec_val=[];
            Atlas_subject_ROI_val=[];
            Atlas_subject_network_val=[];
            for elec_i=1:size(location_data_pair,2)
                channel_inf=location_data_pair{1,elec_i};
                fsavg_coord=[channel_inf.atlases.avg.x channel_inf.atlases.avg.y channel_inf.atlases.avg.z]';
                mni_coord = M*[fsavg_coord ; 1];

                volume_buffer=[mni_coord(1) mni_coord(2) mni_coord(3) 1] / Tmni';
                volume_buffer=round(volume_buffer(1:3));
                atlas_index_buffer=atlas_1mm_image_data(volume_buffer(1),volume_buffer(2),volume_buffer(3));
                atlas_index_network_buffer=network_1mm_image_data(volume_buffer(1),volume_buffer(2),volume_buffer(3));
                if atlas_index_buffer>0
%                     Atlas_all_val{atlas_index_buffer,i}=[Atlas_all_val{atlas_index_buffer,i} elec_i];
                    Atlas_subject_elec_val=[Atlas_subject_elec_val elec_i];
                    Atlas_subject_ROI_val=[Atlas_subject_ROI_val atlas_index_buffer];
                    Atlas_subject_network_val=[Atlas_subject_network_val atlas_index_network_buffer];
                end
            end

            % merge
            session_n=size(RAM_Rec_trials_sessions,2);
            for session_i=1:session_n
                RAM_all_trials_sessions{session_i}=cat(3,RAM_Rec_trials_sessions{session_i},RAM_nonRec_trials_sessions{session_i});
            end
            
            subject_max_voltage=nan(node_n,session_n);
            subject_max_power=nan(node_n,session_n);
            for node_i=1:node_n
                try
                    new_node_i=find(Atlas_subject_elec_val==node_i);
%                     if (ismember(node_i,Atlas_subject_elec_val) & Atlas_subject_network_val(new_node_i))
                    if (ismember(node_i,Atlas_subject_elec_val))
%                         all_n_channel_epochs=[];
                        for session_i=1:session_n
                            try
                                all_trials=[];
                                all_trials=RAM_all_trials_sessions{session_i};
                                channel_epochs=squeeze(all_trials(node_i,:,:));

                                fERP = [];
                                for k = 1:size(channel_epochs,2)
                                    data = squeeze(channel_epochs(:,k));
                                    data_buf = filtfilt(d_filter,data);
                                    fERP(:,k)= data_buf(actual_time,:);
                                end
                                fpower = fERP.^2;
                                
                                n_channel_epochs=[];
                                n_channel_epochs=channel_epochs(actual_time,:);
                                
                                mean_n_channel_epochs=[];
                                mean_n_channel_epochs=nanmean(n_channel_epochs,2);
                                subject_max_voltage(node_i,session_i)=max(abs(mean_n_channel_epochs));
                                
                                mean_fpower=[];
                                mean_fpower=nanmean(fpower,2);
                                subject_max_power(node_i,session_i)=max(abs(mean_fpower));
                                
                                %% plot original voltage&power of all trials
                                TIME = linspace(-0.5,1.5,size(n_channel_epochs, 1));  
                                if session_i==1
                                    vol_max_for_plot=1.3*max(max(abs(n_channel_epochs)));
                                    mean_vol_max_for_plot=1.3*max(max(abs(mean_n_channel_epochs)));
                                    power_max_for_plot=1.3*max(max(abs(fpower)));
                                    mean_power_max_for_plot=1.3*max(max(abs(mean_fpower)));
                                end

                                figure('position',[0 0 4000 8000]);
                                set(gcf, 'color', [1 1 1]);
                                set(gcf,'Visible','off');  

                                subplot(4,1,1)
                                plot(TIME,n_channel_epochs); hold on;
                                xlabel('Time (s)','Fontsize',12);
                                ylabel('Voltage (uV)','Fontsize',12);
                                line([0 0],[-1 1],'color','k');
                                ylim([-vol_max_for_plot,vol_max_for_plot])
                                set(gca,'fontsize',12);

                                subplot(4,1,2)
                                std_upper = NaN(1,size(n_channel_epochs,1));
                                std_lower = NaN(1,size(n_channel_epochs,1));
                                for k = 1:size(n_channel_epochs,1)
                                    std_upper(k) = nanmean(n_channel_epochs(k,:))+nanstd(n_channel_epochs(k,:))/sqrt(size(n_channel_epochs,2));
                                    std_lower(k) = nanmean(n_channel_epochs(k,:))-nanstd(n_channel_epochs(k,:))/sqrt(size(n_channel_epochs,2));
                                end
                                % plot the standard error shaded area
                                fill([TIME,fliplr(TIME)],[std_upper,...
                                    fliplr(std_lower)],'r','edgecolor','none'); hold on;
                                alpha(0.15);
                                % plot the line itself and plot it
                                BC_line = plot(TIME,mean_n_channel_epochs); hold on;

                                % set title and axes for ERP plot
        %                         title([' ERPs after artifact rejection, SNR = ' num2str(SNR_index_all(node_i,1))],'Fontsize',12);
                                xlabel('Time (s)','Fontsize',12);
                                ylabel('Voltage (uV)','Fontsize',12);
                                line([0 0],[-1 1],'color','k');
                                ylim([-mean_vol_max_for_plot,mean_vol_max_for_plot]);
                                set(gca,'fontsize',12);

                                subplot(4,1,3)
                                plot(TIME,fpower); hold on;

                                % set title and axes for ERP plot
        %                         title([' ERPs after artifact rejection, SNR = ' num2str(SNR_index_all(node_i,1))],'Fontsize',12);
                                xlabel('Time (s)','Fontsize',12);
                                ylabel('Gamma Power (40-115Hz)','Fontsize',12);
                                line([0 0],[-1 1],'color','k');
                                ylim([0,power_max_for_plot])
                                set(gca,'fontsize',12);

                                subplot(4,1,4)

                                std_upper = NaN(1,size(fpower,1));
                                std_lower = NaN(1,size(fpower,1));
                                for k = 1:size(fpower,1)
                                    std_upper(k) = nanmean(fpower(k,:))+nanstd(fpower(k,:))/sqrt(size(fpower,2));
                                    std_lower(k) = nanmean(fpower(k,:))-nanstd(fpower(k,:))/sqrt(size(fpower,2));
                                end
                                % plot the standard error shaded area
                                fill([TIME,fliplr(TIME)],[std_upper,...
                                    fliplr(std_lower)],'r','edgecolor','none'); hold on;
                                alpha(0.15);
                                % plot the line itself and plot it
                                BC_line = plot(TIME,mean_fpower); hold on;

                                % set title and axes for ERP plot
        %                         title([' ERPs after artifact rejection, SNR = ' num2str(SNR_index_all(node_i,1))],'Fontsize',12);
                                xlabel('Time (s)','Fontsize',12);
                                ylabel('Gamma Power (40-115Hz)','Fontsize',12);
                                line([0 0],[-1 1],'color','k');
                                ylim([0,mean_power_max_for_plot]);
                                set(gca,'fontsize',12);
                                

                                saveas(gcf,[folder_name '/plot_voltage_elec_'  num2str(node_i) '_session_' num2str(session_i) '.png'])
                                close all;

                            catch
                            end
                        end
%                         %% plot original voltage of all trials
%                         figure('position',[0 0 4000 4000]);
%                         set(gcf, 'color', [1 1 1]);
%                         set(gcf,'Visible','off');  
% 
%                         TIME = linspace(-0.5,1.5,size(all_n_channel_epochs, 1));                       
%                         subplot(2,1,1)
%                         plot(TIME,all_n_channel_epochs); hold on;
% 
%                         % set title and axes for ERP plot
% %                         title([' ERPs after artifact rejection, SNR = ' num2str(SNR_index_all(node_i,1))],'Fontsize',12);
%                         xlabel('Time (s)','Fontsize',12);
%                         ylabel('Voltage (uV)','Fontsize',12);
%                         line([0 0],[-1 1],'color','k');
%                         set(gca,'fontsize',12);
%                         
%                         subplot(2,1,2)
% 
%                         std_upper = NaN(1,size(all_n_channel_epochs,1));
%                         std_lower = NaN(1,size(all_n_channel_epochs,1));
%                         for k = 1:size(all_n_channel_epochs,1)
%                             std_upper(k) = nanmean(all_n_channel_epochs(k,:))+nanstd(all_n_channel_epochs(k,:))/sqrt(size(all_n_channel_epochs,2));
%                             std_lower(k) = nanmean(all_n_channel_epochs(k,:))-nanstd(all_n_channel_epochs(k,:))/sqrt(size(all_n_channel_epochs,2));
%                         end
%                         % plot the standard error shaded area
%                         fill([TIME,fliplr(TIME)],[std_upper,...
%                             fliplr(std_lower)],'r','edgecolor','none'); hold on;
%                         alpha(0.15);
%                         % plot the line itself and plot it
%                         BC_line = plot(TIME,nanmean(all_n_channel_epochs,2)); hold on;
% 
%                         % set title and axes for ERP plot
% %                         title([' ERPs after artifact rejection, SNR = ' num2str(SNR_index_all(node_i,1))],'Fontsize',12);
%                         xlabel('Time (s)','Fontsize',12);
%                         ylabel('Voltage (uV)','Fontsize',12);
%                         line([0 0],[-1 1],'color','k');
%                         set(gca,'fontsize',12);
%                         
%                         saveas(gcf,[folder_name '/plot_voltage_'  num2str(node_i) '.png'])
%                         close all;
                        
                    end
                catch
                 disp(sprintf(' %s subject : no value at %d elec',r_sublist{i,1},node_i));
                end
            end
            all_subjects_max_voltage{i}=subject_max_voltage;
            all_subjects_max_power{i}=subject_max_power;
            disp(sprintf(' %s subject complete',r_sublist{i,1}));
        catch
            
        end
    end
    
end

cd(rootfolder)
save('all_subjects_max_voltage_power.mat','all_subjects_max_voltage','all_subjects_max_power','-v7.3')

close all;







 %% plot original voltage
close all;
folder_name='E:\RAM data set\RAM_Public_Data_all\FR1_final\Maximum_voltage_power_histogram_session';

mkdir(folder_name);

all_buffer=[];
for i=1:251
%     buffer=all_subjects_max_voltage{i};
    buffer=all_subjects_max_power{i};
    all_buffer=[all_buffer ; buffer(:)];
end
figure('position',[0 0 4000 4000]);
set(gcf, 'color', [1 1 1]);
set(gcf,'Visible','on');

edges = [0:50:1000];
histogram(all_buffer,edges)
% set(gca,'xlim',[1 100000])
xlabel('Maximum absolute voltage','Fontsize',12);
saveas(gcf,[folder_name '/all2.png'])
close all

%%
clc;clear;
close all;
folder_name='E:\RAM data set\RAM_Public_Data_all\FR1_final\Maximum_voltage_power_histogram_session';

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder);
load('all_subjects_max_voltage_power.mat')

load r1_all.mat
fid=fopen('Subjects_list_all.txt','r');
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);

for i=1:251
    try
        cd(rootfolder)
        cd('FR1_final')
        cd([num2str(i),'_',r_sublist{i,1}]);
        load('Bad_electrode_categories_index.mat')
        
        buffer_vol=all_subjects_max_voltage{i};
        buffer_power=all_subjects_max_power{i};
        
        buffer_vol(Bad_electrode_categories_index,:)=nan;
        buffer_power(Bad_electrode_categories_index,:)=nan;
        
        session_n = size(buffer_vol,2);

        for session_i=1:session_n
            buffer_session_vol=[];
            buffer_session_vol=buffer_vol(:,session_i);
            buffer_session_vol(isnan(buffer_session_vol))=[];

            buffer_session_power=[];
            buffer_session_power=buffer_power(:,session_i);
            buffer_session_power(isnan(buffer_session_power))=[];

            if (buffer_session_vol)
                figure('position',[0 0 4000 4000]);
                set(gcf, 'color', [1 1 1]);
                set(gcf,'Visible','off'); 

                subplot(2,1,1)
                edges_vol = [0:10:1000];
                histogram(buffer_session_vol(:),edges_vol);
                title('Max voltage histogram','Fontsize',12);

                subplot(2,1,2)
                edges_power = [0:10:1000];
                histogram(buffer_session_power(:),edges_power);
                title('Max gamma power histogram','Fontsize',12);

                saveas(gcf,[folder_name '/' num2str(i) '_' r_sublist{i,1} '_session_' num2str(session_i) '.png'])
                close all;
            end
        end
        disp(sprintf(' %s subject complete',r_sublist{i,1}));
    catch
    end
end



%% test : move figure in Voltage_plot_session at each subject
folder_name='Voltage_plot_session';
for i=1:251
    try
        cd(rootfolder)
        cd('FR1_final')
        cd([num2str(i),'_',r_sublist{i,1}]);
        load('Bad_electrode_categories_index.mat')
        cd(folder_name)
        
        mkdir('Bad_category')
        for bad_i=1:size(Bad_electrode_categories_index,2)
            try
                movefile(['*elec_' num2str(Bad_electrode_categories_index(bad_i)) '*'],'Bad_category')
            catch
            end
        end

        disp(sprintf(' %s subject complete',r_sublist{i,1}));
    catch
    end
end


%% test : move figure in Voltage_plot_session at each subject (visual)


%% test : apply bad session&electrodes at each subject (visual & categoty)
clc;clear;
close all;

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder);

load r1_all.mat
fid=fopen('Subjects_list_all.txt','r');
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);

cd('FR1_final')
Bad_electrodes_visual_all=xlsread('RAM_FR1_demographics_final_artifact_final.csv');
buffer_all_session=[50 51 101 162 179 198];
buffer_all_session_index=[2 2 8 3 1 1];

cnt=1;
for i=1:251
    try
        cd(rootfolder)
        cd('FR1_final')
        cd([num2str(i),'_',r_sublist{i,1}]);
        try
            load('Bad_electrode_categories_index.mat')
        catch
            Bad_electrode_categories_index=[];
        end
        subject_bad_electrodes=[];
        subject_bad_electrodes=Bad_electrodes_visual_all(cnt,:);
        subject_bad_electrodes(isnan(subject_bad_electrodes))=[];
        if isempty(subject_bad_electrodes)
            subject_bad_electrodes=[];
        end
        Bad_electrode_visual_electrodes=floor(subject_bad_electrodes);
        
        Bad_electrode_visual_session=[];
        Bad_electrode_visual_session=floor((subject_bad_electrodes-floor(subject_bad_electrodes))*10);
        
        buffer_session_i=find(buffer_all_session==i);
        if isempty(buffer_session_i)
            Bad_electrode_visual_session_n=[];
        end
        Bad_electrode_visual_session_n=buffer_all_session_index(buffer_session_i);
        
        save('Bad_electrode_categories_index_all.mat','Bad_electrode_categories_index','Bad_electrode_visual_session','Bad_electrode_visual_electrodes','Bad_electrode_visual_session_n','-v7.3')
        disp(sprintf(' %s subject complete',r_sublist{i,1}));
        
        cnt=cnt+1;
    catch
    end
end




