%% artifact rejection & power extraction
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

mask_nifti_path = 'E:\RAM data set\RAM_Public_Data_all\Atlases\yale_broadmann_cortical_dil.nii';
mask_1mm_image = mni2fs_load_nii(mask_nifti_path);
mask_1mm_image_data=mask_1mm_image.img;
actual_time=257:768;
% actual_time=257:512;
frequency_bands=[45 95; 3 8;40 115;13 30];

% design filter
band_i=3;
d_filter= designfilt('bandpassiir','FilterOrder',40, ...
'HalfPowerFrequency1',frequency_bands(band_i,1),'HalfPowerFrequency2',frequency_bands(band_i,2), ...
'SampleRate',fs);

vol_thr=10;
prefix='final_clean_vol_thr_10';
for j=1:1 % FR1 or FR2 or FR3

    for i=1:251%251 % 251 subject
        try
            clearvars RAM_Rec_trials_sessions RAM_nonRec_trials_sessions location_data_pair

            cd(rootfolder)
            cd('FR1_final')
            cd([num2str(i),'_',r_sublist{i,1}]);

            load('Raw_data_extended_4s.mat')
            
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
            
            
            load('Bad_electrode_categories_index_all.mat')
            % remove bad electrodes from category
            if ~isempty(Bad_electrode_categories_index)
                for session_i=1:session_n
                    try
                        all_trials=[];
                        all_trials=RAM_all_trials_sessions{session_i};
                        all_trials(Bad_electrode_categories_index,:,:)=nan;
                        RAM_all_trials_sessions{session_i}=all_trials;
                    catch
                    end
                end
            end
            
            % remove bad electrodes & sessions (visual)
            if ~isempty(Bad_electrode_visual_electrodes)
                for bad_elec_i=1:size(Bad_electrode_visual_electrodes,2)
                    flag=Bad_electrode_visual_session(bad_elec_i);
                    if flag == 0
                        for session_i=1:session_n
                            try
                                all_trials=[];
                                all_trials=RAM_all_trials_sessions{session_i};
                                all_trials(Bad_electrode_visual_electrodes(bad_elec_i),:,:)=nan;
                                RAM_all_trials_sessions{session_i}=all_trials;
                            catch
                            end
                        end
                    else
                        try
                            all_trials=[];
                            all_trials=RAM_all_trials_sessions{flag};
                            all_trials(Bad_electrode_visual_electrodes(bad_elec_i),:,:)=nan;
                            RAM_all_trials_sessions{flag}=all_trials;
                        catch
                        end 
                    end
                    
                end
            end
            
            % remove bad electrodes & sessions (visual) ALL
            if ~isempty(Bad_electrode_visual_session_n)
                all_trials=[];
                all_trials=RAM_all_trials_sessions{Bad_electrode_visual_session_n};
                all_trials(:,:,:)=nan;
                RAM_all_trials_sessions{Bad_electrode_visual_session_n}=all_trials;
            end

            % artifact rejection
            node_n=size(location_data_pair,2);
            noisy_RAM_all_trials_idx=[];
            good_RAM_all_trials_idx=[];
            for session_i=1:session_n
                all_trials=[];
                all_trials=RAM_all_trials_sessions{session_i};
                try
                    for node_i=1:node_n
                        all_trials_node=[];
                        all_trials_node=squeeze(all_trials(node_i,actual_time,:));
                        
%                         buf_index = zscore(var(all_trials_node))>2 | ...
%                             zscore(kurtosis(all_trials_node))>2 | ...
%                             zscore(rms(all_trials_node))>2;
                        buf_index=mean(all_trials_node)<vol_thr & mean(all_trials_node)>-vol_thr;
                        

                        noisy_idx_buffer = find(buf_index==1);
                        good_idx_buffer = find(buf_index==0);
                        noisy_RAM_all_trials_idx{node_i,session_i} = noisy_idx_buffer;
                        good_RAM_all_trials_idx{node_i,session_i} = good_idx_buffer;
                    end
                catch
                end
            end
            save('rejection_idx_extended_all_sessions.mat','noisy_RAM_all_trials_idx','good_RAM_all_trials_idx')

            % power extraction & statistics
            zscore_traces=NaN(size(location_data_pair,2),1,4,39);
%             t_traces=NaN(size(location_data_pair,2),1,4,39);
%             p_traces=NaN(size(location_data_pair,2),1,4,39);

            final_rejection_rate_all=NaN(node_n,1);
            
            for node_i=1:node_n
                try
                    new_node_i=find(Atlas_subject_elec_val==node_i);
                    if (ismember(node_i,Atlas_subject_elec_val))
                        z_trialpowers_window_all_sessions=[];
                        trialpowers_window_all_sessions=[];
                        all_trialpowers_window=[];
                        total_trial_size=[];

                        for session_i=1:session_n
                            try
                                all_trials=[];
                                all_trials=RAM_all_trials_sessions{session_i};
                                channel_epochs=squeeze(all_trials(node_i,:,:));
                                total_trial_size=[total_trial_size size(channel_epochs,2)];

                                % remove bad trials 
                                channel_epochs(:,noisy_RAM_all_trials_idx{node_i,session_i})=[];

    %                             z_channel_epochs = [];
    %                             z_channel_epochs = (channel_epochs(actual_time,:)-nanmean(nanmean(channel_epochs(actual_time,:))))./nanstd(nanmean(channel_epochs(actual_time,:)));
    %                             z_channel_epochs_all = [z_channel_epochs_all z_channel_epochs];

                                fERP = [];
                                for k = 1:size(channel_epochs,2)
                                    data = squeeze(channel_epochs(:,k));
                                    data_buf = filtfilt(d_filter,data);
                                    fERP(:,k)= data_buf(actual_time,:);
                                end
                                fpower = fERP.^2;

    %                             z_channel_epochs_power = [];
    %                             z_channel_epochs_power = (fpower-nanmean(nanmean(fpower)))./nanstd(nanmean(fpower));
    %                             z_channel_epochs_power_all = [z_channel_epochs_power_all z_channel_epochs_power];

                                trialpowers_window=[];
                                for k = 1:size(fpower,2)
                                    fdata_new_trial=fpower(:,k);
                                    for window_i=1:39
                                        trialpowers_window(window_i,k)=nanmean(fdata_new_trial(samples(window_i):samples(window_i)+window_size-1,1));
                                    end
                                end

                                z_trialpowers_window=(trialpowers_window-mean(mean(trialpowers_window)))./std(mean(trialpowers_window));
                                z_trialpowers_window_all_sessions=[z_trialpowers_window_all_sessions z_trialpowers_window];

                                trialpowers_window_all_sessions=[trialpowers_window_all_sessions trialpowers_window];

                            catch
                            end

                        end

    %                     z_trialpowers_window_all_sessions;
    %                     figure;plot(z_trialpowers_window_all_sessions)
    %                     figure;plot(mean(z_trialpowers_window_all_sessions,2))
    %                     
    %                     trialpowers_window_all_sessions;
    %                     figure;plot(trialpowers_window_all_sessions)
    %                     figure;plot(mean(trialpowers_window_all_sessions,2))

                        % baseline z-score
                        mean_trialpowers_window_all_sessions=[];
                        mean_trialpowers_window_all_sessions=nanmean(trialpowers_window_all_sessions,2);
                        std_of_mean=std(mean_trialpowers_window_all_sessions(1:base_sample_idx));
                        mean_of_mean=mean(mean_trialpowers_window_all_sessions(1:base_sample_idx));

                        z_mean_trialpowers_window_all_sessions=[];
                        z_mean_trialpowers_window_all_sessions=(mean_trialpowers_window_all_sessions-mean_of_mean)/std_of_mean;
                        zscore_traces(node_i,1,band_i,:)=z_mean_trialpowers_window_all_sessions';
                        final_rejection_rate_all(node_i,1)=size(trialpowers_window_all_sessions,2)/sum(total_trial_size);
%                         %% plot original voltage
%                         figure('position',[0 0 4000 4000]);
%                         set(gcf, 'color', [1 1 1]);
%                         set(gcf,'Visible','on');  
% 
%                         plot(z_mean_trialpowers_window_all_sessions,'r');
%                         set(gca,'xtick',[1 10 20 30 39]);
%                         set(gca,'xticklabel',[-500 0 500 1000 1500]);
%                         set(gca,'xlim',[1 39])
%                         set(gca,'ylim',[-8 8])
%                         line([1,39],[2,2],'color','k');
%                         line([1,39],[-2,-2],'color','k');
%                         line([10,10],[-8,8],'color','k');
%                         set(gca,'Fontsize',12);
%                         xlabel('Time (ms)');
%                         legend({'Gamma Power z-score'},'Location','northwest');
%                         
%                         saveas(gcf,['test_plot_' Network_name_all{Atlas_subject_network_val(new_node_i)} '_ROI_' num2str(Atlas_subject_ROI_val(new_node_i)) '_elec_' num2str(node_i) '.png'])
%                         close all;

                    end
                catch
                    disp(sprintf(' %s subject : no value at %d elec',r_sublist{i,1},node_i));
                end
            end
            
%             flat_index_ratio=mean(flat_index_ratio_all,2);
            patient_filename = ['./stats_traces_' prefix '.mat'];
            save(patient_filename,'zscore_traces','Atlas_subject_elec_val','Atlas_subject_ROI_val','Atlas_subject_network_val','final_rejection_rate_all','location_data_pair','T')
            disp(sprintf(' %s subject completed!! ',r_sublist{i,1}));
        catch

        end
    end
end





%% plot network
clc;clear;
rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

load r1_all.mat
fid=fopen('Subjects_list_all.txt','r');
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);


% Network_name_all={'MedialFrontal','FrontoParietal','DefaultMode' ... 
%     ,'SubcorticalCerebellum','Motor','Visual1','Visual2','VisualAssociation'};

Network_name_all={'FrontoParietal(FEF)','FrontoParietal(IPL)','FrontoParietal(MTG)' ...
    ,'DefaultMode', 'VisualPrimary','VisualAssociation'};

prefix='final_clean_vol_thr_10';
Image_folder=['E:\RAM data set\RAM_Public_Data_all\FR1_final\Results_network_all_' prefix '_test'];
mkdir(Image_folder);

all_mean_network_zval=nan(251,6,39);
base_sample_idx=9;

for j=1:1 % FR1 or FR2 or FR3

    for i=1:251 % 251 subject
        try
            cd(rootfolder)
            cd('FR1_final')
            cd([num2str(i),'_',r_sublist{i,1}]);
            load(['stats_traces_' prefix '.mat'])  
            
            % threshold
            bad_elec_index = final_rejection_rate_all<0.7; 
            zscore_traces(bad_elec_index,1,3,:)=NaN;

            %z-score
            figure('position',[0 0 8000 4000]);
            set(gcf, 'color', [1 1 1]);
            set(gcf,'Visible','off');  
            
            for network_i=1:6
                subplot(2,3,network_i)
                
                network_electrodes_index=Atlas_subject_elec_val(find(Atlas_subject_network_val==network_i));
                if (network_electrodes_index)
                    network_zval=squeeze(zscore_traces(network_electrodes_index,1,3,:))';
                else
                    network_zval=[];
                end
                plot(network_zval);
                if(size(network_zval,1)>1)
                    mean_network_zval=nanmean(network_zval,2);
                    mean_network_zval=(mean_network_zval-mean(mean_network_zval(1:base_sample_idx)))./std(mean_network_zval(1:base_sample_idx));
                    all_mean_network_zval(i,network_i,:)=mean_network_zval;
                else
                    mean_network_zval=network_zval;
                end
                hold on;
%                 plot(mean_vals,'LineWidth',2,'color','k');
                plot(mean_network_zval,'LineWidth',2,'color','k');
                set(gca,'xtick',[10 20 30]);
                set(gca,'xticklabel',[0 500 1500]);
                set(gca,'xlim',[1 39])
                set(gca,'ylim',[-8 8])
                line([1,39],[2,2],'color','b');
                line([1,39],[-2,-2],'color','b');
                line([10,10],[-8,8],'color','k');
                set(gca,'Fontsize',15);
                title([Network_name_all{network_i}],'Fontsize',12)
                ylabel('Z score');
                xlabel('Time (ms)');

            end
            saveas(gcf,[Image_folder '\all_' r_sublist{i,1} '_zval.png'])
            close all;
            
        catch
        end
    end
end

save('all_mean_network_val.mat','all_mean_network_tval','all_mean_network_zval')
            
%% all plot (ttest)
% cd(rootfolder)
% cd('FR1_final')
% 
% figure('position',[0 0 8000 4000]);
% set(gcf, 'color', [1 1 1]);
% set(gcf,'Visible','on');  
% for network_i=1:6
%     subplot(2,3,network_i)
%     
%     test_all_mean_network_zval=squeeze(all_mean_network_tval(:,network_i,:));
%     test_all_mean_network_zval_base=squeeze(nanmean(all_mean_network_tval(:,network_i,1:9),3));
%     [H,P,CI,STATS]=ttest(test_all_mean_network_zval - test_all_mean_network_zval_base);
%     plot(squeeze(STATS.tstat),'LineWidth',2,'color','k');
%     set(gca,'xtick',[10 20 30]);
%     set(gca,'xticklabel',[0 500 1500]);
%     set(gca,'xlim',[1 39])
%     set(gca,'ylim',[-8 8])
%     line([1,39],[2,2],'color','b');
%     line([1,39],[-2,-2],'color','b');
%     line([10,10],[-8,8],'color','k');
%     set(gca,'Fontsize',15);
%     title([Network_name_all{network_i}],'Fontsize',12)
%     ylabel('T value');
%     xlabel('Time (ms)');
% 
% end
% saveas(gcf,['network_plot_tval_from_tval_all_clean_sum.png'])
% close all;


cd(rootfolder)
cd('FR1_final')

figure('position',[0 0 8000 4000]);
set(gcf, 'color', [1 1 1]);
set(gcf,'Visible','on');  
for network_i=1:6
    subplot(2,3,network_i)
    
    test_all_mean_network_zval=squeeze(all_mean_network_zval(:,network_i,:));
%     test_all_mean_network_zval_base=squeeze(nanmean(all_mean_network_zval(:,network_i,1:9),3));
%     [H,P,CI,STATS]=ttest(test_all_mean_network_zval - test_all_mean_network_zval_base);
%     plot(squeeze(STATS.tstat),'LineWidth',2,'color','k');
    (mean_network_zval-mean(mean_network_zval(1:9)))./std(mean_network_zval(1:9));
    mean_test_all_mean_network_zval=nanmean(test_all_mean_network_zval);
    
    plot(squeeze(),'LineWidth',2,'color','k');
    set(gca,'xtick',[10 20 30]);
    set(gca,'xticklabel',[0 500 1500]);
    set(gca,'xlim',[1 39])
    set(gca,'ylim',[-8 8])
    line([1,39],[2,2],'color','b');
    line([1,39],[-2,-2],'color','b');
    line([10,10],[-8,8],'color','k');
    set(gca,'Fontsize',15);
    title([Network_name_all{network_i}],'Fontsize',12)
%     ylabel('T value');
    ylabel('z score');
    xlabel('Time (ms)');

end
saveas(gcf,['network_plot_tval_from_zval_all_clean_sum.png'])
close all;


