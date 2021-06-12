
% artifact rejection
function RAM_FARNAM(i,subject_name)

addpath(genpath('/home/hk582/project/Codes'));
cd([num2str(i),'_',subject_name]);
load('Raw_data_extended_4s_all.mat')

actual_time=257:768;
node_n=size(location_data_pair,2);
session_n=size(RAM_all_trials_sessions,2);
prefix='Wendy_thr_01';
% prefix='Wendy_thr_02';

fs=256;
load('time_50ms_no_overlap.mat')
window_size=round(fs*0.05); % 50ms
samples=round(T*fs);

% load('time_63ms_no_overlap.mat')
% window_size=round(fs*0.063); % 50ms
% samples=round(T*fs);

base_sample_idx=9;


frequency_bands=[45 95; 3 8;40 115;13 30];
% design filter
band_i=3;
d_filter= designfilt('bandpassiir','FilterOrder',40, ...
'HalfPowerFrequency1',frequency_bands(band_i,1),'HalfPowerFrequency2',frequency_bands(band_i,2), ...
'SampleRate',fs);

% % artifact rejection
% for thr_i=[0.1,0.2]
%     noisy_RAM_all_trials_number=nan(node_n,session_n);
%     RAM_all_trials_total_number=nan(1,session_n);
%     for session_i=1:session_n
%         try
%             all_trials=[];
%             all_trials=RAM_all_trials_sessions{session_i};
%             RAM_all_trials_total_number(1,session_i)=size(all_trials,3);
%             for node_i=1:node_n
%                 try
%                     all_trials_node=[];
%                     all_trials_node=squeeze(all_trials(node_i,actual_time,:));
%                     noisy_idx_buffer=artifact_rejection_all_together(all_trials_node,fs,thr_i);                        
%                     noisy_RAM_all_trials_idx{node_i,session_i} = noisy_idx_buffer;
%                     noisy_RAM_all_trials_number(node_i,session_i) = length(unique(noisy_idx_buffer));
%                     disp(sprintf('session #%d , electrode #%d : artifact rejection is completed! ',session_i,node_i));
%                 catch
%                     disp(sprintf('session #%d , electrode #%d : artifact rejection is failed! ',session_i,node_i));
%                 end
%             end
%         catch
%         end
%     end
%     if thr_i==0.1
%         save(['rejection_idx_all_sessions_' prefix '_thr_01.mat'],'noisy_RAM_all_trials_idx',...
%             'noisy_RAM_all_trials_number','RAM_all_trials_total_number','-v7.3')
%     elseif thr_i==0.2
%         save(['rejection_idx_all_sessions_' prefix '_thr_02.mat'],'noisy_RAM_all_trials_idx',...
%             'noisy_RAM_all_trials_number','RAM_all_trials_total_number','-v7.3')
%     end
% end
load(['rejection_idx_all_sessions_' prefix '.mat'])


zscore_traces=NaN(size(location_data_pair,2),1,4,39);
for node_i=1:node_n
    try
        trialpowers_window_all_sessions=[];
        z_trialpowers_window_all_sessions=[];
        for session_i=1:session_n
            try
                all_trials=[];
                all_trials=RAM_all_trials_sessions{session_i};
                channel_epochs=squeeze(all_trials(node_i,:,:));
                if (noisy_RAM_all_trials_number(node_i,session_i)/RAM_all_trials_total_number(1,session_i)<0.2)
                
                    % remove bad trials 
                    channel_epochs(:,noisy_RAM_all_trials_idx{node_i,session_i})=[];

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
                    trialpowers_window_all_sessions=[trialpowers_window_all_sessions trialpowers_window];

                    %                     % zscoring
%                     mean_trialpowers_window = nanmean(trialpowers_window);                    
%                     z_trialpowers_window=(trialpowers_window-nanmean(mean_trialpowers_window))./nanstd(mean_trialpowers_window);
%                     z_trialpowers_window_all_sessions=[z_trialpowers_window_all_sessions z_trialpowers_window];
   
                end
            catch
            end

        end

        % baseline z-score & remove artifact
        mean_trialpowers_window_all_sessions=[];
        mean_trialpowers_window_all_sessions=nanmean(trialpowers_window_all_sessions,2);
        std_of_mean=std(mean_trialpowers_window_all_sessions);
        mean_of_mean=mean(mean_trialpowers_window_all_sessions);        
        new_thr=40*std_of_mean+mean_of_mean;
        
        bad_power_index=sum(trialpowers_window_all_sessions>new_thr)>0;
        trialpowers_window_all_sessions(:,bad_power_index)=[];
        
%         mean_z_trialpowers_window_all_sessions=[];
%         mean_z_trialpowers_window_all_sessions=nanmean(z_trialpowers_window_all_sessions,2);
%         std_of_mean=std(mean_z_trialpowers_window_all_sessions);
%         mean_of_mean=mean(mean_z_trialpowers_window_all_sessions);
        
        %z-scoring
        mean_trialpowers_window_all_sessions=[];
        mean_trialpowers_window_all_sessions=nanmean(trialpowers_window_all_sessions,2);
        base_std_of_mean=std(mean_trialpowers_window_all_sessions(1:base_sample_idx));
        base_mean_of_mean=mean(mean_trialpowers_window_all_sessions(1:base_sample_idx));
        z_mean_trialpowers_window_all_sessions=[];
        z_mean_trialpowers_window_all_sessions=(mean_trialpowers_window_all_sessions-base_mean_of_mean)/base_std_of_mean;
        
%         set(gcf, 'color', [1 1 1]);
%         set(gcf,'Visible','off');  
%         subplot(2,1,1)
%         plot(trialpowers_window_all_sessions)
%         ylim([0 50])
%         subplot(2,1,2)
%         plot(z_mean_trialpowers_window_all_sessions)
%         ylim([-6 8])
%         saveas(gcf,['./testtesttest_'  num2str(node_i) '_40.png'])
%         close all;
        
        zscore_traces(node_i,1,band_i,:)=z_mean_trialpowers_window_all_sessions';
        disp(sprintf('electrode #%d : power extraction is completed! ',node_i));
    catch
        disp(sprintf('electrode #%d : power extraction is not completed! ',node_i));
    end
end

patient_filename = ['./stats_traces_' prefix '.mat'];
save(patient_filename,'zscore_traces','Atlas_subject_elec_val','location_data_pair','T','-v7.3')



% % Make the Montage
% M=[0.9975 -0.0073 0.0176 -0.0429 ; 
%    0.0146 1.0009 -0.0024 1.5496 ; 
%   -0.0130 -0.0093 0.9971 1.1840]; % https://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems
% 
% load('stats_traces_Wendy.mat')

% masked_elec_n=size(Atlas_subject_elec_val,2);
% left_index=1;
% right_index=1;
% 
% L_MontageMap=[];
% R_MontageMap=[];
% Montage_prefix='final_clean_dil';
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


overlay=4;
inflationstep=5;
file_suffix=prefix;

displayElectrodesInflated_new_FARNAM(overlay, inflationstep ,file_suffix) %, frames_set)

end


