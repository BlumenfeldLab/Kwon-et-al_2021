function GammaPowerExtraction_preprocessing(i,subject_name)
% i=subject_id; subject_name='R1001P';
%% Setting
% Local computer
addpath(genpath('E:\RAM data set\RAM_Public_Data_all/Codes')); 
protocol_path='E:/RAM data set/RAM_Public_Data_all/';
current_path='E:/RAM data set/RAM_Public_Data_all/FR1_FS';

prefix='Final_std_20_fs';

data_path=[current_path '/' num2str(i) '_' subject_name];
% mkdir(data_path)

resampling_rate=500;
Timebuffer=1000; % 
Timebaseline_length=500;
Timetask_length=2500;
baseline_samples=Timebuffer*resampling_rate/1000+1:Timebuffer*resampling_rate/1000+Timebaseline_length*resampling_rate/1000;
task_samples=Timebuffer*resampling_rate/1000+1:(Timebuffer+Timebaseline_length+Timetask_length)*resampling_rate/1000;

%% save RAM data
cd(current_path)
load r1_all.mat
trialtypes = {'RAM_Rec_trials','RAM_nonRec_trials'};
task='FR1';
session={'x0x30_' 'x0x31_' 'x0x32_' 'x0x33_' 'x0x34_' 'x0x35_' 'x0x36_' 'x0x37_' 'x0x38_' 'x0x39_'};
session_folder={'0' '1' '2' '3' '4' '5' '6' '7' '8' '9'}; % session number

RAM_all_trials_session=[];
recInd_session=[];
NonrecInd_session=[];
sampling_rate_session=[];
real_k=1;

for k=1:10 % sessions
    try
        cd(protocol_path)
        eval(sprintf('events_ram=loadjson(r1_list.protocols.r1.subjects.%s.experiments.%s.sessions.%s.all_events);',subject_name,task,session{k}))              
        eval(sprintf('location_pairs=loadjson(r1_list.protocols.r1.subjects.%s.experiments.%s.sessions.%s.pairs);',subject_name,task,session{k}))
        eval(sprintf('location_contacts=loadjson(r1_list.protocols.r1.subjects.%s.experiments.%s.sessions.%s.contacts);',subject_name,task,session{k}))
        
        clearvars events_str
        for l=1:length(events_ram)
            events_str(l)=events_ram{l};
        end

        Events_word = filterStruct(events_str,'strcmp(type,''WORD'')');
        recInd = inStruct(Events_word,"recalled == 1"); %Recalled trials
        NonrecInd = inStruct(Events_word,"recalled == 0"); % Non-recalled trials
        
        recInd_session{real_k}=recInd;
        NonrecInd_session{real_k}=NonrecInd;
        
        subjects_directory = [protocol_path '/protocols/r1/subjects/'];
        
        cd(subjects_directory)
        cd([subject_name '/experiments/' task '/sessions/' session_folder{k} '/ephys/current_processed/']);
        sources_info=loadjson('sources.json');
        eegfile_name=Events_word.eegfile;
        eval(['sampling_rate=sources_info.' eegfile_name '.sample_rate;'])
        sampling_rate_session{real_k}=sampling_rate;
        
        cd(subjects_directory)
        cd([subject_name '/experiments/' task '/sessions/' session_folder{k} '/ephys/current_processed/noreref']);

        eval(sprintf('location_data=struct2dataset(location_contacts.%s.contacts);',subject_name))
        eval(sprintf('location_data_pair=struct2dataset(location_pairs.%s.pairs);',subject_name)) % bipolar
        
        
        %% ERP extraction
        for ii=1:size(location_data_pair,2)
            try

                channel_inf_pair=location_data_pair{1,ii};
                
                [EEG1] = gete_ms_new(channel_inf_pair.channel_1 , Events_word , 5000 ,-1500 , Timebuffer, ...
                    [48 52 ; 58 62], 'stop',1,resampling_rate);
                [EEG2] = gete_ms_new(channel_inf_pair.channel_2 , Events_word , 5000 ,-1500 , Timebuffer, ...
                    [48 52 ; 58 62], 'stop',1,resampling_rate);
                
                EEG = EEG2-EEG1; % difference from 
                
                %baselining option
                MU_EEG_all=nanmean(EEG(:,baseline_samples),2);
                EEG = EEG - MU_EEG_all;
                clear MU_EEG
                
                RAM_all_trials_session{ii,real_k}=EEG';

                disp(['ERP extraction : ' num2str(ii) '/' num2str(size(location_data_pair,2))]);
            catch
                disp(['ERP extraction, failed : ' num2str(ii) '/' num2str(size(location_data_pair,2))]);
            end
        end
        
        disp(sprintf('%s subject, %s Task, %s session completed!!',subject_name,task,session_folder{k}));
        real_k=real_k+1;
    catch
    end
end

% save data
cd(data_path)
[ii,k]=size(RAM_all_trials_session);
RAM_all_trials_together=[];
recInd_together=[];
NonrecInd_together=[];

for new_k=1:k
    RAM_all_trials_together_buffer=[];
    for new_ii=1:ii
        RAM_all_trials_together_buffer(new_ii,:,:)= RAM_all_trials_session{new_ii,new_k};
    end
    RAM_all_trials_together=cat(3,RAM_all_trials_together,RAM_all_trials_together_buffer);
    recInd_together=[recInd_together  recInd_session{new_k}];
    NonrecInd_together=[NonrecInd_together  NonrecInd_session{new_k}];
    subj_session_info(new_k).Rec_total=size(find(recInd_session{new_k}==1),2);
    subj_session_info(new_k).NRec_total=size(find(NonrecInd_session{new_k}==1),2);
    subj_session_info(new_k).All_total=size(RAM_all_trials_together_buffer,3);
    subj_session_info(new_k).sampling_rate=sampling_rate_session{new_k};
end
recInd_together=recInd_together>0;
NonrecInd_together=NonrecInd_together>0;

if ~isempty('RAM_all_trials_together')
    save([data_path '/Raw_data_extended_5s_session.mat'],'RAM_all_trials_together','RAM_all_trials_session','recInd_together','recInd_session','NonrecInd_together','NonrecInd_session','subj_session_info','location_data_pair','location_data','-v7.3')
end

%% artifact rejection
cd(data_path)

[node_n,k]=size(RAM_all_trials_session);
noisy_RAM_all_trials_number=nan(node_n,1);
RAM_all_trials_total_number=size(RAM_all_trials_together,3);
thr_i=0.2;
for node_i=1:node_n
    for new_k=1:k
        try
            all_trials_node=[];
            all_trials_node=squeeze(RAM_all_trials_session{node_i,new_k});
            noisy_idx_buffer=[];
            [~,noisy_idx_buffer]=artifact_rejection_all_together(all_trials_node,resampling_rate,thr_i);                        
            noisy_idx_buffer=unique(noisy_idx_buffer);
    %         histogram(meansqerr,200,'BinLimits',[0,20000],'Normalization', 'probability')
            noisy_RAM_all_trials_idx_session{node_i,new_k} = noisy_idx_buffer;
            disp(sprintf('electrode #%d : artifact rejection is completed! ',node_i));
        catch
            disp(sprintf('electrode #%d : artifact rejection is failed! ',node_i));
        end
    end
end
save(['rejection_idx_' prefix '_session.mat'],'noisy_RAM_all_trials_idx_session','-v7.3')
% load(['rejection_idx_' prefix '_session.mat'])

%% Power extraction from filtering method
cd(data_path)

[node_n,session_n]=size(RAM_all_trials_session);

% design filter
frequency_bands=[45 95; 3 8;40 115;13 30];
band_i=3;
d_filter= designfilt('bandpassiir','FilterOrder',40, ...
'HalfPowerFrequency1',frequency_bands(band_i,1),'HalfPowerFrequency2',frequency_bands(band_i,2), ...
'SampleRate',resampling_rate);

load('time_50ms_overlap.mat')
window_size=round(resampling_rate*0.05); % 50ms
samples=round(T*resampling_rate);
samples(1)=1;
base_sample_idx=19;

zscore_traces_together=NaN(size(location_data_pair,2),1,4,size(T,2));
zscore_traces_rec=NaN(size(location_data_pair,2),1,4,size(T,2));
zscore_traces_Nonrec=NaN(size(location_data_pair,2),1,4,size(T,2));
zscore_traces_diff=NaN(size(location_data_pair,2),1,4,size(T,2));

% noisy_RAM_all_trials=nan(node_n,2);
for node_i=1:node_n
    try
        trialpowers_window=[];
        for new_k=1:session_n
            try
                channel_epochs=squeeze(RAM_all_trials_session{node_i,new_k});
                fpower = [];
                for k = 1:size(channel_epochs,2)
                    data = squeeze(channel_epochs(:,k));
                    data_buf = filtfilt(d_filter,data);
                    % power
                    data_buf=data_buf.^2;
                    fpower(:,k)= data_buf(task_samples,:);
        %                 trialpowers_window(:,k)=nanmean(reshape(data_buf(task_samples,:),[fs*window_size,2*1/window_size]),1);
                end

                % remove bad trials 
                fpower(:,noisy_RAM_all_trials_idx_session{node_i,new_k})=nan;

                % averaging power by window
                trialpowers_window_buf=[];
                for k = 1:size(fpower,2)
                    fdata_new_trial=fpower(:,k);
                    for window_i=1:size(T,2)
                        trialpowers_window_buf(window_i,k)=nanmean(fdata_new_trial(samples(window_i):samples(window_i)+window_size-1,1));
                    end
                end
%                 % additional step
%                 trialpowers_window_buf=trialpowers_window_buf-repmat(nanmean(trialpowers_window_buf(1:19,:)),size(T,2),1);
                
                %z-scroing
                MU_POW = nanmean(nanmean(trialpowers_window_buf(1:base_sample_idx,:),2));
                STD_POW = nanstd(nanmean(trialpowers_window_buf(1:base_sample_idx,:),2));

                trialpowers_window=[trialpowers_window (trialpowers_window_buf-MU_POW)./STD_POW];

                clear MU_POW STD_POW
%                 disp(sprintf('electrode #%d, session #%d : power extraction is completed! ',node_i,new_k));

            catch
%                 disp(sprintf('electrode #%d, session #%d : power extraction is not completed! ',node_i,new_k));
            end
        end

        trialpowers_window_buffer=trialpowers_window;

        %% method 1
        std_num =20;
        mean_mean_all_trials = nanmean(trialpowers_window_buffer(:));
        std_mean_all_trials = nanstd(trialpowers_window_buffer(:));
    %         mean_mean_all_trials = nanmean(nanmean(trialpowers_window_buffer));
    %         std_mean_all_trials = nanstd(nanmean(trialpowers_window_buffer));

        trialpowers_window_buffer(trialpowers_window_buffer > mean_mean_all_trials + std_num*std_mean_all_trials | trialpowers_window_buffer < mean_mean_all_trials - std_num*std_mean_all_trials) = NaN;

        noisy_trials_power_idx=[];
        for h = 1:size(trialpowers_window_buffer,2)
            if ~isempty(find(isnan(trialpowers_window_buffer(:,h)),1))       
                noisy_trials_power_idx = [noisy_trials_power_idx h];       
            end
        end
        trialpowers_window_buffer(:,noisy_trials_power_idx) = nan;

        %% save power
        trialpowers_window=[];
        trialpowers_window=trialpowers_window_buffer;

        recInd_together_noise=recInd_together(noisy_trials_power_idx);
        NonrecInd_together_noise=NonrecInd_together(noisy_trials_power_idx);

        noisy_RAM_all_trials(node_i).All_total=size(trialpowers_window,2);
        noisy_RAM_all_trials(node_i).All_noise=length(noisy_trials_power_idx);
        noisy_RAM_all_trials(node_i).All_ratio=length(noisy_trials_power_idx)/size(trialpowers_window,2);

        noisy_RAM_all_trials(node_i).Rec_total=length(find(recInd_together>0));
        noisy_RAM_all_trials(node_i).Rec_noise=length(find(recInd_together_noise>0));
        noisy_RAM_all_trials(node_i).Rec_ratio=length(find(recInd_together_noise>0))/length(find(recInd_together>0));

        noisy_RAM_all_trials(node_i).NRec_total=length(find(NonrecInd_together>0));
        noisy_RAM_all_trials(node_i).NRec_noise=length(find(NonrecInd_together_noise>0));
        noisy_RAM_all_trials(node_i).NRec_ratio=length(find(NonrecInd_together_noise>0))/length(find(NonrecInd_together>0));

        min_noisy_ratio=min([noisy_RAM_all_trials(node_i).NRec_ratio,noisy_RAM_all_trials(node_i).Rec_ratio, ...
            noisy_RAM_all_trials(node_i).All_ratio]);

        if min_noisy_ratio < 0.5
            zscore_traces_together(node_i,1,band_i,:)=nanmean(trialpowers_window,2);
            zscore_traces_rec(node_i,1,band_i,:)=nanmean(trialpowers_window(:,recInd_together),2);
            zscore_traces_Nonrec(node_i,1,band_i,:)=nanmean(trialpowers_window(:,NonrecInd_together),2);
            zscore_traces_diff(node_i,1,band_i,:)=nanmean(trialpowers_window(:,recInd_together),2) - ...
                nanmean(trialpowers_window(:,NonrecInd_together),2);

            disp(sprintf('electrode #%d : power extraction is completed! ',node_i));
        else
            disp(sprintf('electrode #%d : noise electrode',node_i));
        end
        
    catch
        disp(sprintf('electrode #%d : bad electrode',node_i));
    end
end

all_names={'together','rec','Nonrec','diff'};

for name_i=1:4
    zscore_traces=[];
    zscore_traces=eval(['zscore_traces_' all_names{name_i}]);
    patient_filename = ['./stats_traces_' prefix '_' all_names{name_i} '.mat'];
    save(patient_filename,'zscore_traces','noisy_RAM_all_trials','location_data_pair','T','-v7.3')
    patient_filename_xflip = ['./stats_traces_' prefix '_' all_names{name_i} '_xflip.mat'];
    save(patient_filename_xflip,'zscore_traces','noisy_RAM_all_trials','location_data_pair','T','-v7.3')
end

end