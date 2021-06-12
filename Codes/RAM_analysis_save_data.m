%% save raw
clc;clear;

% r1_list=loadjson('r1.json');
% r1_list=loadjson('r1_fullrelease.json');

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

% count_new=0;
% all_zeeg_rec=[];
% all_zeeg_nonrec=[];
% all_subjects_name_FR=[];

% real_count=0;

sampling_rate=nan(251,1);
%%raw files
for j=1:1 % FR1 or FR2 or FR3

    for i=1:251 % 251 subject
        try
            clearvars -except sampling_rate real_count rootfolder r1_list r_sublist task session session_folder all_subjects_name_FR count_new all_zeeg_rec all_zeeg_nonrec i j fs
            RAM_Rec_trials=[];
            RAM_nonRec_trials=[];

            for k=1:10 % sessions
                try
                    cd(rootfolder)
%                     eval(sprintf('events=loadjson(r1_list.protocols.r1.subjects.%s.experiments.%s.sessions.%s.task_events);',r_sublist{i,1},task{j},session{k}))
                    eval(sprintf('events=loadjson(r1_list.protocols.r1.subjects.%s.experiments.%s.sessions.%s.all_events);',r_sublist{i,1},task{j},session{k}))
%                     
                    eval(sprintf('location_pairs=loadjson(r1_list.protocols.r1.subjects.%s.experiments.%s.sessions.%s.pairs);',r_sublist{i,1},task{j},session{k}))
                    eval(sprintf('location_contacts=loadjson(r1_list.protocols.r1.subjects.%s.experiments.%s.sessions.%s.contacts);',r_sublist{i,1},task{j},session{k}))
% 
                    clearvars events_str;
                    for l=1:length(events)
                        events_str(l)=events{l};
                    end

                    Events_word = filterStruct(events_str,'strcmp(type,''WORD'')');
                    recInd = inStruct(Events_word,"recalled == 1");
                    NonrecInd = inStruct(Events_word,"recalled == 0");

                    %% ERP analysis pair

                    subjects_directory = 'E:\RAM data set\RAM_Public_Data_all\protocols\r1\subjects\';
                    cd(subjects_directory)
                    cd([r_sublist{i,1},'/experiments/',task{j},'/sessions/',session_folder{k},'/ephys/current_processed/']);
                    sources_info=loadjson('sources.json');
                    eegfile_name=events_str.eegfile;
                    eval(['sampling_buff=sources_info.' eegfile_name '.sample_rate;'])
                    sampling_rate(i,1)=sampling_buff;
                    
                    cd(subjects_directory)
                    cd([r_sublist{i,1},'/experiments/',task{j},'/sessions/',session_folder{k},'/ephys/current_processed/noreref']);
                    
                    eval(sprintf('location_data=struct2dataset(location_contacts.%s.contacts);',r_sublist{i,1}))
                    eval(sprintf('location_data_pair=struct2dataset(location_pairs.%s.pairs);',r_sublist{i,1}))

                    for ii=1:size(location_data_pair,2)
                        channel_inf_pair=location_data_pair{1,ii};
                        channel_inf=location_data{1,ii};
                        
    %                     channel_inf_pair.atlases.avg.x
    %                     channel_inf_pair.atlases.avg.y
    %                     channel_inf_pair.atlases.avg.z
    %                     channel_inf_pair.atlases.avg.region


                        %% Power analysis
                        freqs = logspace(log10(2),log10(200),60);
                        [~,pow]=getphasepow_bipolar_new(channel_inf_pair.channel_1,channel_inf_pair.channel_2, Events_word , 3000 ,-800 , 1000, ...
                            'filtfreq',[58 62],'filttype', 'stop','freqs',freqs,'width',6);
                        
                        MU_POW = nanmean(nanmean(pow(:,:,1:800),3));
                        STD_POW = nanstd(nanmean(pow(:,:,1:800),3));
                        
                        [n_ev,n_freq,n_time] = size(pow);
                        pat = [];
                        for f = 1:n_freq
                            all_mu = repmat(MU_POW(f),n_ev,n_time);
                            all_std = repmat(STD_POW(f),n_ev,n_time);
                            pat(:,f,:) = (squeeze(pow(:,f,:))-all_mu)./all_std;
                        end
                        
                        figure;
                        imagesc(flipud(squeeze(nanmean(pat(recInd,:,:)))));
                        F = sort(round(freqs),'descend');
                        set(gca,'ytick',1:10:60);
                        set(gca,'yticklabel',F(1:10:60));
                        caxis([-1 1])
                        
                        figure;
                        imagesc(flipud(squeeze(nanmean(pat(NonrecInd,:,:)))));
                        F = sort(round(freqs),'descend');
                        set(gca,'ytick',1:10:60);
                        set(gca,'yticklabel',F(1:10:60));
                        caxis([-1 1])
                        
%                         figure;
%                         imagesc(flipud(squeeze(nanmean(pat(:,:,:)))));
%                         F = sort(round(freqs),'descend');
%                         set(gca,'ytick',1:10:60);
%                         set(gca,'yticklabel',F(1:10:60));
%                         caxis([-1 1])

                        figure;
                        plot(squeeze(nanmean(nanmean(pat(recInd,freqs>40 & freqs<115,:),2),1)))
                      
                        figure;
                        plot(squeeze(nanmean(nanmean(pat(NonrecInd,freqs>40 & freqs<115,:),2),1)))
                    
                    
                    
                    
                    end


                        
                        
                    %% ERP analysis
                        [EEG1] = gete_ms(channel_inf_pair.channel_1 , Events_word , 4000 , -1500 , 0, ...
                            [58 62], 'stop',1,fs);
                        [EEG2] = gete_ms(channel_inf_pair.channel_2 , Events_word , 4000 , -1500 , 0, ...
                            [55 65], 'stop',1,fs); % original
                        figure;
                        plot(squeeze(nanmean(EEG1(1,:),1)),'r'); hold on;

                        EEG = EEG2-EEG1;

                        MU_EEG = nanmean(nanmean(EEG(:,256:383)));
                        EEG = EEG - MU_EEG;
                        clear MU_EEG

                        RAM_Rec_trials{ii,k}=EEG(recInd,:)';
                        RAM_nonRec_trials{ii,k}=EEG(NonrecInd,:)';

                        figure('position',[0 0 1101 767],'visible','on');
                        set(gcf, 'color', [1 1 1]);
                        
                        subplot(3,1,1)
                        plot(squeeze(nanmean(EEG1(recInd,:),1)),'b'); hold on;
                        plot(squeeze(nanmean(EEG1(NonrecInd,:),1)),'r');
                        set(gca,'xtick',[1 256 384 512 768 1024]);
                        set(gca,'xticklabel',[-1500 -500 0 500 1500 2500]);
                        set(gca,'xlim',[1 1024])
                        set(gca,'Fontsize',15);
                        ylabel('Voltage (uV)');
                        xlabel('Time (ms)');
                        legend({'Recalled','Non recalled'});
                        
                        subplot(3,1,2)
                        plot(squeeze(nanmean(EEG2(recInd,:),1)),'b'); hold on;
                        plot(squeeze(nanmean(EEG2(NonrecInd,:),1)),'r');
                        set(gca,'xtick',[1 256 384 512 768 1024]);
                        set(gca,'xticklabel',[-1500 -500 0 500 1500 2500]);
                        set(gca,'xlim',[1 1024])
                        set(gca,'Fontsize',15);
                        ylabel('Voltage (uV)');
                        xlabel('Time (ms)');
                        legend({'Recalled','Non recalled'});
                        
                        subplot(3,1,3)
                        plot(squeeze(nanmean(EEG(recInd,:),1)),'b'); hold on;
                        plot(squeeze(nanmean(EEG(NonrecInd,:),1)),'r');
                        set(gca,'xtick',[1 256 384 512 768 1024]);
                        set(gca,'xticklabel',[-1500 -500 0 500 1500 2500]);
                        set(gca,'xlim',[1 1024])
                        set(gca,'ylim',[-50 50])
                        set(gca,'Fontsize',15);
                        ylabel('Voltage (uV)');
                        xlabel('Time (ms)');
                        legend({'Recalled','Non recalled'});
                        
                        saveas(gcf,[num2str(ii) '_' channel_inf_pair.code '.png'])
                        close all;
%                     end
                    disp(sprintf('%s subject, %s Task, %s session completed!!',r_sublist{i,1},task{j},session_folder{k}));
                catch
                end
            end

%             [ii,k]=size(RAM_Rec_trials);
%             RAM_Rec_trials_sessions=[];
%             for new_k=1:k
%                 RAM_Rec_trials_final=[];
%                 for new_ii=1:ii
%                     RAM_Rec_trials_final(new_ii,:,:)= RAM_Rec_trials{new_ii,new_k};
%                 end
%                 RAM_Rec_trials_sessions{new_k}=RAM_Rec_trials_final;
%             end
% 
%             [ii,k]=size(RAM_nonRec_trials);
%             RAM_nonRec_trials_sessions=[];
%             for new_k=1:k
%                 RAM_nonRec_trials_final=[];
%                 for new_ii=1:ii
%                     RAM_nonRec_trials_final(new_ii,:,:)= RAM_nonRec_trials{new_ii,new_k};
%                 end
%                 RAM_nonRec_trials_sessions{new_k}=RAM_nonRec_trials_final;
%             end
% 
% 
%             % save values  
%             if ~isempty(RAM_Rec_trials_final)
%                 cd(rootfolder)
%                 cd('FR1_final')
%                 mkdir([num2str(i),'_',r_sublist{i,1}]);
%                 cd([num2str(i),'_',r_sublist{i,1}]);
%                 save('Raw_data_extended_4s.mat','RAM_Rec_trials_sessions','RAM_nonRec_trials_sessions','location_data_pair','location_data','-v7.3')
%             end   
        
        catch
        
        end

    end

end


















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

base_sample_idx=9;


actual_time=257:768;
frequency_bands=[45 95; 3 8;40 115;13 30];

% design filter
band_i=3;
d_filter= designfilt('bandpassiir','FilterOrder',40, ...
'HalfPowerFrequency1',frequency_bands(band_i,1),'HalfPowerFrequency2',frequency_bands(band_i,2), ...
'SampleRate',fs);

for j=1:1 % FR1 or FR2 or FR3

    for i=1:251 % 251 subject
        try
            clearvars RAM_Rec_trials_sessions RAM_nonRec_trials_sessions location_data_pair

            cd(rootfolder)
            cd('FR1_final')
            cd([num2str(i),'_',r_sublist{i,1}]);
%             if exist('rejection_idx_extended_all.mat','file') == 2
%             
%             else
                load('Raw_data_extended_4s.mat')            
                zscore_meanpower_traces=NaN(size(location_data_pair,2),1,4,size(T,2));

                % merge
                session_n=size(RAM_Rec_trials_sessions,2);
                for session_i=1:session_n
                    RAM_all_trials_sessions{session_i}=cat(3,RAM_Rec_trials_sessions{session_i},RAM_nonRec_trials_sessions{session_i});
                end

                % artifact rejection
                node_n=size(location_data_pair,2);
                noisy_RAM_all_trials_idx=[];
                for session_i=1:session_n
                    all_trials=[];
                    all_trials=RAM_all_trials_sessions{session_i};
                    try
                        for node_i=1:node_n
                            all_trials_node=[];
                            all_trials_node=squeeze(all_trials(node_i,actual_time,:));
                            buf_index = abs(zscore(mean(all_trials_node)))>2 | ...
                                abs(zscore(kurtosis(all_trials_node)))>2;
                            noisy_idx_buffer = find(buf_index~=0);
                            noisy_RAM_all_trials_idx{node_i,session_i} = noisy_idx_buffer;
                        end
                    catch
                    end
                end
                save('rejection_idx_extended_all.mat','noisy_RAM_all_trials_idx')

                % power extraction
                for node_i=1:node_n
                    zpower_all_sessions=[];
                    zpower_mean_all_sessions=[];
                    zpower_mean_all_sessions2=[];
                    for session_i=1:session_n
                        try
                            all_trials=[];
                            all_trials=RAM_all_trials_sessions{session_i};
                            channel_epochs=squeeze(all_trials(node_i,:,:));
                            channel_epochs(:,noisy_RAM_all_trials_idx{node_i,session_i})=[];

                            fERP = [];
                            for k = 1:size(channel_epochs,2)
                                data = squeeze(channel_epochs(:,k));
                                data_buf = filtfilt(d_filter,data);
                                fERP(:,k)= data_buf(actual_time,:);
                            end
                            fpower = fERP.^2;

                            trialpowers_window=[];
                            for k = 1:size(fpower,2)
                                fdata_new_trial=fpower(:,k);
                                for window_i=1:size(T,2)
            %                              trialpowers_window(window_i,k)=nanmean(fdata_new_trial(samples(window_i)-window_size/2+1:samples(window_i)+window_size/2,1));
                                     trialpowers_window(window_i,k)=nanmean(fdata_new_trial(samples(window_i):samples(window_i)+window_size-1,1));
                                end
                            end
                            session_mean=nanmean(nanmean(trialpowers_window(1:base_sample_idx,:),1)); %% check edge effect
                            session_std=nanstd(nanmean(trialpowers_window(1:base_sample_idx,:),1));   %% check
                            [n_time n_ev]=size(trialpowers_window);
                            all_session_mean = repmat(session_mean,n_time,n_ev);
                            all_session_std = repmat(session_std,n_time,n_ev);

                            zpower_session=[];
                            zpower_session=(trialpowers_window-all_session_mean)./all_session_std;
                            zpower_all_sessions=[zpower_all_sessions  zpower_session];
        %                     all_mean=nanmean(trialpowers_window(1:base_sample_idx,:)); %% check edge effect
        %                     all_std=nanstd(trialpowers_window(1:base_sample_idx,:));   %% check
        %                     zpower=(trialpowers_window-all_mean)./all_std;
        %                     figure;
        %                     plot(mean(zpower'))
        %                     
        %                     zscore_meanpower_traces(node_i,1,band_i,:)=zpower;
        % 
        % 
                            mean_trialpowers_window=nanmean(trialpowers_window,2);
                            all_mean_mean=nanmean(mean_trialpowers_window(1:base_sample_idx,1)); %% check edge effect
                            all_mean_std=nanstd(mean_trialpowers_window(1:base_sample_idx,1));   %% check
                            zpower_mean_session=[];
                            zpower_mean_session=(mean_trialpowers_window-all_mean_mean)/all_mean_std;
                            zpower_mean_all_sessions=[zpower_mean_all_sessions  zpower_mean_session];
    %                     
    %                     zscore_meanpower_traces(node_i,1,band_i,:)=zpower;
                        catch
                        end

                    end
                    zscore_meanpower_traces(node_i,1,band_i,:)=nanmean(zpower_all_sessions,2);
                    zscore_mean_meanpower_traces(node_i,1,band_i,:)=nanmean(zpower_mean_all_sessions,2);
                end
                patient_filename = ['./meanpower_traces_zscore_filtering.mat'];
                save(patient_filename,'zscore_meanpower_traces','zscore_mean_meanpower_traces','T')
                disp(sprintf(' %s subject completed!! ',r_sublist{i,1}));
%             end
        catch
            
        end
    end
end


%% artifact rejection & power extraction
clc;clear;

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

fid=fopen('Subjects_list_all.txt','r');
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);

mni_nifti_path = 'C:\yale\bioimagesuite30\images\MNI_T1_1mm_stripped.nii';
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
Tmni = mnit1_1mm_stripped.transform;

electrode_volume=zeros(181,217,181);
subn='fs2mni';
% subn_atlas='broadmann';
subn_atlas='shen';

M=[0.9975 -0.0073 0.0176 -0.0429 ; 
   0.0146 1.0009 -0.0024 1.5496 ; 
  -0.0130 -0.0093 0.9971 1.1840]; % https://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems

% mask_nifti_path = 'E:\RAM data set\RAM_Public_Data_all\Atlases\yale_broadmann_cortical.nii';
atlas_nifti_path = 'E:\RAM data set\RAM_Public_Data_all\Atlases\shen_1mm_268_parcellation_cortical.nii';
atlas_1mm_image = mni2fs_load_nii(atlas_nifti_path);
atlas_1mm_image_data=atlas_1mm_image.img;
Atlas_index=unique(atlas_1mm_image_data);
Atlas_index(1,:)=[];
test=find(Atlas_index>0);
Atlas_index_n=size(Atlas_index,1);


for test_time_point=10:13
% Atlas_all_val=cell(Atlas_index_n,251);
% for i=1:251
%     try
%         clear zscore_mean_meanpower_traces zscore_meanpower_traces
%         cd(rootfolder)
%         cd('FR1_final')
%         cd([num2str(i),'_',r_sublist{i,1}]);
%         load('meanpower_traces_zscore_filtering.mat')
%         load('Raw_data_extended_4s.mat')
%         
%         for elec_i=1:size(location_data_pair,2)
%             channel_inf=location_data_pair{1,elec_i};
%             fsavg_coord=[channel_inf.atlases.avg.x channel_inf.atlases.avg.y channel_inf.atlases.avg.z]';
%             mni_coord = M*[fsavg_coord ; 1];
% 
%             volume_buffer=[mni_coord(1) mni_coord(2) mni_coord(3) 1] / Tmni';
%             volume_buffer=round(volume_buffer(1:3));
%             atlas_index_buffer=atlas_1mm_image_data(volume_buffer(1),volume_buffer(2),volume_buffer(3));
%             if atlas_index_buffer>0
%                 atlas_index_buffer_final=find(Atlas_index==atlas_index_buffer);
%                 Atlas_all_val{atlas_index_buffer_final,i}=[Atlas_all_val{atlas_index_buffer_final,i} zscore_mean_meanpower_traces(elec_i,1,3,test_time_point)];
%             end
%         end
% 
%     catch
%         i
%     end
%     
% end
cd(rootfolder)
load(['shen_val_all_time' num2str(test_time_point) '.mat'])   
% save(['shen_val_all_time' num2str(test_time_point) '.mat'],'Atlas_all_val','-v7.3')        

% Atlas_all_val_max=nan(188,251);
Atlas_all_val_mean=nan(188,251);
for i=1:251
    for j=1:Atlas_index_n
%         Atlas_single=Atlas_all_val{j,i};
%         [~ ,max_index]=max(abs(Atlas_single));
%         if max_index
%             Atlas_all_val_max(j,i)=Atlas_single(max_index);
%         end
        % mean
        Atlas_all_val_mean(j,i)=mean(Atlas_all_val{j,i});

    end
end

%% test
% [a b c d]=ttest(Atlas_all_val_max');
[a b c d]=ttest(Atlas_all_val_mean');


test_ii=find(b<0.05);
for k=1:size(test_ii,2)
    electrode_volume(atlas_1mm_image_data==Atlas_index(test_ii(k)))=Atlas_index(test_ii(k));
end

mnit1_1mm_stripped_new=fmris_read_nifti(mni_nifti_path);
mnit1_1mm_stripped_new.file_name=['test_shen_time' num2str(test_time_point) '_mean.nii'];
mnit1_1mm_stripped_new.precision='float32';
mnit1_1mm_stripped_new.byte=4;
mnit1_1mm_stripped_new.data=electrode_volume;
fmris_write_nifti(mnit1_1mm_stripped_new);

end




% Atlas_all_val_mean=nan(188,251);
Atlas_all_val_max=nan(188,251);
for i=1:251
    for j=1:Atlas_index_n
        Atlas_single=Atlas_all_val{j,i};
        [~ ,max_index]=max(abs(Atlas_single));
        if max_index
            Atlas_all_val_max(j,i)=Atlas_single(max_index);
        end
%         % mean
%         Atlas_all_val_mean(j,i)=mean(Atlas_all_val{j,i});

    end
end

%% test
[a b c d]=ttest(Atlas_all_val_max');
d.tstat
b(find(b<0.05))




% ans =
% 
% ROI index = 4    38    47    48    84    85    98   102   116   136   180
% Atlas index = 4    38    47    48    84    85   138   142   156   176   221
% pval : 0.0262    0.0334    0.0485    0.0389    0.0209    0.0230    0.0162    0.0262    0.0231    0.0482    0.0067
% tval : -2.4489    2.5091    2.0427    2.2630    2.5612    2.6057   -2.7091    2.4007   -2.3856   -2.2840   -3.3316

Atlas_index_roi=47;
roi_index=find(Atlas_index==Atlas_index_roi);


buffer=Atlas_all_val_max(roi_index,:);
buffer(~isnan(buffer))


nanmean(Atlas_all_val_max,2)

cd(rootfolder)
cd('Atlases')

for jj=1:size(Atlas_all_val_max,1)
    buffer=Atlas_all_val_max(jj,:);
    hist(buffer(~isnan(buffer)))
    xlim([-5 5])
    saveas(gcf,['E:\RAM data set\RAM_Public_Data_all\FR1_final\z_hist\' num2str(jj) '_z_histogram.png'])
    close all;
end



test_ii=find(b<0.05);
Atlas_index(test_ii(k))
for k=1:size(test_ii,2)
    electrode_volume(atlas_1mm_image_data==Atlas_index(test_ii(k)))=Atlas_index(test_ii(k));
end

mnit1_1mm_stripped_new=fmris_read_nifti(mni_nifti_path);
mnit1_1mm_stripped_new.file_name=['test_shen_time10.nii'];
mnit1_1mm_stripped_new.precision='float32';
mnit1_1mm_stripped_new.byte=4;
mnit1_1mm_stripped_new.data=electrode_volume;
fmris_write_nifti(mnit1_1mm_stripped_new);



