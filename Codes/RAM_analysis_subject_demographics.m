clc;clear;
close all;

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

load r1_all.mat
fid=fopen('Subjects_list_all.txt','r');
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);

data_location='E:\RAM data set\RAM_Public_Data_all\FR1_FS';
cd(data_location)

cnt=1;
% subject_elec_trial_inf
Number_of_trials=[];
for subject_i=1:251
    try
        cd(data_location)
        cd([num2str(subject_i),'_',r_sublist{subject_i,1}]);
        
        load('rejection_idx_Final_std_20_fs.mat')
        Number_of_trials(cnt)=RAM_all_trials_total_number;
        
        load('stats_traces_Final_std_20_fs_together.mat')
        load('L_MontageMap_final_roi_fs.mat')
        load('R_MontageMap_final_roi_fs.mat')
        
        subject_elec_trial_inf(cnt).Subject_name=r_sublist{subject_i,1};
        subject_elec_trial_inf(cnt).all_e=size(location_data_pair,2);
        subject_elec_trial_inf(cnt).rem_e=size(L_MontageMap,1)+size(R_MontageMap,1);
        subject_elec_trial_inf(cnt).rem_e_left=size(L_MontageMap,1);
        subject_elec_trial_inf(cnt).rem_e_right=size(R_MontageMap,1);
        
        % left
        trial_good_all_left=0;
        trial_good_Rec_left=0;
        trial_good_NRec_left=0;
        trial_all_left=0;
        additional_noise_count=0;
        if size(L_MontageMap,1) > 0
            rem_e_all_index_left=L_MontageMap(:,1)';
            for ii=rem_e_all_index_left
                if noisy_RAM_all_trials(ii).All_ratio < 0.5
                    trial_buf=noisy_RAM_all_trials(ii).All_total-noisy_RAM_all_trials(ii).All_noise;
                    trial_buf_Rec=noisy_RAM_all_trials(ii).Rec_total-noisy_RAM_all_trials(ii).Rec_noise;
                    trial_buf_NRec=noisy_RAM_all_trials(ii).NRec_total-noisy_RAM_all_trials(ii).NRec_noise;
                    trial_good_all_left=trial_good_all_left+trial_buf;
                    trial_good_Rec_left=trial_good_Rec_left+trial_buf_Rec;
                    trial_good_NRec_left=trial_good_NRec_left+trial_buf_NRec;
                    trial_all_left=trial_all_left+noisy_RAM_all_trials(ii).All_total;
                else
                    additional_noise_count=additional_noise_count+1;
                end
            end
            subject_elec_trial_inf(cnt).bad_e_left=additional_noise_count;
        else
            subject_elec_trial_inf(cnt).bad_e_left=additional_noise_count;
        end
        
        % right
        trial_good_all_right=0;
        trial_good_Rec_right=0;
        trial_good_NRec_right=0;
        trial_all_right=0;
        additional_noise_count=0;
        if size(R_MontageMap,1) > 0
            rem_e_all_index_right=R_MontageMap(:,1)';
            for ii=rem_e_all_index_right
                if noisy_RAM_all_trials(ii).All_ratio < 0.5
                    trial_buf=noisy_RAM_all_trials(ii).All_total-noisy_RAM_all_trials(ii).All_noise;
                    trial_buf_Rec=noisy_RAM_all_trials(ii).Rec_total-noisy_RAM_all_trials(ii).Rec_noise;
                    trial_buf_NRec=noisy_RAM_all_trials(ii).NRec_total-noisy_RAM_all_trials(ii).NRec_noise;
                    trial_good_all_right=trial_good_all_right+trial_buf;
                    trial_good_Rec_right=trial_good_Rec_right+trial_buf_Rec;
                    trial_good_NRec_right=trial_good_NRec_right+trial_buf_NRec;
                    trial_all_right=trial_all_right+noisy_RAM_all_trials(ii).All_total;
                else
                    additional_noise_count=additional_noise_count+1;
                end
            end
            subject_elec_trial_inf(cnt).bad_e_right=additional_noise_count;
        else
            subject_elec_trial_inf(cnt).bad_e_right=additional_noise_count;
        end
        
        trial_all=trial_all_left+trial_all_right;
        trial_good_all=trial_good_all_left+trial_good_all_right;
        subject_elec_trial_inf(cnt).all_t=trial_all;
        subject_elec_trial_inf(cnt).rem_t=trial_good_all;
        subject_elec_trial_inf(cnt).ratio_t=(trial_all-trial_good_all)/trial_all;
        
        trial_good_Rec=trial_good_Rec_left+trial_good_Rec_right;
        trial_good_NRec=trial_good_NRec_left+trial_good_NRec_right;

        subject_elec_trial_inf(cnt).rem_t_Rec=trial_good_Rec;
        subject_elec_trial_inf(cnt).rem_t_NRec=trial_good_NRec;
        subject_elec_trial_inf(cnt).avg_rem_t=trial_good_all/(size(L_MontageMap,1)+size(R_MontageMap,1));
        subject_elec_trial_inf(cnt).avg_rem_t_Rec=trial_good_Rec/(size(L_MontageMap,1)+size(R_MontageMap,1));
        subject_elec_trial_inf(cnt).avg_rem_t_NRec=trial_good_NRec/(size(L_MontageMap,1)+size(R_MontageMap,1));
        
        cnt=cnt+1;
        disp(cnt)
    catch
    end
end

mean([subject_elec_trial_inf.avg_rem_t])
std([subject_elec_trial_inf.avg_rem_t])
mean([subject_elec_trial_inf.avg_rem_t_Rec])
mean([subject_elec_trial_inf.avg_rem_t_NRec])

std([subject_elec_trial_inf.avg_rem_t_Rec]./[subject_elec_trial_inf.avg_rem_t])



mean([subject_elec_trial_inf.ratio_t])
sum([subject_elec_trial_inf.all_e])
sum([subject_elec_trial_inf.rem_e_left])
sum([subject_elec_trial_inf.rem_e_right])
sum([subject_elec_trial_inf.bad_e_left])
sum([subject_elec_trial_inf.bad_e_right])

%% extract subject information
rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

[Subject_name Gender Age Handedness] = csvimport('RAM_subject_demographics_all.csv', 'columns', {'Subject Number', 'Gender', 'Age', 'Handedness'});

load r1_all.mat
fid=fopen('Subjects_list_all.txt','r');
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);

data_location='E:\RAM data set\RAM_Public_Data_all\FR1_FS';
cd(data_location)

age_all=nan(158,1);
sex_all=nan(158,1);
hand_all=nan(158,1);
cnt=1;
All_subject_info=[Subject_name Gender Age Handedness];
subject_index_final=[];
for subject_i=1:251
    try
        cd(data_location)
        cd([num2str(subject_i),'_',r_sublist{subject_i,1}]);
        if strcmpi(Gender{subject_i,1},'Male') == 1
            sex_all(cnt,1)=1;
        elseif strcmpi(Gender{subject_i,1},'Female') == 1
            sex_all(cnt,1)=2;
        end
%         if strcmpi(Handedness{subject_i,1},'Male') == 1
%             hand_all(cnt,1)=1;
%         elseif strcmpi(Handedness{subject_i,1},'Female') == 1
%             hand_all(cnt,1)=2;
%         end
        
        age_all(cnt,1)=str2double(Age{subject_i,1});
        subject_all{cnt,1}=r_sublist{subject_i,1};
        cnt=cnt+1;
        disp(cnt)
    catch
        subject_index_final=[subject_index_final subject_i];
    end
end
All_subject_info(subject_index_final,:)=[];

cd(rootfolder)
save('subject_demographics_final.mat', 'All_subject_info', 'subject_elec_trial_inf');


