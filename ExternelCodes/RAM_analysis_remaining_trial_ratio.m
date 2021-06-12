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

% prefix='Wendy_thr_01';
prefix='Wendy_thr_02';
prefix2='Wendy';
cnt=1;
exel_list=[];
artifact_ratio_all=[];
for i=1:251 % 251 subject
    try
        cd(rootfolder)
        cd('FR1_final')
        cd([num2str(i),'_',r_sublist{i,1}]);
        
        %buffer
        load('rejection_idx_extended_all_sessions.mat')
        session_n=size(good_RAM_all_trials_idx,2);
        
        load(['./rejection_idx_all_sessions_' prefix '.mat'])
        load(['./stats_traces_' prefix2 '.mat'])
        
        a=size(location_data_pair,2);
        b=size(Atlas_subject_elec_val,2);
        
        artifact_ratio=[];
        artifact_ratio=noisy_RAM_all_trials_number(Atlas_subject_elec_val,1:session_n)./RAM_all_trials_total_number(1:session_n);
        artifact_ratio_all=[artifact_ratio_all ; artifact_ratio(:)];
        
        artifact_ratio(artifact_ratio>0.25)=NaN;
        c=length(find(~isnan(nanmean(artifact_ratio,2))>0));
        d=c/b;
        
        e=a*session_n;
        f=b*session_n;
        g=length(find(~isnan(artifact_ratio)>0));
        h=g/f;

        exel_list(cnt,:)=[a b c d e f g h];
        
        
        
        cnt=cnt+1;
        disp(sprintf(' %s subject complete',r_sublist{i,1}));
    catch

    end
end

figure;histogram(artifact_ratio_all,100)

