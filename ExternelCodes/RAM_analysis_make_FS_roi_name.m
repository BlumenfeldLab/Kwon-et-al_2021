%% Freesurfer roi name
clc;clear;
close all;

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

load r1_all.mat
fid=fopen('Subjects_list_all.txt','r');
cnt=1;
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);

% data_location='E:\RAM data set\RAM_Public_Data_all\FR1_final';
data_location='E:\RAM data set\RAM_Public_Data_all\FR1_FARNAM';
cd(data_location)

load('fs_atlas_roi.mat')
Total_elec_roi_number=zeros(36,156);
subject_name=[];
cnt=1;
for i=1:251 % 251 subject
    try
        cd(rootfolder)
        cd('FR1_FARNAM')
        cd([num2str(i),'_',r_sublist{i,1}]);
        load('stats_traces_Wendy_hunki_std_5_20.mat', 'location_data_pair')
        subject_name{1,cnt}=[num2str(i),'_',r_sublist{i,1}];
%         roi_name=[];
        for elec_i=1:size(location_data_pair,2)
            roi_name_buff=[];
            roi_name_index=[];
            channel_inf=location_data_pair{1,elec_i};
            roi_name_buff=channel_inf.atlases.avg.region;
            if ~isempty(roi_name_buff)
                roi_name_index=find(strcmp(all_roi_name,roi_name_buff));
                Total_elec_roi_number(roi_name_index,cnt)=Total_elec_roi_number(roi_name_index,cnt)+1;
            else
                Total_elec_roi_number(1,cnt)=Total_elec_roi_number(1,cnt)+1;
            end
        end
%         [Au,~,idx2]=uniquecell(roi_name);
        cnt=cnt+1
    catch
    end
end
Total_subject_roi_number=logical(Total_elec_roi_number);

