
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

data_location='E:\RAM data set\RAM_Public_Data_all\FR1_FARNAM';

% file_suffix='Wendy_ftest_fdr';
% file_suffix='Wendy';
file_suffix='Wendy_fdr';
% file_suffix='final_clean';
% file_suffix='Wendy_thr_02';
% file_suffix='Wendy_thr_01';
% file_suffix='Wendy_merge';

fs=256;
load('time_31ms_no_overlap.mat') % 31.125ms
base_sample_idx=16;
window_size=0.03125;

% load('time_63ms_no_overlap.mat') % 63.25ms
% base_sample_idx=8;
% window_size=0.0626;

for i=1:251
    try
        cd(data_location)
        cd([num2str(i),'_',r_sublist{i,1}]);
        mkdir(file_suffix)
        
        for side_index = 1:2   % Generate frames for L and R
            if side_index == 1
                side = 'L';
            else
                side = 'R';
            end
            
            load([side '_vertex_values_' file_suffix '.mat'])
            for time_i=1:size(T,2)
                vertex_values_time=[];
                vertex_values_time=vertex_values(:,time_i);
                save([file_suffix '/' side '_vertex_values_' num2str(time_i) '.mat'],'vertex_values_time');  
            end
        end
        disp(sprintf('%s is completed! ',r_sublist{i,1}));
    catch
    end
end



