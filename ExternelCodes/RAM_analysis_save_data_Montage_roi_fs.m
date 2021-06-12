clc;clear;

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

load r1_all.mat

% load fs_atlas_roi.mat
% all_roi_name(1,:)=[];
% all_roi_name(4,:)=[];
% all_roi_name_cortical=all_roi_name;
% save('fs_atlas_roi_cortical.mat','all_roi_name_cortical')

% load fs_DK_atlas.mat
% DK_atlas_names(1,:)=[];
% DK_atlas_names(4,:)=[];
load fs_DK_atlas_cortical.mat

fid=fopen('Subjects_list_all.txt','r');
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);

Montage_prefix='final_roi_fs';
cnt=0;
for i=1:251 % 251 subject
    try
        cd(rootfolder)
        cd('FR1_FS')
        cd([num2str(i),'_',r_sublist{i,1}]);
        load('Bad_electrode_categories_index.mat')
        load('Raw_data_extended_5s.mat', 'location_data_pair')
        
        elec_n=size(location_data_pair,2);

        left_index=1;
        right_index=1;

        L_MontageMap=[];
        R_MontageMap=[];
        L_MontageMap_roi=[];
        R_MontageMap_roi=[];
        
        good_electrodes=1:elec_n;
        good_electrodes(Bad_electrode_categories_index)=[];
        for elec_i=good_electrodes
            region_name=[];
            channel_inf=location_data_pair{1,elec_i};
            try
                region_name=channel_inf.atlases.avg.region;
            catch
            end
            if(~isempty(region_name))
                if(ismember(region_name,DK_atlas_names_cortical))
                    fsavg_coord=[channel_inf.atlases.avg.x channel_inf.atlases.avg.y channel_inf.atlases.avg.z]';
                    if(fsavg_coord(1)<0)
                        L_MontageMap(left_index,:)=[elec_i fsavg_coord(1) fsavg_coord(2) fsavg_coord(3)];
                        L_MontageMap_roi{left_index,1}=region_name;
                        left_index=left_index+1;
                    else
                        R_MontageMap(right_index,:)=[elec_i fsavg_coord(1) fsavg_coord(2) fsavg_coord(3)]; 
                        R_MontageMap_roi{right_index,1}=region_name;
                        right_index=right_index+1;
                    end
                end
            end
        end
        save(['R_MontageMap_' Montage_prefix '.mat'], 'R_MontageMap', 'R_MontageMap_roi');
        save(['L_MontageMap_' Montage_prefix '.mat'], 'L_MontageMap', 'L_MontageMap_roi');

        cd(rootfolder)
        cd('FR1_copy')
        mkdir([num2str(i),'_',r_sublist{i,1}]);
        cd([num2str(i),'_',r_sublist{i,1}]);
        save(['R_MontageMap_' Montage_prefix '.mat'], 'R_MontageMap', 'R_MontageMap_roi');
        save(['L_MontageMap_' Montage_prefix '.mat'], 'L_MontageMap', 'L_MontageMap_roi');

        disp(sprintf(' %s subject completed!! ',r_sublist{i,1}));
        cnt=cnt+1;
    catch
    end
end


%% x flip
clc;clear;

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

load r1_all.mat
load fs_DK_atlas_cortical.mat

fid=fopen('Subjects_list_all.txt','r');
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);

Montage_prefix='final_roi_fs';
cnt=0;
for i=1:251 % 251 subject
    try
        cd(rootfolder)
        cd('FR1_FS')
        cd([num2str(i),'_',r_sublist{i,1}]);
        load('Bad_electrode_categories_index.mat')
        load('Raw_data_extended_5s.mat', 'location_data_pair')

        elec_n=size(location_data_pair,2);

        left_index=1;
        right_index=1;

        L_MontageMap=[];
        R_MontageMap=[];
        L_MontageMap_roi=[];
        R_MontageMap_roi=[];
        
        good_electrodes=1:elec_n;
        good_electrodes(Bad_electrode_categories_index)=[];
        for elec_i=good_electrodes
            region_name=[];
            channel_inf=location_data_pair{1,elec_i};
            try
                region_name=channel_inf.atlases.avg.region;
            catch
            end
            if(~isempty(region_name))
                if(ismember(region_name,DK_atlas_names_cortical))
                    fsavg_coord=[channel_inf.atlases.avg.x channel_inf.atlases.avg.y channel_inf.atlases.avg.z]';
                    if(fsavg_coord(1)<0)
                        L_MontageMap(left_index,:)=[elec_i fsavg_coord(1) fsavg_coord(2) fsavg_coord(3)];
                        L_MontageMap_roi{left_index,1}=region_name;
                        left_index=left_index+1;
                    else
                        L_MontageMap(left_index,:)=[elec_i -fsavg_coord(1) fsavg_coord(2) fsavg_coord(3)];
                        L_MontageMap_roi{left_index,1}=region_name;
                        left_index=left_index+1;
                    end
                end
            end
        end
        R_MontageMap=L_MontageMap;
        R_MontageMap(:,2)=abs(R_MontageMap(:,2));
        
        R_MontageMap_roi=L_MontageMap_roi;
        
        save(['L_MontageMap_' Montage_prefix '_xflip.mat'], 'L_MontageMap', 'L_MontageMap_roi');
        save(['R_MontageMap_' Montage_prefix '_xflip.mat'], 'R_MontageMap', 'R_MontageMap_roi');
        
        cd(rootfolder)
        cd('FR1_copy')
        cd([num2str(i),'_',r_sublist{i,1}]);
        
        save(['L_MontageMap_' Montage_prefix '_xflip.mat'], 'L_MontageMap', 'L_MontageMap_roi');
        save(['R_MontageMap_' Montage_prefix '_xflip.mat'], 'R_MontageMap', 'R_MontageMap_roi');
        disp(sprintf(' %s subject completed!! ',r_sublist{i,1}));
        cnt=cnt+1;
    catch
    end
end

