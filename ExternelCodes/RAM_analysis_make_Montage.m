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

mni_nifti_path = 'C:\yale\bioimagesuite30\images\MNI_T1_1mm_stripped.nii';
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
Tmni = mnit1_1mm_stripped.transform;

M=[0.9975 -0.0073 0.0176 -0.0429 ; 
   0.0146 1.0009 -0.0024 1.5496 ; 
  -0.0130 -0.0093 0.9971 1.1840]; % https://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems

prefix='Wendy';
Montage_prefix='final_clean_dil';

for i=1:251 % 251 subject
    try
        clearvars RAM_Rec_trials_sessions RAM_nonRec_trials_sessions location_data_pair

        cd(rootfolder)
        cd('FR1_final')
        cd([num2str(i),'_',r_sublist{i,1}]);

        load(['stats_traces_' prefix '.mat'])

        % Make the Montage

        masked_elec_n=size(Atlas_subject_elec_val,2);

        left_index=1;
        right_index=1;

        L_MontageMap=[];
        R_MontageMap=[];

        for elec_i=1:masked_elec_n
            masked_elec=Atlas_subject_elec_val(elec_i);
            channel_inf=location_data_pair{1,masked_elec};
            fsavg_coord=[channel_inf.atlases.avg.x channel_inf.atlases.avg.y channel_inf.atlases.avg.z]';
            mni_coord = M*[fsavg_coord ; 1];

            if(mni_coord(1)<0)
                L_MontageMap(left_index,:)=[masked_elec mni_coord(1) mni_coord(2) mni_coord(3)];
                left_index=left_index+1;
            else
                R_MontageMap(right_index,:)=[masked_elec mni_coord(1) mni_coord(2) mni_coord(3)]; 
                right_index=right_index+1;
            end
        end
        save(['R_MontageMap_' Montage_prefix '.mat'], 'R_MontageMap');
        save(['L_MontageMap_' Montage_prefix '.mat'], 'L_MontageMap');

        disp(sprintf(' %s subject completed!! ',r_sublist{i,1}));
    catch

    end
end

