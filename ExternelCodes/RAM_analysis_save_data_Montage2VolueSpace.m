%% cortical area || final!!!!
%% save Montage on Volume space 
clc;clear;

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)
load('subjects.mat')

mni_nifti_path = 'C:\yale\bioimagesuite30\images\MNI_T1_1mm_stripped.nii';
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
Tmni = mnit1_1mm_stripped.transform;

electrode_volume=zeros(181,217,181);
subn='fs2mni';
% subn_atlas='broadmann';
subn_atlas='shen';

% mask_nifti_path = 'E:\RAM data set\RAM_Public_Data_all\Atlases\yale_broadmann_cortical.nii';
mask_nifti_path = 'E:\RAM data set\RAM_Public_Data_all\Atlases\shen_1mm_268_parcellation_cortical.nii';
mask_1mm_image = mni2fs_load_nii(mask_nifti_path);
mask_1mm_image_data=mask_1mm_image.img;
for i=1:142 % 142 subject
    try
        clear R_MontageMap L_MontageMap R_MontageMap_ROI_name L_MontageMap_ROI_name R_MontageMap_ROI_hem L_MontageMap_ROI_hem          
        cd(rootfolder)
        cd('FR1')
        cd(patients{i});
        load(['R_MontageMap_' subn '.mat']);
        load(['L_MontageMap_' subn '.mat']);
        L_electrode_size=size(L_MontageMap,1);
        R_electrode_size=size(R_MontageMap,1);
        
        wm_index=[];
        if (L_electrode_size>0)
            L_MontageMap_volume = [L_MontageMap(:,2:4), ones(L_electrode_size, 1)] / Tmni';
            L_MontageMap_volume = round(L_MontageMap_volume(:,1:3));
            for elec_i=1:size(L_MontageMap_volume,1)
                buf_mask_data=mask_1mm_image_data(L_MontageMap_volume(elec_i,1),L_MontageMap_volume(elec_i,2),L_MontageMap_volume(elec_i,3));
                if(buf_mask_data>0)
                    electrode_volume(L_MontageMap_volume(elec_i,1),L_MontageMap_volume(elec_i,2),L_MontageMap_volume(elec_i,3))=i;
                else
                    wm_index=[wm_index elec_i];
                end
            end
            L_MontageMap(wm_index,:)=[];
        end
        
        wm_index=[];
        if (R_electrode_size>0)
            R_MontageMap_volume = [R_MontageMap(:,2:4), ones(R_electrode_size, 1)] / Tmni';
            R_MontageMap_volume = round(R_MontageMap_volume(:,1:3));
            for elec_i=1:size(R_MontageMap_volume,1)
                buf_mask_data=mask_1mm_image_data(R_MontageMap_volume(elec_i,1),R_MontageMap_volume(elec_i,2),R_MontageMap_volume(elec_i,3));
                if(buf_mask_data>0)
                    electrode_volume(R_MontageMap_volume(elec_i,1),R_MontageMap_volume(elec_i,2),R_MontageMap_volume(elec_i,3))=i;
                else
                    wm_index=[wm_index elec_i];
                end
            end
            R_MontageMap(wm_index,:)=[];
        end
               
        save(['R_MontageMap_fs2mni_cortical_' subn_atlas '.mat'], 'R_MontageMap');
        save(['L_MontageMap_fs2mni_cortical_' subn_atlas '.mat'], 'L_MontageMap');
        
    catch
        i
    end
    
end

cd(rootfolder)
cd('Atlases')

mnit1_1mm_stripped_new=fmris_read_nifti(mni_nifti_path);
mnit1_1mm_stripped_new.file_name=['electrodes_distribution_' subn '_cortical_' subn_atlas '.nii'];
mnit1_1mm_stripped_new.precision='float32';
mnit1_1mm_stripped_new.byte=4;
mnit1_1mm_stripped_new.data=electrode_volume;
fmris_write_nifti(mnit1_1mm_stripped_new);


%% each network

clc;clear;

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)
load('subjects.mat')

Network_name_all={'MedialFrontal','FrontoParietal','DefaultMode' ... 
    ,'SubcorticalCerebellum','Motor','Visual1','Visual2','VisualAssociation'};

network_nifti_path = 'E:\RAM data set\RAM_Public_Data_all\Atlases\shen_1mm_268_parcellation_network_cortical.nii';
network_1mm_image = mni2fs_load_nii(network_nifti_path);
network_1mm_image_data=network_1mm_image.img;
Tmni=network_1mm_image.transform;

subn='fs2mni';
subn_atlas='shen';

for network_index=1:8
    electrode_volume=zeros(181,217,181);
    Network_name=Network_name_all{network_index};
    for i=1:142 % 142 subject
        try
            clear R_MontageMap L_MontageMap          
            cd(rootfolder)
            cd('FR1')
            cd(patients{i});

            load(['R_MontageMap_' subn '_cortical_shen.mat']);
            load(['L_MontageMap_' subn '_cortical_shen.mat']);
            
            
            L_electrode_size=size(L_MontageMap,1);
            R_electrode_size=size(R_MontageMap,1);
            
            wm_index=[];
            if (L_electrode_size>0)
                L_MontageMap_volume = [L_MontageMap(:,2:4), ones(L_electrode_size, 1)] / Tmni';
                L_MontageMap_volume = round(L_MontageMap_volume(:,1:3));
                for elec_i=1:size(L_MontageMap_volume,1)
                    buf_network_data=network_1mm_image_data(L_MontageMap_volume(elec_i,1),L_MontageMap_volume(elec_i,2),L_MontageMap_volume(elec_i,3));
                    if(buf_network_data == network_index)
                        electrode_volume(L_MontageMap_volume(elec_i,1),L_MontageMap_volume(elec_i,2),L_MontageMap_volume(elec_i,3))=i;
                    else
                        wm_index=[wm_index elec_i];
                    end
                end
                L_MontageMap(wm_index,:)=[];
            end
            
            wm_index=[];
            if (R_electrode_size>0)
                R_MontageMap_volume = [R_MontageMap(:,2:4), ones(R_electrode_size, 1)] / Tmni';
                R_MontageMap_volume = round(R_MontageMap_volume(:,1:3));
                for elec_i=1:size(R_MontageMap_volume,1)
                    buf_network_data=network_1mm_image_data(R_MontageMap_volume(elec_i,1),R_MontageMap_volume(elec_i,2),R_MontageMap_volume(elec_i,3));
                    if(buf_network_data == network_index)
                        electrode_volume(R_MontageMap_volume(elec_i,1),R_MontageMap_volume(elec_i,2),R_MontageMap_volume(elec_i,3))=i;
                    else
                        wm_index=[wm_index elec_i];
                    end
                end
                R_MontageMap(wm_index,:)=[];
            end
            save(['R_MontageMap_fs2mni_cortical_' subn_atlas '_' num2str(network_index) '_' Network_name '.mat'], 'R_MontageMap');
            save(['L_MontageMap_fs2mni_cortical_' subn_atlas '_' num2str(network_index) '_' Network_name '.mat'], 'L_MontageMap');
            
        catch
            i
        end
    end
    cd(rootfolder)
	cd('Atlases/shen_1mm_268_networks')

    shen_rois_new=fmris_read_nifti(network_nifti_path);
    shen_rois_new.file_name=['electrodes_distribution_' num2str(network_index) '_' Network_name '.nii'];
    shen_rois_new.precision='float32';
    shen_rois_new.byte=4;
    shen_rois_new.data=electrode_volume;
    fmris_write_nifti(shen_rois_new);
    
    network_index
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% save Montage on Volume space
clc;clear;

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)
load('subjects.mat')

%% all
mni_nifti_path = 'C:\yale\bioimagesuite30\images\MNI_T1_1mm_stripped.nii';
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
Tmni = mnit1_1mm_stripped.transform;

electrode_volume=zeros(181,217,181);
subn='fs2mni_blank';

for i=1:142 % 142 subject
    try
    clear R_MontageMap L_MontageMap R_MontageMap_ROI_name L_MontageMap_ROI_name R_MontageMap_ROI_hem L_MontageMap_ROI_hem          
    cd(rootfolder)
    cd('FR1')
    cd(patients{i});
    load(['R_MontageMap_' subn '.mat']);
    load(['L_MontageMap_' subn '.mat']);
    L_electrode_size=size(L_MontageMap,1);
    R_electrode_size=size(R_MontageMap,1);
    if (L_electrode_size>0)
        L_MontageMap_volume = [L_MontageMap(:,2:4), ones(L_electrode_size, 1)] / Tmni';
        L_MontageMap_volume = round(L_MontageMap_volume(:,1:3));
        for elec_i=1:size(L_MontageMap_volume,1)
            electrode_volume(L_MontageMap_volume(elec_i,1),L_MontageMap_volume(elec_i,2),L_MontageMap_volume(elec_i,3))=i;
        end
    end
    
    if (R_electrode_size>0)
        R_MontageMap_volume = [R_MontageMap(:,2:4), ones(R_electrode_size, 1)] / Tmni';
        R_MontageMap_volume = round(R_MontageMap_volume(:,1:3));
        for elec_i=1:size(R_MontageMap_volume,1)
            electrode_volume(R_MontageMap_volume(elec_i,1),R_MontageMap_volume(elec_i,2),R_MontageMap_volume(elec_i,3))=i;
        end
    end
    catch
        i
    end
    
end

cd(rootfolder)
cd('Atlases')

mnit1_1mm_stripped_new=fmris_read_nifti(mni_nifti_path);
mnit1_1mm_stripped_new.file_name=['electrodes_distribution_' subn '.nii'];
mnit1_1mm_stripped_new.precision='float32';
mnit1_1mm_stripped_new.byte=4;
mnit1_1mm_stripped_new.data=electrode_volume;
fmris_write_nifti(mnit1_1mm_stripped_new);

%% compare between mni and fs2mni
clc;clear;

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)
load('subjects.mat')

mni_nifti_path = 'C:\yale\bioimagesuite30\images\MNI_T1_1mm_stripped.nii';
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
Tmni = mnit1_1mm_stripped.transform;

electrode_volume_mni=zeros(181,217,181);
electrode_volume_fs2mni=zeros(181,217,181);
subn_mni='mni';
subn_fs2mni='fs2mni';

for i=1:142 % 142 subject
    try
    clear R_MontageMap L_MontageMap R_MontageMap_ROI_name L_MontageMap_ROI_name R_MontageMap_ROI_hem L_MontageMap_ROI_hem          
    cd(rootfolder)
    cd('FR1')
    cd(patients{i});
    load(['R_MontageMap_' subn_mni '.mat']);
    load(['L_MontageMap_' subn_mni '.mat']);
    L_electrode_size=size(L_MontageMap,1);
    R_electrode_size=size(R_MontageMap,1);
    if (L_electrode_size>0)
        L_MontageMap_volume = [L_MontageMap(:,2:4), ones(L_electrode_size, 1)] / Tmni';
        L_MontageMap_volume = round(L_MontageMap_volume(:,1:3));
        for elec_i=1:size(L_MontageMap_volume,1)
            electrode_volume_mni(L_MontageMap_volume(elec_i,1),L_MontageMap_volume(elec_i,2),L_MontageMap_volume(elec_i,3))=i;
        end
    end
    
    if (R_electrode_size>0)
        R_MontageMap_volume = [R_MontageMap(:,2:4), ones(R_electrode_size, 1)] / Tmni';
        R_MontageMap_volume = round(R_MontageMap_volume(:,1:3));
        for elec_i=1:size(R_MontageMap_volume,1)
            electrode_volume_mni(R_MontageMap_volume(elec_i,1),R_MontageMap_volume(elec_i,2),R_MontageMap_volume(elec_i,3))=i;
        end
    end
    
    clear R_MontageMap L_MontageMap R_MontageMap_ROI_name L_MontageMap_ROI_name R_MontageMap_ROI_hem L_MontageMap_ROI_hem          
    
    load(['R_MontageMap_' subn_fs2mni '.mat']);
    load(['L_MontageMap_' subn_fs2mni '.mat']);
    L_electrode_size=size(L_MontageMap,1);
    R_electrode_size=size(R_MontageMap,1);
    if (L_electrode_size>0)
        L_MontageMap_volume = [L_MontageMap(:,2:4), ones(L_electrode_size, 1)] / Tmni';
        L_MontageMap_volume = round(L_MontageMap_volume(:,1:3));
        for elec_i=1:size(L_MontageMap_volume,1)
            electrode_volume_fs2mni(L_MontageMap_volume(elec_i,1),L_MontageMap_volume(elec_i,2),L_MontageMap_volume(elec_i,3))=i;
        end
    end
    
    if (R_electrode_size>0)
        R_MontageMap_volume = [R_MontageMap(:,2:4), ones(R_electrode_size, 1)] / Tmni';
        R_MontageMap_volume = round(R_MontageMap_volume(:,1:3));
        for elec_i=1:size(R_MontageMap_volume,1)
            electrode_volume_fs2mni(R_MontageMap_volume(elec_i,1),R_MontageMap_volume(elec_i,2),R_MontageMap_volume(elec_i,3))=i;
        end
    end
    
    
    catch
        i
    end
    
end

cd(rootfolder)
cd('Atlases')

mnit1_1mm_stripped_new=fmris_read_nifti(mni_nifti_path);
mnit1_1mm_stripped_new.file_name=['electrodes_distribution_' subn_mni '.nii'];
mnit1_1mm_stripped_new.precision='float32';
mnit1_1mm_stripped_new.byte=4;
mnit1_1mm_stripped_new.data=electrode_volume_mni;
fmris_write_nifti(mnit1_1mm_stripped_new);

mnit1_1mm_stripped_new=fmris_read_nifti(mni_nifti_path);
mnit1_1mm_stripped_new.file_name=['electrodes_distribution_' subn_fs2mni '_only_mni.nii'];
mnit1_1mm_stripped_new.precision='float32';
mnit1_1mm_stripped_new.byte=4;
mnit1_1mm_stripped_new.data=electrode_volume_fs2mni;
fmris_write_nifti(mnit1_1mm_stripped_new);





%% each roi

clc;clear;

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)
load('subjects.mat')
load('all_ROI_names.mat')

mni_nifti_path = 'C:\yale\bioimagesuite30\images\MNI_T1_1mm_stripped.nii';
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
Tmni = mnit1_1mm_stripped.transform;

subn='mni';

for roi_index=2:36
    electrode_volume=zeros(181,217,181);
    for i=1:142 % 142 subject
        try
            clear R_MontageMap L_MontageMap R_MontageMap_ROI_name L_MontageMap_ROI_name R_MontageMap_ROI_hem L_MontageMap_ROI_hem          
            cd(rootfolder)
            cd('FR1')
            cd(patients{i});


            load(['R_MontageMap_' subn '.mat']);
            load(['L_MontageMap_' subn '.mat']);
            ROI_name=all_ROI_names{roi_index};
            if ~isempty(L_MontageMap_ROI_name)
                L_WM_index=find(~strcmp(L_MontageMap_ROI_name,ROI_name));
            else
                L_WM_index=[];
            end

            if ~isempty(R_MontageMap_ROI_name)
                R_WM_index=find(~strcmp(R_MontageMap_ROI_name,ROI_name));
            else
                R_WM_index=[];
            end

            L_MontageMap(L_WM_index,:)=[];
            R_MontageMap(R_WM_index,:)=[];

            L_electrode_size=size(L_MontageMap,1);
            R_electrode_size=size(R_MontageMap,1);
            if (L_electrode_size>0)
                L_MontageMap_volume = [L_MontageMap(:,2:4), ones(L_electrode_size, 1)] / Tmni';
                L_MontageMap_volume = round(L_MontageMap_volume(:,1:3));
                for elec_i=1:size(L_MontageMap_volume,1)
                    electrode_volume(L_MontageMap_volume(elec_i,1),L_MontageMap_volume(elec_i,2),L_MontageMap_volume(elec_i,3))=i;
                end
            end

            if (R_electrode_size>0)
                R_MontageMap_volume = [R_MontageMap(:,2:4), ones(R_electrode_size, 1)] / Tmni';
                R_MontageMap_volume = round(R_MontageMap_volume(:,1:3));
                for elec_i=1:size(R_MontageMap_volume,1)
                    electrode_volume(R_MontageMap_volume(elec_i,1),R_MontageMap_volume(elec_i,2),R_MontageMap_volume(elec_i,3))=i;
                end
            end

        catch
        i
        end
        
        
    end
    cd(rootfolder)
	cd('Atlases/Electrodes_volume')

    mnit1_1mm_stripped_new=fmris_read_nifti(mni_nifti_path);
    mnit1_1mm_stripped_new.file_name=['electrodes_' ROI_name '_' subn '_volume.nii'];
    mnit1_1mm_stripped_new.precision='float32';
    mnit1_1mm_stripped_new.byte=4;
    mnit1_1mm_stripped_new.data=electrode_volume;
    fmris_write_nifti(mnit1_1mm_stripped_new);
    roi_index
end









%% each subject
clc;clear;

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)
load('subjects.mat')
load('all_ROI_names.mat')

mni_nifti_path = 'C:\yale\bioimagesuite30\images\MNI_T1_1mm_stripped.nii';
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
Tmni = mnit1_1mm_stripped.transform;

subn='mni';

for i=1:142 % 142 subject
    electrode_volume=zeros(181,217,181);
    try
        clear R_MontageMap L_MontageMap R_MontageMap_ROI_name L_MontageMap_ROI_name R_MontageMap_ROI_hem L_MontageMap_ROI_hem          
        cd(rootfolder)
        cd('FR1')
        cd(patients{i});


        load(['L_MontageMap_' subn '.mat']);
        load(['R_MontageMap_' subn '.mat']);

        L_electrode_size=size(L_MontageMap,1);
        R_electrode_size=size(R_MontageMap,1);
        if (L_electrode_size>0)
            L_MontageMap_volume = [L_MontageMap(:,2:4), ones(L_electrode_size, 1)] / Tmni';
            L_MontageMap_volume = round(L_MontageMap_volume(:,1:3));
            for elec_i=1:size(L_MontageMap_volume,1)
                electrode_volume(L_MontageMap_volume(elec_i,1),L_MontageMap_volume(elec_i,2),L_MontageMap_volume(elec_i,3))=1;
            end
        end

        if (R_electrode_size>0)
            R_MontageMap_volume = [R_MontageMap(:,2:4), ones(R_electrode_size, 1)] / Tmni';
            R_MontageMap_volume = round(R_MontageMap_volume(:,1:3));
            for elec_i=1:size(R_MontageMap_volume,1)
                electrode_volume(R_MontageMap_volume(elec_i,1),R_MontageMap_volume(elec_i,2),R_MontageMap_volume(elec_i,3))=1;
            end
        end
        cd(rootfolder)
        cd('Atlases/Electrodes_volume_patients')

        mnit1_1mm_stripped_new=fmris_read_nifti(mni_nifti_path);
        mnit1_1mm_stripped_new.file_name=[patients{i} '_electrodes_mni_volume_' subn '.nii'];
        mnit1_1mm_stripped_new.precision='float32';
        mnit1_1mm_stripped_new.byte=4;
        mnit1_1mm_stripped_new.data=electrode_volume;
        fmris_write_nifti(mnit1_1mm_stripped_new);
        
        
    catch
    i
    end


end


















