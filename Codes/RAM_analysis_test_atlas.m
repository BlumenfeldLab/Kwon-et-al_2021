%% reorient shen 268 roi template
clc;clear;
rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

mni_nifti_path = 'C:\yale\bioimagesuite30\images\MNI_T1_1mm_stripped.nii';
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
Tmni_yale = mnit1_1mm_stripped.transform;
transformed_volume_shen=zeros(size(mnit1_1mm_stripped.img));

mni_nifti_path_shen = 'E:\RAM data set\RAM_Public_Data_all\Atlases\shen_1mm_268_parcellation.nii';
shen_rois = mni2fs_load_nii(mni_nifti_path_shen);
Tmni_shen=shen_rois.transform;
electrode_volume_shen=shen_rois.img;

all_index_shen=find(electrode_volume_shen>0);
[a b c]=ind2sub(size(electrode_volume_shen),all_index_shen);
volume_coordinates_shen=[a b c ones(size(all_index_shen,1),1)];
mni_coordinates_shen=Tmni_shen*volume_coordinates_shen';
new_volume_coordinates_shen=mni_coordinates_shen'/Tmni_yale';
new_volume_coordinates_shen=new_volume_coordinates_shen(:,1:3);

ori_val=electrode_volume_shen(all_index_shen);
for t_i=1:size(ori_val,1)
    aa=new_volume_coordinates_shen(t_i,1);
    bb=new_volume_coordinates_shen(t_i,2);
    cc=new_volume_coordinates_shen(t_i,3);
    transformed_volume_shen(aa,bb+2,cc)=ori_val(t_i);
end

cd('E:\RAM data set\RAM_Public_Data_all\Atlases')
broadmann_rois_new=fmris_read_nifti(mni_nifti_path);
broadmann_rois_new.file_name=['shen_1mm_268_parcellation_yale.nii'];
broadmann_rois_new.precision='float32';
broadmann_rois_new.byte=4;
broadmann_rois_new.data=transformed_volume_shen;
fmris_write_nifti(broadmann_rois_new);



%% shen 268 rois --> networks
%% all
clc;clear;

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

load('shen_268_rois.mat')

cd('E:\RAM data set\RAM_Public_Data_all\Atlases')
mni_nifti_path = 'E:\RAM data set\RAM_Public_Data_all\Atlases\shen_1mm_268_parcellation_yale.nii';
shen_rois = mni2fs_load_nii(mni_nifti_path);
electrode_volume_data=shen_rois.img;

electrode_volume=zeros(181,217,181);
electrode_volume(electrode_volume_data>0)=roi_labels(electrode_volume_data(electrode_volume_data>0),2);

shen_rois_new=fmris_read_nifti(mni_nifti_path);
shen_rois_new.file_name=['shen_1mm_268_parcellation_network.nii'];
shen_rois_new.precision='float32';
shen_rois_new.byte=4;
shen_rois_new.data=electrode_volume;
fmris_write_nifti(shen_rois_new);

%% each networks
clc;clear;

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

load('shen_268_rois.mat')

cd('E:\RAM data set\RAM_Public_Data_all\Atlases')
mni_nifti_path = 'E:\RAM data set\RAM_Public_Data_all\Atlases\shen_1mm_268_parcellation_yale.nii';
shen_rois = mni2fs_load_nii(mni_nifti_path);
electrode_volume_data=shen_rois.img;

Network_name={'MedialFrontal','FrontoParietal','DefaultMode' ... 
    ,'SubcorticalCerebellum','Motor','Visual1','Visual2','VisualAssociation'};

electrode_volume=zeros(181,217,181);
electrode_volume(electrode_volume_data>0)=roi_labels(electrode_volume_data(electrode_volume_data>0),2);

for network_i=1:8
    cd('E:\RAM data set\RAM_Public_Data_all\Atlases')
    cd('shen_1mm_268_networks')
    
    electrode_volume_network=zeros(181,217,181);
    electrode_volume_network(find(electrode_volume==network_i))=network_i;
    
    shen_rois_new=fmris_read_nifti(mni_nifti_path);
    shen_rois_new.file_name=[num2str(network_i) '_' Network_name{network_i} '.nii'];
    shen_rois_new.precision='float32';
    shen_rois_new.byte=4;
    shen_rois_new.data=electrode_volume_network;
    fmris_write_nifti(shen_rois_new);
end


%% Broadmann atlas
clc;clear;

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

cd('E:\RAM data set\RAM_Public_Data_all\Atlases')
mni_nifti_path = 'E:\RAM data set\RAM_Public_Data_all\Atlases\yale_broadmann.nii';

broadmann_rois = mni2fs_load_nii(mni_nifti_path);
electrode_volume_data=broadmann_rois.img;

electrode_volume_data(find(electrode_volume_data>100))=electrode_volume_data(find(electrode_volume_data>100))-100;

electrode_volume=zeros(181,217,181);
cortical_index=(electrode_volume_data<48 & electrode_volume_data>0);
electrode_volume(cortical_index)=1;

broadmann_rois_new=fmris_read_nifti(mni_nifti_path);
broadmann_rois_new.file_name=['yale_broadmann_cortical.nii'];
broadmann_rois_new.precision='float32';
broadmann_rois_new.byte=4;
broadmann_rois_new.data=electrode_volume;
fmris_write_nifti(broadmann_rois_new);


%% Shen 268 rois with Broadmann atlas
clc;clear;

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

cd('E:\RAM data set\RAM_Public_Data_all\Atlases')
mni_nifti_path_broadmann = 'E:\RAM data set\RAM_Public_Data_all\Atlases\yale_broadmann.nii';
broadmann_rois = mni2fs_load_nii(mni_nifti_path_broadmann);
electrode_volume_broadmann=broadmann_rois.img;
electrode_volume_broadmann(find(electrode_volume_broadmann>100))=electrode_volume_broadmann(find(electrode_volume_broadmann>100))-100;

electrode_volume_broadmann_mask=zeros(181,217,181);
cortical_index=(electrode_volume_broadmann<48 & electrode_volume_broadmann>0);
electrode_volume_broadmann_mask(cortical_index)=1;

mni_nifti_path_shen = 'E:\RAM data set\RAM_Public_Data_all\Atlases\shen_1mm_268_parcellation_yale.nii';
shen_rois = mni2fs_load_nii(mni_nifti_path_shen);
electrode_volume_shen=shen_rois.img;
masked_electrode_volume_shen = electrode_volume_shen.*electrode_volume_broadmann_mask;
masked_index_all=unique(masked_electrode_volume_shen);
masked_index_all(1)=[];

masked_counts=countmember(masked_index_all,masked_electrode_volume_shen);
ori_counts=countmember(masked_index_all,electrode_volume_shen)/2;

masked_index_all(masked_counts < ori_counts)=[];

electrode_volume_shen_cortical=zeros(181,217,181);
for masked_index_i=1:size(masked_index_all,1)
    cortical_index_shen=electrode_volume_shen==masked_index_all(masked_index_i);
    electrode_volume_shen_cortical(cortical_index_shen)=masked_index_all(masked_index_i);
end

broadmann_rois_new=fmris_read_nifti(mni_nifti_path_broadmann);
broadmann_rois_new.file_name=['shen_1mm_268_parcellation_cortical.nii'];
broadmann_rois_new.precision='float32';
broadmann_rois_new.byte=4;
broadmann_rois_new.data=electrode_volume_shen_cortical;
fmris_write_nifti(broadmann_rois_new);




%% shen 268 rois --> networks (cortical)
%% all
clc;clear;

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

load('shen_268_rois.mat')
load('shen_268_rois_cortical.mat')

cd('E:\RAM data set\RAM_Public_Data_all\Atlases')
mni_nifti_path = 'E:\RAM data set\RAM_Public_Data_all\Atlases\shen_1mm_268_parcellation_cortical.nii';
shen_rois = mni2fs_load_nii(mni_nifti_path);
electrode_volume_data=shen_rois.img;

electrode_volume=zeros(181,217,181);
electrode_volume(electrode_volume_data>0)=roi_labels(electrode_volume_data(electrode_volume_data>0),2);

shen_rois_new=fmris_read_nifti(mni_nifti_path);
shen_rois_new.file_name=['shen_1mm_268_parcellation_network_cortical.nii'];
shen_rois_new.precision='float32';
shen_rois_new.byte=4;
shen_rois_new.data=electrode_volume;
fmris_write_nifti(shen_rois_new);

%% each networks
clc;clear;

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

load('shen_268_rois.mat')

cd('E:\RAM data set\RAM_Public_Data_all\Atlases')
mni_nifti_path = 'E:\RAM data set\RAM_Public_Data_all\Atlases\shen_1mm_268_parcellation_cortical.nii';
shen_rois = mni2fs_load_nii(mni_nifti_path);
electrode_volume_data=shen_rois.img;

Network_name={'MedialFrontal','FrontoParietal','DefaultMode' ... 
    ,'SubcorticalCerebellum','Motor','Visual1','Visual2','VisualAssociation'};

electrode_volume=zeros(181,217,181);
electrode_volume(electrode_volume_data>0)=roi_labels(electrode_volume_data(electrode_volume_data>0),2);

for network_i=1:8
    cd('E:\RAM data set\RAM_Public_Data_all\Atlases')
    cd('shen_1mm_268_networks')
    
    electrode_volume_network=zeros(181,217,181);
    electrode_volume_network(find(electrode_volume==network_i))=network_i;
    
    shen_rois_new=fmris_read_nifti(mni_nifti_path);
    shen_rois_new.file_name=[num2str(network_i) '_' Network_name{network_i} '_cortical.nii'];
    shen_rois_new.precision='float32';
    shen_rois_new.byte=4;
    shen_rois_new.data=electrode_volume_network;
    fmris_write_nifti(shen_rois_new);
end






