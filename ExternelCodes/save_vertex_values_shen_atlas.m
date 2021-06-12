clc;clear;

cd('E:\RAM data set\RAM_Public_Data_all')
load('fs_surf.mat')

atlas_nifti_path = 'E:\RAM data set\RAM_Public_Data_all\Atlases\shen_1mm_268_parcellation_cortical.nii';
atlas_1mm_image = mni2fs_load_nii(atlas_nifti_path);
atlas_1mm_image_data=atlas_1mm_image.img;
Atlas_index=unique(atlas_1mm_image_data);
Atlas_index(1,:)=[];
k=2;
all_k=(k^2+1)^3;

for side_index=1:2
    vertices=fs_surf(side_index).vertices;
    vertices_volume=round(vertices);
    size_vertices=size(fs_surf(side_index).vertices,1);
    vertex_values_shen_atlas=zeros(size_vertices,1);
    for atlas_i=1:size_vertices
        cnt=1;
        vertex_values_shen_atlas_buffer=zeros(all_k,1);
        for a=-k:k
            for b=-k:k
                for c=-k:k
                    try
                        new_x=vertices_volume(atlas_i,1)+a;
                        new_y=vertices_volume(atlas_i,2)+b;
                        new_z=vertices_volume(atlas_i,3)+c;
                        vertex_values_shen_atlas_buffer(cnt)=atlas_1mm_image_data(new_x,new_y,new_z);
                        cnt=cnt+1;
                    catch
                        cnt=cnt+1;
                    end
                end
            end
        end
        vertex_values_shen_atlas_buffer(vertex_values_shen_atlas_buffer==0)=[];
        
        if side_index == 1
            vertex_values_shen_atlas_buffer(vertex_values_shen_atlas_buffer<100)=[];
        else
            vertex_values_shen_atlas_buffer(vertex_values_shen_atlas_buffer>100)=[];
        end
        
        vertex_values_shen_atlas(atlas_i,1)=mode(vertex_values_shen_atlas_buffer);
    end
    
    if side_index == 1
        side = 'L';
    else
        side = 'R';
    end
    save([side '_vertex_values_shen_atlas.mat'],'vertex_values_shen_atlas');

end


%% plot
clc;clear;

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

for side_index=1:2
    if side_index == 1
        side = 'L';
        load([side '_vertex_values_shen_atlas.mat']);
        L_vertex_values_frame=vertex_values_shen_atlas;
    else
        side = 'R';
        load([side '_vertex_values_shen_atlas.mat']);    
        R_vertex_values_frame=vertex_values_shen_atlas;
    end
end
save(['Shen_atlas.mat'], 'R_vertex_values_frame','L_vertex_values_frame');

%% all plot
image_storage_folder='E:\RAM data set\RAM_Public_Data_all\FR1_final\Atlas';
overlay=4;
laterality=0;
views=4;
inflationstep=5;
color_range=[1 235];
L_vertex_values_frame_buf=L_vertex_values_frame;
R_vertex_values_frame_buf=R_vertex_values_frame;

displayElectrodesInflated_single(image_storage_folder, overlay,laterality, views, ... 
    inflationstep,L_vertex_values_frame_buf,R_vertex_values_frame_buf)








%% individual plot
image_storage_folder='E:\RAM data set\RAM_Public_Data_all\FR1_final\atlas';
overlay=4;
laterality=0;
views=4;
inflationstep=5;
% color_range=[1 235];
for buf_i=1:188
    L_vertex_values_frame_buf=L_vertex_values_frame;
    L_vertex_values_frame_buf(L_vertex_values_frame~=Atlas_index(buf_i))=0;
    
    R_vertex_values_frame_buf=R_vertex_values_frame;
    R_vertex_values_frame_buf(R_vertex_values_frame~=Atlas_index(buf_i))=0;
    
    displayElectrodesInflated_single(image_storage_folder, overlay,laterality, views, ... 
        inflationstep,L_vertex_values_frame_buf,R_vertex_values_frame_buf)
    close all;
    
    copyfile('combined_views_full_1.tiff',['Individual ROI/ROI_' num2str(Atlas_index(buf_i)) '.tiff'])
end





