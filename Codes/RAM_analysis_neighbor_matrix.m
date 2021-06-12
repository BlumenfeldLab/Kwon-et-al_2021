%% find nearest ROIs
clc;clear
start_side = 1;
sides_num = 2;
neighbor_step=2;
atlas_index='400'; % '1000' '200','DK','400'

rootfolder='E:\RAM data set\RAM_Public_Data_all\Codes';
cd(rootfolder)

% L and R separately
for side_index = start_side:sides_num    % Generate frames for L and R
    
    % Determine what side we are working on
    if side_index == 1
        side = 'L';
        switch atlas_index
            case '1000'
                [vertices, label, colortable]=read_annotation('lh.Schaefer2018_1000Parcels_7Networks_order.annot');
            case '400'
                [vertices, label, colortable]=read_annotation('lh.Schaefer2018_400Parcels_7Networks_order.annot');
            case '200'
                [vertices, label, colortable]=read_annotation('lh.Schaefer2018_200Parcels_7Networks_order.annot');
            case '100'
                [vertices, label, colortable]=read_annotation('lh.Schaefer2018_100Parcels_7Networks_order.annot');
            case 'DK'
                [vertices, label, colortable] = read_annotation('lh.aparc.annot');
            otherwise
                disp('Error, no such atlas.')
        end
    else
        side = 'R';
        switch atlas_index
            case '1000'
                [vertices, label, colortable]=read_annotation('rh.Schaefer2018_1000Parcels_7Networks_order.annot');
            case '400'
                [vertices, label, colortable]=read_annotation('rh.Schaefer2018_400Parcels_7Networks_order.annot');
            case '200'
                [vertices, label, colortable]=read_annotation('rh.Schaefer2018_200Parcels_7Networks_order.annot');
            case '100'
                [vertices, label, colortable]=read_annotation('rh.Schaefer2018_100Parcels_7Networks_order.annot');
            case 'DK'
                [vertices, label, colortable] = read_annotation('rh.aparc.annot');
            otherwise
                disp('Error, no such atlas.')
        end
    end

    load('FS_ori_avg_surface.mat')
    faces_buff=fs_surf(start_side).faces;
    fs_roi_index_all=colortable.table(:,5);
    roi_neighbor_matrix=zeros(size(fs_roi_index_all,1),size(fs_roi_index_all,1));

    for atlas_i=2:size(colortable.table,1)
        roi_vertex_index=[];
        roi_neighbor_vertex_index=[];
        roi_neighbor_index=[];

        roi_vertex_index=find(label == colortable.table(atlas_i,5));
        roi_neighbor_vertex_index=find_neighbor_vertices(faces_buff,roi_vertex_index,neighbor_step);

        roi_neighbor_index=unique(label(roi_neighbor_vertex_index,1));
        for neighbor_i=1:size(roi_neighbor_index,1)
            roi_neighbor_index_buff=find(fs_roi_index_all==roi_neighbor_index(neighbor_i));
            if roi_neighbor_index_buff ~= 1 % Remove background roi, 1=background
                roi_neighbor_matrix(atlas_i,roi_neighbor_index_buff)=1;
                roi_neighbor_matrix(roi_neighbor_index_buff,atlas_i)=1;
            end
        end

        disp([side ', ROI : ' num2str(atlas_i)])

    end

%     rootfolder='E:\RAM data set\RAM_Public_Data_all\FR1_FS\Atlas_info';
    cd(rootfolder)

    save([side '_roi_neighbor_matrix_' atlas_index '_' num2str(neighbor_step) '.mat'],'roi_neighbor_matrix');
end

% figure('position',[0 0 8000 8000]);
% set(gcf, 'color', [1 1 1]);
% set(gcf,'Visible','on');
% 
% imagesc(roi_neighbor_matrix)
% title(['ROI neighbor matrix'])
% set(gca,'ydir','normal');
% set(gca,'Fontsize',15);
% ylabel('ROI index');
% xlabel('ROI index');


%% both
atlas_index='400'; % '1000' '200','DK','400'
for side_index = start_side:sides_num    % Generate frames for L and R
    
    if side_index == 1
        side = 'L';
        roi_neighbor_matrix=[];
        load([side '_roi_neighbor_matrix_' atlas_index '_2.mat' ]);
        L_roi_neighbor_matrix=roi_neighbor_matrix;
    else
        side = 'R';
        roi_neighbor_matrix=[];
        load([side '_roi_neighbor_matrix_' atlas_index '_2.mat' ]);
        R_roi_neighbor_matrix=roi_neighbor_matrix;
    end
end
buff_neighbor_matrix=zeros(size(L_roi_neighbor_matrix,1),size(L_roi_neighbor_matrix,2));
Both_roi_neighbor_matrix = [L_roi_neighbor_matrix buff_neighbor_matrix ;buff_neighbor_matrix R_roi_neighbor_matrix];
roi_neighbor_matrix=[];
roi_neighbor_matrix=Both_roi_neighbor_matrix;
save(['Both_roi_neighbor_matrix_' atlas_index '_' num2str(neighbor_step) '.mat'],'roi_neighbor_matrix');


figure('position',[0 0 8000 8000]);
set(gcf, 'color', [1 1 1]);
set(gcf,'Visible','on');

imagesc(roi_neighbor_matrix)
title(['ROI neighbor matrix'])
set(gca,'ydir','normal');
set(gca,'Fontsize',15);
ylabel('ROI index');
xlabel('ROI index');

