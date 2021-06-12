function RAM_FARNAM(i,subject_name)

addpath(genpath('/home/hk582/project/Codes'));
cd([num2str(i),'_',subject_name]);

% M=[0.9975 -0.0073 0.0176 -0.0429 ; 
%    0.0146 1.0009 -0.0024 1.5496 ; 
%   -0.0130 -0.0093 0.9971 1.1840]; % https://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems
% 
% load('stats_traces_Wendy.mat')
% % Make the Montage
% masked_elec_n=size(Atlas_subject_elec_val,2);
% left_index=1;
% right_index=1;
% 
% L_MontageMap=[];
% R_MontageMap=[];
% Montage_prefix='final_clean_dil';
% for elec_i=1:masked_elec_n
%     masked_elec=Atlas_subject_elec_val(elec_i);
%     channel_inf=location_data_pair{1,masked_elec};
%     fsavg_coord=[channel_inf.atlases.avg.x channel_inf.atlases.avg.y channel_inf.atlases.avg.z]';
%     mni_coord = M*[fsavg_coord ; 1];
% 
%     if(mni_coord(1)<0)
%         L_MontageMap(left_index,:)=[masked_elec mni_coord(1) mni_coord(2) mni_coord(3)];
%         left_index=left_index+1;
%     else
%         R_MontageMap(right_index,:)=[masked_elec mni_coord(1) mni_coord(2) mni_coord(3)]; 
%         right_index=right_index+1;
%     end
% end
% save(['R_MontageMap_' Montage_prefix '.mat'], 'R_MontageMap');
% save(['L_MontageMap_' Montage_prefix '.mat'], 'L_MontageMap');


overlay=4;
inflationstep=5;

% file_suffix='Wendy_thr_02';
file_suffix='Wendy_thr_01';
% file_suffix='Wendy_merge';
% file_suffix='Wendy';

displayElectrodesInflated_new_FARNAM(overlay, inflationstep ,file_suffix) %, frames_set)

end