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

fs=256;
trialtypes = {'RAM_Rec_trials','RAM_nonRec_trials'};
task={'FR1' 'FR2' 'FR3'};
session={'x0x30_' 'x0x31_' 'x0x32_' 'x0x33_' 'x0x34_' 'x0x35_' 'x0x36_' 'x0x37_' 'x0x38_' 'x0x39_'};
session_folder={'0' '1' '2' '3' '4' '5' '6' '7' '8' '9'};

load('time_50ms_no_overlap.mat')
prefix='50ms_no_overlap_together';
window_size=round(fs*0.05); % 50ms
samples=round(T*fs);

% load('time_63ms_no_overlap.mat')
% window_size=round(fs*0.063); % 50ms
% samples=round(T*fs);

base_sample_idx=9;

% atlas

mni_nifti_path = 'C:\yale\bioimagesuite30\images\MNI_T1_1mm_stripped.nii';
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
Tmni = mnit1_1mm_stripped.transform;

M=[0.9975 -0.0073 0.0176 -0.0429 ; 
   0.0146 1.0009 -0.0024 1.5496 ; 
  -0.0130 -0.0093 0.9971 1.1840]; % https://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems

% mask_nifti_path = 'E:\RAM data set\RAM_Public_Data_all\Atlases\yale_broadmann_cortical.nii';
network_nifti_path = 'E:\RAM data set\RAM_Public_Data_all\Atlases\shen_1mm_268_parcellation_network_cortical.nii';
network_1mm_image = mni2fs_load_nii(network_nifti_path);
network_1mm_image_data=network_1mm_image.img;

atlas_nifti_path = 'E:\RAM data set\RAM_Public_Data_all\Atlases\shen_1mm_268_parcellation_cortical.nii';
atlas_1mm_image = mni2fs_load_nii(atlas_nifti_path);
atlas_1mm_image_data=atlas_1mm_image.img;
Atlas_index=unique(atlas_1mm_image_data);
Atlas_index(1,:)=[];
test=find(Atlas_index>0);
Atlas_all_val=cell(max(Atlas_index),251);

Network_name_all={'MedialFrontal','FrontoParietal','DefaultMode' ... 
    ,'SubcorticalCerebellum','Motor','Visual1','Visual2','VisualAssociation'};


actual_time=257:768;
% actual_time=257:512;
frequency_bands=[45 95; 3 8;40 115;13 30];

% design filter
band_i=3;
d_filter= designfilt('bandpassiir','FilterOrder',40, ...
'HalfPowerFrequency1',frequency_bands(band_i,1),'HalfPowerFrequency2',frequency_bands(band_i,2), ...
'SampleRate',fs);

prefix='final_clean';
for j=1:1 % FR1 or FR2 or FR3

    for i=1:251 % 251 subject
        try
            clearvars RAM_Rec_trials_sessions RAM_nonRec_trials_sessions location_data_pair

            cd(rootfolder)
            cd('FR1_final')
            cd([num2str(i),'_',r_sublist{i,1}]);

            load(['stats_traces_' prefix '.mat'])
            
            % Make the Montage
            
            masked_elec_n=size(Atlas_subject_elec_val,2);
            
            L_MontageMap=[];
            R_MontageMap=[];

            for elec_i=1:masked_elec_n
                masked_elec=Atlas_subject_elec_val(elec_i);
                channel_inf=location_data_pair{1,masked_elec};
                fsavg_coord=[channel_inf.atlases.avg.x channel_inf.atlases.avg.y channel_inf.atlases.avg.z]';
                mni_coord = M*[fsavg_coord ; 1];

                if(mni_coord(1)<0)
                    L_MontageMap(left_index,:)=[elec_i mni_coord(1) mni_coord(2) mni_coord(3)];
                    left_index=left_index+1;
                else
                    R_MontageMap(right_index,:)=[elec_i mni_coord(1) mni_coord(2) mni_coord(3)]; 
                    right_index=right_index+1;
                end
            end
            save('R_MontageMap_fs2mni.mat', 'R_MontageMap', 'R_MontageMap_ROI_name','R_MontageMap_ROI_hem');
            save('L_MontageMap_fs2mni.mat', 'L_MontageMap', 'L_MontageMap_ROI_name','L_MontageMap_ROI_hem');
        
            
            
            
            
            for elec_i=1:masked_elec_n
                masked_elec=Atlas_subject_elec_val(elec_i);
                channel_inf=location_data_pair{1,masked_elec};
                fsavg_coord=[channel_inf.atlases.avg.x channel_inf.atlases.avg.y channel_inf.atlases.avg.z]';
                mni_coord = M*[fsavg_coord ; 1];

                volume_buffer=[mni_coord(1) mni_coord(2) mni_coord(3) 1] / Tmni';
                volume_buffer(1:3)
                volume_buffer=round();
                
            end
            
            
            
            disp(sprintf(' %s subject completed!! ',r_sublist{i,1}));
        catch

        end
    end
end





%% plot network
clc;clear;
rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

load r1_all.mat
fid=fopen('Subjects_list_all.txt','r');
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);


% Network_name_all={'MedialFrontal','FrontoParietal','DefaultMode' ... 
%     ,'SubcorticalCerebellum','Motor','Visual1','Visual2','VisualAssociation'};

Network_name_all={'FrontoParietal(FEF)','FrontoParietal(IPL)','FrontoParietal(MTG)' ...
    ,'DefaultMode', 'VisualPrimary','VisualAssociation'};

prefix='final_clean_vol_thr_10';
Image_folder=['E:\RAM data set\RAM_Public_Data_all\FR1_final\Results_network_all_' prefix '_test'];
mkdir(Image_folder);

all_mean_network_zval=nan(251,6,39);
base_sample_idx=9;

for j=1:1 % FR1 or FR2 or FR3

    for i=1:251 % 251 subject
        try
            cd(rootfolder)
            cd('FR1_final')
            cd([num2str(i),'_',r_sublist{i,1}]);
            load(['stats_traces_' prefix '.mat'])  
            
            % threshold
            bad_elec_index = final_rejection_rate_all<0.7; 
            zscore_traces(bad_elec_index,1,3,:)=NaN;

            %z-score
            figure('position',[0 0 8000 4000]);
            set(gcf, 'color', [1 1 1]);
            set(gcf,'Visible','off');  
            
            for network_i=1:6
                subplot(2,3,network_i)
                
                network_electrodes_index=Atlas_subject_elec_val(find(Atlas_subject_network_val==network_i));
                if (network_electrodes_index)
                    network_zval=squeeze(zscore_traces(network_electrodes_index,1,3,:))';
                else
                    network_zval=[];
                end
                plot(network_zval);
                if(size(network_zval,1)>1)
                    mean_network_zval=nanmean(network_zval,2);
                    mean_network_zval=(mean_network_zval-mean(mean_network_zval(1:base_sample_idx)))./std(mean_network_zval(1:base_sample_idx));
                    all_mean_network_zval(i,network_i,:)=mean_network_zval;
                else
                    mean_network_zval=network_zval;
                end
                hold on;
%                 plot(mean_vals,'LineWidth',2,'color','k');
                plot(mean_network_zval,'LineWidth',2,'color','k');
                set(gca,'xtick',[10 20 30]);
                set(gca,'xticklabel',[0 500 1500]);
                set(gca,'xlim',[1 39])
                set(gca,'ylim',[-8 8])
                line([1,39],[2,2],'color','b');
                line([1,39],[-2,-2],'color','b');
                line([10,10],[-8,8],'color','k');
                set(gca,'Fontsize',15);
                title([Network_name_all{network_i}],'Fontsize',12)
                ylabel('Z score');
                xlabel('Time (ms)');

            end
            saveas(gcf,[Image_folder '\all_' r_sublist{i,1} '_zval.png'])
            close all;
            
        catch
        end
    end
end

save('all_mean_network_val.mat','all_mean_network_tval','all_mean_network_zval')
            
%% all plot (ttest)
% cd(rootfolder)
% cd('FR1_final')
% 
% figure('position',[0 0 8000 4000]);
% set(gcf, 'color', [1 1 1]);
% set(gcf,'Visible','on');  
% for network_i=1:6
%     subplot(2,3,network_i)
%     
%     test_all_mean_network_zval=squeeze(all_mean_network_tval(:,network_i,:));
%     test_all_mean_network_zval_base=squeeze(nanmean(all_mean_network_tval(:,network_i,1:9),3));
%     [H,P,CI,STATS]=ttest(test_all_mean_network_zval - test_all_mean_network_zval_base);
%     plot(squeeze(STATS.tstat),'LineWidth',2,'color','k');
%     set(gca,'xtick',[10 20 30]);
%     set(gca,'xticklabel',[0 500 1500]);
%     set(gca,'xlim',[1 39])
%     set(gca,'ylim',[-8 8])
%     line([1,39],[2,2],'color','b');
%     line([1,39],[-2,-2],'color','b');
%     line([10,10],[-8,8],'color','k');
%     set(gca,'Fontsize',15);
%     title([Network_name_all{network_i}],'Fontsize',12)
%     ylabel('T value');
%     xlabel('Time (ms)');
% 
% end
% saveas(gcf,['network_plot_tval_from_tval_all_clean_sum.png'])
% close all;


cd(rootfolder)
cd('FR1_final')

figure('position',[0 0 8000 4000]);
set(gcf, 'color', [1 1 1]);
set(gcf,'Visible','on');  
for network_i=1:6
    subplot(2,3,network_i)
    
    test_all_mean_network_zval=squeeze(all_mean_network_zval(:,network_i,:));
%     test_all_mean_network_zval_base=squeeze(nanmean(all_mean_network_zval(:,network_i,1:9),3));
%     [H,P,CI,STATS]=ttest(test_all_mean_network_zval - test_all_mean_network_zval_base);
%     plot(squeeze(STATS.tstat),'LineWidth',2,'color','k');
    (mean_network_zval-mean(mean_network_zval(1:9)))./std(mean_network_zval(1:9));
    mean_test_all_mean_network_zval=nanmean(test_all_mean_network_zval);
    
    plot(squeeze(),'LineWidth',2,'color','k');
    set(gca,'xtick',[10 20 30]);
    set(gca,'xticklabel',[0 500 1500]);
    set(gca,'xlim',[1 39])
    set(gca,'ylim',[-8 8])
    line([1,39],[2,2],'color','b');
    line([1,39],[-2,-2],'color','b');
    line([10,10],[-8,8],'color','k');
    set(gca,'Fontsize',15);
    title([Network_name_all{network_i}],'Fontsize',12)
%     ylabel('T value');
    ylabel('z score');
    xlabel('Time (ms)');

end
saveas(gcf,['network_plot_tval_from_zval_all_clean_sum.png'])
close all;


