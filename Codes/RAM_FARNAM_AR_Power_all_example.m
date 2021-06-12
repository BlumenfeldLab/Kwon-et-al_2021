
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
data_location_final='E:\RAM data set\RAM_Public_Data_all\FR1_final';

for subject_i=1:10
    try
        cd(data_location)
        RAM_FARNAM_AR_Power_all(subject_i,r_sublist{subject_i,1})
        close all;
    catch
    end
end


%% check bad electrodes
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
data_location_final='E:\RAM data set\RAM_Public_Data_all\FR1_final';

all_subject_noisy_elec_ratio=nan(251,1);
all_subject_noisy_elec_ratio_x=[1:1:251];
for subject_i=1:251
    try
        cd(data_location_final)
        cd([num2str(subject_i),'_',r_sublist{subject_i,1}]);
        load('L_MontageMap_final_clean_dil.mat')
        load('R_MontageMap_final_clean_dil.mat')
        
        cd(data_location)
        cd([num2str(subject_i),'_',r_sublist{subject_i,1}]);
        load('stats_traces_Wendy_hunki_2.mat')
        
        noisy_RAM_all_trials_plot=[];
        
        if ~isempty(L_MontageMap)
            plot_buffer=noisy_RAM_all_trials(L_MontageMap(:,1),:);
            noisy_RAM_all_trials_plot=[noisy_RAM_all_trials_plot;plot_buffer];
        end
        
        if ~isempty(R_MontageMap)
            plot_buffer=noisy_RAM_all_trials(R_MontageMap(:,1),:);
            noisy_RAM_all_trials_plot=[noisy_RAM_all_trials_plot;plot_buffer];
        end
        figure('position',[0 0 2000 2000]);
        set(gcf, 'color', [1 1 1]);
        set(gcf,'Visible','off');  
        scatter(noisy_RAM_all_trials_plot(:,1)',noisy_RAM_all_trials_plot(:,2)','filled')
        line([1,max(noisy_RAM_all_trials_plot(:,1))],[0.2,0.2],'color','b');
        line([100,100],[0,1],'color','b');
        set(gca,'Fontsize',15);
        noise_elecs=noisy_RAM_all_trials_plot(:,1)<100 | noisy_RAM_all_trials_plot(:,2)<0.2 | noisy_RAM_all_trials_plot(:,2)==NaN;
        noise_elecs_size=size(find(noise_elecs>0),1);
        title(['Total : ' num2str(size(noisy_RAM_all_trials_plot,1)) ',    Number of artifact electrodes : ' num2str(noise_elecs_size) ],'Fontsize',12)
        ylabel('Ratio of remaining trials');
        xlabel('Number of remaining trials');
        
        cd('E:\RAM data set\RAM_Public_Data_all\FR1_FARNAM\noisy_electrode_plot')
        saveas(gcf,['noisy_electrode_plot_' num2str(subject_i),'_',r_sublist{subject_i,1} '.png'])
        close all;
        
        all_subject_noisy_elec_ratio(subject_i,1)=(size(noisy_RAM_all_trials_plot,1)-noise_elecs_size)/size(noisy_RAM_all_trials_plot,1);
        
    catch
    end
end

figure('position',[0 0 2000 2000]);
set(gcf, 'color', [1 1 1]);
set(gcf,'Visible','on');  
scatter(all_subject_noisy_elec_ratio_x,all_subject_noisy_elec_ratio,'filled')
set(gca,'xlim',[1 220])
line([1,251],[0.8,0.8],'color','k');
set(gca,'Fontsize',15);
% title(' : ','Fontsize',12)
ylabel('Ratio of remaining electrodes');
xlabel('Subject index');

saveas(gcf,'all_subject_noisy_electrode_plot.png')
close all





%% check zscore (electrode level) // define threshold 
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
data_location_final='E:\RAM data set\RAM_Public_Data_all\FR1_final';

prefix='Wendy_hunki_std_5_20';

zscore_elec_all=[];
for subject_i=1:251
    try
        cd(data_location_final)
        cd([num2str(subject_i),'_',r_sublist{subject_i,1}]);
        load('L_MontageMap_final_clean_dil.mat')
        load('R_MontageMap_final_clean_dil.mat')
        
        cd(data_location)
        cd([num2str(subject_i),'_',r_sublist{subject_i,1}]);
        load(['stats_traces_' prefix '.mat'])

        if ~isempty(L_MontageMap)
            noisy_RAM_all_trials_buffer=[];
            noisy_RAM_all_trials_buffer=noisy_RAM_all_trials(L_MontageMap(:,1),:);
            
            noise_elec_index=[];
            noise_elec_index=noisy_RAM_all_trials_buffer(:,1)<100 | noisy_RAM_all_trials_buffer(:,2)<0.2 | noisy_RAM_all_trials_buffer(:,2)==NaN;
            noise_elec_index=find(noise_elec_index>0);
            L_MontageMap(noise_elec_index,:)=[];
            
            zscore_elec_buffer=[];
            zscore_elec_buffer=squeeze(zscore_traces(L_MontageMap(:,1),1,3,:));
            zscore_elec_all=[zscore_elec_all;zscore_elec_buffer];
        end
        
        if ~isempty(R_MontageMap)
            noisy_RAM_all_trials_buffer=[];
            noisy_RAM_all_trials_buffer=noisy_RAM_all_trials(R_MontageMap(:,1),:);
            
            noise_elec_index=[];
            noise_elec_index=noisy_RAM_all_trials_buffer(:,1)<100 | noisy_RAM_all_trials_buffer(:,2)<0.2 | noisy_RAM_all_trials_buffer(:,2)==NaN;
            noise_elec_index=find(noise_elec_index>0);
            R_MontageMap(noise_elec_index,:)=[];
            
            zscore_elec_buffer=[];
            zscore_elec_buffer=squeeze(zscore_traces(R_MontageMap(:,1),1,3,:));
            zscore_elec_all=[zscore_elec_all;zscore_elec_buffer];
        end
        subject_i
%         save(['R_MontageMap_' prefix '.mat'], 'R_MontageMap');
%         save(['L_MontageMap_' prefix '.mat'], 'L_MontageMap');
    catch
    end
end


rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

figure('position',[0 0 2000 2000]);
set(gcf, 'color', [1 1 1]);
set(gcf,'Visible','on');

% Plot histogram of the vertex values
bins = 200;
subplot(2,1,1)
H = histogram([zscore_elec_all(:)], bins,...
    'Normalization', 'probability','BinLimits',[-20,120]);
% H = histogram([log(abs(vertex_values_sumL(:))); log(abs(vertex_values_sumR(:)))], bins,...
%     'Normalization', 'probability');
title('Original')
xlabel('z-score gamma power','Fontsize',12)
ylabel('Normalized counts','Fontsize',12)

subplot(2,1,2)
TIME = T-0.5;
plot(TIME,zscore_elec_all')
ylim([-20,120])
xlabel('Time(ms)','Fontsize',12)
ylabel('z-score gamma power','Fontsize',12)

saveas(gcf,'zscore_histogram_elec.png')
close all;

%% Compute the outer boundaries to retain some percentage of the data in the
% visualization
[~, zero_init] = min(abs(H.BinEdges));
alpha = .05;
bounds = [zero_init - 2, zero_init + 2];
loss = sum(H.Values([1 : bounds(1), bounds(2) : end]));
while loss > alpha
    
    % Compute half-probability distribution partitions
    tail_loss_upper = sum(H.Values(bounds(2) - 1 : end));
    trunk_loss_upper = sum(H.Values(zero_init + 1 : bounds(2) - 1));
    tail_loss_lower = sum(H.Values(1 : bounds(1) + 1));
    trunk_loss_lower = sum(H.Values(bounds(1) + 1 : zero_init - 1));
    
    % Compute the partition ratios for each bound
    upper_ratio = trunk_loss_upper / tail_loss_upper;
    lower_ratio = trunk_loss_lower / tail_loss_lower;
    
    if upper_ratio > lower_ratio
        bounds(1) = bounds(1) - 1;
    else
        bounds(2) = bounds(2) + 1;
    end
    loss = sum(H.Values([1 : bounds(1), bounds(2) : end]));
end
color_bar_outer_boundary = H.BinEdges([bounds(1), bounds(2) + 1]); % a and b
tail_percentages = [sum(H.Values(1 : bounds(1))), sum(H.Values(bounds(2) : end))];









%% zscore to vertex
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
cd(data_location)

for subject_i=1:251
    try
        cd(data_location)
        RAM_FARNAM_AR_Power_all_vertex(subject_i,r_sublist{subject_i,1})
    catch
        
    end
end




