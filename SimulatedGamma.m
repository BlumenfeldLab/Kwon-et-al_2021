clc;clear; close all;

% Sampling
fs = 500;     % Sampling rate [Hz]
Ts = 1/fs;     % Sampling period [s]
fNy = fs / 2;  % Nyquist frequency [Hz]
duration = 1; % Duration [s]
t = 0 : Ts : duration-Ts; % Time vector
t_50 = 0.5 : Ts : 0.55-Ts; % Time vector
t_150 = 0.6 : Ts : 0.65-Ts; % Time vector
noSamples = length(t);    % Number of samples

d_filter= designfilt('bandpassiir','FilterOrder',40, ...
'HalfPowerFrequency1',40,'HalfPowerFrequency2',115, ...
'SampleRate',fs);

for num=1:100

% simulation data
simulation_data=zeros(1,fs*duration);
for f=1:140
    simulation_data = simulation_data + (rand-0.5)*sin((2*pi*f*t));
end

% gamma power (0-50ms)
simulation_data_gamma=zeros(1,fs*duration);
% simulation_data_gamma_buff=0;
% for f=40:115
%     simulation_data_gamma_buff = simulation_data_gamma_buff+ (rand-0.5)*10*sin((2*pi*f*t_50));
% end
% simulation_data_gamma(251:275)=simulation_data_gamma_buff;

% gamma power (100-150ms)
simulation_data_gamma2=zeros(1,fs*duration);
simulation_data_gamma_buff=0;
for f=40:115
    simulation_data_gamma_buff = simulation_data_gamma_buff+ (rand-0.5)*10*sin((2*pi*f*t_150));
end
simulation_data_gamma2(301:325)=simulation_data_gamma_buff;

simulation_data_final=simulation_data+simulation_data_gamma+simulation_data_gamma2;

figure('position',[0 0 500 1500],'visible','on');
set(gcf,'color', [1 1 1]);

subplot(4,1,1)
plot(simulation_data_final)
xlim([126 375])
xticks([126 150 175 200 225 250 275  300  325  350 375])
xticklabels({'-250','200','-150','-100','-50','0','50','100','150','200','250'})
xlabel('Time(ms)')
ylabel('Magnitude')
xline(250,'-r','LineWidth',2)
y1=min(get(gca,'ylim')); y2=max(get(gca,'ylim'));
ylim([y1*1.1 y2*1.1])
title('Simulated data, fs=500')



subplot(4,1,2)
simulation_data_final_filtered = filtfilt(d_filter,simulation_data_final);
plot(simulation_data_final_filtered)
xlim([126 375])
xticks([126 150 175 200 225 250 275  300  325  350 375])
xticklabels({'-250','200','-150','-100','-50','0','50','100','150','200','250'})
xlabel('Time(ms)')
ylabel('Magnitude')
xline(250,'-r','LineWidth',2)
y1=min(get(gca,'ylim')); y2=max(get(gca,'ylim'));
ylim([y1*1.1 y2*1.1])
title('Filtered (40-115Hz)')

subplot(4,1,3)
simulation_data_final_gamma=simulation_data_final_filtered.^2;
plot(simulation_data_final_gamma)
xlim([126 375])
xticks([126 150 175 200 225 250 275  300  325  350 375])
xticklabels({'-250','200','-150','-100','-50','0','50','100','150','200','250'})
xlabel('Time(ms)')
ylabel('Power(Magnitude^2)')
xline(250,'-r','LineWidth',2)
y1=min(get(gca,'ylim')); y2=max(get(gca,'ylim'));
ylim([y1*1.1 y2*1.1])
title('Gamma power')

subplot(4,1,4)
average_window=[];
cnt=1;
window_size_overlap=12.5;
time_x=[125:window_size_overlap:250-2*window_size_overlap 250:window_size_overlap:375-2*window_size_overlap];
for i=time_x
    average_window(cnt)=nanmean(simulation_data_final_gamma(i:i+round(window_size_overlap)*2-1));
    cnt=cnt+1;
end
plot(time_x+12.5,average_window,'.k','MarkerSize',20)
xlim([125 375])
xticks([time_x]+12.5)
xticklabels({'-250, -200','-225, -175','-200 150','-175, -125','-150,-100','-125, -75','-100, -50','-75 -25','-50, 0', ...
    '0, 50','25, 75','50, 100','75, 125','100, 150','125, 175','150, 200','175, 225','200, 250'})
xtickangle(45)
xlabel('Time(ms)')
ylabel('Power(Magnitude^2)')
y1=min(get(gca,'ylim')); y2=max(get(gca,'ylim'));
xline(250,'-r','LineWidth',2)
ylim([y1*1.1 y2*1.1])
title('Averaged power within 50ms windows with 25ms overlap')
print('-dtiff','-r300',['SimulationForFilter_100_150_' num2str(num)])
close all
end
