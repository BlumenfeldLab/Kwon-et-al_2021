
clc;clear;

rootfolder='E:\RAM data set\RAM_Public_Data_all\';
cd(rootfolder)

load r1_all.mat
fid=fopen('Subjects_list_all.txt','r');
for i=1:251
    r_sublist{i,1}=fgetl(fid);
end
fclose(fid);

all_Categories={'Bad Electrodes','Bad electrodes',...
    'Brain Lesions','Brain lesions','Broken leads','Broken Leads',...
    'Brain lesions:','Brain Lesions:','Broken leads:','Broken Leads:'};

for i=1:251 % 251 subject
    try
        cd(rootfolder)
        cd('FR1_FS')
        cd([num2str(i),'_',r_sublist{i,1}]);

        folder = 'E:\RAM data set\RAM_Public_Data_all\electrode_categories_all\final';
        fileList = dir(fullfile(folder, ['*' r_sublist{i,1} '*']));
        fid=fopen([folder '\' fileList.name],'r');
        tline = fgetl(fid);
        Bad_electrode_categories=[];
        cnt=1;
        try
        while ischar(tline)
                tline = fgetl(fid);    
                if ~isempty(tline) & sum(ismember(all_Categories,tline)) & ischar(tline)
                    tline=fgetl(fid);
                    if isempty(tline) & ischar(tline)
                        tline=fgetl(fid);
                        while ~isempty(tline) & ischar(tline)
                            Bad_electrode_categories{cnt,1} = tline;
                            cnt=cnt+1;
                            tline = fgetl(fid);
                        end
                    else
                        while ~isempty(tline) & ischar(tline)
                            Bad_electrode_categories{cnt,1} = tline;
                            cnt=cnt+1;
                            tline = fgetl(fid);
                        end
                    end
                end
        end
        catch
        end
        fclose(fid);
        save('Bad_electrode_categories.mat','Bad_electrode_categories','-v7.3')
        disp(sprintf(' %s subject complete',r_sublist{i,1}));
    catch
    end
end



for i=1:251 % 251 subject
    try
        cd(rootfolder)
        cd('FR1_FS')
        cd([num2str(i),'_',r_sublist{i,1}]);
%             load('stats_traces_final.mat')  
        load('Bad_electrode_categories.mat')  
        load('Raw_data_extended_5s.mat','location_data','location_data_pair') 

        Bad_electrode_categories_buffer=[];

        for ii=1:size(location_data,2)
            channel_inf=[];
            channel_inf=location_data{1,ii};
            flag=[];
            flag=sum(ismember(Bad_electrode_categories,channel_inf.code));
            if flag>0
                Bad_electrode_categories_buffer=[Bad_electrode_categories_buffer channel_inf.channel];
            end
        end

        Bad_electrode_categories_index=[];
        for ii=1:size(location_data_pair,2)
            channel_inf=[];
            channel_inf=location_data_pair{1,ii};
            flag=[];
            flag=sum(ismember(Bad_electrode_categories_buffer,[channel_inf.channel_1 channel_inf.channel_2]));

            if flag>0
                Bad_electrode_categories_index=[Bad_electrode_categories_index ii];
            end
        end

        save('Bad_electrode_categories_index.mat','Bad_electrode_categories_index','-v7.3')
        disp(sprintf(' %s subject complete',r_sublist{i,1}));
    catch
    end
end





