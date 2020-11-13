clear all
close all
clc

file='/Volumes/LaCie/DATOS_OLFATO_HELENA_AMATOMIA/PRUEBAS_PIPELINE/DATOS_OLFATO_HELENA_AMATOMIA_PRUEBAS_PIPELINE_datafiles.txt';
fid = fopen(file,'r');
Text = fread(fid,'*char');
Text=Text';
Info=strsplit(Text,' ');
count=1;
for i=1:3:length(Info)
    grupo1(count).file=Info{i};
    grupo1(count).age=str2num(Info{i+1});
    grupo1(count).sex=Info{i+2};
    count=count+1;
end 