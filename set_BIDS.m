%% Date : 22 - Feb - 2019
% Author : Mario Gil 
%
% PROYECTO AGES-CM II
% Preparing of BIDS 
% In this scrip the main folder names will be treated as variables.
%
%% Initialization
clear all
clc
%% Main Folder paths:

%   This variables has to be changed for every compueter where it is executed
%   or any other study.
main = '/Users/administrator/Documents/SENECA/SENECA-PICASO';

%% -------------------------------------------------------------------------- 
% Variables declaration:

anat     = strcat(filesep,'anat');
T1w      = strcat(filesep,'T1w');
T2w      = strcat(filesep,'T2w');
func     = strcat(filesep,'func');
%En este estudio solamente tenemos imágene sfuncionales de estado de
%reposo: restng state fMRI.
rest     = strcat(filesep,'task-rest_bold');
dwi      = strcat (filesep,'dwi');
fieldmap = strcat(filesep,'fmap');
% Para el proyecto AGES, el fieldmap es de tipo: Two TE with
% Real-Imaginary images.
fmap_dwi_ri = strcat(filesep, 'acq-dwi_fmap');
fmap_rest_ri = strcat(filesep, 'acq-rest_fmap');
%Existen otro tipo de fieldmaps: Oposite Encoding directions & Two TE with
%Phase-Magnitude images.
% fmap_dwi_epi = strcat(filesep,'acq-dwi_dir-PA_epi');
% fmap_func_mag = strcat(filesep,'acq-rest_magnitude');
% fmap_func_phs = strcat(filesep,'acq-rest_phasediff');


%% 
for num=10:25
    % num=13
    numero= sprintf('sub-sp%d*',num)
    Subj=dir(strcat(main,filesep,numero));

    for i = 1:length(Subj)
        fprintf('\nStarting Subject: %s\n',Subj(i).name)
        Ses = dir(fullfile(main,Subj(i).name,'ses-*')); % take number of files (sesions) in the directory
        for p = 1:length(Ses)
    %     %if not(isempty(dir(strcat(main,filesep,Subj(i).name,filesep,'Ages_Ii*'))))        
    %         fprintf('\nRe-Structuring folder system into BIDS standard\n')    
    % 
    %         %Find the name of folder containing the sequences (It changes for each subject) 
    %         fileA = dir(strcat(main,filesep,Subj(i).name,filesep,'Ages_Ii*'));
    %         
    %         %Find the name of the folder contaning the images (it can be
    %         %different at the end for each subject)        
    %         fileB = dir(strcat(main,filesep,Subj(i).name,filesep,fileA.name,filesep,'SAG_IRSPGR*'));
    %         if not(isdir(strcat(main,filesep,Subj(i).name,anat)))
    %             mkdir(strcat(main,filesep,Subj(i).name,anat))
    %         end;
    %         if not(isempty(fileB))
    %             fprintf('\nSAG_IRSPGR founded. Converting into T1w')
    %             %First we need to rename the folder in the same parent folder
    %             movefile(strcat(main,filesep,Subj(i).name,filesep,fileA.name,filesep,fileB.name),strcat(main,filesep,Subj(i).name,filesep,fileA.name,T1w));
    %             %Then we can move that folder to another folder
    %             movefile(strcat(main,filesep,Subj(i).name,filesep,fileA.name,T1w),strcat(main,filesep,Subj(i).name,anat,T1w));
    %         end;
    %         clear fileB
    %         
    %         %Repeat for each sequence.
    %         fileB = dir(strcat(main,filesep,Subj(i).name,filesep,fileA.name,filesep,'T2W_FSE*'));
    %         if not(isempty(fileB))
    %             fprintf('\nT2W_FSE founded. Converting into T2w')            
    %             movefile(strcat(main,filesep,Subj(i).name,filesep,fileA.name,filesep,fileB.name),strcat(main,filesep,Subj(i).name,filesep,fileA.name,T2w));
    %             movefile(strcat(main,filesep,Subj(i).name,filesep,fileA.name,T2w),strcat(main,filesep,Subj(i).name,anat,T2w));
    %         end;
    %         clear fileB
    %         if not(isdir(strcat(main,filesep,Subj(i).name,func)))
    %             mkdir(strcat(main,filesep,Subj(i).name,func))
    %         end;
    %         fileB = dir(strcat(main,filesep,Subj(i).name,filesep,fileA.name,filesep,'rs_fMRI*'));
    %         if not(isempty(fileB))
    %             fprintf('\nrs_fMRI founded. Converting into task-rest_bold')
    %             movefile(strcat(main,filesep,Subj(i).name,filesep,fileA.name,filesep,fileB.name),strcat(main,filesep,Subj(i).name,filesep,fileA.name,rest));
    %             movefile(strcat(main,filesep,Subj(i).name,filesep,fileA.name,rest),strcat(main,filesep,Subj(i).name,func,rest));        
    %         end;
    %         clear fileB
    %         if not(isdir(strcat(main,filesep,Subj(i).name,dwi)))
    %             mkdir(strcat(main,filesep,Subj(i).name,dwi))
    %         end;
    %         fileB = dir(strcat(main,filesep,Subj(i).name,filesep,fileA.name,filesep,'DTI_64_dir*'));
    %         if not(isempty(fileB))
    %             fprintf('\nDTI_64_dir founded. Converting into dir-AP_dwi')
    %             movefile(strcat(main,filesep,Subj(i).name,filesep,fileA.name,filesep,fileB.name),strcat(main,filesep,Subj(i).name,filesep,fileA.name,filesep,'dir-AP_dwi'));
    %             movefile(strcat(main,filesep,Subj(i).name,filesep,fileA.name,filesep,'dir-AP_dwi'),strcat(main,filesep,Subj(i).name,dwi,filesep,'dir-AP_dwi'));
    %         end;
    %         clear fileB
    %         fileB = dir(strcat(main,filesep,Subj(i).name,filesep,fileA.name,filesep,'Fractional_Ansio*'));
    %         if not(isempty(fileB))
    %             fprintf('\nFractional_Ansio founded. Converting into dir-AP_FA')
    %             movefile(strcat(main,filesep,Subj(i).name,filesep,fileA.name,filesep,fileB.name),strcat(main,filesep,Subj(i).name,filesep,fileA.name,filesep,'dir-AP_FA'));
    %             movefile(strcat(main,filesep,Subj(i).name,filesep,fileA.name,filesep,'dir-AP_FA'),strcat(main,filesep,Subj(i).name,dwi,filesep,'dir-AP_FA'));
    %         end;
    %         clear fileB
    %         if not(isdir(strcat(main,filesep,Subj(i).name,fieldmap)))
    %             mkdir(strcat(main,filesep,Subj(i).name,fieldmap))
    %         end;
    %         fileB = dir(strcat(main,filesep,Subj(i).name,filesep,Ses(p).name,filesep,'fmap',filesep,'Field_map_RL*'));
    %         if not(isempty(fileB))
    %             fprintf('\nFieldmap for Resting founded. Converting into acq-rest_fmap')
    %             movefile(strcat(main,filesep,Subj(i).name,filesep,Ses(p).name,filesep,'fmap',filesep,fileB.name),strcat(main,filesep,Subj(i).name,filesep,Ses(p).name,filesep,'fmap',fmap_rest_ri));
    %             movefile(strcat(main,filesep,Subj(i).name,filesep,Ses(p).name,filesep,'fmap',fmap_rest_ri),strcat(main,filesep,Subj(i).name,filesep,Ses(p).name,fieldmap,fmap_rest_ri));
    %         end;
    %         clear fileB
    %         fileB = dir(strcat(main,filesep,Subj(i).name,filesep,fileA.name,filesep,'Field_map_tenssor*'));
    %         if not(isempty(fileB))
    %             fprintf('\nFieldmap dor DTI founded. Converting into acq-dwi_fmap')
    %             movefile(strcat(main,filesep,Subj(i).name,filesep,fileA.name,filesep,fileB.name),strcat(main,filesep,Subj(i).name,filesep,fileA.name,fmap_dwi_ri));
    %             movefile(strcat(main,filesep,Subj(i).name,filesep,fileA.name,fmap_dwi_ri),strcat(main,filesep,Subj(i).name,fieldmap,fmap_dwi_ri));
    %         end;
    %         clear fileB 
    %         %rmdir(strcat(main,filesep,Subj(i).name,filesep,fileA.name),'s')
    %         clear fileA
    %     end;
     %% 
     %En esta parte del codigo se reordenan las carpetas de fieldmap ya que 
     %contienen todas las imágenes juntas. Se separan por Tiempo de Eco (TE1 y TE2) y por tipo de imagen (real, imaginaria, magnitud)

%         if isdir ([main,filesep,Subj(i).name,filesep,Ses(p).name,fieldmap,fmap_rest_ri])
        if isdir ([main,filesep,Subj(i).name,filesep,Ses(p).name,fieldmap,fmap_rest_ri])
            fprintf('\nReal-Imaginary fieldmap founded for functional acquisition.\nSorting into Real, Imaginary and Magnitude Images\n')
                fmap_real = strcat(filesep,'acq-rest_real');
                fmap_imaginary = strcat(filesep,'acq-rest_imag');
                fmap_magnitude = strcat(filesep,'acq-rest_magnitude');
                if isdir([main,filesep,Subj(i).name,filesep,Ses(p).name,fieldmap,fmap_rest_ri]) && not(isdir([main,filesep,Subj(i).name,filesep,Ses(p).name,fieldmap,fmap_real,fmap_rest_ri]))
                    % CHECKED that there IS a directory to the images and 
                    mkdir([main filesep Subj(i).name filesep Ses(p).name fieldmap fmap_real]);
                    mkdir([main filesep Subj(i).name filesep Ses(p).name fieldmap fmap_imaginary]);
                    mkdir([main filesep Subj(i).name filesep Ses(p).name fieldmap fmap_magnitude]);
                    fmap_img = dir(strcat(main,filesep,Subj(i).name,filesep, Ses(p).name,fieldmap,filesep,'acq-rest_fmap',filesep,'*.dcm'));
                    for k = 1 : size(fmap_img,1)
                        tmp = strcat(main,filesep,Subj(i).name,filesep,Ses(p).name,fieldmap,filesep,'acq-rest_fmap',filesep,fmap_img(k).name);
                        info = dicominfo(tmp);
    %                   Multi-echo (e.g. QSM, R2*) data from GE scanner
    %                   Campo DICOM: (0043,102F) 
    %                   (0=Mag, 1=Phase, 2=Real, 3=Imag).
                        if info.Private_0043_102f(1) == 0
                            movefile(strcat(main,filesep,Subj(i).name,filesep,Ses(p).name,fieldmap,filesep,'acq-rest_fmap',filesep,fmap_img(k).name),strcat(main,filesep,Subj(i).name,filesep,Ses(p).name,fieldmap,filesep,fmap_magnitude,filesep,fmap_img(k).name));
                        end;
                        if info.Private_0043_102f(1) == 2
                            movefile(strcat(main,filesep,Subj(i).name,filesep,Ses(p).name,fieldmap,filesep,'acq-rest_fmap',filesep,fmap_img(k).name),strcat(main,filesep,Subj(i).name,filesep,Ses(p).name,fieldmap,filesep,fmap_real,filesep,fmap_img(k).name));                    
                        end;
                        if info.Private_0043_102f(1) == 3
                            movefile(strcat(main,filesep,Subj(i).name,filesep,Ses(p).name,fieldmap,filesep,'acq-rest_fmap',filesep,fmap_img(k).name),strcat(main,filesep,Subj(i).name,filesep,Ses(p).name,fieldmap,filesep,fmap_imaginary,filesep,fmap_img(k).name));                    
                        end;
                    end;
                    rmdir([main,filesep,Subj(i).name,filesep,Ses(p).name,fieldmap,fmap_rest_ri])
                end;
        end;
        clear fmap_real fmap_imaginary fmap_magnitude
        if isdir ([main,filesep,Subj(i).name,fieldmap,fmap_dwi_ri])
            fprintf('\nReal-Imaginary fieldmap founded for DTI acquisition.\nSorting into Real, Imaginary and Magnitude Images\n')
                fmap_real = strcat(filesep,'acq-dwi_real');
                fmap_imaginary = strcat(filesep,'acq-dwi_imag');
                fmap_magnitude = strcat(filesep,'acq-dwi_magn');
                if isdir([main,filesep,Subj(i).name,fieldmap,fmap_dwi_ri]) && not(isdir([main,filesep,Subj(i).name,fieldmap,fmap_real]))
                    mkdir([main filesep Subj(i).name fieldmap fmap_real]);
                    mkdir([main filesep Subj(i).name fieldmap fmap_imaginary]);
                    mkdir([main filesep Subj(i).name fieldmap fmap_magnitude]);
                    fmap_img = dir(strcat(main,filesep,Subj(i).name,fieldmap,filesep,'acq-dwi_fmap',filesep,'*.dcm'));
                    for k = 1 : size(fmap_img,1)
                        tmp = strcat(main,filesep,Subj(i).name,fieldmap,filesep,'acq-dwi_fmap',filesep,fmap_img(k).name);
                        info = dicominfo(tmp);
    %                   Multi-echo (e.g. QSM, R2*) data from GE scanner
    %                   [GEMS_PARM_01] (0043,102F) 
    %                   (0=Mag, 1=Phase, 2=Real, 3=Imag).
                        if info.Private_0043_102f(1) == 0
                            movefile(strcat(main,filesep,Subj(i).name,fieldmap,filesep,'acq-dwi_fmap',filesep,fmap_img(k).name),strcat(main,filesep,Subj(i).name,fieldmap,filesep,fmap_magnitude,filesep,fmap_img(k).name));
                        end;
                        if info.Private_0043_102f(1) == 2
                            movefile(strcat(main,filesep,Subj(i).name,fieldmap,filesep,'acq-dwi_fmap',filesep,fmap_img(k).name),strcat(main,filesep,Subj(i).name,fieldmap,filesep,fmap_real,filesep,fmap_img(k).name));                    
                        end;
                        if info.Private_0043_102f(1) == 3
                            movefile(strcat(main,filesep,Subj(i).name,fieldmap,filesep,'acq-dwi_fmap',filesep,fmap_img(k).name),strcat(main,filesep,Subj(i).name,fieldmap,filesep,fmap_imaginary,filesep,fmap_img(k).name));                    
                        end;
                    end;
                    rmdir([main,filesep,Subj(i).name,fieldmap,fmap_dwi_ri])
                end;
        end;
        end;
        fprintf('\nSubject %s DONE!\n',Subj(i).name)
    end;
end    