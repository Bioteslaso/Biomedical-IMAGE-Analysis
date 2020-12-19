%% Date : 18 - Oct - 2016
% Author : Mario Gilipichis 
%
% In this scrip the main folder names will be treated as variables.
%
%% The inputs needed are:
%   - Origin directory of subjects with their corresponding dicoms.
%       Note: The directory must have the BIDS structure.
%       >FolderX
%           >sub-idxxxx
%               >anat
%                   >T1w
%                       dicom images
%               >func
%                   >task-rest
%                       dicom images
%               >modality
%                   >   
%
%   - Output directory
%   - Folder where the dcm2nii method is stored.
%   - Folder where the STMP text files are stored.
%
%% Initialization
clear all
clc
%% Main Folder paths:

%   This variables has to be changed for every compueter where it is executed
%   or any other study.
main = '/Users/administrator/Documents/SENECA/SENECA-PICASO';
mainOut = '/Users/administrator/Documents/SENECA/SP_nifti';
mainDCM2NII= '/Applications/MRIcroGL.app/Contents/Resources';


% Variables declaration:

anat     = strcat(filesep,'anat');
T1w      = strcat(filesep,'T1w');
T2w      = strcat(filesep,'T2w');
func     = strcat(filesep,'func');
fieldmap = strcat(filesep,'fmap');
fmap_dwi = strcat(filesep,'acq-dwi_dir-PA_epi');
fmap_func_mag = strcat(filesep,'acq-rest_magnitude');
fmap_func_phs = strcat(filesep,'acq-rest_phasediff');
rest     = strcat(filesep,'task-rest_bold');
dwi      = strcat (filesep,'dwi');
asl      = strcat (filesep,'asl');
relCBF   = strcat (filesep,'relCBF');

DCM2NII = strcat (mainDCM2NII,filesep,'dcm2niix');
Subj=dir(strcat(main,filesep,'sub-sp*'));

%--------------------------------------------------------------------------
% Dicom to NifTi conversion

%% IDEACA
%Dirs = dir(strcat(main,filesep,'sub*/**'));
%dirFlags = [Dirs.isdir] & ~strcmp({Dirs.name},'.') & ~strcmp({Dirs.name},'..');
% 
for i = 1:length(Subj)
    
    %check for the existance of different sessions within a subject.
    Ses = dir(fullfile(main,Subj(i).name,'ses-*'));
    Ses_flag = isempty(Ses);
    
     %if not(isdir(strcat(mainOut,filesep,Subj(i).name))) %Check if the participant has been already converted from DICOM to Nifti.
         
        
         
         
         
%% ------------- CAT12 ----------------------------------------------------
%         if Ses_flag % Si no hay sesiones dentro de cada sujeto
%             command = strcat(DCM2NII,{' -f "'},Subj(i).name,{'_%f" -b n -z n -x y -o '},mainOut,{' '},main,filesep,Subj(i).name,anat,T1w);
%             system(command{1});
%             fileA = dir(strcat(mainOut,filesep,'*Crop*'));
%             command = strcat({'mv '},mainOut,filesep,fileA.name,{' '},mainOut,filesep,Subj(i).name,'_acq-crop_T1w.nii');
%             system(command{1});
%             clear fileA
%         else % Si hay una o mas sesiones dentro de cada sujeto
%             for j = 1:length(Ses)
%                 command = strcat(DCM2NII,{' -f "'},Subj(i).name,{'_%f" -b n -z n -x y -o '},mainOut,{' '},main,filesep,Subj(i).name,filesep,Ses(j).name,anat,T1w);
%                 system(command{1});
%                 fileA = dir(strcat(mainOut,filesep,'*Crop*'));
%                 command = strcat({'mv '},mainOut,filesep,fileA.name,{' '},mainOut,filesep,Subj(i).name,'_acq-crop_T1w.nii');
%                 system(command{1});
%                 clear fileA
%             end;
%         end;

%% ------------- ANAT -----------------------------------------------------
        if Ses_flag % Si no hay sesiones dentro de cada sujeto
            if isdir([main,filesep,Subj(i).name,anat]) % Comprueba que haya adquisiciones anatómicas.
                if not(isdir ([mainOut filesep Subj(i).name anat]))
                    mkdir([mainOut filesep Subj(i).name anat]);
                    command = strcat(DCM2NII,{' -f "'},Subj(i).name,{'_%f" -b n -z n -x y -o '},mainOut,filesep,Subj(i).name,anat,{' '},main,filesep,Subj(i).name,anat,T1w);
                    system(command{1});
                    fileA = dir(strcat(mainOut,filesep,Subj(i).name,anat,filesep,'*Crop*'));
                    command = strcat({'mv '},mainOut,filesep,Subj(i).name,anat,filesep,fileA.name,{' '},mainOut,filesep,Subj(i).name,anat,filesep,Subj(i).name,'_acq-crop_T1w.nii');
                    system(command{1});
                    clear fileA
                end;
            end;
        else % Si hay una o mas sesiones dentro de cada sujeto
            for j = 1:length(Ses)
                if not(isdir (fullfile(mainOut,Subj(i).name,Ses(j).name,'anat'))) && isdir(fullfile(main,Subj(i).name,Ses(j).name,'anat','T1w'))
                    mkdir(fullfile(mainOut,Subj(i).name,Ses(j).name,'anat'));
                    command = strcat(DCM2NII,{' -f "'},Subj(i).name,{'_'},Ses(j).name,{'_%f" -b n -z n -x y -o '},fullfile(mainOut,Subj(i).name,Ses(j).name,'anat'),{' '},fullfile(main,Subj(i).name,Ses(j).name,'anat','T1w'));
                    system(command{1});
                    fileA = dir(fullfile(mainOut,Subj(i).name,Ses(j).name,'anat','*Crop*'));
                    if length(fileA)>1
                        disp('warning:',fileA)
                        for k = 1:length(fileA)
                            movefile(fullfile(mainOut,Subj(i).name,Ses(j).name,'anat',fileA(k).name),strcat(mainOut,filesep,Subj(i).name,filesep,Ses(j).name,anat,filesep,Subj(i).name,'_',Ses(j).name,'_acq-crop_T1w.nii'))
                        end
                    
                    else
                        movefile(fullfile(mainOut,Subj(i).name,Ses(j).name,'anat',fileA.name),strcat(mainOut,filesep,Subj(i).name,filesep,Ses(j).name,anat,filesep,Subj(i).name,'_',Ses(j).name,'_acq-crop_T1w.nii'))
                    end    
                        %                     command = strcat({'mv '},mainOut,filesep,Subj(i).name,filesep,Ses(j).name,anat,filesep,fileA.name,{' '},mainOut,filesep,Subj(i).name,filesep,Ses(j).name,anat,filesep,Subj(i).name,'_',Ses(j).name,'_acq-crop_T1w.nii');
%                     system(command{1});
                    clear fileA
                end;
            end;
        end;
        
% %% ------------- FUNC -----------------------------------------------------
        if Ses_flag %Si no hay sesiones dentro de cada sujeto

                if isdir([main,filesep,Subj(i).name,func]) && not(isempty(dir(fullfile(main,Subj(i).name,'func','task-*_bold'))))
                    if not(isdir ([mainOut filesep Subj(i).name func]))
                        mkdir([mainOut filesep Subj(i).name func]);
                    end;
                    func_acq=dir(strcat(main,filesep,Subj(i).name,func,filesep,'task*bold'));
                    for l = 1:length(func_acq)
                        command = strcat(DCM2NII,{' -f "'},Subj(i).name,{'_%f" -z n -o '},mainOut,filesep,Subj(i).name,func,{' '},main,filesep,Subj(i).name,func,filesep,func_acq(l).name);
                        system(command{1}); 
                    end;
                end;
        else
            for j = 1:length(Ses)
                if not(isdir(fullfile(mainOut,Subj(i).name,Ses(j).name,'func'))) && not(isempty(dir(fullfile(main,Subj(i).name, Ses(j).name,'func','task-*_bold'))))  
                    mkdir(fullfile(mainOut,Subj(i).name, Ses(j).name,'func'));
                    func_acq=dir(fullfile(main,Subj(i).name,Ses(j).name,'func','task*bold'));
                    for l = 1:length(func_acq)
                        command = strcat(DCM2NII,{' -f "'},Subj(i).name,{'_'},Ses(j).name,{'_%f" -z n -o '},fullfile(mainOut,Subj(i).name,Ses(j).name,'func'),{' '},fullfile(main,Subj(i).name,Ses(j).name,'func',func_acq(l).name));
                        system(command{1}); 
                    end;
                end;
            end;
        end;

% %% ------------- DWI ------------------------------------------------------
        if Ses_flag %Si no hay sesiones dentro de cada sujeto

                if isdir([main,filesep,Subj(i).name,dwi]) && not(isempty(dir(fullfile(main,Subj(i).name,'dwi','dir-AP*'))))
                    if not(isdir ([mainOut filesep Subj(i).name dwi]))
                        mkdir([mainOut filesep Subj(i).name dwi]);
                    end;
                    dwi_acq=dir(strcat(main,filesep,Subj(i).name,dwi,filesep,'dir-AP*'));
                    for l = 1:length(dwi_acq)
                        command = strcat(DCM2NII,{' -f "'},Subj(i).name,{'_%f" -z n -o '},mainOut,filesep,Subj(i).name,dwi,{' '},main,filesep,Subj(i).name,dwi,filesep,dwi_acq(l).name);
                        system(command{1}); 
                    end;
                end;
        else
            for j = 1:length(Ses) %% checked by Don Laso and working
                if not(isdir(fullfile(mainOut,Subj(i).name,Ses(j).name,'dwi'))) && not(isempty(dir(fullfile(main,Subj(i).name, Ses(j).name,'dwi','dir-AP*'))))  
                    % output folder does not exit && input folder is not
                    % create output folder to store secrets
                    mkdir(fullfile(mainOut,Subj(i).name, Ses(j).name,'dwi'));
                    % take all files that fulfill...:
                    dwi_acq=dir(fullfile(main,Subj(i).name,Ses(j).name,'dwi','dir-AP*'));
                    for l = 1:length(dwi_acq) % and convert them to niifty uajajaja
                        command = strcat(DCM2NII,{' -f "'},Subj(i).name,{'_'},Ses(j).name,{'_%f" -z n -o '},fullfile(mainOut,Subj(i).name,Ses(j).name,'dwi'),{' '},fullfile(main,Subj(i).name,Ses(j).name,'dwi',dwi_acq(l).name));
                        system(command{1}); 
                    end;
                end;
            end;
        end;
        %%%%%
%         if isdir ([main,filesep,Subj(i).name,dwi])  
%             if not(isdir ([mainOut filesep Subj(i).name dwi]))
%                 mkdir([mainOut filesep Subj(i).name dwi]);
%             end;
%             dwi_acq=dir(strcat(main,filesep,Subj(i).name,dwi,filesep,'dir*'));
%             for j = 1:length(dwi_acq)
%                 command2=strcat(DCM2NII,{' -f "'},Subj(i).name,{'_%f" -b y -z n -o '},mainOut,filesep,Subj(i).name,dwi,{' '},main,filesep,Subj(i).name,dwi,filesep,dwi_acq(j).name);
%                 system(command2{1});
%             end;                  
%         end;
%             
% %% ------------ FMAP ------------------------------------------------------           

        if Ses_flag %Si no hay sesiones dentro de cada sujeto

                if isdir([main,filesep,Subj(i).name,fieldmap]) && not(isempty(dir(fullfile(main,Subj(i).name,'fmap','acq-*'))))
                    if not(isdir ([mainOut filesep Subj(i).name fieldmap]))
                        mkdir([mainOut filesep Subj(i).name fieldmap]);
                    end;
                    fmap_acq=dir(strcat(main,filesep,Subj(i).name,fieldmap,filesep,'acq-*'));
                    for l = 1:length(fmap_acq)
                        command = strcat(DCM2NII,{' -f "'},Subj(i).name,{'_%f" -z n -o '},mainOut,filesep,Subj(i).name,fieldmap,{' '},main,filesep,Subj(i).name,fieldmap,filesep,fmap_acq(l).name);
                        system(command{1}); 
                    end;
                end;
        else
            for j = 1:length(Ses) %% checked by Don Laso and working
                if not(isdir(fullfile(mainOut,Subj(i).name,Ses(j).name,'fmap'))) && not(isempty(dir(fullfile(main,Subj(i).name, Ses(j).name,'fmap','acq-*'))))  
                    % output folder does not exit && input folder is not
                    % create output folder to store secrets
                    mkdir(fullfile(mainOut,Subj(i).name, Ses(j).name,'fmap'));
                    % take all files that fulfill...:
                    fmap_acq=dir(fullfile(main,Subj(i).name,Ses(j).name,'fmap','acq-*'));
                    for l = 1:length(fmap_acq) % and convert them to niifty uajajaja
                        command = strcat(DCM2NII,{' -f "'},Subj(i).name,{'_'},Ses(j).name,{'_%f" -z n -o '},fullfile(mainOut,Subj(i).name,Ses(j).name,'fmap'),{' '},fullfile(main,Subj(i).name,Ses(j).name,'fmap',fmap_acq(l).name));
                        system(command{1}); 
                    end;
                end;
            end;
        end;
%         if isdir ([main,filesep,Subj(i).name,fieldmap]) 
%             if not(isdir ([mainOut filesep Subj(i).name fieldmap]))
%                 mkdir([mainOut filesep Subj(i).name fieldmap]);
%             end;
%             fmap_acq = dir(strcat(main,filesep,Subj(i).name,fieldmap,filesep,'acq-rest*'));
%             for j = 1:length(fmap_acq)
%                 command2=strcat(DCM2NII,{' -f "'},Subj(i).name,{'_%f" -b y -z n -o '},mainOut,filesep,Subj(i).name,fieldmap,{' '},main,filesep,Subj(i).name,fieldmap,filesep,fmap_acq(j).name);
%                 system(command2{1});
%             end;
%         end;
% %% ------------- ASL ------------------------------------------------------
%         if isdir ([main,filesep,Subj(i).name,func,asl])
%             if not(isdir ([mainOut filesep Subj(i).name func]))
%                 mkdir([mainOut filesep Subj(i).name func]);
%             end;
%             command2 = strcat(DCM2NII,{' -f "'},Subj(i).name,{'_%f" -b y -z n -o '},mainOut,filesep,Subj(i).name,func,{' '},main,filesep,Subj(i).name,func,asl);
%             system(command2{1});
%         end;
%         
% %% ------------- relCBF ---------------------------------------------------
%         if isdir ([main,filesep,Subj(i).name,func,relCBF])
%             if not(isdir ([mainOut filesep Subj(i).name func]))
%                 mkdir([mainOut filesep Subj(i).name func]);
%             end;
%             command2 = strcat(DCM2NII,{' -f "'},Subj(i).name,{'_%f" -b y -z n -o '},mainOut,filesep,Subj(i).name,func,{' '},main,filesep,Subj(i).name,func,relCBF);
%             system(command2{1});
%         end;
        

end;
%--------------------------------------------------------------------------


fprintf('\nDONE!\n')
