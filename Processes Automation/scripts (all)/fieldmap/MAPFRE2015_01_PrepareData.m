%% Date : 18 - Oct - 2016
% Author : Mario Gil 
%
% PROYECTO MAPFRE 2015
% Preparing of Olfaction fMRI data
% In this scrip the main folder names will be treated as variables.
%
%% The inputs needed are:
%   - Origin directory of subjects with their corresponding dicoms.
%       Note: The directory must have the following structure.
%       >FolderX
%           >>idxxxx
%               >>>ANATOMICO
%                   >>>>SAG_3D
%                       dicom images
%
%               >>>GE_EPI
%                   >>>>OLFATO
%                       dicom images
%                   >>>>RESTING
%                       dicom images
%
%               >>>GE_EPI_fieldmal
%                   dicom images
%
%   - Output directory
%   - Folder where the dcm2nii method is stored.
%   - Folder where the STMP text files are stored.
%
%% Steps of the script
% Step 1: Dicom to NifTi conversion
% Step 2: STMP reading and EV files writing
%% Initialization
clear all
clc
%% Main Folder paths:

%   This variables has to be changed for every compueter where it is executed
%   or any other study.
main = '/Volumes/LaCie/ProyectoMAPFRE2015/ProyectoMAPFRE_DICOM';
mainOut = '/Volumes/LaCie/ProyectoMAPFRE2015/ProyectoMAPFRE_Sujetos';
mainDCM2NII= '/Users/frs16g/Downloads/mricronmac';
mainSTMP = '/Volumes/LaCie/ProyectoMAPFRE2015/ProyectoMAPFRE_Olfatometro/STMP_files';

%% Step 1. Dicom to NifTi conversion
%
%   This scipt converts all dicom files of every subject into nifti files
%   without changin the archiving structure. Which is the following: 
%       idxxxx/ANATOMICO/SAG_3D/*.dcm
%       idxxxx/GE_EPI/OLFATO/*.dcm
%       idxxxx/GE_EPI/RESTING/*.dcm
%       idxxxx/GE_EPI_fieldmap/*.dcm
%   The resulting NifTi files will be also separated into different grops
%   which are:
%       -Controls
%       -Pre
%       -Post
%-------------------------------------------------------------------------- 
% Variables declaration:

anat     = strcat(filesep,'ANATOMICO',filesep,'SAG_3D');
ge_epi   = strcat(filesep,'GE_EPI');
fieldmap = strcat(filesep,'GE_EPI_fieldmap');
rest     = strcat(filesep,'RESTING');
olfat    = strcat(filesep,'OLFATO');
outCont  = strcat(filesep,'controles');
outPre   = strcat(filesep,'pacientes_pre');
outPost  = strcat(filesep,'pacientes_post');

DCM2NII = strcat (mainDCM2NII,filesep,'dcm2nii');
PreSubj=dir(strcat(main,filesep,'id0*'));
PostSubj=dir(strcat(main,filesep,'id1*'));
ContSubj=dir(strcat(main,filesep,'id2*'));

%--------------------------------------------------------------------------
% Dicom to NifTi conversion in Pre Subjects

 command = {strcat('ln -s ',main,' /Users/frs16g/Desktop/test')};
 system(command{1});

for i = 1:length(PreSubj)
    
   if not(isdir ([mainOut outPre filesep PreSubj(i).name anat]))
       mkdir([mainOut outPre filesep PreSubj(i).name anat]);
   end;
   command = strcat({'dcm2nii -d N -e N -f N -g N -p N -r N -x Y '},main,filesep,PreSubj(i).name,anat);
   system(command{1});
   fileA = dir(strcat(main,filesep,PreSubj(i).name,anat,filesep,'*.nii'));
   command = strcat({'mv '},main,filesep,PreSubj(i).name,anat,filesep,fileA.name,{' '},mainOut,outPre,filesep,PreSubj(i).name,anat,filesep,PreSubj(i).name,'_anat.nii');
   system(command{1});
   
   
   if not(isdir ([mainOut outPre filesep PreSubj(i).name ge_epi rest]))
       mkdir([mainOut outPre filesep PreSubj(i).name ge_epi rest]);
   end;
   command = strcat(DCM2NII,{' -d N -e N -f N -g N -p N -r N '},main,filesep,PreSubj(i).name,ge_epi,rest);
   system(command{1});
   fileB = dir(strcat(main,filesep,PreSubj(i).name,ge_epi,rest,filesep,'*.nii'));
   command = strcat({'mv '},main,filesep,PreSubj(i).name,ge_epi,rest,filesep,fileB.name,{' '},mainOut,outPre,filesep,PreSubj(i).name,ge_epi,rest,filesep,PreSubj(i).name,'_rest.nii');
   system(command{1});
   
   
   if not(isdir ([mainOut outPre filesep PreSubj(i).name ge_epi olfat]))
       mkdir([mainOut outPre filesep PreSubj(i).name ge_epi olfat]);
   end;
   command = strcat(DCM2NII,{' -d N -e N -f N -g N -p N -r N '},main,filesep,PreSubj(i).name,ge_epi,olfat);
   system(command{1});
   fileC = dir(strcat(main,filesep,PreSubj(i).name,ge_epi,olfat,filesep,'*.nii'));
   command = strcat({'mv '},main,filesep,PreSubj(i).name,ge_epi,olfat,filesep,fileC.name,{' '},mainOut,outPre,filesep,PreSubj(i).name,ge_epi,olfat,filesep,PreSubj(i).name,'_olfat.nii');
   system(command{1});
   
     
   if not(isdir ([mainOut outPre filesep PreSubj(i).name fieldmap]))
       mkdir([mainOut outPre filesep PreSubj(i).name fieldmap]);
   end;

   fieldmap2nii(strcat(main,filesep,PreSubj(i).name,fieldmap),mainDCM2NII)
   movefile(strcat(main,filesep,PreSubj(i).name,fieldmap,filesep,'TE*'),strcat(mainOut,outPre,filesep,PreSubj(i).name,fieldmap))

%    fileD = dir(strcat(main,filesep,PreSubj(i).name,fieldmap,filesep,'TE*'));
%    command = strcat({'cp -R '},main,filesep,PreSubj(i).name,fieldmap,filesep,fileD(1).name,{' '},mainOut,outPre,filesep,PreSubj(i).name,fieldmap,filesep);
%    system(command{1});
%    command = strcat({'cp -R '},main,filesep,PreSubj(i).name,fieldmap,filesep,fileD(2).name,{' '},mainOut,outPre,filesep,PreSubj(i).name,fieldmap,filesep);
%    system(command{1});  
end;
%--------------------------------------------------------------------------
% Dicom to NifTi conversion in Post Subjects

for i = 1:length(PostSubj)
    
    if not(isdir ([mainOut outPost '/' PostSubj(i).name anat]))
        mkdir([mainOut outPost '/' PostSubj(i).name anat]);
    end;
    command = strcat(DCM2NII,{' -d N -e N -f N -g N -p N -r N '},main,filesep,PostSubj(i).name,anat);
    system(command{1});
    fileE = dir(strcat(main,filesep,PostSubj(i).name,anat,filesep,'*.nii'));
    command = strcat({'mv '},main,filesep,PostSubj(i).name,anat,filesep,fileE.name,{' '},mainOut,outPost,filesep,PostSubj(i).name,anat,filesep,PostSubj(i).name,'_anat.nii');
    system(command{1});
    
    
    if not(isdir ([mainOut outPost filesep PostSubj(i).name ge_epi rest]))
        mkdir([mainOut outPost filesep PostSubj(i).name ge_epi rest]);
    end;
    command = strcat(DCM2NII,{' -d N -e N -f N -g N -p N -r N '},main,filesep,PostSubj(i).name,ge_epi,rest);
    system(command{1});
    fileF = dir(strcat(main,filesep,PostSubj(i).name,ge_epi,rest,filesep,'*.nii'));
    command = strcat({'mv '},main,filesep,PostSubj(i).name,ge_epi,rest,filesep,fileF.name,{' '},mainOut,outPost,filesep,PostSubj(i).name,ge_epi,rest,filesep,PostSubj(i).name,'_rest.nii');
    system(command{1});
    
    
    if not(isdir ([mainOut outPost filesep PostSubj(i).name ge_epi olfat]))
        mkdir([mainOut outPost filesep PostSubj(i).name ge_epi olfat]);
    end;
    command = strcat(DCM2NII,{' -d N -e N -f N -g N -p N -r N '},main,filesep,PostSubj(i).name,ge_epi,olfat);
    system(command{1});
    fileG = dir(strcat(main,filesep,PostSubj(i).name,ge_epi,olfat,filesep,'*.nii'));
    command = strcat({'mv '},main,filesep,PostSubj(i).name,ge_epi,olfat,filesep,fileG.name,{' '},mainOut,outPost,filesep,PostSubj(i).name,ge_epi,olfat,filesep,PostSubj(i).name,'_olfat.nii');
    system(command{1});
    
    
    if not(isdir ([mainOut outPost filesep PostSubj(i).name fieldmap]))
        mkdir([mainOut outPost filesep PostSubj(i).name fieldmap]);
    end;
    
    fieldmap2nii(strcat(main,filesep,PostSubj(i).name,fieldmap),mainDCM2NII)  
    movefile(strcat(main,filesep,PostSubj(i).name,fieldmap,filesep,'TE*'),strcat(mainOut,outPost,filesep,PostSubj(i).name,fieldmap))

%     fileH = dir(strcat(main,filesep,PostSubj(i).name,fieldmap,filesep,'TE*'));
%     command = strcat({'cp -R '},main,filesep,PostSubj(i).name,fieldmap,filesep,fileH(1).name,{' '},mainOut,outPost,filesep,PostSubj(i).name,fieldmap,filesep);
%     system(command{1});
%     command = strcat({'cp -R '},main,filesep,PostSubj(i).name,fieldmap,filesep,fileH(2).name,{' '},mainOut,outPost,filesep,PostSubj(i).name,fieldmap,filesep);
%     system(command{1});   
end;

%--------------------------------------------------------------------------
% Dicom to NifTi conversion in Control Subjects

for i = 1:length(ContSubj) 
    %i=length(ContSubj) ;
    if not(isdir ([mainOut outCont filesep ContSubj(i).name anat]))
        mkdir([mainOut outCont filesep ContSubj(i).name anat]);
    end;
    command = strcat(DCM2NII,{' -d N -e N -f N -g N -p N -r N '},main,filesep,ContSubj(i).name,anat);
    system(command{1});
    fileI = dir(strcat(main,filesep,ContSubj(i).name,anat,filesep,'*.nii'));
    command = strcat({'mv '},main,filesep,ContSubj(i).name,anat,filesep,fileI.name,{' '},mainOut,outCont,filesep,ContSubj(i).name,anat,filesep,ContSubj(i).name,'_anat.nii');
    system(command{1});
    
    if not(isdir ([mainOut outCont '/' ContSubj(i).name ge_epi rest]))
        mkdir([mainOut outCont '/' ContSubj(i).name ge_epi rest]);
    end;
    command = strcat(DCM2NII,{' -d N -e N -f N -g N -p N -r N '},main,filesep,ContSubj(i).name,ge_epi,rest);
    system(command{1});
    fileJ = dir(strcat(main,filesep,ContSubj(i).name,ge_epi,rest,filesep,'*.nii'));
    command = strcat({'mv '},main,filesep,ContSubj(i).name,ge_epi,rest,filesep,fileJ.name,{' '},mainOut,outCont,filesep,ContSubj(i).name,ge_epi,rest,filesep,ContSubj(i).name,'_rest.nii');
    system(command{1});
    
    if not(isdir ([mainOut outCont filesep ContSubj(i).name ge_epi olfat]))
        mkdir([mainOut outCont filesep ContSubj(i).name ge_epi olfat]);
    end;
    command = strcat(DCM2NII,{' -d N -e N -f N -g N -p N -r N '},main,filesep,ContSubj(i).name,ge_epi,olfat);
    system(command{1});
    fileK = dir(strcat(main,filesep,ContSubj(i).name,ge_epi,olfat,filesep,'*.nii'));
    command = strcat({'mv '},main,filesep,ContSubj(i).name,ge_epi,olfat,filesep,fileK.name,{' '},mainOut,outCont,filesep,ContSubj(i).name,ge_epi,olfat,filesep,ContSubj(i).name,'_olfat.nii');
    system(command{1});
   
    
    if not(isdir ([mainOut outCont filesep ContSubj(i).name fieldmap]))
        mkdir([mainOut outCont filesep ContSubj(i).name fieldmap]);
    end;
    
    fieldmap2nii(strcat(main,filesep,ContSubj(i).name,fieldmap),mainDCM2NII)
    movefile(strcat(main,filesep,ContSubj(i).name,fieldmap,filesep,'TE*'),strcat(mainOut,outCont,filesep,ContSubj(i).name,fieldmap))
  
%     fileL = dir(strcat(main,filesep,ContSubj(i).name,fieldmap,filesep,'TE*'));
%     command = strcat({'cp -R '},main,filesep,ContSubj(i).name,fieldmap,filesep,fileL(1).name,{' '},mainOut,outCont,filesep,ContSubj(i).name,fieldmap,filesep);
%     system(command{1});
%     command = strcat({'cp -R '},main,filesep,ContSubj(i).name,fieldmap,filesep,fileL(2).name,{' '},mainOut,outCont,filesep,ContSubj(i).name,fieldmap,filesep);
%     system(command{1});
end;
%% STEP 2. STMP files reading and EV files Writing and relocation
%Lectura de STMP y Escritura de EVs
%experiment sera la forma del pulso del estimulo
%maximo_testimulo es el tiempo maximo de estimulos, para el hpf
%archivo.txt
%TR en segundos
files=dir(strcat(mainSTMP,filesep,'*.txt'));
for k = 1:length(files)
    archivo = strcat(mainSTMP,filesep,files(k).name);
    fid = fopen(archivo, 'r');
    if fid==-1 
        fprintf('Error archivo') 
    end
    F = fread(fid,'*char')';
    F = regexprep(F,'(?:?)','a');
    fclose(fid);
    C2 = strsplit(F,'vula 2, Instante: ');
    tiempo_apertura_v2=zeros(1,(size(C2,2)-1)/2);

    C6 = strsplit(F,'vula 6, Instante: ');
    tiempo_apertura_v6=zeros(1,(size(C6,2)-1)/2);

    contador=0;
    for i=2:2:size(C2,2)
    
        contador=contador+1;
        cadena = strsplit(C2{i},',');
    
        tiempo_apertura_v2(contador) = str2double(cadena{1});
 %   tiempo_apertura_v2(contador)=tiempo_apertura_v2(contador)-Tdummy;
    end

    contador=0;
    for i=2:2:size(C6,2)
    
        contador=contador+1;
        cadena = strsplit(C6{i},',');
    
        tiempo_apertura_v6(contador) = str2double(cadena{1});
  %  tiempo_apertura_v6(contador)=tiempo_apertura_v6(contador)-Tdummy;
    end

    %Crear el fichero EV

    amplitud=1;
    duracion=12; 

    %Tiempo_inicio Duracion Amplitud
%     name_fsl = strsplit(archivo,filesep);
%     name_fsl=name_fsl{length(name_fsl)};
    name_fsl=strsplit(archivo,'.txt');
    
    name_archivo=[name_fsl{1},'_MD_fsl_v2.txt'];
    abs_path = strsplit(name_archivo,strcat(filesep,'STMP_files'));
    if not(isdir(strcat(abs_path{1},filesep,'EV_files')))
        mkdir(abs_path{1},'EV_files');
    end;

    name_archivo = strcat(abs_path{1},filesep,'EV_files',abs_path{2});
    fid = fopen(name_archivo, 'w');
    for i=1:length(tiempo_apertura_v2)
        %cadena=[num2str(tiempo_apertura(i)),'  ',num2str(duracion), '  ', num2str(amplitud), '\n'];
        %F = fwrite(fid,cadena,'*char');
        fprintf(fid,'%d %d %d\n',tiempo_apertura_v2(i),duracion,amplitud);
    end
    fclose(fid);

    name_archivo=[name_fsl{1},'_MD_fsl_v6.txt'];
    abs_path = strsplit(name_archivo,strcat(filesep,'STMP_files'));
    name_archivo = strcat(abs_path{1},filesep,'EV_files',abs_path{2});
    fid = fopen(name_archivo, 'w');
    for i=1:length(tiempo_apertura_v6)
        %cadena=[num2str(tiempo_apertura(i)),'  ',num2str(duracion), '  ', num2str(amplitud), '\n'];
        %F = fwrite(fid,cadena,'*char');
        fprintf(fid,'%d %d %d\n',tiempo_apertura_v6(i),duracion,amplitud);
    end

    fclose(fid);

    clear archivo
end;

%--------------------------------------------------------------------------
% Relocation of EV files
    
mainEVs = strcat(abs_path{1},filesep,'EV_files');
    
PreSubjEV=dir(strcat(mainEVs,filesep,'MAPFRE_0*'));
PostSubjEV=dir(strcat(mainEVs,filesep,'MAPFRE_POST*'));
ContSubjEV=dir(strcat(mainEVs,filesep,'MAPFRE_CONTROL*'));

% PreSubj=dir(strcat(mainOut,filesep,'id0*'));
% PostSubj=dir(strcat(mainOut,filesep,'id1*'));
% ContSubj=dir(strcat(mainOut,filesep,'id2*'));

for i = 1:2:length(PreSubjEV)-1
    
   if not(isdir (strcat(mainOut,outPre,filesep,PreSubj(round(i/2)).name,ge_epi,olfat,filesep,'onsets')))
       mkdir(strcat(mainOut,outPre,filesep,PreSubj(round(i/2)).name,ge_epi,olfat,filesep,'onsets'));
   end;
   file = dir(strcat(mainEVs,filesep,PreSubjEV(i).name));
   command = strcat({'cp '},mainEVs,filesep,PreSubjEV(i).name,{' '},mainOut,outPre,filesep,PreSubj(round(i/2)).name,ge_epi,olfat,filesep,'onsets');
   system(command{1});
   command = strcat({'cp '},mainEVs,filesep,PreSubjEV(i+1).name,{' '},mainOut,outPre,filesep,PreSubj(round(i/2)).name,ge_epi,olfat,filesep,'onsets');
   system(command{1});
   
end;

for i = 1:2:length(PostSubjEV)-1
    
   if not(isdir (strcat(mainOut,outPost,filesep,PostSubj(round(i/2)).name,ge_epi,olfat,filesep,'onsets')))
       mkdir(strcat(mainOut,outPost,filesep,PostSubj(round(i/2)).name,ge_epi,olfat,filesep,'onsets'));
   end;
   file = dir(strcat(mainEVs,filesep,PostSubjEV(i).name));
   command = strcat({'cp '},mainEVs,filesep,PostSubjEV(i).name,{' '},mainOut,outPost,filesep,PostSubj(round(i/2)).name,ge_epi,olfat,filesep,'onsets',filesep,PostSubjEV(i).name);
   system(command{1});
   command = strcat({'cp '},mainEVs,filesep,PostSubjEV(i+1).name,{' '},mainOut,outPost,filesep,PostSubj(round(i/2)).name,ge_epi,olfat,filesep,'onsets',filesep,PostSubjEV(i+1).name);
   system(command{1});
    
end;

for i = 1:2:length(ContSubjEV) - 1
        
   if not(isdir (strcat(mainOut,outCont,filesep,ContSubj(round(i/2)).name,ge_epi,olfat,filesep,'onsets')))
       mkdir(strcat(mainOut,outCont,filesep,ContSubj(round(i/2)).name,ge_epi,olfat,filesep,'onsets'));
   end;
   file = dir(strcat(mainEVs,filesep,ContSubjEV(i).name));
   command = strcat({'cp '},mainEVs,filesep,ContSubjEV(i).name,{' '},mainOut,outCont,filesep,ContSubj(round(i/2)).name,ge_epi,olfat,filesep,'onsets',filesep,ContSubjEV(i).name);
   system(command{1});
   command = strcat({'cp '},mainEVs,filesep,ContSubjEV(i+1).name,{' '},mainOut,outCont,filesep,ContSubj(round(i/2)).name,ge_epi,olfat,filesep,'onsets',filesep,ContSubjEV(i+1).name);
   system(command{1});
end;
