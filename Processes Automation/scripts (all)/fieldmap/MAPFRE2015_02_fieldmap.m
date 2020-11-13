% Date : 04 - Nov - 2016
% Author : Mario Gil 
%
%PROYECTO MAPFRE 2015
% Preprocessing of Olfaction fMRI data
%
% Step 1. Calculation of the fieldma.p
%
%   This scrip extracts the short and long the Echo Time for the fieldmap
%   series, the Effective Eco Spacing (0043,102c) from the funcional
%   series,the Echo Number (0018,0086) from the funtional and the asset
%   (1043,1083) from the functional. 
% 
%   It calculates the Total Readout Time according to the following formula. 
%   ReadOut = EchoNumber * EchoSpacing * Asset
%
%   The resulting values will be saved in a txt file where you can find the
%   subject identifier and the corresponding extracted values.
%  
%% Ibicializar el entorno 
clc
clear all
close all
%% Variables declaration:

anat     = strcat(filesep,'ANATOMICO',filesep,'SAG_3D');
ge_epi   = strcat(filesep,'GE_EPI');
fieldmap = strcat(filesep,'GE_EPI_fieldmap');
rest     = strcat(filesep,'RESTING');
olfat    = strcat(filesep,'OLFATO');
outCont  = strcat(filesep,'controles');
outPre   = strcat(filesep,'pacientes_pre');
outPost  = strcat(filesep,'pacientes_post');

% Main Folder path:
%   This variable has to be changed for every compueter where it is executed
%   or any other study.

main_raw = '/Volumes/PassportMarioGilCorrea/URJC/PRACTICAS_TFG/Datos_TFG/MAPFRE_fMRI_Sujetos_raw';
mainOut ='/Volumes/PassportMarioGilCorrea/URJC/PRACTICAS_TFG/Datos_TFG/MAPFRE_fMRI_Sujetos';

all_Subj=dir(fullfile(mainOut,'id*'));

%%

for i = 1:length(all_Subj)
   
 
   fieldmaps = dir(strcat(main_raw,filesep,all_Subj(i).name,fieldmap,filesep,'*.dcm'));
   funcionals = dir(strcat(main_raw,filesep,all_Subj(i).name,ge_epi,olfat,filesep,'*.dcm'));

   fieldmapIfo.short = dicominfo(strcat(main_raw,filesep,all_Subj(i).name,fieldmap,filesep,fieldmaps(1).name)); 
   fieldmapIfo.long = dicominfo(strcat(main_raw,filesep,all_Subj(i).name,fieldmap,filesep,fieldmaps(length(fieldmaps)).name)); 
   functionalInfo = dicominfo (strcat(main_raw,filesep,all_Subj(i).name,ge_epi,olfat,filesep,funcionals(1).name)); 
   
   TE1 = double(fieldmapIfo.short.EchoTime);
   TE2 = double(fieldmapIfo.long.EchoTime);
   EchoSpacing = double(functionalInfo.Private_0043_102c);
   EchoNumber = double(functionalInfo.AcquisitionMatrix(1));
   asset = double(functionalInfo.Private_0043_1083(1));
   ReadOut = EchoSpacing* 10^(-3) * EchoNumber * asset ;
   
   fid=fopen(strcat(mainOut,filesep,all_Subj(i).name,fieldmap,filesep,all_Subj(i).name,'_FieldmapCalc.txt'),'w');

   fprintf(fid,'%f %f %f %f',TE1,TE2,EchoSpacing,ReadOut); 
   
   fclose(fid);
   clear fid
   
   fieldmapIfo = [];
   functionalInfo = []; 
   fieldmaps = [];
   funcionals = [];
   TE1 = [];
   TE2 = [];
   EchoSpacing = []; 
   EchoNumber = [];
   asset = []; 
   ReadOut = []; 

end;

%%
%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%----------------------------------------------------------------------


%% Inicializacion del entorno de spm

spm('defaults','fmri');
spm_jobman('initcfg')


for i=1:length(all_Subj)
        
        paciente_actual = fullfile(mainOut,all_Subj(i).name);
        paciente_ge_epi_butanol = [paciente_actual, '/GE_EPI/OLFATO'];
        paciente_ge_fieldmap = [paciente_actual, '/GE_EPI_fieldmap/'];
        
        shortreal=dir([paciente_ge_fieldmap, '/TE1/TE1_r/*.nii']);
        shortimag=dir([paciente_ge_fieldmap, '/TE1/TE1_i/*.nii']);
        longreal=dir([paciente_ge_fieldmap, '/TE2/TE2_r/*.nii']);
        longimag=dir([paciente_ge_fieldmap, '/TE2/TE2_i/*.nii']);
        
        FieldmapFile = dir(fullfile(paciente_ge_fieldmap,'*.txt'));
        fid=fopen(strcat(paciente_ge_fieldmap,filesep,FieldmapFile.name),'r');
        FieldmapData=fscanf(fid,'%f');
        
        fprintf('Preproceso funcional: %s\n',paciente_actual)
        
        ruta_dir_func_nii=strcat(paciente_ge_epi_butanol,filesep,'*.nii');
        list_func_dir=dir(ruta_dir_func_nii);
        data_r_u=cell(196,1);
        
        %comprobamos que no haya mas .nii de otras veces
        for n=1:length(list_func_dir)
            temp_func=list_func_dir(n).name;
            
            %Estamos suponiendo que los ficheros empiezan por I
            if temp_func(1)=='I'
                list_temp=temp_func;
            end
        end
        for k=1:196
            data_r_u{k}=strcat(paciente_ge_epi_butanol,filesep,temp_func,',',int2str(k));
        end

matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.realimag.shortreal = {strcat(paciente_ge_fieldmap,'/TE1/TE1_r/',shortreal.name,',1')};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.realimag.shortimag = {strcat(paciente_ge_fieldmap,'/TE1/TE1_i/',shortimag.name,',1')};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.realimag.longreal = {strcat(paciente_ge_fieldmap,'/TE2/TE2_r/',longreal.name,',1')};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.realimag.longimag = {strcat(paciente_ge_fieldmap,'/TE2/TE2_i/',longimag.name,',1')};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [FieldmapData(1) FieldmapData(2)];
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 1;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = -1;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = FieldmapData(4);
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {'/Volumes/PassportMarioGilCorrea/Programs/SPM/spm12/toolbox/FieldMap/T1.nii'};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session.epi = {'/Volumes/PassportMarioGilCorrea/Programs/SPM/spm12/toolbox/OldNorm/EPI.nii,1'};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 0;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat = '';
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;
%%
matlabbatch{2}.spm.tools.fieldmap.applyvdm.data.scans = data_r_u;
%%
matlabbatch{2}.spm.tools.fieldmap.applyvdm.data.vdmfile(1) = cfg_dep('Calculate VDM: Voxel displacement map (Subj 1, Session 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','vdmfile', '{}',{1}));
matlabbatch{2}.spm.tools.fieldmap.applyvdm.roptions.pedir = 2;
matlabbatch{2}.spm.tools.fieldmap.applyvdm.roptions.which = [2 1];
matlabbatch{2}.spm.tools.fieldmap.applyvdm.roptions.rinterp = 4;
matlabbatch{2}.spm.tools.fieldmap.applyvdm.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.tools.fieldmap.applyvdm.roptions.mask = 1;
matlabbatch{2}.spm.tools.fieldmap.applyvdm.roptions.prefix = 'B0';
        
output_list_3 =spm_jobman('run',matlabbatch);
        
clear matlabbatch
        
    %end
    
end

fprintf('/nFieldmap calculated and applied for all the founded subjects./n')

