% Date : 04 - Nov - 2016
% Author : Mario Gil 
%
%
%% Ibicializar el entorno 
clc
clear all
close all
%% Variables declaration:

% Main Folder path:
%   This variable has to be changed for every compueter where it is executed
%   or any other study.

BIDS_dir = '/Users/administrator/Documents/SENECA-PICASSO/SP_nifti';
Subj=dir(fullfile(BIDS_dir,'sub-sp01*'));
fprintf('\nStarting fieldmap reconstruction with a pair of real and imaginary images.\nSubjects founded %d',length(Subj))
spms(12);
spm('defaults','fmri');
spm_jobman('initcfg');
%% Main
for i = 1:length(Subj)
   Ses = dir(fullfile(BIDS_dir,Subj(i).name,'ses-01*'));
   
   for j = 1:length(Ses)
       
       fprintf('\nFieldmap preprocessing: %s_%s\n', Subj(i).name, Ses(j).name)
       
       fmap_echo1 = dir(fullfile(BIDS_dir,Subj(i).name,Ses(j).name,'fmap','sub-*_real_e1.json'));
       fmap_echo2 = dir(fullfile(BIDS_dir,Subj(i).name,Ses(j).name,'fmap','sub-*_real_e2.json'));
       bold = dir(fullfile(BIDS_dir,Subj(i).name,Ses(j).name,'func','sub-*_bold.json'));
       
       fmap_echo1_json = loadjson(fullfile(fmap_echo1.folder,fmap_echo1.name));
       fmap_echo2_json = loadjson(fullfile(fmap_echo2.folder,fmap_echo2.name));
       bold_json = loadjson(fullfile(bold.folder,bold.name));
       
       TE1 = fmap_echo1_json.EchoTime * 1000;
       TE2 = fmap_echo2_json.EchoTime * 1000;
%        EPI_ReadOut = bold_json.TotalReadoutTime * 1000;
       EPI_ReadOut = bold_json.AcquisitionMatrixPE*bold_json.EffectiveEchoSpacing*1000*0.5; % hola
       
       clear fmap_echo1_json fmap_echo2_json bold_json
      
        
       shortreal = dir( fullfile( BIDS_dir, Subj(i).name, Ses(j).name, 'fmap', 'sub-*_real_e1.nii'));
       shortimag = dir( fullfile( BIDS_dir, Subj(i).name, Ses(j).name, 'fmap', 'sub-*_imag_e1.nii'));
       longreal = dir( fullfile( BIDS_dir, Subj(i).name, Ses(j).name, 'fmap', 'sub-*_real_e2.nii'));
       longimag = dir( fullfile( BIDS_dir, Subj(i).name, Ses(j).name, 'fmap', 'sub-*_imag_e2.nii'));
       bold_nii = dir( fullfile( BIDS_dir, Subj(i).name, Ses(j).name, 'func', 'sub-*_bold.nii'));
       
       bold_data=cell(120,1);
       for k = 1:120
           bold_data{k}=strcat(bold_nii.folder,filesep,bold_nii.name,',',int2str(k));
       end

matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.realimag.shortreal = {strcat(shortreal.folder,filesep,shortreal.name,',1')};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.realimag.shortimag = {strcat(shortimag.folder,filesep,shortimag.name,',1')};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.realimag.longreal = {strcat(longreal.folder,filesep,longreal.name,',1')};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.realimag.longimag = {strcat(longimag.folder,filesep,longimag.name,',1')};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [TE1 TE2];
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 1;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = -1;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = EPI_ReadOut;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {'/Users/administrator/Documents/MATLAB/spm12/toolbox/FieldMap/T1.nii'};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session.epi = {strcat(bold_nii.folder,filesep,bold_nii.name,',1')};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 0;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat = '';
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;
%%
matlabbatch{2}.spm.tools.fieldmap.applyvdm.data.scans = bold_data;
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
       
       
       
   end
end;



fprintf('\nFieldmap calculated and applied for all the founded subjects.\n')

