clear all
close all
clc

SubjDir='/Volumes/LaCie/2015_ProyectoMAPFRE/ProyectoMAPFRE_Resultados/Estudio_Anatomico/imagenes';
GroupPre=dir([SubjDir filesep 'id0*.nii']);
GroupPost=dir([SubjDir filesep 'id1*.nii']);
tpmfile='/Users/frs16g/Documents/spm12/tpm/TPM.nii';
templatePrefix='TemplateMAPFRE';
prefiles = [GroupPre;GroupPost];
files={};
for i=1:length(prefiles)
    path=strcat(SubjDir,filesep,prefiles(i).name);
    files=[files;path];
end
segmented = NewSegment12(files, tpmfile);

gm = {};
wm = {};
csf = {};
for i=1:length(segmented)
    gm{i} = segmented{i}{1};
    wm{i} = segmented{i}{2};
    csf{i} = segmented{i}{3};
end;

volumenes = {gm' wm'};
templates = CreateDartelTemplate12(volumenes, [templatePrefix 'GMWM']);
deformation_fields = ExistingDartelTemplate12(volumenes, templates');
GMWM_files = Normalise12(volumenes,deformation_fields', templates{6});

volumenescsf = {csf'};
templatescsf = CreateDartelTemplate12(volumenescsf, [templatePrefix 'CSF']);
deformation_fields_csf = ExistingDartelTemplate12(volumenescsf, templatescsf);
CSF_files = Normalise12(volumenescsf,deformation_fields_csf, templatescsf{6});