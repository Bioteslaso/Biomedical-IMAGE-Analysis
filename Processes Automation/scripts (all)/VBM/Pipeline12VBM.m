%function [final_files, sgm,swm,scsf] = PipelineVBM(group1, group2,tpmfile,templatePrefix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PipelineVBM
%--------------------------------------------------------------------------
% 
% Run a pipeline to perform the VBM DARTEL. The steps are: segmentation,
% create a template with DARTEL, deformate to this template and
% normalisation.
%
%__________________________________________________________________________
% INPUTS:
%  - group1: vector with all subject(file, age, sex) of the first group.
%  - group2: vector with all subject(file, age, sex) of the second group.
%  - tpmfile: path to TPM.nii file included in SPM software.
%  - templatePrefix: a string to name templates files.
%     
% OUTPUTS:
%  - final_files: a cell array with the path to all deformated and
%  normalised subjects.
%     
%__________________________________________________________________________
% %tpmfile =
% '/Users/lauranunez/Downloads/spm8/spm8_mcr/spm8/toolbox/Seg/TPM.nii';
%
%__________________________________________________________________________
%
% Author: Laura Nunez Gonzalez
% URJC - 08 / March / 2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Group1=PrepareData12('/Volumes/LaCie/DATOS_OLFATO_HELENA_AMATOMIA/HOMBRES','/Volumes/LaCie/DATOS_OLFATO_HELENA_AMATOMIA/HOMBRES_PREPARED');
Group2=PrepareData12('/Volumes/LaCie/DATOS_OLFATO_HELENA_AMATOMIA/MUJERES','/Volumes/LaCie/DATOS_OLFATO_HELENA_AMATOMIA/MUJERES_PREPARED');
file1=input('Please, introduce the complete path to the .txt file generated in the previous step for GROUP 1\nGROUP 1:','s');
file2=input('Please, introduce the complete path to the .txt file generated in the previous step for GROUP 2:\nGroup 2:','s');
tpmfile='/Users/lauranunez/Documents/MATLAB/spm12/tpm/TPM.nii';
templatePrefix='TrialTemplate';


group1={};
[group1.file group1.age group1.sex]=textread(file1,'%s %d %s');
group2={};
[group2.file group2.age group2.sex]=textread(file2,'%s %d %s');

prefiles = [group1;group2];
segmented = NewSegment12({prefiles.file}, tpmfile);
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

ng1 =length(group1);
for j=1:ng1
    group1(j).file = GMWM_files{1}{j};
end;
for k=1:length(group2)
        group2(k).file = GMWM_files{1}{ng1 + k};
end;

dirout = '/Volumes/LaCie/DATOS_OLFATO_HELENA_AMATOMIA/RESULTADOS_PRUEBAS_PIPELINE/ResultadosEstadistica';
sgm = Statistic12(group1, group2, [dirout '/GM']);
gr_wm1 = group1;
gr_wm2 = group2;
for j=1:ng1
    gr_wm1(j).file = GMWM_files{2}{j};
end;
for k=1:length(gr_wm2)
        gr_wm2(k).file = GMWM_files{2}{ng1 + k};
end;
swm = Statistic12(gr_wm1, gr_wm2, [dirout '/WM']);

volumenescsf = {csf'};
templatescsf = CreateDartelTemplate12(volumenescsf, [templatePrefix 'CSF']);
deformation_fields_csf = ExistingDartelTemplate12(volumenescsf, templatescsf);
CSF_files = Normalise12(volumenescsf,deformation_fields_csf, templatescsf{6});

gr_csf1 = group1;
gr_csf2 = group2;
for j=1:ng1
    gr_csf1(j).file = CSF_files{1}{j};
end;
for k=1:length(gr_csf2)
        gr_csf2(k).file = CSF_files{1}{ng1 + k};
end;
scsf = Statistic12(gr_csf1, gr_csf2, [dirout '/CSF']);

final_files = [GMWM_files CSF_files];