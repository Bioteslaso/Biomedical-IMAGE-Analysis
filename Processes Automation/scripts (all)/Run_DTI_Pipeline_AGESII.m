function Run_DTI_Pipeline_AGESII(Dir)
%
% Syntax :
%   Run_DTI_Pipeline_AGESII(Dir)
%
% This f'n runs the suggested DTI pipeline for AGES-II. It includes Motion  
% Correction and Tracula. It assumes that FS has already been run on the 
% subjects.
%
% Input Parameters:
%     Dir           : Directory with:
%                      1.- FS Outputs
%                      2.- subfolder called 'subjects_Original_Data', which contains the
%                           original '.nii' images.
%                      3.- Tracula Configuration File 'Tracula_template_AGESII.txt'
%                      4.- Motion Correction Configuration file 
%                          'DWI_Correction_Pipeline_AGESII'
%                       
%
%
%
% Requirements:
%     Freesurfer installed in the system
%     spm12 Toolbox added to the MATLAB path
%     ANTS installed in the system
%__________________________________________________
% Authors: Javi Santonja, UADO, HUGGM, 2018
tic
if nargin < 1
    Dir = '/Users/administrator/Documents/AGESII/AgesII_DTI';
end

MotionCorrectionBool = 1; % Perform Motion Correction
setenv('PATH','/Users/administrator/Applications/ants/bin:/usr/local/fsl/bin:/Applications/freesurfer/bin:/Applications/freesurfer/fsfast/bin:/Applications/freesurfer/tktools:/Applications/freesurfer/mni/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:')
%'/Users/frs16g/antsbin/bin:/Applications/freesurfer/bin:/Applications/freesurfer/fsfast/bin:/Applications/freesurfer/tktools:/Applications/fsl/bin:/Applications/freesurfer/mni/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/Users/frs16g/abin')
setenv('FREESURFER_HOME','/Applications/freesurfer')
setenv('FSFAST_HOME','/Applications/freesurfer/fsfast')
setenv('ANTSPATH','/Users/administrator/Applications/ants/bin/')
setenv('FSLDIR','/usr/local/fsl')
setenv('FSLOUTPUTTYPE','NIFTI_GZ')
setenv('FSF_OUTPUT_FORMAT', 'nii.gz')
setenv('SUBJECTS_DIR','/Applications/freesurfer/subjects')
setenv('MNI_DIR','/Applications/freesurfer/mni')
setenv('FSL_DIR','/usr/local/fsl')
spms(12)
%% ========= INITIALIZATION (# OF SUBJECTS AND MODIFYING CONF FILES) =========
startTime = clock();
% List directories to find the number of subjects
files= dir([Dir filesep 'Subjects_Original_Data' filesep 'sub-ag2*']);
dirFlags=[files.isdir];
subFolders= files(dirFlags);
cellnames = {subFolders.name}; 
%cellnames(1:2) = [];

% Loop through each subject
for s = 1:length(cellnames)
    
    sub = deblank(cellnames{s});
    
    % Modify the Conf files: Motion Correction
    fin = fopen([Dir filesep 'Configuration_File_DWI_Template.txt']);
    fout = fopen([Dir filesep 'Configuration_File_DWI.txt'],'w+');
    while ~feof(fin)
        s = fgetl(fin);
        s = strrep(s, 'HEREGOESTHEDIRECTORY', Dir); % Add the directory
        s = strrep(s, 'HEREGOESTHESUBJECTSNAME', sub); % Add the subject
        fprintf(fout,'%s\n',s);
    end
    fclose(fin);
    fclose(fout);
    
    % Modify the Conf files: Tracula
    fin = fopen([Dir filesep 'Tracula_AGESII_Template.txt']);
    fout = fopen([Dir filesep 'Tracula_AGESII.txt'],'w+');
    while ~feof(fin)
        s = fgetl(fin);
        s = strrep(s, 'HEREGOESTHEDIRECTORY', Dir); % Add the directory
        s = strrep(s, 'HEREGOESTHESUBJECTSNAME', sub); % Add the subject
        fprintf(fout,'%s\n',s);
    end
    fclose(fin);
    fclose(fout);

    %% ========= RUN MOTION CORRECTION =========
    MoCoStarTime = clock();
    toc
    if MotionCorrectionBool 
    
        opts = DWI_Correction_Pipeline_AGESII([Dir filesep 'Configuration_File_DWI.txt']);
    
        TracDir = [Dir filesep sub];
        Connectomesub = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId];
        
        %% Estimate FA and others
        Image=[Connectomesub filesep 'preproc' filesep sub '_DC.nii'];
        Mask=[Connectomesub filesep 'preproc' filesep sub '_BDMask.nii'];
        BVEC=[Connectomesub filesep 'preproc' filesep sub '_DC.bvec'];
        BVAL=[Connectomesub filesep 'preproc' filesep sub '_DC.bval'];
        Outtemp = [Connectomesub filesep 'fsl' filesep sub]; mkdir([Connectomesub filesep 'fsl']);
        cad = ['dtifit --data=' Image ' --out=' Outtemp ' --mask=' Mask ' --bvecs=' BVEC ' --bvals=' BVAL];
        system(cad);
        
        system(['gzip ' Connectomesub filesep '*' filesep '*.nii']);
        
        system(['mkdir ' TracDir filesep 'dmri']);
        system(['mkdir ' TracDir filesep 'dlabel']);
        system(['mkdir ' TracDir filesep 'dlabel' filesep 'diff']);
        
        disp('......copying files......');
        
        system(['cp ' Connectomesub filesep 'preproc' filesep sub '*.bval ' TracDir filesep 'dmri/bvals']);
        system(['cp ' Connectomesub filesep 'preproc' filesep sub '*.bvec ' TracDir filesep 'dmri/bvecs']);
%        system(['mv ' TracDir filesep 'dmri/*_EC.bvec ' TracDir filesep 'dmri/bvecs']);
        system(['cp ' TracDir filesep 'dmri/bvecs ' TracDir filesep 'dmri/bvecs.norot']);
        
        system(['cp ' Connectomesub filesep 'preproc' filesep sub '*_DC.nii.gz ' TracDir filesep 'dmri/dwi.nii.gz']);
        system(['ln -sf ' TracDir filesep 'dmri/dwi.nii.gz ' TracDir filesep 'dmri/data.nii.gz']);
        system(['cp ' Connectomesub filesep 'fsl' filesep sub '*_FA.nii.gz ' TracDir filesep 'dmri/dtifit_FA.nii.gz']);
        system(['cp ' Connectomesub filesep 'fsl' filesep sub '*_L1.nii.gz ' TracDir filesep 'dmri/dtifit_L1.nii.gz']);
        system(['cp ' Connectomesub filesep 'fsl' filesep sub '*_L2.nii.gz ' TracDir filesep 'dmri/dtifit_L2.nii.gz']);
        system(['cp ' Connectomesub filesep 'fsl' filesep sub '*_L3.nii.gz ' TracDir filesep 'dmri/dtifit_L3.nii.gz']);
        system(['cp ' Connectomesub filesep 'fsl' filesep sub '*_MD.nii.gz ' TracDir filesep 'dmri/dtifit_MD.nii.gz']);
        system(['cp ' Connectomesub filesep 'fsl' filesep sub '*_MO.nii.gz ' TracDir filesep 'dmri/dtifit_MO.nii.gz']);
        system(['cp ' Connectomesub filesep 'fsl' filesep sub '*_S0.nii.gz ' TracDir filesep 'dmri/dtifit_S0.nii.gz']);
        system(['cp ' Connectomesub filesep 'fsl' filesep sub '*_V1.nii.gz ' TracDir filesep 'dmri/dtifit_V1.nii.gz']);
        system(['cp ' Connectomesub filesep 'fsl' filesep sub '*_V2.nii.gz ' TracDir filesep 'dmri/dtifit_V2.nii.gz']);
        system(['cp ' Connectomesub filesep 'fsl' filesep sub '*_V3.nii.gz ' TracDir filesep 'dmri/dtifit_V3.nii.gz']);
        system(['cp ' Dir filesep 'Subjects_Original_Data' filesep sub filesep 'dwi' filesep sub '_dir-AP_dwi.nii ' TracDir filesep 'dmri/dwi_orig.nii.gz']);
        system(['cp ' Connectomesub filesep 'preproc' filesep sub '*_b0.nii.gz ' TracDir filesep 'dmri/lowb.nii.gz']);
        system(['cp ' Connectomesub filesep 'preproc' filesep sub '*_BDMask.nii.gz ' TracDir filesep 'dlabel/diff/lowb_brain_mask.nii.gz']);
        
        system(['fslmaths ' TracDir filesep 'dmri/lowb.nii.gz ' '-mul ' TracDir filesep 'dlabel/diff/lowb_brain_mask.nii.gz ' TracDir filesep 'dmri/lowb_brain.nii.gz']);
        
        disp('......end of copying files......');
    
    end
    MoCoEndTime = clock();
    toc
    %% ========= RUN TRACULA =========
    traculaStartTime = clock();
    toc
    disp('TRACULA.......Part 1: Preprocessing.........');
    system(['trac-all -prep -c ' Dir filesep 'Tracula_AGESII.txt -nocorr -notensor']);
    disp('TRACULA.......Part 2: Bedpost.........');
    system(['mv ' Dir filesep sub filesep 'dmri/bvecs.norot ' Dir filesep sub filesep 'dmri/bvecs']);
    system(['trac-all -bedp -c ' Dir filesep 'Tracula_AGESII.txt']);
    
    disp('TRACULA.......Part 3: Pathway Reconstruction.........');
    system(['trac-all -path -c ' Dir filesep 'Tracula_AGESII.txt']);
    traculaEndTIme = clock();
    toc

end
return