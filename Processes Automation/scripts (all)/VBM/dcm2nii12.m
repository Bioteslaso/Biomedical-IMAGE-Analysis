function [niifile, patientage, patientsex] = dcm2nii12(dirdcm,outdir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dcm2nii
%--------------------------------------------------------------------------
% 
% Executed the 'DICOM Import' task in SPM8.
%
%__________________________________________________________________________
% INPUTS:
%  - dirdcm: path to folder with dicom files.
%  - outdir: path to folder where data will be:'outdir/niifile/[new file]'
%     
% OUTPUTS:
%  - niifile: full path of the new .nii file.
%  - patientage: subject age.
%  - patientsex: subject sex.
%__________________________________________________________________________
%
% Author: Laura Nunez Gonzalez
% URJC - 08 / March / 2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try,
    
    
        dirTrabajo = what;
        cd(dirdcm);
        nombreFicheros = dir('*.dcm');
        cd(dirTrabajo.path);
        [fi, co] = size(nombreFicheros);
        dicomfiles = {};
        for c=1:co,
            for f=1:fi,
                nombre = nombreFicheros(f,c).name;
                dicomfiles{f} = [dirdcm '/' nombre];
            end;
        end;
        info = dicominfo(dicomfiles{1});
        patientsex = info.PatientSex == 'M';
        patientage = str2num(info.PatientAge(1:3));
        
        outdir2 = [outdir '/niifile'];
        mkdir(outdir2)
        

        if isempty(which('spm')),
             throw(MException('SPMCheck:NotFound', 'SPM not in matlab path'));
        end
        [name, version] = spm('ver');
        fprintf('SPM version: %s Release: %s\n',name, version);
        fprintf('SPM path: %s\n', which('spm'));
        spm('defaults', 'FMRI');


        if strcmp(name, 'SPM8') || strcmp(name(1:5), 'SPM12'),
           spm_jobman('initcfg');
           spm_get_defaults('cmdline', 1);
        end

%%      Tarea de spm a ejecutar       

        matlabbatch{1}.spm.util.import.dicom.data = dicomfiles';
        matlabbatch{1}.spm.util.import.dicom.root = 'flat';
        matlabbatch{1}.spm.util.import.dicom.outdir = {outdir2};
        matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
        matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
        matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;


        spm_jobman('run', matlabbatch);

        
        if strcmp(name, 'SPM8') || strcmp(name(1:5), 'SPM12'),
            close('all', 'force');
        end;
        
%%      Recoleccion de los datos            
        f = dir([outdir2 '/*.nii']);
        niifile = [outdir2 '/' f(1).name];
        
,catch ME,
fprintf(2,'MATLAB code threw an exception:\n');
fprintf(2,'%s\n',ME.message);
if length(ME.stack) ~= 0, fprintf(2,'File:%s\nName:%s\nLine:%d\n',ME.stack.file,ME.stack.name,ME.stack.line);, end;
end;
