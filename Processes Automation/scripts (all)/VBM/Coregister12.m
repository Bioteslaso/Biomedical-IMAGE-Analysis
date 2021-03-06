function Coregister12(volume, ref)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coregister
%--------------------------------------------------------------------------
% 
% Executed the 'Coregister' task in SPM8.
%
%__________________________________________________________________________
% INPUTS:
%  - volume: image that is jiggled about to best match the reference. 
%  - ref: image  that  is  assumed  to  remain stationary while the
%         source image is moved to match it.
%     
%
%__________________________________________________________________________
%
% Author: Laura Nunez Gonzalez URJC - 08 / March / 2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try,
        %% Generated by nipype.interfaces.spm
        if isempty(which('spm')),
             throw(MException('SPMCheck:NotFound', 'SPM not in matlab path'));
        end
        [name, version] = spm('ver');
        fprintf('SPM version: %s Release: %s\n',name, version);
        fprintf('SPM path: %s\n', which('spm'));
        spm('Defaults','fMRI');

        if strcmp(name, 'SPM8') || strcmp(name(1:5), 'SPM12'),
           spm_jobman('initcfg');
           spm_get_defaults('cmdline', 1);
        end
        
%%      Tarea de spm a ejecutar

        matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[ref ',1']};
        matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[volume ',1']};
        matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 1;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';


        spm_jobman('run', matlabbatch);

        
        if strcmp(name, 'SPM8') || strcmp(name(1:5), 'SPM12'),
            close('all', 'force');
        end;
            
,catch ME,
fprintf(2,'MATLAB code threw an exception:\n');
fprintf(2,'%s\n',ME.message);
if length(ME.stack) ~= 0, fprintf(2,'File:%s\nName:%s\nLine:%d\n',ME.stack.file,ME.stack.name,ME.stack.line);, end;
end;
