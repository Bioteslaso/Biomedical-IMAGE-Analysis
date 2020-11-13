function resultado = Statistic12(group1, group2, dirout)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistic
%--------------------------------------------------------------------------
% 
% Run t-test of group1 vs group2 using SPM. Include age as covariate. Save
% the results in dirout.
%
%__________________________________________________________________________
% INPUTS:
%  - group1: vector of subjects with following struct: file, age, sex.
%  - group2: vector of subjects with following struct: file, age, sex.
%  - dirout: string with an existing path to save the results.
%     
% OUTPUTS:
%  - resultado: path to SPM.mat file.
%
%__________________________________________________________________________
%
% Author: Laura Nunez Gonzalez
% URJC - 10 / March / 2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try,
        %% Model
        
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
        
        
        matlabbatch{1}.spm.stats.factorial_design.dir = {dirout};
        matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = {group1.file}';
        matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = {group2.file}';
        matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
        matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
%%
        matlabbatch{1}.spm.stats.factorial_design.cov(1).c = [group1.age group2.age]';
%%
        matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'edad';
        matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
        matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
%%
        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;


        spm_jobman('run', matlabbatch);

        
        if strcmp(name, 'SPM8') || strcmp(name(1:5), 'SPM12'),
            close('all', 'force');
        end;
            
,catch ME,
fprintf(2,'MATLAB code threw an exception:\n');
fprintf(2,'%s\n',ME.message);
if length(ME.stack) ~= 0, fprintf(2,'File:%s\nName:%s\nLine:%d\n',ME.stack.file,ME.stack.name,ME.stack.line);, end;
end;

try,
        %% Estimate
        
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


        dirmap = dir([dirout '/*.mat']);
        namemap = dirmap.name;
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {[dirout '/' namemap]};
        matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

        spm_jobman('run', matlabbatch);

        
        if strcmp(name, 'SPM8') || strcmp(name(1:5), 'SPM12'),
            close('all', 'force');
        end;
            
,catch ME,
fprintf(2,'MATLAB code threw an exception:\n');
fprintf(2,'%s\n',ME.message);
if length(ME.stack) ~= 0, fprintf(2,'File:%s\nName:%s\nLine:%d\n',ME.stack.file,ME.stack.name,ME.stack.line);, end;
end;


try,
        %% Contrast
        
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


        matlabbatch{2}.spm.stats.con.spmmat(1) = {[dirout '/' namemap]};
        matlabbatch{2}.spm.stats.con.consess{1}.tcon.name = 'H > M';
        matlabbatch{2}.spm.stats.con.consess{1}.tcon.weights = [1 -1 0];
        matlabbatch{2}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{2}.spm.stats.con.consess{2}.tcon.name = 'M > H';
        matlabbatch{2}.spm.stats.con.consess{2}.tcon.weights = [-1 1 0];
        matlabbatch{2}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        matlabbatch{2}.spm.stats.con.delete = 0;
        
        spm_jobman('run', matlabbatch);

        
        if strcmp(name, 'SPM8') || strcmp(name(1:5), 'SPM12'),
            close('all', 'force');
        end;
            
,catch ME,
fprintf(2,'MATLAB code threw an exception:\n');
fprintf(2,'%s\n',ME.message);
if length(ME.stack) ~= 0, fprintf(2,'File:%s\nName:%s\nLine:%d\n',ME.stack.file,ME.stack.name,ME.stack.line);, end;
end;

resultado = [dirout '/' namemap];

