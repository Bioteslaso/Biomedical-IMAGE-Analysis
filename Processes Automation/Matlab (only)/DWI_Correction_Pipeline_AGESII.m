function opts = DWI_Correction_Pipeline_AGESII(Dir)
%
% Syntax :
%   opts = DWI_Correction_Pipeline(opts);
%
% In order to detect and correct artifacts in DTI data, the following pipeline 
% was implemented. (1) Slice-wise intensity related artifacts were detected by 
% computing the Normalized Correlation between successive slices and averaging 
% it over all the diffusion volumes (NCmean) (Liu et al.). A large deviation from the 
% NCmean anywhere in the DTI data indicates a dramatic intensity change due to 
% subject motion or an intensity artifact. The outliers were defined as NC values
% smaller than NCmean - α·σ, where σ is the standard deviation of the computed NC 
% values and α is a constant parameter to be adjusted. In line with Liu et al. (year?),
% the current study used α = 3.5. Diffusion volumes containing one or more outliers
% were excluded from further analysis. (2) Eddy-current and head motion artifact 
% correction was performed using the eddy_correct programme implemented in FSL,
% using default settings. In this step all diffusion volumes were affinely registered
% to the T2-b0 image. The corresponding diffusion gradient vectors were properly 
% reoriented using the resulting affine transformations (Leemans and Jones, year). 
% (3) Machine-related spatial distortions (i.e. main field B0 inhomogeneity) 
% and subject-related spatial distortions (i.e., susceptibility and chemical 
% shift artifacts) were corrected following the image-based technique proposed 
% by Huang et al. (year). This technique consists of warping a subject's T2-b0 
% image to the anatomical T2-weighted image of the same individual. The technique
% produces (a) a warp-field, which was applied to all the subject's diffusion volumes 
% and (b) a Jacobian map of the warp-field, which was multiplied with the subject's 
% warped diffusion volumes in order to restore true image intensity after warping. 
% To achieve high dimensional and robust warping a large-deformation diffeomorphic
% mapping was computed using the symmetric normalization (SyN) technique as implemented 
% in the Advance Normalization Tools (ANTS) software package (Avants et al. year). 
% In a recent evaluation of 14 algorithms for high dimensional warping of anatomical 
% brain images, SyN was found to be the top-rank performer, producing the highest 
% accuracy across subjects and providing the best results according to overlap 
% and distance measures (Klein et al. year). In addition ANTS-SyN also allows for
% exclusively applying the warp-field to the anterior-posterior axis of the 
% subject, i.e. the phase-encoding direction (y-coordinate), reducing the geometric 
% distortion present only along that axis while preserving the signal in the other axes.
%
% Input Parameters:
%     opts            :  Input Options
%
% Output Parameters:
%     opts            :  Output Options
%
% Requirements:
%     spm12 Toolbox added to the MATLAB path
%     ANTS installed in the system
%__________________________________________________
% Authors: Yasser Iturria-Medina and Erick J Canales-Rodriguez
% LIM, HUGGM
% March 22th 2012
% Version $1.0
% 
% Edited by: Javi Santonja, UADO, HUGGM, 2018


if nargin ==0
    [Dir] = '/media/jjanssen/DATA/projects/Motion_Correction/Configuration_File_DWI.txt';
end
try
    [opts] = Reading_Pipe_Configuration_File(Dir);
catch
    error('Please enter a correct configuration file');
    return;
end
%=========================================================================%
%% ============== DATA VERIFICATION ==========%
% Verifying subject Id
if isempty(opts.pipe.subjId)
    error('ID field is missing in configuration file');
    return;
end
% Verifying Output directory
if isempty(opts.pipe.outdir)
    error('Output directory field is missing in configuration file');
    return;
elseif ~exist(opts.pipe.outdir,'dir')&~isempty(opts.pipe.outdir)
    try
        mkdir(opts.pipe.outdir);
    catch
        error('Output directory do not exist');
        return;
    end
end
% Verifying Anatomical File (T1 Image);
if isempty(opts.anat.ndata.t1)
    error('T1 image field is missing in configuration file');
    return;
end
if ~exist(opts.anat.ndata.t1,'file')
    error([opts.anat.ndata.t1 ':  Do not exist']);
    return;
end
% Verifying Anatomical File (T2 Image);
if ~exist(opts.anat.ndata.t2,'file')&~isempty(opts.anat.ndata.t2)
    error([opts.anat.ndata.t2 ':  Do not exist']);
    return;
end

% Verifying Diffusion Files (Diffusion Image)
if isempty(opts.diff.ndata.diff)&isempty(opts.diff.pdata.diff)
    error('Diffusion weighted image field is missing in configuration file');
    return;
end

% Verifying Diffusion Files (Gradient Table)
if isempty(opts.diff.ndata.bvec)&isempty(opts.diff.pdata.bvec)
    error('Gradient table field is missing in configuration file');
    return;
end
% Verifying Diffusion Files (B-Value Table)
if isempty(opts.diff.ndata.bval)&isempty(opts.diff.pdata.bval)
    error('B-values table is missing in configuration file');
    return;
end

if ~exist(opts.diff.pdata.diff,'file')
    opts.diff.pdata.diff = opts.diff.ndata.diff;
end
if ~exist(opts.diff.pdata.bvec,'file')
    opts.diff.pdata.bvec = opts.diff.ndata.bvec;
end
if ~exist(opts.diff.pdata.bval,'file')
    opts.diff.pdata.bval = opts.diff.ndata.bval;
end

% Selecting freesurfer output directory
if ~isempty(opts.pipe.freesdir)&exist(opts.pipe.freesdir,'dir')
    setenv('SUBJECTS_DIR',opts.pipe.freesdir);
else
    [a,temp] = system('echo $SUBJECTS_DIR');
    opts.pipe.freesdir = temp;
end

% indfsep = strfind(opts.anat.ndata.t1,filesep);  % Detecting Paths
IFiles = logical([exist(opts.anat.ndata.t1,'file') exist(opts.anat.ndata.t2,'file') exist(opts.diff.ndata.diff,'file') exist(opts.diff.ndata.bvec,'file') exist(opts.diff.ndata.bval,'file')]);

%% ====== END OF DATA VERIFICATION ==========%

%% ========= ANATOMICAL PROCESSING (FreeSurfer Autorecon 1) ==============%

if IFiles(1)==1
    [Outboolean] = FreeSurfer_verif(opts.pipe.freesdir,opts.pipe.subjId,'autorecon1'); %% Verifing freesurfer autorecon1 result files
    [pth_T2, nm_T2, ext_T2] = fileparts(opts.anat.ndata.t2);
    if Outboolean %
        if ~exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'brainmask.auto.mgz'],'file');  % If there is not brain mask compute it!.
            if strcmp(ext_T2,'.gz')
                system(['gunzip ' pth_T2 filesep nm_T2 ext_T2]);
                opts.anat.ndata.t2 = [pth_T2 filesep nm_T2];
                gzipbool = 1;
            end
            cad = ['mri_convert -i ' opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'T1.mgz' ' -o ' opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1.nii'];
            system(cad);
            outputs = flip_to_axial_nii([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1.nii'],0);
            delete([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1.nii']);
            movefile(outputs,[opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1.nii']);
            % Creating p0 Brain Mask
            disp(['T1 Brainmask Computation.....................................']);tic;
            [opts] = T1_Brainmask_Computation(opts);
            if gzipbool; system(['gzip ' pth_T2 filesep nm_T2 ext_T2]); end;
            toc;
        end
    else % If there is no 'autorecon1' files, run 'autorecon1'
        gzipbool = 0; 
        if strcmp(ext_T2,'.gz') 
            system(['gunzip ' pth_T2 filesep nm_T2 ext_T2]); 
            opts.anat.ndata.t2 = [pth_T2 filesep nm_T2]; 
            gzipbool = 1; 
        end
        mkdir(opts.pipe.freesdir);
        disp(['    Running FreeSurfer ==> Autorecon 1...................']);
        freefilename = [opts.pipe.input filesep opts.pipe.subjId filesep opts.pipe.subjId '.mgz'];
        cad = ['mri_convert -i ' opts.anat.ndata.t1 ' -o ' freefilename];
        system(cad);
        cad = ['recon-all -i ' freefilename  ' -subjid ' opts.pipe.subjId];
        system(cad);
        delete(freefilename);
        %cad = ['recon-all -autorecon1 -subjid ' opts.pipe.subjId];
        cad = ['recon-all -all -subjid ' opts.pipe.subjId];
        system(cad);
        cad = ['mri_convert -i ' opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'T1.mgz' ' -o ' opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1.nii'];
        system(cad);
        outputs = flip_to_axial_nii([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1.nii'],0);
        delete([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1.nii']);
        movefile(outputs,[opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1.nii']);
        disp(['T1 Brainmask Computation.....................................']);tic;
        [opts] = T1_Brainmask_Computation(opts);
        if gzipbool; system(['gzip ' pth_T2 filesep nm_T2 ext_T2]); end;
        toc;
    end
    if ~exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_binBMask.nii'],'file')
        if exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_BMask.nii'],'file') % If there is a Manual segmentation mask, save it in FS's 'temp' folder as "Subjid_Bmask.nii"
            Vtemp = spm_vol([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_BMask.nii']);
            Itemp = logical(spm_read_vols(Vtemp));
            Vtemp.fname = [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_binBMask.nii'];
            Vtemp.descrip = [Vtemp.descrip '/ Binary Mask created from Manual segmentation (BMask)'];
            Vtemp.dt = [2 0];
            spm_write_vol(Vtemp,Itemp); clear Itemp;
            opts.anat.ndata.bmask = [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_binBMask.nii'];
        elseif ~exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_BMask.nii'],'file')&(exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'brainmask.auto.mgz'],'file'))
            if ~exist([opts.pipe.outdir filesep opts.pipe.subjId filesep 'mri' filesep '_p0.nii'],'file')
                gzipbool = 0; if strcmp(ext_T2,'.gz'); system(['gunzip ' pth_T2 filesep nm_T2 ext_T2]); opts.anat.ndata.t2 = [pth_T2 filesep nm_T2]; gzipbool = 1; end
                cad = ['mri_convert -i ' opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'T1.mgz' ' -o ' opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1.nii'];
                system(cad);
                outputs = flip_to_axial_nii([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1.nii'],0);
                delete([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1.nii']);
                movefile(outputs,[opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1.nii']);
                disp(['T1 Brainmask Computation.....................................']);tic;
                [opts] = T1_Brainmask_Computation(opts);
                if gzipbool; system(['gzip ' pth_T2 filesep nm_T2 ext_T2]); end;
            end
            Vtemp = spm_vol([opts.pipe.outdir filesep opts.pipe.subjId filesep 'mri' filesep '_p0.nii']);
            Itemp = logical(spm_read_vols(Vtemp));
            Vtemp.fname = [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_binBMask.nii'];
            Vtemp.descrip = [Vtemp.descrip '/ Binary Mask created from p0 File'];
            Vtemp.dt = [2 0];
            spm_write_vol(Vtemp,Itemp); clear Itemp;
            opts.anat.ndata.bmask = [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_binBMask.nii'];
        end
    else
        opts.anat.ndata.bmask = [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_binBMask.nii'];
    end
    disp('                                                                    ')
end
%% ====================END OF ANATOMICAL PROCESSING(Autorecon1) ==========%

%% ====================DIFFUSION PREPROCESSING STEPS ======================%
%====       1------Slice-wise intensity related artifacts checking    ====%
if sum(IFiles(3:end))==3
    disp(['Correcting Bad Volumes ............................................']);tic;
    gzipbool2 = 0;
    [pthd, nmd, extd] = fileparts(opts.diff.ndata.diff);
    if strcmp(extd,'.gz'); 
        system(['gunzip ' pthd filesep nmd extd]); 
        opts.diff.ndata.diff = [pthd filesep nmd]; 
        opts.diff.pdata.diff = [pthd filesep nmd]; 
        gzipbool2 = 1; 
    end
    try
        mkdir([opts.pipe.outdir filesep 'connectome'],opts.pipe.subjId);
    end
    try
        mkdir([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId],'preproc');
        mkdir([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId],'tmp');
    end
    if opts.diff.preproc.icorr.bool
        switch opts.diff.preproc.icorr.method
            case 'sinten'
                CorrectFiles = correct_diff_images(opts.diff.pdata.diff, opts.diff.pdata.bvec, opts.diff.pdata.bval,[opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_IC.nii']);
                opts.diff.pdata.diff = deblank(CorrectFiles(1,:));
                opts.diff.pdata.bvec = deblank(CorrectFiles(2,:));
                opts.diff.pdata.bval = deblank(CorrectFiles(3,:));
        end
    end
    toc;
    disp('                                                                    ');
    %=====================================================================%
    %====      2----- Extracting B0    ===================================%
    disp(['Separating B0 Image ............................................']);tic;
    if ~exist(opts.diff.ndata.b0,'file')
        Vdif = spm_vol(opts.diff.ndata.diff);
        grads = load(opts.diff.ndata.bvec);
        if size(grads,2)>size(grads,1)
            grads = grads';
        end
        bvals = load(opts.diff.ndata.bval);
        if size(bvals,2)>size(bvals,1)
            bvals = bvals';
        end
        ind = find(sum([grads bvals]') ==0);
        Vb0s = Vdif(ind);
        Vb0 = Vdif(1);
        It = zeros(Vb0.dim(1:3));
        for k = 1:length(ind)
            It = It+spm_read_vols(Vb0s(ind(k)));
        end
        It = It/k;
        Vb0.fname = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_b0.nii'];
        spm_write_vol(Vb0,It); clear It;
        opts.diff.ndata.b0 = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_b0.nii'];
    end
    opts.diff.pdata.b0 = opts.diff.ndata.b0;
    if gzipbool2; system(['gzip ' pthd filesep nmd]); end
    toc;
    disp('                                                                    ');
    %=====================================================================%
    %=====       3------Eddy current correction            ===============%
    if opts.diff.preproc.eddyc.bool
        disp(['Eddy Correction Step ............................................']);tic;
        grads = load(opts.diff.pdata.bvec);
        if size(grads,2)>size(grads,1)
            grads = grads';
        end
        bvals = load(opts.diff.pdata.bval);
        if size(bvals,2)>size(bvals,1)
            bvals = bvals';
        end
        indb0 = find(sum([grads bvals]') ==0);
        switch opts.diff.preproc.eddyc.method
            case 'eddy_correct'
                %  Eddy current correction using FSL %
                if exist([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_EC.ecclog'],'file')
                    try
                        delete([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_EC.ecclog']);
                    end
                end
                cad =['eddy_correct ' opts.diff.pdata.diff ' ' opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_EC.nii' ' ' num2str(indb0(1)-1)];
                system(cad);
                gunzip([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_EC.nii.gz']);
                delete([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_EC.nii.gz']);
                switch opts.diff.preproc.eddyc.gradrot
                    case 'saad'
                        %=========== Gradient Rotation 1 =================%
                        cad = ['bash fdt_rotate_bvecs.sh ' opts.diff.pdata.bvec ' ' opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_EC.bvec ' opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_EC.ecclog'];
                        system(cad);
                    case 'our'
                        %=========== Gradient Rotation 2 =================%
                        fid = fopen([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_EC.ecclog']);
                        cont = 0; lines = '';
                        while 1
                            cont = cont + 1;
                            line = fgetl(fid);
                            if ~ischar(line),   break,   end
                            lines = strvcat(lines,line);
                        end
                        fclose(fid);
                        clear tline cont;
                        ind = find(sum(ismember(lower(lines),'processing')')>0);
                        ind2 = find(sum(ismember(lower(lines),'result')')>0);
                        ind = unique([ind; ind2]);
                        lines(ind,:) = [];
                        temp = str2num(lines);
                        for i = 0:size(lines,1)/4-1
                            rgrad(i+1,:) = temp(4*i+1:4*i+3,1:3)*[grads(i+1,:)]';
                        end
                        rgrad = rgrad';
                        save([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_EC.bvec'],'rgrad','-ascii');
                        %=================================================================%
                end
            case 'ants'
                %  Eddy current correction using ANTS %
                Vdif = spm_vol(opts.diff.pdata.diff);
                rvols=find(ismember([1:length(Vdif)],indb0(1))==0);
                try
                    copyfile(opts.diff.pdata.diff,[opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_EC.nii']);
                end
                Vtemp = spm_vol([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_EC.nii']);
                rgrad = zeros(length(Vdif),3);
                Nv = length(rvols);
                for i = 1:Nv
                    disp(['Correcting volume ' num2str(i) ' of ' num2str(Nv)]);
                    Vt = Vtemp(rvols(i));
                    It = spm_read_vols(Vtemp(rvols(i)));
                    Vt.fname = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_temp.nii'];
                    Vt = rmfield(Vt,'private');Vt.n = [1 1];
                    Vt = spm_create_vol(Vt);
                    spm_write_vol(Vt,squeeze(It));
                    cad = ['ANTS 3 -m MI[' opts.diff.pdata.b0 ',' Vt.fname ',1,32] -i 0 -o ' [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_EddyCurrent.nii']];
                    system(cad);
                    cad = ['WarpImageMultiTransform 3 ' Vt.fname ' ' [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_EddyCurrent.nii'] ' -R ' Vt.fname ' ' opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_EddyCurrentAffine.txt  --use-BSpline'];
                    system(cad);
                    Vt = spm_vol([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_EddyCurrent.nii']);
                    Ic = spm_read_vols(Vt);
                    spm_write_vol(Vtemp(rvols(i)),Ic);
                    fidr = fopen([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_EddyCurrentAffine.txt'],'rt');
                    cont = 0;
                    while 1
                        cont = cont+1;
                        line = fgetl(fidr);
                        if ~ischar(line),   break,   end
                        if cont ==4
                            params = strread(line,'%s','delimiter',' ');
                            params = char(params);params(1,:) = [];params = str2num(params);
                            Mat = [reshape(params(1:9,:),[3 3]) [params(10) params(11) params(12)]'];  %CHECK PARAMETER ORGANIZATION
                            Mat = [Mat;0 0 0 1];
                        elseif cont == 5
                            fparams = strread(line,'%s','delimiter',' ');
                            fparams = char(fparams);fparams(1,:) = [];fparams = str2num(fparams); %CHECK PARAMETER ORGANIZATION
                        end
                        
                    end
                    fclose(fidr);
                    ngrad = Mat*[grads(rvols(i),:)-[fparams(1) fparams(2) fparams(3)] 1]';
                    rgrad(rvols(i),:) = ngrad(1:3)'+[fparams(1) fparams(2) fparams(3)];
                end
                rgrad(rvols,:) = rgrad(rvols,:)./repmat((sqrt(sum((rgrad(rvols,:).^2)'))'),[1 3]);
                rgrad = rgrad';
                save([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_EC.bvec'],'rgrad','-ascii');
                delete([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_EddyCurrentAffine.txt']);
                delete([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_EddyCurrent.nii']);
                delete([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_temp.nii']);
                clear Ic It;
        end
        Vdif = spm_vol([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_EC.nii']);
        for i =1:length(Vdif);
            Vdif(i).descrip = [Vdif(i).descrip '/===> Proc: Eddy Current Correction ('  datestr(now) ')' ];
        end
        Vdif = spm_create_vol(Vdif);
        try
            copyfile(opts.diff.pdata.bval,[opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_EC.bval']);
        end
        if ~strcmp(opts.diff.pdata.diff,opts.diff.ndata.diff)
            delete(opts.diff.pdata.diff);
            [pth,nm,ext] = fileparts(opts.diff.pdata.diff);
            if exist([pth filesep nm '.mat'],'file')
                delete([pth filesep nm '.mat']);
            end
        end
        if ~strcmp(opts.diff.pdata.bvec,opts.diff.ndata.bvec)
            delete(opts.diff.pdata.bvec);
        end
        if ~strcmp(opts.diff.pdata.bval,opts.diff.ndata.bval)
            delete(opts.diff.pdata.bval);
        end
        opts.diff.pdata.diff = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_EC.nii'];
        opts.diff.pdata.bvec = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_EC.bvec'];
        opts.diff.pdata.bval = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_EC.bval'];
    end
    toc;
    disp('                                                                   ');
    %=====================================================================%
    % -- Masking B0
    cad = ['bet2 ' opts.diff.pdata.b0 ' ' opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep 'B0_temp.nii -m -f 0.2'];
    system(cad);
    gunzip([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep 'B0_temp.nii.gz']);delete([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep 'B0_temp.nii.gz'])
    gunzip([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep 'B0_temp_mask.nii.gz']);delete([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep 'B0_temp_mask.nii.gz'])
    B0Temp = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep 'B0_temp.nii'];
    % -- End of Masking B0
    %=====   4------EPI Distortions correction using registration  =======%
    if sum(IFiles(1:2)) ==2 % Verifing the existence of both T2 and T1 image
        mkdir([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId],'tmp');
        mkdir([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId],'transforms');
        T2bet = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_temp.nii'];
        % Affine registration - from T2 to T1 space
        disp('"Affine registration - from T2 to T1: Mutual information"');
        if ~exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1.nii'],'file')
            cad = ['mri_convert -i ' opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'T1.mgz' ' -o ' opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1.nii'];
            system(cad);
        end
        if  ~exist(opts.transforms.t2_2_t1.affine,'file')
            cad = ['ANTS 3 -m MI[' opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1.nii,' opts.anat.ndata.t2 ',1,32] -i 0 -o ' [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_T2_2_T1.nii']];
            system(cad);
            opts.transforms.t2_2_t1.affine = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_T2_2_T1Affine.txt'];
        end
        T2T1Res = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_T2_2_T1.nii'];
        if ~exist(T2T1Res,'file')
            cad = ['WarpImageMultiTransform 3 ' opts.anat.ndata.t2 ' ' T2T1Res ' -R ' [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1.nii'] ' ' opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_T2_2_T1Affine.txt  --use-BSpline'];
            system(cad);
        end
        %-----------------------------------------------------------------%
        % T2 Masking: Brain Extraction
        VT2 = spm_vol(T2T1Res);
        outnames = change_spacen(opts.anat.ndata.bmask,T2T1Res,0,[opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp']);
        It2 = spm_read_vols(VT2);
        VM = spm_vol(outnames);
        IM = spm_read_vols(VM);
        It2 = It2.*logical(IM);
        Vol = VT2;Vol.fname = T2bet;
        spm_write_vol(Vol,It2);clear It2 IM;
        %-----------------------------------------------------------------%
        % Affine registration - from T2 betted to diffusion space
        disp('"Affine registration - from T2 to S0(T2): Mutual information"')
        if ~exist(opts.transforms.t1_2_diff.affine,'file')
            cad = ['ANTS 3 -m MI[' B0Temp ',' T2bet ',1,32] -i 0 -o ' opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_T2Bet_2_diff.nii.gz'];
            system(cad);
        end
        T2DifRes = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_T2Bet_2_diff.nii'];
        if ~exist(T2DifRes,'file')
            cad = ['WarpImageMultiTransform 3 ' T2bet ' ' T2DifRes ' -R ' opts.diff.pdata.b0 ' ' opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_T2Bet_2_diffAffine.txt  --use-BSpline'];
            system(cad);
        end
        opts.transforms.t1_2_diff.affine = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_T2Bet_2_diffAffine.txt'];
        if opts.diff.preproc.distc.bool
            disp(['EPI Distortion Correction Step ............................................']);tic;
            switch opts.diff.preproc.distc.method
                case 'regist' % Distortion correction using registration
                    %-----------------------------------------------------------------%
                    % Non Linear registration - from T2 Betted to diffusion space:
                    % Deformations are restricted to Antero-Posterior direction
                    defFile = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_Diff_2_T2BetNL.nii.gz'];
                    if ~exist([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_Diff_2_T2BetNLWarpzvec.nii.gz'],'file')
                        cad = [ 'ANTS 3 -m CC[' T2DifRes ',' B0Temp ',1,4] --number-of-iterations 100x100x100x20 -o ' defFile ' --Restrict-Deformation 0x1x0 -t SyN[0.1] -r Gauss[3,0] --number-of-affine-iterations 10000x10000x10000'];
                        system(cad);
                    end
                    [pth,nm,ext] = fileparts(opts.diff.pdata.diff);
                    defFilej = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_Diff_2_T2BetNLWarp.nii.gz'];
               
                    cad = ['CreateJacobianDeterminantImage 3 ' defFilej ' ' [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_jacobian.nii.gz'] ' 0'];
                    system(cad);
                    if exist([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_jacobian.nii.gz'],'file')
                        gunzip([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_jacobian.nii.gz']);
                    end
                    Vjac = spm_vol([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_jacobian.nii']);
                    Ij = spm_read_vols(Vjac);
                    delete([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_jacobian.nii']);
                    Vdif = spm_vol(opts.diff.pdata.diff);
                    Nv = length(Vdif);
                    VEC = Vdif;
                    for i =1:Nv
                        VEC(i).fname = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_DC.nii'];
                        VEC(i).descrip = [VEC(i).descrip '/===> Proc: EPI Distortions Correction ('  datestr(now) ')' ];
                    end
                    VEC = spm_create_vol(VEC);
                    H = waitbar(0,['Applying Non Linear Transformation  to ' num2str(Nv) '  volumes' ],'Resize','on','Position',[233.25 237.75 273 50.25],'Resize','off');
                    for i = 1:Nv
                        waitbar(i/Nv,H,['Volume  ' num2str(i) ' of ' num2str(Nv)]);
                        Idif = spm_read_vols(Vdif(i));
                        Ic = Idif*0;
                        Vtemp = Vdif(1);
                        Vtemp.fname = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_temp.nii'];
                        spm_write_vol(Vtemp,Idif);
                        cad = ['WarpImageMultiTransform 3 '  Vtemp.fname ' ' [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_EPICORR.nii'] ' -R ' Vtemp.fname ' ' [defFile(1:end-7) 'Warp.nii.gz  --use-BSpline'] ];
                        system(cad);
                        Vtot = spm_vol([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_EPICORR.nii']);
                        if opts.diff.preproc.distc.modul
                            Ic = spm_read_vols(Vtot).*Ij;
                        else
                            Ic = spm_read_vols(Vtot);
                        end
                        spm_write_vol(VEC(i),Ic);
                        delete([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_EPICORR.nii']);
                    end
                    delete([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_temp.nii']);
                    clear Idif Ic;
                    close(H);
                    if ~strcmp(opts.diff.pdata.diff,opts.diff.ndata.diff)
                        delete(opts.diff.pdata.diff);
                        [pth,nm,ext] = fileparts(opts.diff.pdata.diff);
                        if exist([pth filesep nm '.mat'],'file')
                            delete([pth filesep nm '.mat']);
                        end
                    end
                    try
                        movefile(opts.diff.pdata.bvec,[opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_DC.bvec']);
                        movefile(opts.diff.pdata.bval,[opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_DC.bval']);
                    end
                    opts.diff.pdata.diff =  VEC(1).fname;
                    opts.diff.pdata.bval=[opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_DC.bval'];
                    opts.diff.pdata.bvec=[opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_DC.bvec'];
                    %-----------------------------------------------------------------%
                    %
%                 case 'fieldmap' % Field Map correction
            end
            disp('                                                                   ');
            toc;
        end
        % ---Moving Brain mask from T1 space to Diffusion space ----------%
        opts.diff.dbmask.filename = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_BDMask.nii'];
        if ~exist(opts.diff.dbmask.filename,'file')
            cad = ['WarpImageMultiTransform 3 ' opts.anat.ndata.bmask ' ' opts.diff.dbmask.filename ' -R ' opts.diff.pdata.b0 ' ' opts.transforms.t1_2_diff.affine ' --use-NN'];
            system(cad);
        end
        Vtemp = spm_vol(opts.diff.dbmask.filename);
        Itemp = logical(spm_read_vols(Vtemp));
        spm_write_vol(Vtemp,Itemp);
        
        %-----------------------------------------------------------------%
    elseif (IFiles(1) ==1)&(IFiles(2) ==0) % There is T1 image and there is not T2 image
        mkdir([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId],'tmp');
        mkdir([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId],'transforms');
        %% Pruebas
        % -- Masking
        Vdif = spm_vol(opts.diff.pdata.diff);
        VtempM = spm_vol([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep 'B0_temp_mask.nii']);
        IM = spm_read_vols(VtempM);
        Nv = length(Vdif);
        Vtemp = Vdif;
        if exist([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_Diff_masked.nii'],'file');
            delete([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_Diff_masked.nii']);
        end
        if exist([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_Diff_masked.mat'],'file')
            delete([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_Diff_masked.mat']);
        end
        for i = 1:Nv
            Vtemp(i).fname = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_Diff_masked.nii'];
            spm_write_vol(Vtemp(i),spm_read_vols(Vdif(i)).*IM);
        end
        % -- End of masking
        % --- Estimating FA
        Outtemp = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId];
        cad = ['dtifit --data=' [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_Diff_masked.nii'] ' --out=' Outtemp ' --mask=' [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep 'B0_temp_mask.nii'] ' --bvecs=' opts.diff.pdata.bvec ' --bvals=' opts.diff.pdata.bval];
        system(cad);
        List = dir([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep '*.gz']);
        for k = 1:size(List)
            gunzip([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep List(k).name]);
            delete([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep List(k).name]);
        end
        delete([Outtemp '_L1.nii']);
        delete([Outtemp '_L2.nii']);
        delete([Outtemp '_L3.nii']);
        delete([Outtemp '_V1.nii']);
        delete([Outtemp '_V2.nii']);
        delete([Outtemp '_V3.nii']);
        delete([Outtemp '_MD.nii']);
        delete([Outtemp '_MO.nii']);
        delete([Outtemp '_S0.nii']);
        % --- End of Estimating FA
        disp('"Affine registration - from T1 to S0(T1): Mutual information"')
        if ~exist([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_WM_2_diffAffine.txt'],'file')
            cad = ['ANTS 3 -m MI[' [Outtemp '_FA.nii'] ',' opts.anat.ndata.wm ',1,32] -i 0 -o ' opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_WM_2_diff.nii.gz'];
            system(cad);
        end
        WMDifRes = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_WM_2_diff.nii'];
        if ~exist(WMDifRes,'file')
            cad = ['WarpImageMultiTransform 3 ' opts.anat.ndata.wm ' ' WMDifRes ' -R ' opts.diff.pdata.b0 ' ' opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_WM_2_diffAffine.txt  --use-BSpline'];
            system(cad);
        end
        opts.transforms.t1_2_diff.affine = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_WM_2_diffAffine.txt'];
        if opts.diff.preproc.distc.bool
            disp(['EPI Distortion Correction Step ............................................']);tic;
            switch opts.diff.preproc.distc.method
                case 'regist'
                    %-----------------------------------------------------------------%
                    % Non Linear registration - from T1 Betted to diffusion space:
                    % Deformations are restricted to Antero-Posterior direction
                    defFile = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_Diff_2_WMNL.nii.gz'];
                    if ~exist([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_Diff_2_WMNLWarpzvec.nii.gz'],'file')
                        cad = [ 'ANTS 3 -m MI[' WMDifRes ','  [Outtemp '_FA.nii'] ',1,32] --number-of-iterations 100x100x100x20 -o ' defFile ' --Restrict-Deformation 0x1x0 -t SyN[0.1] -r Gauss[3,0] --number-of-affine-iterations 10000x10000x10000'];
                        system(cad);
                    end
                    [pth,nm,ext] = fileparts(opts.diff.pdata.diff);
                    defFilej = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_Diff_2_WMNLWarp.nii.gz'];
                    cad = ['CreateJacobianDeterminantImage 3 ' defFilej ' ' [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_jacobian.nii.gz'] ' 0'];
                    system(cad);
                    if exist([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_jacobian.nii.gz'],'file')
                        gunzip([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_jacobian.nii.gz']);
                    end
                    Vjac = spm_vol([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_jacobian.nii']);
                    Ij = spm_read_vols(Vjac);
                    delete([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'transforms' filesep opts.pipe.subjId '_jacobian.nii']);
                    Vdif = spm_vol(opts.diff.pdata.diff);
                    Nv = length(Vdif);
                    VEC = Vdif;
                    for i =1:Nv
                        VEC(i).fname = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_DC.nii'];
                        VEC(i).descrip = [VEC(i).descrip '/===> Proc: EPI Distortions Correction ('  datestr(now) ')' ];
                    end
                    VEC = spm_create_vol(VEC);
                    H = waitbar(0,['Applying Non Linear Transformation  to ' num2str(Nv) '  volumes' ],'Resize','on','Position',[233.25 237.75 273 50.25],'Resize','off');
                    for i = 1:Nv
                        waitbar(i/Nv,H,['Volume  ' num2str(i) ' of ' num2str(Nv)]);
                        Idif = spm_read_vols(Vdif(i)); Ic = Idif*0;
                        Vtemp = Vdif(1);
                        Vtemp.fname = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_temp.nii'];
                        spm_write_vol(Vtemp,Idif);
                        cad = ['WarpImageMultiTransform 3 '  Vtemp.fname ' ' [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_EPICORR.nii'] ' -R ' Vtemp.fname ' ' [defFile(1:end-7) 'Warp.nii.gz  --use-BSpline'] ];
                        system(cad);
                        Vtot = spm_vol([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_EPICORR.nii']);
                        if opts.diff.preproc.distc.modul
                            Ic = spm_read_vols(Vtot).*Ij;
                        else
                            Ic = spm_read_vols(Vtot);
                        end
                        spm_write_vol(VEC(i),Ic);
                        delete([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_EPICORR.nii']);
                    end
                    delete([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep opts.pipe.subjId '_temp.nii']);
                    clear Idif Ic;
                    
                    % Non Linear deformation applied to BET Mask
                    cad = ['WarpImageMultiTransform 3 '  [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep 'B0_temp_mask.nii'] ' ' [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_BDMask.nii'] ' -R ' opts.diff.pdata.b0 ' ' [defFile(1:end-7) 'Warp.nii.gz  --use-NN'] ];
                    system(cad);
                    
                    close(H);
                    if ~strcmp(opts.diff.pdata.diff,opts.diff.ndata.diff)
                        delete(opts.diff.pdata.diff);
                        [pth,nm,ext] = fileparts(opts.diff.pdata.diff);
                        if exist([pth filesep nm '.mat'],'file')
                            delete([pth filesep nm '.mat']);
                        end
                    end
                    try
                        movefile(opts.diff.pdata.bvec,[opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_DC.bvec']);
                        movefile(opts.diff.pdata.bval,[opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_DC.bval']);
                    end
                    opts.diff.pdata.diff =  VEC(1).fname;
                    opts.diff.pdata.bval=[opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_DC.bval'];
                    opts.diff.pdata.bvec=[opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_DC.bvec'];
                    %-----------------------------------------------------------------%
                    %
                case 'fieldmap' % Field Map correction
            end
        end
        disp('                                                                   ');
        toc;
        %=====   6------ Creating Brain Mask in diffusion space  =============%
        % ---Moving Brain mask from T1 space to Diffusion space ----------%
        opts.diff.dbmask.filename = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_BDMask.nii'];
        try
            copyfile([opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep 'B0_temp_mask.nii'] , opts.diff.dbmask.filename);
        end
        if ~exist(opts.diff.dbmask.filename,'file')
            cad = ['WarpImageMultiTransform 3 ' opts.anat.ndata.bmask ' ' opts.diff.dbmask.filename ' -R ' opts.diff.pdata.b0 ' ' opts.transforms.t1_2_diff.affine ' --use-NN'];
            %system(cad);
        end
        Vtemp = spm_vol(opts.diff.dbmask.filename);
        Itemp = logical(spm_read_vols(Vtemp));
        spm_write_vol(Vtemp,Itemp);
        %-----------------------------------------------------------------%
        % To be included---- Erick's T1 aproximation
    end
    %
    %===========      Extracting B0    ===================================%
    disp(['Separating Corrected B0 Image ............................................']);tic;
    Vdif = spm_vol(opts.diff.pdata.diff);
    grads = load(opts.diff.pdata.bvec);
    if size(grads,2)>size(grads,1)
        grads = grads';
    end
    bvals = load(opts.diff.pdata.bval);
    if size(bvals,2)>size(bvals,1)
        bvals = bvals';
    end
    ind = find(sum([grads bvals]') ==0);
    Vb0s = Vdif(ind);
    Vb0 = Vdif(1);
    It = zeros(Vb0.dim(1:3));
    for k = 1:length(ind)
        It = It+spm_read_vols(Vb0s(ind(k)));
    end
    It = It/k;
    Vb0.fname = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_b0C.nii'];
    spm_write_vol(Vb0,It); clear It;
    opts.diff.pdata.b0 = [opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_b0C.nii'];
    toc;
    disp('                                                                    ');
    %=====================================================================%
    %================ END OF PREPROCESSING STEPS =========================%
else
    error('Some diffusion information is missing. Please enter the correct diffusion files. Diffusion Image, Gradient table, and B values');
    return;
end
% == Saving text file for TRACULA
[pth,nm,ext] = fileparts(opts.diff.pdata.diff);
fid  = fopen([pth filesep opts.pipe.subjId '.txt'],'wt');
fprintf(fid,'%s\n',['Diffusion_Image:=====>' opts.diff.pdata.diff]);
fprintf(fid,'%s\n',['Gradient_Table:=====>' opts.diff.pdata.bvec]);
fprintf(fid,'%s\n',['B-values_Table:=====>' opts.diff.pdata.bval]);
fprintf(fid,'%s\n',['Diffusion_B0:=====>' opts.diff.pdata.b0]);
fprintf(fid,'%s\n',['Diffusion_BrainMask:=====>' opts.diff.dbmask.filename]);
fclose(fid);
% == End of Saving text file for TRACULA               
% 
% === Deleting Temporal Files  ===========================================%
cad = ['rm ' opts.pipe.outdir filesep 'connectome' filesep opts.pipe.subjId filesep 'tmp' filesep '*'];
system(cad);
cad = ['rm -r ' opts.pipe.outdir filesep '5-freesurfer_processing' filesep opts.pipe.subjId filesep 'tmp' filesep '*'];
system(cad);
% === End of Deleting Temporal Files  ====================================%
return


%   =================== Internal scripts ==============================   %

function [Outboolean] =FreeSurfer_verif(InputDir,Id,cad);
%
% Syntax :
% [Outboolean] =FreeSurfer_verif(InputDir,Id,cad);
%
% This function computes  Brain masks for MRI data. It creates also the GM,
% WM and CSF segmentation files previously obtained by a segmentation process.
%
% Input Parameters:
%  InputDir       : Freesurfer output directory
%  Id             : Subject ID
%  cad            : Freesurfer stage (ie. autorecon1, autorecon2,
%  autorecon3);
%
% Output Parameter:
%   Outboolean    : Boolean variable to express if the freesurfer stage is
%   completed and all the results are saved inside the freesurfer
%   directory.
%
% See also: i_spm_segment
%__________________________________________________
% Authors: Yasser Aleman Gomez
%LIM
% February 1th 2012
% Version $1.0

%Autorecon1 No skull
filenamesa1 = strvcat([ 'mri' filesep 'orig' filesep '001.mgz'],[ 'mri' filesep 'rawavg.mgz'],[ 'mri' filesep 'orig.mgz'],[ 'mri' filesep 'nu.mgz'],[ 'mri' filesep 'T1.mgz'],[ 'mri' filesep 'transforms' filesep 'talairach.auto.xfm']);
%Autorecon1
filenamesa1 = strvcat([ 'mri' filesep 'orig' filesep '001.mgz'],[ 'mri' filesep 'rawavg.mgz'],[ 'mri' filesep 'orig.mgz'],[ 'mri' filesep 'nu.mgz'],[ 'mri' filesep 'T1.mgz'],[ 'mri' filesep 'brainmask.mgz'],[ 'mri' filesep 'brainmask.auto.mgz'],[ 'mri' filesep 'transforms' filesep 'talairach.auto.xfm']);

%Autorecon2
filenamesa2 = '';
filenamesa2 = strvcat(filenamesa2,[ 'mri' filesep 'transforms' filesep 'talairach.lta'],[ 'mri' filesep 'norm.mgz'],[ 'mri' filesep 'transforms' filesep 'talairach.m3z']);
filenamesa2 = strvcat(filenamesa2,[ 'mri' filesep 'transforms' filesep 'talairach.m3z.inv.x.mgz'],[ 'mri' filesep 'transforms' filesep 'talairach.m3z.inv.y.mgz'],[ 'mri' filesep 'transforms' filesep 'talairach.m3z.inv.z.mgz']);
filenamesa2 = strvcat(filenamesa2,[ 'mri' filesep 'nu_noneck.mgz'],[ 'mri' filesep 'transforms' filesep 'talairach_with_skull.lta'],[ 'mri' filesep 'aseg.auto_noCCseg.mgz'],[ 'mri' filesep 'aseg.auto.mgz'],[ 'mri' filesep 'aseg.mgz']);
filenamesa2 = strvcat(filenamesa2,[ 'mri' filesep 'brain.mgz'],[ 'mri' filesep 'brain.finalsurfs.mgz'],[ 'mri' filesep 'wm.seg.mgz'],[ 'mri' filesep 'wm.asegedit.mgz'],[ 'mri' filesep 'wm.mgz'],[ 'mri' filesep 'filled.mgz'],['scripts' filesep 'ponscc.cut.log']);
filenamesa2 = strvcat(filenamesa2,['surf' filesep 'lh.orig.nofix'],['surf' filesep 'lh.orig.nofix'],['surf' filesep 'lh.smoothwm.nofix']);
filenamesa2 = strvcat(filenamesa2,['surf' filesep 'rh.smoothwm.nofix'],['surf' filesep 'lh.inflated.nofix'],['surf' filesep 'rh.inflated.nofix'],['surf' filesep 'lh.qsphere.nofix'],['surf' filesep 'rh.qsphere.nofix']);
filenamesa2 = strvcat(filenamesa2,['surf' filesep 'lh.orig'],['surf' filesep 'rh.orig'],['surf' filesep 'lh.inflated'],['surf' filesep 'rh.inflated'],['surf' filesep 'lh.white'],['surf' filesep 'rh.white'],['surf' filesep 'lh.pial'],['surf' filesep 'rh.pial']);
filenamesa2 = strvcat(filenamesa2,['surf' filesep 'lh.thickness'],['surf' filesep 'rh.thickness'],['surf' filesep 'lh.curv'],['surf' filesep 'rh.curv'],['surf' filesep 'lh.area'],['surf' filesep 'rh.area'],['label' filesep 'lh.cortex.label'],['label' filesep 'rh.cortex.label']);
filenamesa2 = strvcat(filenamesa2,['surf' filesep 'lh.area.mid'],['surf' filesep 'rh.area.mid'],['surf' filesep 'lh.volume'],['surf' filesep 'rh.volume'],['surf' filesep 'lh.smoothwm'],['surf' filesep 'rh.smoothwm'],['surf' filesep 'lh.sulc'],['surf' filesep 'rh.sulc']);
filenamesa2 = strvcat(filenamesa2,['surf' filesep 'lh.inflated.H'],['surf' filesep 'rh.inflated.H'],['surf' filesep 'lh.inflated.K'],['surf' filesep 'rh.inflated.K']);


%Autorecon3
filenamesa3 = '';
filenamesa3 = strvcat(filenamesa3,['surf' filesep 'lh.sphere'],['surf' filesep 'rh.sphere'],['surf' filesep 'lh.sphere.reg'],['surf' filesep 'rh.sphere.reg'],['surf' filesep 'lh.jacobian_white'],['surf' filesep 'rh.jacobian_white'],['surf' filesep 'lh.avg_curv'],['surf' filesep 'rh.avg_curv']);
filenamesa3 = strvcat(filenamesa3,['label' filesep 'lh.aparc.annot'],['label' filesep 'rh.aparc.annot'],['stats' filesep 'lh.aparc.stats'],['stats' filesep 'rh.aparc.stats'],['label' filesep 'aparc.annot.ctab'],['label' filesep 'lh.aparc.a2009s.annot'],['label' filesep 'rh.aparc.a2009s.annot']);
filenamesa3 = strvcat(filenamesa3,['stats' filesep 'lh.aparc.a2009s.stats'],['stats' filesep 'rh.aparc.a2009s.stats'],['label' filesep 'aparc.annot.a2009s.ctab'],['stats' filesep 'aseg.stats'],[ 'mri' filesep 'lh.ribbon.mgz'],[ 'mri' filesep 'rh.ribbon.mgz'],[ 'mri' filesep 'aparc+aseg.mgz'],[ 'mri' filesep 'aparc.a2009s+aseg.mgz']);
filenamesa3 = strvcat(filenamesa3,[ 'mri' filesep 'wmparc.mgz'],['stats' filesep 'wmparc.stats']);

%TRACULA
%filentrabedpostx

if strcmp(cad,'autorecon1')
    filenames = filenamesa1;
elseif strcmp(cad,'autorecon2')
    filenames = strvcat(filenamesa1,filenamesa2);
elseif strcmp(cad,'autorecon3')
    filenames = strvcat(filenamesa1,filenamesa2,filenamesa3);
elseif strcmp(cad,'autorecon23')
    filenames = strvcat(filenamesa2,filenamesa3);
elseif strcmp(cad,'all')
    filenames = strvcat(filenamesa1,filenamesa2,filenamesa3);
end
A = zeros(1,size(filenames,1));
names = [repmat([InputDir filesep Id filesep],[size(filenames,1) 1]) filenames];
for j = 1:size(filenames,1)
    
    ind=  exist(deblank(names(j,:)));
    if ind>0
        A(j) = 1;
    end
end
Outboolean = prod(A);
return;

function outnames = change_spacen(Ii0,If0,order,Outdir);
if size(If0,1) == 1, 
    If0 = If0(ones(size(Ii0,1),1),:);
elseif ~(size(Ii0,1)==size(If0,1)), 
    error('it must be the same number of images')
end
outnames = '';
for image = 1:size(Ii0),
    Ii = deblank(Ii0(image,:));
    If = deblank(If0(image,:));
    Vi = spm_vol(Ii); n = length(Vi); name = Vi.fname;
    Vf = spm_vol(If);
    [a1,a2,a3] = fileparts(Ii); [b1,b2,b3] = fileparts(If);
    counter = 0; 
    [X,Y] = ndgrid(1:Vf.dim(1),1:Vf.dim(2));
    for i = 1:n,
        Vn = Vf; c = spm_bsplinc(Vi(i),[order*ones(1,3) 0 0 0]);
        Vn.pinfo = Vi(i).pinfo; Vn.dt = Vi(i).dt; Vn.n = Vi(i).n;
        [pathstr,name,ext] = fileparts(Vi(i).fname);
        if nargin<4
            Vn.fname = fullfile(pathstr,['s' name ext]);
        elseif nargin == 4
            Vn.fname = fullfile(Outdir,['s' name ext]);
        end
        try, Vn.userdata.g = Vf.userdata.g; Vn.userdata.mat = Vn.mat; end
        Vn = spm_create_vol(Vn); M = Vi(i).mat\Vn.mat;
        for z = 1:Vn.dim(3),
            counter = counter+1;
            A = spm_bsplins(c,M(1,1)*X + M(1,2)*Y + M(1,3)*z + M(1,4),...
                M(2,1)*X + M(2,2)*Y + M(2,3)*z + M(2,4),...
                M(3,1)*X + M(3,2)*Y + M(3,3)*z + M(3,4),...
                [order*ones(1,3) 0 0 0]);
            spm_write_plane(Vn,A,z);
        end
        outnames = strvcat(outnames,Vn.fname);
    end
end
outnames = unique(outnames,'rows');
fclose all;
return;

function OutFilenames = correct_diff_images(diff_image_name, bvecs_name, bvals_name, OutFilename);
%
% Syntax :
% OutFilename = correct_diff_images(diff_image_name, bvecs_name, bvals_name, OutFilename);
%
% This script performs slice-wise intensity related artifacts checking
% as recommened by Zhexing Liua et al. "Quality Control of Diffusion Weighted Images"
%
% Input Parameters:
%   diff_image_name        : Diffusion weighted image
%   bvecs_name             : Gradient table
%   bvals_name             : B values
%   OutFilename            : Corrected diffusion weighted image.
%
% Output Parameters:
%   OutFilenames           : Corrected diffusion weighted image, Corrected
%                            gradient table and corrected b-values table
%
% See also: 
%__________________________________________________
% Authors: Yasser Iturria-Medina and Erick J Canales-Rodriguez
% LIM, HUGGM
% March 22th 2012
% Version $1.0
spms(8)
diff_gradients = load(bvecs_name);
bvalues        = load(bvals_name);

if size(diff_gradients,1) > size(diff_gradients,2)
    diff_gradients = diff_gradients';
end
if size(bvalues,1) > size(bvalues,2)
    bvalues = bvalues';
end

Vdiff = spm_vol(diff_image_name);
%figure; title('Correlation vs slice number');
for gradient = 1:size(diff_gradients,2)
    disp(['Testing gradient ' num2str(gradient)])
    volumen_data = spm_read_vols(Vdiff(gradient));
    for slice = 1:Vdiff(gradient).dim(3)-1
        vector_1 = volumen_data(:,:,slice); vector_1 = vector_1(:);
        vector_2 = volumen_data(:,:,slice+1); vector_2 = vector_2(:);
        % correlation_values(gradient,slice) = (vector_1'*vector_2)/sqrt(sum(vector_1)*sum(vector_2)); % normalized correlation
        correlation_values(gradient,slice) = corr(vector_1,vector_2);
    end
end
ind_Nob0s = find(bvalues ~= 0);
NC_mean = mean(correlation_values(ind_Nob0s,:),1);
NC_std = std(correlation_values(ind_Nob0s,:),0,1);

alfa = 3.5; % The default is a number between 3 and 4

indices_bad_volumens = [];
% neglecting the three first and last slices: not contain brain information
num_slices_neglect = 3;
for slice = (1 + num_slices_neglect):(Vdiff(gradient).dim(3) - 1 - num_slices_neglect)
    indices_bad_volumens = [indices_bad_volumens; find(correlation_values(:,slice) < NC_mean(slice)-alfa*NC_std(slice) & bvalues' ~= 0)];
end
indices_bad_volumens = unique(indices_bad_volumens);
disp(['Volumens with intensity artifacts -> ' num2str(indices_bad_volumens')])

if ~isempty(indices_bad_volumens)
    % creating corrected image volume
    New_Number_gradients = size(diff_gradients,2)-length(indices_bad_volumens);
    for i = 1:New_Number_gradients
        V(i) = Vdiff(1);
        % V(i).fname = fullfile(path_m,[nam_m '_intensity_corrected' ext_m]);
        V(i).fname = OutFilename;
        V(i).pinfo = [1 0 0]';
        V(i).descrip = [V(i).descrip '===> Proc: Bad Volumes Correction('  datestr(now) ')' ];
        if strcmp(spm('ver'),'SPM2')
            V(i).dim(4) = spm_type('float');
            V(i).n = i;
        elseif strcmp(spm('ver'),'SPM5') || strcmp(spm('ver'),'SPM8')
            V(i).dt = [16 0];
            V(i).n = [i 1];
        end
    end
    V = spm_create_vol(V);
    indices_god_volumens = setdiff(1:size(diff_gradients,2),indices_bad_volumens');
    % writting new gradient volumes
    I = spm_read_vols(Vdiff);I = I(:,:,:,indices_god_volumens);
    for j =1:size(I,4);
        spm_write_vol(V(j),squeeze(I(:,:,:,j)));
    end
    new_gradients = diff_gradients(:,indices_god_volumens);
    new_bvalues = bvalues(indices_god_volumens);
    [pthf,nmf,ext] = fileparts(OutFilename);
    newbvec = [pthf filesep nmf '.bvec'];
    newbval = [pthf filesep nmf '.bval'];
    newl = [pthf filesep nmf '_Bad_volumes.txt'];
    if sum(new_bvalues(:) ~= 0)&sum(new_gradients(:) ~= 0)
        %new_gradients = new_gradients';
        save(newbvec,'new_gradients','-ascii');
        save(newbval,'new_bvalues','-ascii');
    end
    save(newl,'indices_bad_volumens','-ascii');
    status = fclose('all');
    close all;
    OutFilenames = strvcat(OutFilename,newbvec,newbval);
else
    OutFilenames = strvcat(diff_image_name,bvecs_name,bvals_name);
end
spms(12)
return;

function outputs = flip_to_axial_nii(filenames,opt)

%FLIP_TO_AXIAL flips sagittal or coronal images to axial view.
%
% USAGE: FLIP_TO_AXIAL_NII(FILENAME,OPT)
% 
% FILENAME: name of the file to be flipped.
%
% OPT: 1 erases original image.
%
% Note: Coordinate system where image is embbeded has to be with z axis 
% toward the top of the head. This information has to be in the corresponding 
% mat file. In case of no mat file present image must be flipped before
% (e.g. using MRIcro). Using an appriopiate Dicom to Analyze converter a mat file
% with the correct orientation of the image will be quarantized.
% 
% Author: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% January 13 2006
% Version $4.0

if ~nargin || isempty(filenames),
    filenames = spm_get(Inf,'*.nii','Select images to flip');
end
V = spm_vol(filenames);

for i = 1:length(V),
    [pathstr,name] = fileparts(deblank(filenames(i,:)));
    Vn = V(i); 
    Vn.fname = fullfile(pathstr,['f' name '.nii']);
    R = getrot(V(i).mat(1:3,1:3));
    if abs(R) == eye(3), outputs{i} = deblank(filenames(i,:)); continue, end %#ok
    ind = find(R); if length(ind) ~= 3 && abs(det(R)) == 1,
        error('Imposible to get axial plane normal direction:'); end
    [r,c] = find(R<0); R = [R zeros(3,1);zeros(1,3) 1];
    for j = 1:length(c), R(r(j),4) = V(i).dim(c(j))+1; end
    Vn.dim(1:3) = abs(R(1:3,1:3)*V(i).dim(1:3)')';
    Vn.mat = V(i).mat*inv(R); if det(Vn.mat)>0,
        Vn.mat = Vn.mat*inv([diag([-1 1 1]) (Vn.dim(1)+1)*eye(3,1); zeros(1,3) 1]); end
    Vn = spm_create_vol(Vn);
    for z = 1:Vn.dim(3),
        A = spm_slice_vol(V(i),V(i).mat\Vn.mat*spm_matrix([0 0 z]),Vn.dim(1:2),0);
        spm_write_plane(Vn,A,z);
    end
    fclose all;
    outputs{i} = fullfile(pathstr,['f' name '.nii']); %#ok
    if nargin && opt,
        delete(fullfile(pathstr,[name '.nii']));
        movefile(fullfile(pathstr,['f' name '.nii']),fullfile(pathstr,[name '.nii']));
        outputs{i} = fullfile(pathstr,[name '.nii']); %#ok
    end
    try %#ok
        usrdata = dti_get_dtidata(outputs{i});
        usrdata.mat = Vn.mat;
        dti_get_dtidata(outputs{i},usrdata);
    end
end
outputs = strvcat(outputs{:}); %#ok
return
% ----------------------------------------------------------
function R = getrot(A)

M = diag(sqrt(sum(A.^2))); A = A*inv(M); [s,ind] = max(abs(A)); %#ok
R(3,3) = 0; for i = 1:3, R(ind(i),i) = sign(A(ind(i),i)); end
return
% ----------------------------------------------------------
function [opts] = T1_Brainmask_Computation(opts)

Outfiles = Extract_brain([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1.nii'],opts.anat.ndata.t2); % NEEDS spm12
system(['gzip ' opts.anat.ndata.t2]);
Vbm = spm_vol(deblank(Outfiles(1,:))); Vi = spm_vol([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1.nii']); %% Masking T1 to create brainmask.auto.mgz
Im = logical(spm_read_vols(Vbm));      Ii = spm_read_vols(Vi);
Ii = Im.*Ii; Vi.fname = [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1_Masked.nii']; spm_write_vol(Vi,Ii); clear Ii Im;
cad = ['mri_convert -c -i ' Vi.fname ' -o ' [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1_Masked.mgz']];   % Converting T1 Masked to MGZ
system(cad);
movefile([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'brainmask.mgz'],     [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'brainmask.FS.mgz']); % Moving freesurfer skullstripped T1
movefile([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'brainmask.auto.mgz'],[opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'brainmask.auto.FS.mgz']);
movefile([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1_Masked.mgz'],     [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'brainmask.mgz']); % Moving P0brains to skullstripped T1
try
    copyfile([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'brainmask.mgz'],     [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'brainmask.auto.mgz']);
end
delete([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1_Masked.nii']);
delete([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1_Masked.mgz']);
movefile(deblank(Outfiles(1,:)),[opts.pipe.outdir filesep opts.pipe.subjId filesep 'mri' filesep '_p0.nii']);
movefile(deblank(Outfiles(2,:)),[opts.pipe.outdir filesep opts.pipe.subjId filesep 'mri' filesep '_p1.nii']);
movefile(deblank(Outfiles(3,:)),[opts.pipe.outdir filesep opts.pipe.subjId filesep 'mri' filesep '_p2.nii']);
movefile(deblank(Outfiles(4,:)),[opts.pipe.outdir filesep opts.pipe.subjId filesep 'mri' filesep '_p3.nii']);
movefile(deblank(Outfiles(5,:)),[opts.pipe.outdir filesep opts.pipe.subjId filesep 'mri']);
movefile(deblank(Outfiles(6,:)),[opts.pipe.outdir filesep opts.pipe.subjId filesep 'mri']);
opts.anat.ndata.p0bmask = [opts.pipe.outdir filesep opts.pipe.subjId filesep 'mri' filesep '_p0.nii']; % BrainMask
opts.anat.ndata.gm = [opts.pipe.outdir filesep opts.pipe.subjId filesep 'mri' filesep '_p1.nii'];
opts.anat.ndata.wm = [opts.pipe.outdir filesep opts.pipe.subjId filesep 'mri' filesep '_p2.nii'];
opts.anat.ndata.csf = [opts.pipe.outdir filesep opts.pipe.subjId filesep 'mri' filesep '_p3.nii'];
return
% ----------------------------------------------------------
function Outfiles = Extract_brain(T1Images,T2Images);
%
% Syntax :
% Outfiles = Extract_brain(T1Images,T2Images);
%
% This function computes  Brain masks for MRI data. It creates also the GM,
% WM and CSF segmentation files previously obtained by a segmentation process.
%
% Input Parameters:
%  T1Images       : T1 Weighted Image.
%  T2Images       : T2 Weighted Image.
%
%
% Output Parameter:
%   Outfiles       : Brain mask, GM, WM and CSF segmentation files
%
% See also: i_spm_segment
%__________________________________________________
% Authors: Yasser Aleman Gomez
%LIM
% February 1th 2012
% Version $1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tempd = which('spm_preproc8.m');
Outfiles = '';
if ~isempty(tempd)
    V = spm_vol(T1Images);
    Vs = spm_vol(T2Images);
    Ns = size(V,1);
    for i = 1:Ns
        T1 = deblank(T1Images(i,:));
        %outnames = change_spacen(Vs(i).fname,V(i).fname,7);Vt = spm_vol(outnames);
        [path,name,ext] = fileparts(T1);
        job.data{1} = T1;
        job.channel(1).vols{1} = V.fname;
        job.channel(1).biasreg = 1.0000e-04;
        job.channel(1).biasfwhm = 60;
        job.channel(1).write = [0 0];
        
%         job.channel(2).vols{1} = Vs.fname;
%         job.channel(2).biasreg = 1.0000e-04;
%         job.channel(2).biasfwhm = 60;
%         job.channel(2).write = [0 0];
        
        job.opts.tpm = {[which('TPM.nii')]};
        job.opts.ngaus = [2 2 2 3 4 2];
        job.opts.biasreg = 1.0000e-04;
        job.opts.biasfwhm = 60;
        job.opts.affreg = 'mni';
        job.opts.warpreg = 4;
        job.opts.samp = 3;
        job.extopts.dartelwarp = 1;
        job.extopts.sanlm = 1;
        job.extopts.mrf = 0.1500;
        job.extopts.cleanup = 1;
        job.extopts.print = 1;
        job.output.GM.native = 1;
        job.output.GM.warped = 0;
        job.output.GM.modulated = 0;
        job.output.GM.dartel = 0;
        job.output.WM.native = 1;
        job.output.WM.warped = 0;
        job.output.WM.modulated = 0;
        job.output.WM.dartel = 0;
        job.output.CSF.native = 1;
        job.output.CSF.warped = 0;
        job.output.CSF.modulated = 0;
        job.output.CSF.dartel = 0;
        job.output.bias.native = 0;
        job.output.bias.warped = 0;
        job.output.bias.affine = 0;
        job.output.label.native = 1;
        job.output.label.warped = 0;
        job.output.label.dartel = 0;
        job.output.jacobian.warped = 0;
        job.output.warps = [0 0];
        job.bias = [0 0 0];
        job.label = [1 0 0 0];
        job.jacobian = 0;
        job.biasreg = 1.0000e-04;
        job.biasfwhm = 60;
        job.warp.affreg = 'mni';
        job.warp.samp = 3;
        job.warp.reg = 4;
        job.warp.write = [0 0];
        job.warp.sanlm = 1;
        job.warp.mrf = 0.1500;
        job.warp.print = 1;
        job.warp.cleanup = 1;
        job.warp.dartelwarp = 1;
        job.warps = [0 0];
        native = [1 0 0;1 0 0;1 0 0;0 0 0;0 0 0;0 0 0];
        for j = 1:6
            job.tissue(j).tpm = {[which('TPM.nii') ',' num2str(i)]};
            job.tissue(j).ngaus = job.opts.ngaus(j);
            job.tissue(j).native = native(j,:);
            job.tissue(j).warped = [0 0 0];
        end
        cg_vbm8_run(job);
        Outfiles = strvcat(Outfiles,[repmat([path filesep],[4 1]) strvcat('p0','p1','p2','p3') repmat([name ext(1:4)],[4 1])]);
        Outfiles = strvcat(Outfiles,[path filesep name '_seg8.mat'],[path filesep 'p' name '_seg8.txt']);
%         v1 = spm_vol([path filesep 'T1.nii']);
%         v2 = spm_vol([path filesep 'p0T1.nii']);
%         v1_1 = spm_read_vols(v1);
%         v2_1 = spm_read_vols(v2);
%         v2_1 = logical(v2_1).*v1_1;
%         spm_write_vol(v1, v2_1);
    end
else
    global defaults
    spm_defaults;
    dseg				= defaults.segment;
    V = spm_vol(T1Images);Vs = spm_vol(T2Images);
    Ns = size(V,1);
    VG0 = which('T1.nii');
    for i =1:Ns
        disp(['Extracting    ' num2str(i) ' off ' num2str(Ns) ': ======>>   ' deblank(T1Images(i,:))]);
        outnames = change_spacen(Vs(i).fname,V(i).fname,7);Vt = spm_vol(outnames);
        i_spm_segment(V(i),VG0,dseg,Vt);
        [pth,nm,ext] = fileparts(V(i).fname);
        [pth1,nm1,ext1] = fileparts(Vt.fname);
        GMImages = [pth filesep 'p1' nm ext(1:4)];
        WMImages = [pth filesep 'p2' nm ext(1:4)];
        CSFImages = [pth filesep 'p3' nm ext(1:4)];
        Create_Brain_Masks(GMImages, WMImages, CSFImages);
        delete(outnames);
        %delete([pth filesep 'p1' nm ext(1:4)]);
        %delete([pth filesep 'p2' nm ext(1:4)]);
        %delete([pth filesep 'p3' nm ext(1:4)]);
        delete([pth filesep 'p4' nm ext(1:4)]);
        delete([pth filesep 'p5' nm ext(1:4)]);
        delete([pth filesep 'p6' nm ext(1:4)]);
        [pth1, nm1,ext1] = fileparts(GMImages);
        Outfiles = strvcat(Outfiles,[pth1 filesep nm1(3:end) '_BM' ext1(1:4)],GMImages,WMImages,CSFImages,[pth filesep nm '_seg8.mat'],[pth1 filesep 'y_' nm ext(1:4)]);
        %delete([pth filesep nm '_seg8.mat']);
        %delete([pth1 filesep 'y_' nm1 ext1(1:4)]);
    end
end
spms(12)
return
% ----------------------------------------------------------
function [opts] = Reading_Pipe_Configuration_File(Conf_File);
%
% Syntax :
% [opts] = Reading_Pipe_Configuration_File(Conf_File);
%
% Reading configuration files
%
% Input Parameters:
%   Conf_File         :  Configuration file
%
% Output Parameters:
%
%       opts          : Pipeline options
%
% Related references:
%
%
% See also: 
% 
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% April 27th 2012
% Version $1.0

warning off;

%% ====================== Reading Configuration File =====================%
fio = fopen(Conf_File,'rt');lines = '';cont = 0;
conts = 0;
while 1
    cont = cont + 1;
    line = fgetl(fio);
    if ~ischar(line),   break,   end
    line = [deblank(line) ';'];
    eval(line);
end
fclose(fio);
return;