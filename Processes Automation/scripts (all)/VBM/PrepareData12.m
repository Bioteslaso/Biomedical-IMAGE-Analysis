function datafiles = PrepareData12(fldr,outdir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PrepareData
%--------------------------------------------------------------------------
% 
% This module take a folder with subfolders where there are dicom files and  
% convert those to a NIfTI file(.nii), coregister and put in 'outdir'. It
% also produces a txt file witht all the information of the data: file
% path, age and gender. 
%
%__________________________________________________________________________
% INPUTS:
%  - fldr: path to folder with dicoms subfolders.
%  - outdir: output directory.
% OUTPUTS:
%  - datafiles: vector with full paths of the .nii files, age and sex.
%
%__________________________________________________________________________
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ref = '/Users/lauranunez/Documents/MATLAB/spm12/canonical/avg152T1.nii';

fls = dir([fldr '/*']);
subjects = {fls.name};
inx = strncmpi('.',subjects,1);
SubID=subjects(~inx);
subjects = strcat(fldr,'/',subjects(~inx));

sub.file ='';
sub.age = 0;
sub.sex = 0;
datafiles = [sub];

Folders=strsplit(fldr,'/');
File=strcat('/',Folders{length(Folders)-1},'_',Folders{length(Folders)},'_datafiles.txt');
fid=fopen(strcat(outdir,File),'w');

for i=1:length(subjects)
    [niifile, age, sex] = dcm2nii12(subjects{i},subjects{i});
    auto_reorient(niifile)
    Coregister12(niifile, ref)
    dirvol = dir(niifile);
    name = dirvol.name;
    movefile(niifile,strcat(outdir,'/',SubID{i},'.nii'));
    movefile([subjects{i} '/niifile/r' name],strcat(outdir,'/r',SubID{i},'.nii'));
    rmdir([subjects{i} '/niifile'], 's');
    sub.file = [outdir '/r' SubID{i} '.nii'];
    sub.age = age;
    sub.sex = sex;
    datafiles(i) = sub;
    if sub.sex
        Sex='M';
    else
        Sex='F';
    end;
    fprintf(fid,'%s %d %s\n',sub.file,sub.age,Sex); 
end;
fclose(fid);
