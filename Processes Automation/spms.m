function spms(spmver)
%SPM switcher
%You can select to run eithr SPM8 or SPM12
%Usage; spms will bring up a dialogue to select which SPM you want to run
% spms(8) will run SPM8 and spms(12) will run SPM12 directory.

%%Specification of the parent directory of spm
%Please specify the location of your spm installation.
%This script assumes that your spm12 or spm8 is installed under this
%directory.
%For example
%Windows: c:\spm
%MacOSX: /Users/your_username/spm
%Linux: /home/your_username/spm
spm_path='/Users/laimbio/Documents/MATLAB';

if nargin == 0
    SelectVer = questdlg('Which SPM do you want to run?',...
        'Select SPM ver.','spm12','spm8','spm12');
else
    SelectVer = ['spm' num2str(spmver)];
end

%remove spm paths
while true
    try
        evalc('spm_rmpath;');
        catch break;
    end
end

%add paths for your selected spm
path=genpath(fullfile(spm_path,SelectVer));
addpath(path);
%clear
clear classes

% run spm
% spm;
disp('done')        
return;
