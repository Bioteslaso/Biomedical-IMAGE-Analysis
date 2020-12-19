function [VBMrel, VBMver]  = cg_vbm8_version
% check for VBM revision
%
% FORMAT [VBMrel, VBMver] = cg_vbm8_version
% 
% This function will retrieve the VBM release and version and is a
% modified version of spm_version.m
%_______________________________________________________________________
% Christian Gaser
% $Id: cg_vbm8_version.m 445 2015-12-17 13:26:56Z gaser $

persistent VBM_VER;
v = VBM_VER;
if isempty(VBM_VER)
    v = struct('Name','','Version','','Release','','Date','');
    try
        vfile = fullfile(spm('Dir'),'toolbox','vbm8','Contents.txt');
        fid = fopen(vfile,'rt');
        if fid == -1, error('Can''t open %s.',vfile); end
        l1 = fgetl(fid); l2 = fgetl(fid);
        fclose(fid);
        l1 = strtrim(l1(2:end)); l2 = strtrim(l2(2:end));
        t  = textscan(l2,'%s','delimiter',' '); t = t{1};
        v.Name = l1; v.Date = t{4};
        v.Version = t{2}; v.Release = t{3}(2:end-1);
    catch
        error('Can''t obtain VBM Revision information.');
    end
    VBM_VER = v;
end
VBMrel = v.Release;
VBMver = v.Version;
