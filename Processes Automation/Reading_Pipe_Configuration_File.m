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
    cont = cont + 1;/Users/laimbio/Documents/INVESTIGACION/AGES_II/agesII_DTI
    line = fgetl(fio);
    if ~ischar(line),   break,   end
    line = [deblank(line) ';'];
    eval(line);
end
fclose(fio);
return;
