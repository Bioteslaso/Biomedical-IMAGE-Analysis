function []=fieldmap(DirIn,mainDCM2NII)

%% Inicializacion del entorno de spm

%spm('defaults','fmri');
%spm_jobman('initcfg')

%% Seleccion de la ruta a procesar

% path_inicio=pwd;

%Elegimos directorio donde estan las imagenes
path_imagenes=DirIn;
%path_imagenes_fieldmap = [path_imagenes '/fieldmap.nii'];

DCM2NII=strcat(mainDCM2NII,filesep,'dcm2nii');

temp = [path_imagenes '/*.dcm'];
imagenes_fieldmap = dir(temp);

clear temp
%% Extraemos informacion de las cabeceras

path_TE1 = [path_imagenes '/TE1'];
path_TE2 = [path_imagenes '/TE2'];
mkdir(path_TE1)
mkdir(path_TE2)

path_TE1_real = [path_TE1 '/TE1_r'];
path_TE1_imaginaria = [path_TE1 '/TE1_i'];
path_TE1_magnitud = [path_TE1 '/TE1_magnitud'];

path_TE2_real = [path_TE2 '/TE2_r'];
path_TE2_imaginaria = [path_TE2 '/TE2_i'];
path_TE2_magnitud = [path_TE2 '/TE2_magnitud'];

mkdir(path_TE1_real)
mkdir(path_TE1_imaginaria)
mkdir(path_TE1_magnitud)

mkdir(path_TE2_real)
mkdir(path_TE2_imaginaria)
mkdir(path_TE2_magnitud)


temp = [path_imagenes '/' imagenes_fieldmap(1).name];   
info = dicominfo(temp);    
TE1 = info.EchoTime;  
fprintf('TE1 = %s', TE1);
%cmd = ['cp '  temp ' ' path_TE1]
%system(cmd);

cmd = ['cp '  temp ' ' path_TE1_magnitud];
system(cmd);

for i=2:size(imagenes_fieldmap,1)
    
    
    temp = [path_imagenes '/' imagenes_fieldmap(i).name];
    info = dicominfo(temp);
    TE = info.EchoTime;
   
    
    if (i==2)
        contador_real = i;
        
    elseif (i==3)
        contador_imaginaria=i;
    end
    
    if (TE==TE1)
        %cmd = ['cp '  temp ' ' path_TE1];
        %system(cmd);
        
        if (contador_real==i)
            
            cmd = ['cp '  temp ' ' path_TE1_real];
            system(cmd);
            contador_real=contador_real+3;
        elseif contador_imaginaria==i
            
            cmd = ['cp '  temp ' ' path_TE1_imaginaria];
            system(cmd);
            contador_imaginaria=contador_imaginaria+3;
            
        else
            
          cmd = ['cp '  temp ' ' path_TE1_magnitud];
          system(cmd);  
        end
        
        
    else
        TE2=TE;
        %cmd = ['cp '  temp ' ' path_TE2];
        %system(cmd);
        if (contador_real==i)
            
            cmd = ['cp '  temp ' ' path_TE2_real];
            system(cmd);
            contador_real=contador_real+3;
        elseif contador_imaginaria==i
            
            cmd = ['cp '  temp ' ' path_TE2_imaginaria];
            system(cmd);
            contador_imaginaria=contador_imaginaria+3;
            
        else
            
          cmd = ['cp '  temp ' ' path_TE2_magnitud];
          system(cmd);  
        
        end
        
    end
end

fprintf('TE2 = %s', TE2)

%% Convertimos las imagenes de DICOM a NIFTI

imagenes_list = dir([path_TE2_real '/*.dcm']);
cmd = strcat(DCM2NII,{' -a y -d n -e n -f y -g y -i n -g n -r n -x n -o '}, path_TE2_real, {' '}, path_TE2_real, '/', imagenes_list(1).name);
system(cmd{1})

imagenes_list = dir([path_TE2_imaginaria '/*.dcm']);
cmd =strcat(DCM2NII,{' -a y -d n -e n -f y -g y -i n -g n -r n -x n -o '}, path_TE2_imaginaria, {' '},  path_TE2_imaginaria, '/', imagenes_list(1).name);
system(cmd{1})

imagenes_list = dir([path_TE1_real '/*.dcm']);
cmd = strcat(DCM2NII,{' -a y -d n -e n -f y -g y -i n -g n -r n -x n -o '}, path_TE1_real, {' '}, path_TE1_real, '/', imagenes_list(1).name);
system(cmd{1})

imagenes_list = dir([path_TE1_imaginaria '/*.dcm']);
cmd = strcat(DCM2NII,{' -a y -d n -e n -f y -g y -i n -g n -r n -x n -o '}, path_TE1_imaginaria, {' '}, path_TE1_imaginaria, '/', imagenes_list(1).name);
system(cmd{1})



















