clear all
clc
while 1 
        try
        anw= menu('Choose function to run:','Sort','DCM2NII','Fieldmap', 'Laso action','EXIT');

        if anw==5
            msgbox('See you soon.','BYE')
            break
        elseif anw==1
            direc= string(inputdlg('Enter directory where files are stored.','Directory'));        
            strcat(set_BIDS(char(direc(1))))
            msgbox('Images were converted to .nii succesfully.','DONE')
        elseif anw==2
            direc1= string(inputdlg('Enter directory where files are mow stored.','Directory'));   
            direc2= string(inputdlg('Enter directory where files are to be saved.','Directory'));  
            direc3= string(inputdlg('Enter directory where dcm2niix is found.','Directory'));  
            strcat(BIDS_dcm2nii(char(direc1(1)),char(direc2(1)),char(direc3(1))))
            msgbox('Images were sorted succesfully.','DONE')
        elseif anw==3
            direc= string(inputdlg('Enter directory where files are stored.','Directory'));        
            strcat(SP_fieldmap(char(direc(1))))
        elseif anw==4
            msgbox('All dicom images were deleted successfully.','Congratulations')
        else
            msgbox('Invalid option. Try again.', 'Error');
        end
        catch
            msgbox('An unpredicted error ocurred. Please, restart the Menu.','ERROR unkwon')
        end
end
   