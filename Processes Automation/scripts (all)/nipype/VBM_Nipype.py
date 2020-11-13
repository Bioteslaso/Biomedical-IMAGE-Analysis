# -*- coding: utf-8 -*-
"""
@author: lauranunez
@eddited by: Mario Gil

"""
#Import of needed interfaces---------------------------------------------------

import nipype.interfaces.io as nio           # nipype i/o routines
import nipype.interfaces.spm as spm          # spm
import nipype.pipeline.engine as pe          # pypeline engine
import nipype.interfaces.utility as util     # utility
import normalise as norm
import glob
import os
import existingTemplate as et
import groupanalysis as stst
import nipype.interfaces.matlab as mlab      # how to run matlab
import DARTELinterface as dartelI

#Definition of needed functions------------------------------------------------
def getRCclasses(rc_folder):
    rc1_files=glob.glob(rc_folder+'/rc1*.nii')
    rc2_files=glob.glob(rc_folder+'/rc2*.nii')
    rc3_files=glob.glob(rc_folder+'/rc3*.nii')

    return [[rc1_files,rc2_files], [rc3_files]]
   
def getIterationParameters(templates):
    iterationParameters = []
    iterationParameters.append((3, (4, 2, 0.000001), 1, templates[0]))
    iterationParameters.append((3, (2, 1, 0.000001), 1, templates[1]))
    iterationParameters.append((3, (1, 0.5, 0.000001), 2, templates[2]))
    iterationParameters.append((3, (0.5, 0.25, 0.000001), 4, templates[3]))
    iterationParameters.append((3, (0.25, 0.125, 0.000001), 16, templates[4]))
    iterationParameters.append((3, (0.25, 0.125, 0.000001), 64, templates[5]))

    return iterationParameters
    

matlab_cmd = '/Volumes/PassportMarioGilCorrea/Programs/SPM/spm8_runtime/run_spm8.sh /Applications/MATLAB/MATLAB_Compiler_Runtime/v84/ script'
spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)


#Import of rc_files------------------------------------------------------------

rc_folderG1 = '/Volumes/PassportMarioGilCorrea/Programs/ScriptsVBM/Segmentation_Results/rc_filesG1' #Variable para la ruta al directorio
rc_folderG2 = '/Volumes/PassportMarioGilCorrea/Programs/ScriptsVBM/Segmentation_Results/rc_filesG2'
[rc1_2_G1, rc3G1]=getRCclasses(rc_folderG1)
[rc1_2_G2, rc3G2]=getRCclasses(rc_folderG2)
rc1_images = rc1_2_G1[0]
rc1_images.extend (rc1_2_G2[0])
rc2_images = rc1_2_G1[1]
rc2_images.extend (rc1_2_G2[1])
rc1_2_images = [rc1_images,rc2_images]
rc3_images = rc3G1[0]
rc3_images.extend(rc3G2[0])
rc3_images=[rc3_images]
#Build different nodes up------------------------------------------------------ 

dartelGMWM = pe.Node(interface=dartelI.DARTEL(), name="dartelGMWM")
dartelGMWM.inputs.iteration_parameters = [(3, (4, 2, 0.000001), 1, 16) , (3, (2, 1, 0.000001), 1, 8), (3, (1, 0.5, 0.000001), 2, 4), (3, (0.5, 0.25, 0.000001), 4, 2), (3, (0.25, 0.125, 0.000001), 16, 1), (3, (0.25, 0.125, 0.000001), 64, 0.5) ]
dartelGMWM.inputs.optimization_parameters = (0.01, 3, 3)
dartelGMWM.inputs.regularization_form = 'Linear'
dartelGMWM.inputs.template_prefix = 'Template_GM_WM'
dartelGMWM.inputs.image_files = rc1_2_images

dartelCSF = pe.Node(interface=dartelI.DARTEL(), name="dartelCSF")
dartelCSF.inputs.iteration_parameters = [(3, (4, 2, 0.000001), 1, 16) , (3, (2, 1, 0.000001), 1, 8), (3, (1, 0.5, 0.000001), 2, 4), (3, (0.5, 0.25, 0.000001), 4, 2), (3, (0.25, 0.125, 0.000001), 16, 1), (3, (0.25, 0.125, 0.000001), 64, 0.5) ]
dartelCSF.inputs.optimization_parameters = (0.01, 3, 3)
dartelCSF.inputs.regularization_form = 'Linear'
dartelCSF.inputs.template_prefix = 'Template_CSF'
dartelCSF.inputs.image_files = rc3_images

datasink = pe.Node(interface=nio.DataSink(), name="datasink")
datasink.inputs.base_directory = os.path.abspath('./Results')

darteltemplateGMWM = pe.Node(interface=et.DARTELExisting(), name="darteltemplateGMWM")
darteltemplateGMWM.inputs.optimization_parameters = (0.01, 3, 3)
darteltemplateGMWM.inputs.regularization_form = 'Linear'
darteltemplateGMWM.inputs.image_files = rc1_2_images

darteltemplateCSF = pe.Node(interface=et.DARTELExisting(), name="darteltemplateCSF")
darteltemplateCSF.inputs.optimization_parameters = (0.01, 3, 3)
darteltemplateCSF.inputs.regularization_form = 'Linear'
darteltemplateCSF.inputs.image_files = rc3_images
          
nmGMWM = pe.Node(interface=norm.DARTELNorm2MNI(), name="nmGMWM")
nmGMWM.inputs.modulate = True
nmGMWM.inputs.fwhm = [8, 8, 8]
nmGMWM.inputs.voxel_size = (1.5, 1.5, 1.5)
nmGMWM.inputs.apply_to_files = rc1_2_images
    
nmCSF = pe.Node(interface=norm.DARTELNorm2MNI(), name="nmCSF")
nmCSF.inputs.modulate = True
nmCSF.inputs.fwhm = [8, 8, 8]
nmCSF.inputs.voxel_size = (1.5, 1.5, 1.5)
nmCSF.inputs.apply_to_files = rc3_images

#Build the workflow up---------------------------------------------------------

create_template = pe.Workflow(name="create_template")

#GMWM nodes conection----------------------------------------------------------
create_template.connect([(dartelGMWM, darteltemplateGMWM,[(('template_files', getIterationParameters),'iteration_parameters')])])
create_template.connect([(dartelGMWM,nmGMWM,[('final_template_file','template_file')])])
create_template.connect([(darteltemplateGMWM,nmGMWM,[('dartel_flow_fields','flowfield_files')])])

#CSF nodes conection-----------------------------------------------------------
create_template.connect([(dartelCSF, darteltemplateCSF,[(('template_files', getIterationParameters),'iteration_parameters')])])
create_template.connect([(dartelCSF,nmCSF,[('final_template_file','template_file')])])
create_template.connect([(darteltemplateCSF,nmCSF,[('dartel_flow_fields','flowfield_files')])])

#Results saving----------------------------------------------------------------
create_template.connect([(dartelGMWM,datasink,[('dartel_flow_fields','flow_fields_GMWM'), ('template_files','template_files_GMWM')])])
create_template.connect([(dartelCSF,datasink,[('dartel_flow_fields','flow_fields_CSF'), ('template_files','template_files_CSF')])])
create_template.connect([(darteltemplateGMWM,datasink,[('dartel_flow_fields','ff_GMWM_existingTempl')])])
create_template.connect([(darteltemplateCSF,datasink,[('dartel_flow_fields','ff_CSF_existingTempl')])])  
create_template.connect([(nmGMWM,datasink,[('normalized_files','normalizedFiles_WMGM')])])
create_template.connect([(nmCSF,datasink,[('normalized_files','normalizedFiles_CSF')])])

""""""""""""
#Estatistic--------------------------------------------------------------------
""""""""""""

#def getTipoClasses(dartel_files, tipo):
    
#    resultado = []
#    n = len(dartel_files)/2
#    if tipo==0:
#        resultado = dartel_files[0:n]
#    elif tipo==1:
#        resultado = dartel_files[-n:]
#    return resultado

        

#staticsGM = stst.statisticWF()

#staticsGM.inputs.datasink.base_directory = os.path.abspath('./ResultsGM')
#staticsGM.inputs.conestimate.contrasts = stst.inputContrast()
#staticsGM.inputs.onesamplettestdes.covariates = stst.inputCovariates(ngrupo1, ngrupo2)
##staticsGM.inputs.onesamplettestdes.in_files = gray

#staticsWM = staticsGM.clone(name='staticsWM')
##staticsWM.inputs.onesamplettestdes.in_files = white
#staticsWM.inputs.datasink.base_directory = os.path.abspath('./ResultsWM')

#staticsCSF = staticsGM.clone(name='staticsCSF')
##staticsCSF.inputs.onesamplettestdes.in_files = csf
#staticsCSF.inputs.datasink.base_directory = os.path.abspath('./ResultsCSF')


#create_template.connect([(nm,staticsGM,[(('normalized_files', getTipoClasses, 0),'onesamplettestdes.in_files')])])
#create_template.connect([(nm,staticsWM,[(('normalized_files', getTipoClasses, 1),'onesamplettestdes.in_files')])])
#create_template.connect([(nmCSF,staticsCSF,[('normalized_files','onesamplettestdes.in_files')])])

""""""""""""
""""""""""""

#Execution of the workflow-----------------------------------------------------

try:
    try:
        create_template.write_graph(graph2use='colored')
    except Exception as e:
        print 'Error al crear el gráfico'
    try:
        create_template.write_graph(graph2use='flat')
    except Exception as e:
        print 'Error al crear el gráfico'
        
    create_template.run(plugin='Linear')
    filepath = datasink.inputs.base_directory

    outs1 = glob.glob(filepath + "/*.nii.gz")
    
 #   print str(outs1)


        
except Exception as e:
    raise e
