# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 13:10:10 2016

@author: Laura Nuñez
@eddited by: Mario Gil

"""
#Import needed interfaces------------------------------------------------------
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

#Define nedded functions-------------------------------------------------------
def dicom2nifti(input, method):

	import subprocess
	import nipype.interfaces.freesurfer as fs
	
	if method == "dcm2nii":
	
		mkdir = "mkdir -p ./tmp"
		out, err = subprocess.Popen(mkdir, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
		if err != '':
			raise Exception (out,err)
		
		#Dicom to nifti conversion
		dicom_files = glob.glob(input + "/*")[0]
		dicom2nii = "dcm2nii -a y -d n -e n -f y -g n -i n -r n -x n -o " + "./tmp" + "/" + " " + dicom_files
		out, err = subprocess.Popen(dicom2nii, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
		print out
		if err != '' or "Unsupported Transfer Syntax" in out or "Error" in out:
			raise Exception (out, err)
		
		nifti_file = glob.glob("./tmp/*.nii")[0]
		chname =  "mv " + nifti_file + " " + input + ".nii"
		out, err = subprocess.Popen(chname, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
		if err != '':
			raise Exception (out,err)

	elif method == "FS":
		
		# Dicom to nifti conversion
		mc = fs.MRIConvert()
		mc.inputs.out_file = input + ".nii"
		mc.inputs.out_type = 'nii'
		mc.inputs.in_file = glob.glob(input + "/*")[0]
		mc.inputs.terminal_output = 'none'
		mc.run()
			
	return input + ".nii"


def orientationACPC(input):
	
	import subprocess

	# Orientantion in AC-PC plane
	acpcd = "acpcdetect -i " + input + " -o " + input
	out, err = subprocess.Popen(acpcd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if err != '':
		raise Exception (out, err)
	
	return "Done"		

def coregistration(input):
    
    import os
    import nipype.interfaces.spm as spm
    
    cor = spm.Coregister()
    cor.inputs.target = os.path.abspath('./TemplateCoregister/16032005TemplateCoregisterT1.nii')
    cor.inputs.source = str(input)
    cor.inputs.cost_function = 'nmi'
    cor.inputs.fwhm = [7,7]
    cor.inputs.tolerance =[0.02, 0.02, 0.02, 0.001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001]
    cor.inputs.separation = [4,2]
    cor.inputs.out_prefix = ''
    
    outs = cor.run()
    return outs

def add_extension_gz (input):

	import subprocess
	
	addextension = "gzip -q " + input		
	out, err = subprocess.Popen(addextension, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if err != '':
		raise Exception (out, err)	
	
	return input + ".gz"


#File settings----------------------------------------------------------------- 
originalPath = os.getcwd()

#spm.SPMCommand.set_mlab_paths(matlab_cmd=None,use_mcr=True)    
#spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)

#matlab_cmd = '/Users/Alex/Desktop/spm8/run_spm8.sh /Users/Alex/Documents/MATLAB/v713/ script'
#matlab_cmd = '/Users/lauranunez/Documents/MATLAB/spm12_old/spm12_S/run_spm12.sh /Applications/MATLAB/MATLAB_Compiler_Runtime/v717/ script'
matlab_cmd = '/Volumes/PassportMarioGilCorrea/Programs/SPM/spm8_runtime/run_spm8.sh /Applications/MATLAB/MATLAB_Compiler_Runtime/v84/ script'
spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)
#mlab = matlab.MatlabCommand()
#mlab.inputs.script = "which('who')"
#mlab.set_default_paths('/Users/lauranunez/Documents/MATLAB/spm12')
#mlab.set_default_matlab_cmd("/Applications/MATLAB_R2012a.app/bin/matlab -nodesktop -nosplash")
#out = mlab.run() 


#vol_grp_t1 = glob.glob('/Users/Helena/Desktop/OLFATO/OLFATO_ANATOMICO/DICOM/HOMBRES/*')
#vol_grp_t2 = glob.glob('/Users/Helena/Desktop/OLFATO/OLFATO_ANATOMICO/DICOM/MUJERES/*')
vol_grp_t1 = ['/Volumes/PassportMarioGilCorrea/Programs/ProyectoMAPFRE_fMRI_preprocesados/LAGUELA_CARBALLOSA_CONCEPCION_PRE/id0014/ANATOMICO/SAG_3D']
vol_grp_t2 = ['/Volumes/PassportMarioGilCorrea/Programs/ProyectoMAPFRE_fMRI_preprocesados/LAGUELA_CARBALLOSA_CONCEPCION_PRE/id0014/ANATOMICO/SAG_3D']

#vol_grp_t1 = glob.glob('./HOMBRES/*')
#vol_grp_t2 = glob.glob('./MUJERES/*')
#
nifti_T1_group1 = []
for vol in vol_grp_t1:
    try:
        
        nifti_T1 = dicom2nifti (vol, "dcm2nii")
        out = coregistration(nifti_T1)
        nifti_T1_group1.append(os.path.abspath(nifti_T1))
	
    except Exception as e:
        try:
            nifti_T1 = dicom2nifti (vol, "FS")
            coregistration(nifti_T1)
            nifti_T1_group1.append(os.path.abspath(nifti_T1))
        except Exception as e:
            print 'Error: ' + str(e)
            

nifti_T1_group2 = []
for vol_t1 in vol_grp_t2:
    try:
        nifti_T1 = dicom2nifti (vol_t1, "dcm2nii")
        coregistration(nifti_T1)
        nifti_T1_group2.append(os.path.abspath(nifti_T1))
    	
    except Exception as e:
        try:
            nifti_T1 = dicom2nifti (vol_t1, "FS")
            coregistration(nifti_T1)
            nifti_T1_group2.append(os.path.abspath(nifti_T1))
        except Exception as e:
            print 'Error: ' + str(e)
      
os.chdir(originalPath)


#nifti_T1_group1 = glob.glob(os.path.abspath('./HOMBRES/id*.nii'))
#nifti_T1_group2 = glob.glob(os.path.abspath('./MUJERES/id*.nii'))

#nifti_T1_group1 = ['/Users/lauranunez/Desktop/LAURA/Balsalobre_Munoz_Jose_Vicente/Olfato-5584/SAG_3D_IR_4.nii','/Users/lauranunez/Desktop/LAURA/Sagues_Monente_Pedro/Olfato-5561/SAG_3D_IR_4_3.nii']
#nifti_T1_group2 = ['/Users/lauranunez/Desktop/LAURA/Heredia_Nieto_Enrique/Olfato-4871/SAG_3D_IR_4_2.nii','/Users/lauranunez/Desktop/LAURA/Mov14/6124960351256576_.nii']
#nifti_T1_group1 = ['/Volumes/PassportMarioGilCorrea/Programs/ScriptsVBM/Pruebas/SAG_3D_IR_4.nii','/Volumes/PassportMarioGilCorrea/Programs/ScriptsVBM/Pruebas/SAG_3D_IR_4_3.nii']
#nifti_T1_group2 = ['/Volumes/PassportMarioGilCorrea/Programs/ScriptsVBM/Pruebas/id2861_SAG_3D.nii','/Volumes/PassportMarioGilCorrea/Programs/ScriptsVBM/Pruebas/id3387_SAG_3D.nii']


print nifti_T1_group1
print nifti_T1_group2

#Building up of Segmentation node ---------------------------------------------

seg = pe.Node(interface=spm.NewSegment(), name="seg")
seg.inputs.channel_files = nifti_T1_group1
seg.inputs.channel_info = (0.0001,60,(False,False))
path_tpm = '/Volumes/PassportMarioGilCorrea/Programs/SPM/spm8_updates_r6313/toolbox/Seg/TPM.nii'
tissue1 = ((os.path.abspath(path_tpm), 1), 2, (True,True), (False, False))
tissue2 = ((os.path.abspath(path_tpm), 2), 2, (True,True), (False, False))
tissue3 = ((os.path.abspath(path_tpm), 3), 2, (True,True), (False, False))
tissue4 = ((os.path.abspath(path_tpm), 4), 3, (False,False), (False, False))
tissue5 = ((os.path.abspath(path_tpm), 5), 4, (False,False), (False, False))
tissue6 = ((os.path.abspath(path_tpm), 6), 2, (False,False), (False, False))
seg.inputs.tissues = [tissue1, tissue2, tissue3, tissue4, tissue5, tissue6]
seg.inputs.warping_regularization = 4
seg.inputs.sampling_distance = 3
seg.inputs.write_deformation_fields = [False, False]
seg.inputs.affine_regularization = 'none'

seg2=seg.clone(name='seg2')
seg2.inputs.channel_files = nifti_T1_group2
#Building up of datasink node -------------------------------------------------
datasink = pe.Node(interface=nio.DataSink(), name="datasink")
datasink.inputs.base_directory = os.path.abspath('./Segmentation_Results')

#Build workflow up------------------------------------------------------------- 
Segmentation = pe.Workflow(name="Segmentation")

Segmentation.connect([(seg,datasink,[('native_class_images','segmentationsG1')])])
Segmentation.connect([(seg,datasink,[('dartel_input_images','rc_filesG1')])])

Segmentation.connect([(seg2,datasink,[('native_class_images','segmentationsG2')])])
Segmentation.connect([(seg2,datasink,[('dartel_input_images','rc_filesG2')])])

#Execution of the workflow-----------------------------------------------------
try:
    try:
        Segmentation.write_graph(graph2use='colored')
    except Exception as e:
        print 'Error al crear el gráfico'
    try:
        Segmentation.write_graph(graph2use='flat')
    except Exception as e:
        print 'Error al crear el gráfico'
        
    Segmentation.run(plugin='Linear')
    filepath = datasink.inputs.base_directory

    outs1 = glob.glob(filepath + "/*.nii.gz")
    
 #   print str(outs1)


        
except Exception as e:
    raise e
