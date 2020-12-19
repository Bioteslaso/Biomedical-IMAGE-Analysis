from nipype.interfaces import fsl
from nipype.interfaces.fsl import Eddy

'''
BET (Brain Extraction Tool) deletes non-brain tissue from an image
of the whole head. It can also estimate the inner and the oter skull surfaces, and outer 
scalp surface. 

We will also create the brain mask in this step. 


Andy's process: andy has a folder with acqparams.txt, bvecs, bvals, dwidata.nii.gz, nodif_PA.nii.gz

fslroi dwidata.nii.gz nodif 0 1 (get the image in the other direction)

fslflirt in order to correct dimensions

fslmerge -t AP_PA_b0 nodif nodif_PA (get the two images merged in a single one AP_PA_b0)

topup --imain= AP_PA_b0 --datain= acqparams.txt --config=b02b0.cnf --out= topup_AP_PA_b0

applytopup --imain=nodif, nodif_PA --topup=topup_AP_PA_b0 --datain=acqparams.txt --inindex=1,2 --out=hifi_nodif

Now we must run Eddy. In order to do Eddy, first we must to the brain extraction and then the brain mask

'''

#bet
btr = fsl.BET()

#Input file, which is the output from the applytopup (hifi)
btr.inputs.in_file = ''

#Out
btr.inputs.out_file ='input_brain.nii' #Output file, which will be the brain extraction file

#-f: fractional intensity threshol. 0.5 is the default
btr.inputs.frac = 0.5

#-m: generate the brainmask
btr.inputs.mask = True

res = btr.run()


'''
Now we can apply Eddy, a tool which corrects Eddy distortions. 

'''
#Eddy function
eddy = Eddy()

#--imain= File containing all the images to estimate distortions for flag (dwidata)
eddy.inputs.in_file = 'epi.nii'

#--mask= File containing the brain mask corresponding to flag
eddy.inputs.in_mask  = 'epi_mask.nii'

#--index= file containing the indices for all volumes in --imain into --acqp and --topup
eddy.inputs.in_index = 'epi_index.txt'

#--acqp= File containing the acquisition parameters
eddy.inputs.in_acqp  = 'epi_acqp.txt'

#--bvecs= File containing the b-vectors for all volumes in --imain
eddy.inputs.in_bvec  = 'bvecs.scheme'

#--bvals= File containing the b-values for all volumes in --imain
eddy.inputs.in_bval  = 'bvals.scheme'

#Run Eddy
res = eddy.run()