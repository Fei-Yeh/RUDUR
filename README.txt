

--> The folder 'RUDUR_Algo/' contains the RUDUR algorithm's MATLAB files. 
--> The folder 'Demo/' containes the demo files.
--> The folder 'Dataset/' contains several dataset with the groundtruth

-------------------------------------------------------------------------
---------------------------------- DEMO ---------------------------------

All results found in the thesis "Séparation de sources en imagerie nucléaire"
can be reproduced with the demo files in the folder 'Demo/'. User must be
in the folder 'Demo/' in order to execute these scripts. Only results on 
the clinical renography dataset cannot be reproduced, because the dataset 
cannot be redistributed. Please visit http://www.dynamicrenalstudy.org to
download this dataset.
	
	Section 4.1 : - /Demo/Demo_SimpleDataset.m
                  - /Demo/Demo_SimpleDataset100.m
    
    Section 4.2 : - /Demo/Demo_SynteticRenography.m
                  - /Demo/Demo_SyntheticRenography6.m

    Section 4.4.1 :  - /Demo/Demo_RobustParam.m
    Section 4.4.2 :  - /Demo/Demo_Alpha.m
    Section 4.4.3 :  - /Demo/Demo_Beta.m
    Section 4.4.4 :  - /Demo/Demo_RobustROI.m

    Section 4.5 : - /Demo/Demo_WeakSource.m

    Section 4.6.1 : - /Demo/Demo_6DIG.m
    Section 4.6.3 : - /Demo/Demo_6DIGsynth.m


-------------------------------------------------------------------------
--------------------------------- RUDUR ---------------------------------

Folder 'RUDUR_Algo/' : Contains the RUDUR algorithm's MATLAB files.
    - /RUDUR_Algo/rudur.m : Main function.

    - /RUDUR_Algo/buildPrior.m : Build matrix 'D' of distance to ROIs.
    - /RUDUR_Algo/buildWeight.m : Build matrix of weights 'W'.

    - /RUDUR_Algo/optimA_ConjGrad.m : Optimize according to 'A' with 'F' fixed.
    - /RUDUR_Algo/optimF_ConjGrad.m : Optimize according to 'F' with 'A' fixed.

    - /RUDUR_Algo/conjugateGrad.m : Compute conjugate gradient.
    - /RUDUR_Algo/descentA.m : Update matrix 'A'.
    - /RUDUR_Algo/descentF.m : Update matrix 'F'.
    - /RUDUR_Algo/grad_fObj.m : Compute gradient of objective function 'fRUDUR'.
        - /RUDUR_Algo/grad_fROI.m : Compute gradient of 'fPrior' criterion.
        - /RUDUR_Algo/grad_fTik.m : Compute gradient of 'fTik' criterion.
        - /RUDUR_Algo/grad_fWLS.m : Compute gradient of 'fROI' criterion.

    - /RUDUR_Algo/compute_fObj.m : Compute objective function 'fRUDUR'.
    - /RUDUR_Algo/compute_fObj_bis.m : Compute objective function 'fRUDUR'.
        - /RUDUR_Algo/compute_fROI.m : Compute 'fROI' criterion .
        - /RUDUR_Algo/compute_fTik.m : Compute  'fTik' criterion.
        - /RUDUR_Algo/compute_fQWLS.m : Compute 'fQWLS' criterion.

    - /RUDUR_Algo/CorrectAndNormalize.m : Delete negative values and normalize
		
