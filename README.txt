


	--> The folder 'RUDUR_Algo/' contains the RUDUR algorithm's MATLAB files. 


---------------------------------------------------------
-------------------------- DEMO -------------------------

	--> A Demo is available in the script 'Demo_toyExample.m' to test RUDUR on the toy example.
	--> You can modify parameters of the dataset (background/no background, Noise)
	--> You can also modify parameters of RUDUR (alpha, beta, gamma, mu)
	
	

---------------------------------------------------------
------------ Description of the MATLAB files ------------
		
		- /dataAF.mat : Contains ground truth matrix A and F to genere data
		- /Demo_toyExample.m : Demo of RUDUR algorithm on the toy example
		-	/normalize.m : Normalize matrix A and F to compare estimation with ground truth
		-	/sortSources.m : Sort estimate sources in the same order than ground truth
		
		
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
			- /RUDUR_Algo/compute_fWLS.m : Compute 'fWLS' criterion.
			
		- /RUDUR_Algo/CorrectAndNormalize.m : Delete negative values and normalize
		

