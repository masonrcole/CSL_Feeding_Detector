# CSL_Feeding_Detector
The California sea lion accelerometer-based feeding detector (MatLab code) described in Cole et al., 2021 (in Animal Biotelemetry)

Notes about the code 'CSL_Feeding_Detector.m'

1. The function CSL_Feeding_Detector.m calls the script CSL_Feeding.m

2. Changes (i.e. inputs) should only be made to CSL_Feeding_Detector.m, 
	unless one wishes to make changes to the function of the code
	
3. the freely-available function 'movingmean.m' is needed for the code 
	to function.  That function is available here: https://www.mathworks.com/matlabcentral/fileexchange/41859-moving-average-function
	
4. Directions for use are given at the beginning of the function 'CSL_Feeding_Detector.m'

5. Outputs are found in the workspace 


Notes about proper use of this code (as of 2021) 

1. See Cole et al. (2021) in Animal Biotelemetry for the effect of sampling rate on true positive and false positive detection rates 
	
2. As of mid 2021, this code has only been validated on adult california sea lions in a captive setting. 
	Use of this code on wild California sea lions requires either animal-borne video validation or 
	several major assumptions. It is likely that true positive and false positive detection rates 
	will differ from published values in a wild foraging context. Species-specific validation should occur before or with
	the use of this code on other species. 
	
3. It is possible that the 'stock' surge-axis acceleration thresholds included in this code will need to be adjusted 
	in wild studies or to better detect feeding in some individuals. Currently, this can be accomplished
	by adjusting values set in lines 59 and 60 of the 'internal' CSL_Feeding.m script. 
	

Please contact Mason Cole at masonrcole@gmail.com for help, questions, or comments regarding this code. This was some 
of the first code that Mason wrote, so please forgive the inefficient style with which some of it is written. 

Cheers,
Mason 
