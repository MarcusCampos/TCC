This code is a MATLAB implementation of the tracking algorithm in CVIU2013 paper 
        "Segmentation of the left ventricle in cardiac cine MRI using shape-constrained snake model.
 Computer Vision and Image Understanding(CVIU), 117(9): 990-1003, 2013." 
               by Yuwei Wu, Yuanquan Wang, and Yunde Jia.
-------------------------------------------------------------------------------------

This source code is shared  for research purposes only. Commercial uses are not allowed 
without our permission. This is the first version of the distribution.
Please use them at your own risk, and we appreciate any comments/suggestions. 
For more quetions, please contact us at wuyuwei@bit.edu.cn,yqwang@bit.edu.cn

The main function is "Run_Segmentation.m"
	
Note: Since the snake model is sensitive to the contour initialization, different image should
be initialized by different contour. Gerenrally the initializaion is completed mannually. In addition, 
The results are slightly sensitive to some parameters. You may obtain better results by adjusting the parameters. 