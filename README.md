# ls-svm-project 
v1.0
The first full version of the app that allows you to build a regression model based on least squares support vector machine (LS SVM) to set the parameters: 
- form of the functional dependence of and the number of factors;
- sample size; 
- varying noise; 
- varying intervals factors; 
- the type and parameters of the kernel function.  

To set up automatic selection of optimal parameters is used the learning algorithm with the teacher (leave-one-out CV).
The result of the program is to assess the quality of the restored functional dependence according the given parameters, which include:
- graphical representation of the produced resumptive model;
- the value of the mean square error (MSE); 
- the value of the quality criterion model (loo cv); 
- the input value and the resulting best estimaten of response variable.

This application allows you to find the optimal customizable parameters of the algorithm depending on a number of conditions of the experiment, and carry out analysis of their choice.

v2.0
In developing. Software implemented the practical part of the modifications. Currently being finalized application interface.
Added:
 - set of valid parameter values sigma^2 RBF-kernel replaced the generation of the grid;
 - weighted procedure LS SVM (robust LS SVM);
 - auto-tuning procedure is the boundary of the parameter for the weighted algorithm;
 - formulated and added generic regularity criterion;
 - formulated and added to the weighted criteria loo CV;
 - graphics rendering results added procedures.
 
 
