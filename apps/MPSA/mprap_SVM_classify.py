###PROGRAM running svm classify

import sys,re,os
import operator
#from string import *
#from shutil import *
#from math import *

############################################################################
#input file
inputSVM=sys.argv[1] #svm vector from script create SVM_table
#output file
outSVMclassify=sys.argv[2]


############################################################################
# directory where SVM_light (svm_classify) is installed
path2=sys.argv[3]
# directory containing the SVM model file
path3=sys.argv[4]


######################
#run svm_classify
print path2+'svm_classify'+' '+inputSVM+' '+path3+'SVM_model_mpRAP_12_12'+' '+outSVMclassify

argsystem = os.system(path2+'svm_classify'+' '+inputSVM+' '+path3+'SVM_model_mpRAP_12_12'+' '+outSVMclassify)


