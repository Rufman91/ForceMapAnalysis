import sys
import os
from deconvolute_by_mathematical_morphology import deconvolute_by_mathematical_morphology

def run_deconvolution(python_file_dir, InImage, ErodingGeometry):
    sys.path.insert(0, python_file_dir)
    
    return deconvolute_by_mathematical_morphology(InImage, ErodingGeometry)
