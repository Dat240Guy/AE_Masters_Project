import sys
import os
import numpy as np
import pandas as pd

sys.path.append(r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Lib")

from DatFileParser import DatFileParsing
from StressStrainCalc import StressStrainCalc, rThetaStress
from FringePlotter import Contour
import DispCalc

def resultsDir(file):
    baseName = os.path.splitext(os.path.basename(file))[0]
    cwd = os.getcwd()
    newFolder = os.path.join(r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Results", baseName)
    os.makedirs(newFolder, exist_ok=True)
    return newFolder

if __name__ == "__main__":
    print("welcome")
    
    '''Hole in Plate Models'''
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\HoleInPlate\HoleInPlate-0001.dat"
    # File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\Final_Models\Hole_In_Square_Plate\MyFEM\My_Hole_In_Square_Plate.dat"
    
    '''Compounding Tension Strip Models'''
    #PyDFEM Model
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\!FinalModels\CompoundTensionStrip\PyDFEM\CompoundTensionStripPyDFEM.dat"
    File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\!FinalModels\CompoundTensionStrip\PyDFEM\CompoundTensionStripPyDFEM_Refined_GoodBoundary.dat"
    
    ''' Creating Results Dir '''
    outDir = resultsDir(File)
    
    ''' Parsing Nastran Dat File '''
    Parsed = DatFileParsing(File)
    dfNodes, dfEle4, dfEle3, dfEle6, dfEle7, dfEle8, dfMatProps, dfForces, dfConstraints = Parsed
    
    ''' Calculating gloabl displacements '''
    dfDisp = DispCalc.DispCalc(dfNodes, dfEle4, dfEle3, dfEle8, dfEle7, dfForces, dfConstraints, dfMatProps)
    dfDisp.to_csv(outDir + "\\01_gloabal_displacements.csv")    
    
    ''' Calculating Strain and Stresses '''
    dfTemp = [dfEle4, dfEle8, dfEle7, dfEle6, dfEle3]
    eTemp = ["CQ4", "CQ8", "CQ7", "CQ6", "CTRIA3"]
    dfEles = [x for x in dfTemp if isinstance(x, pd.DataFrame)]
    eTypes = [x for x, y in zip(eTemp, dfTemp) if isinstance(y, pd.DataFrame)]
    
    dfStress, dfStrain = StressStrainCalc(dfEles, eTypes, dfDisp, dfNodes, "PlaneStress", dfMatProps)
    dfStrain.to_csv(outDir + "\\02_gloabal_strains.csv")
    dfStress.to_csv(outDir + "\\03_gloabal_stresses.csv")
    
    
    '''Plotting'''
    dfDefNodes = DispCalc.build_deformed_nodes(dfNodes, dfDisp, scale = 50)
    DispCalc.plot_deformed_mesh(dfDefNodes, dfEles)
    Contour(dfNodes, dfDisp, ["U1", "U2"], dfEles)
    # # Contour(dfNodes, dfStrain, ["E1", "E_max"], dfEles)
    # Contour(dfNodes, dfStress, ["S1", "S2", "S12", "S_max"], dfEles)
    Contour(dfNodes, dfStress, ["S1", "S2", "S12", "S_max"], dfEles, Averaging="Nodal")
        
    print("stop")