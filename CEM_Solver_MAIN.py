import sys
import os
import pandas as pd

sys.path.append(r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Lib")

from DatFileParser import DatFileParsing
from DispCalc import DispCalc
from StressStrainCalc import StressStrainCalc
from FringePlotter import Contour

def resultsDir(file):
    baseName = os.path.splitext(os.path.basename(file))[0]
    cwd = os.getcwd()
    newFolder = os.path.join(r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Results", baseName)
    os.makedirs(newFolder, exist_ok=True)
    return newFolder

if __name__ == "__main__":
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Sinlge_Ele_Natural_Sys\Single_Ele_Natural_SYS-0000.dat"
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Sinlge_Ele_Natural_Sys\Single_Ele_Natural_SYS-0001.dat"
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Sinlge_Ele_Natural_Sys\Single_Ele_Natural_SYS-0002.dat"
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Sinlge_Ele_Natural_Sys\Single_Ele_Natural_SYS-0003.dat"
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Sinlge_Ele_Natural_Sys\Single_Ele_Natural_SYS-0004.dat"
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Sinlge_Ele_Natural_Sys\Single_Ele_Natural_SYS-0005.dat"
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Sinlge_Ele_Origin\SingleEle_Translated-0000.dat"
    File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Sinlge_Ele_Origin\SingleEle_Translated_Rotated-0000.dat"
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Sinlge_Ele_Origin\SingleEle_Translated_Rotated_180-0000.dat"
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Sinlge_Ele_Origin\SingleEle_Translated_Rotated_270-0000.dat"
    
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Sinlge_Ele_Origin\SingleEle_Translated_Skewed-0000.dat"
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Sinlge_Ele_Origin\SingleEle_Translated_Skewed_4ELe-0001.dat"
    
    ''' Creating Results Dir '''
    outDir = resultsDir(File)
    
    ''' Parsing Nastran Dat File '''
    Parsed = DatFileParsing(File)
    dfNodes, dfEle4, dfEle6, dfEle7, dfEle8, dfMatProps, dfForces, dfConstraints = Parsed
    
    ''' Calculating gloabl displacements '''
    dfDisp = DispCalc(dfNodes, dfEle4, dfEle8, dfEle7, dfForces, dfConstraints, dfMatProps)
    dfDisp.to_csv(outDir + "\\01_gloabal_displacements.csv")    
    
    dfTemp = [dfEle4, dfEle8, dfEle7, dfEle6]
    eTemp = ["CQ4", "CQ8", "CQ7", "CQ6"]
    dfEles = [x for x in dfTemp if isinstance(x, pd.DataFrame)]
    eTypes = [x for x, y in zip(eTemp, dfTemp) if isinstance(y, pd.DataFrame)]
    
    dfStress, dfStrain = StressStrainCalc(dfEles, eTypes, dfDisp, dfNodes, "PlaneStress", dfMatProps)

    Contour(dfNodes, dfDisp, ["U1", "U2"])
    Contour(dfNodes, dfStress, ["S1", "S2", "S12"])
    print("stop")