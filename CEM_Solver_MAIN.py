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
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Sinlge_Ele_Natural_Sys\Single_Ele_Natural_SYS-0000.dat"
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Sinlge_Ele_Natural_Sys\Single_Ele_Natural_SYS-0001.dat"
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Sinlge_Ele_Natural_Sys\Single_Ele_Natural_SYS-0002.dat"
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Sinlge_Ele_Natural_Sys\Single_Ele_Natural_SYS-0003.dat"
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Sinlge_Ele_Natural_Sys\Single_Ele_Natural_SYS-0004.dat"
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Sinlge_Ele_Natural_Sys\Single_Ele_Natural_SYS-0005.dat"
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Sinlge_Ele_Origin\SingleEle_Translated-0000.dat"
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Sinlge_Ele_Origin\SingleEle_Translated_Rotated-0000.dat"
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Sinlge_Ele_Origin\SingleEle_Translated_Rotated_180-0000.dat"
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Sinlge_Ele_Origin\SingleEle_Translated_Rotated_270-0000.dat"
    
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Single_45_Element\Pos45-0001.dat"
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Single_45_Element\Neg45-0000.dat"
    
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Sinlge_Ele_Origin\SingleEle_Translated_Skewed-0000.dat"
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\Sinlge_Ele_Origin\SingleEle_Translated_Skewed_4ELe-0001.dat"
    
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\DistortedPlate\DistortedPlate-0000.dat"
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\DistortedPlate\DistortedPlate-0001.dat"
    
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\SimplePlateModel\SimplePlateModel-0000.dat"
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\SimplePlateModel\SimplePlateModel_CQ8-0000.dat"
    
    # File = r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Nastran\HoleInPlate\HoleInPlate-0001.dat"
    File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\Final_Models\Hole_In_Square_Plate\MyFEM\My_Hole_In_Square_Plate.dat"
    
    ''' Creating Results Dir '''
    outDir = resultsDir(File)
    
    ''' Parsing Nastran Dat File '''
    Parsed = DatFileParsing(File)
    dfNodes, dfEle4, dfEle6, dfEle7, dfEle8, dfMatProps, dfForces, dfConstraints = Parsed
    
    ''' Calculating gloabl displacements '''
    dfDisp = DispCalc.DispCalc(dfNodes, dfEle4, dfEle8, dfEle7, dfForces, dfConstraints, dfMatProps)
    dfDisp.to_csv(outDir + "\\01_gloabal_displacements.csv")    
    
    ''' Calculating Strain and Stresses '''
    dfTemp = [dfEle4, dfEle8, dfEle7, dfEle6]
    eTemp = ["CQ4", "CQ8", "CQ7", "CQ6"]
    dfEles = [x for x in dfTemp if isinstance(x, pd.DataFrame)]
    eTypes = [x for x, y in zip(eTemp, dfTemp) if isinstance(y, pd.DataFrame)]
    
    dfStress, dfStrain = StressStrainCalc(dfEles, eTypes, dfDisp, dfNodes, "PlaneStress", dfMatProps)
    dfStrain.to_csv(outDir + "\\02_gloabal_strains.csv")
    dfStress.to_csv(outDir + "\\03_gloabal_stresses.csv")
    
    
    '''Plotting'''
    # dfDefNodes = DispCalc.build_deformed_nodes(dfNodes, dfDisp, scale = 50)
    # DispCalc.plot_deformed_mesh(dfDefNodes, dfEles)
    Contour(dfNodes, dfDisp, ["U1", "U2"], dfEles)
    # # Contour(dfNodes, dfStrain, ["E1", "E_max"], dfEles)
    # Contour(dfNodes, dfStress, ["S1", "S2", "S12", "S_max"], dfEles)
    Contour(dfNodes, dfStress, ["S1", "S2", "S12", "S_max"], dfEles, Averaging="Nodal")
    
    # holeNodeIDs = [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 37, 41, 44, 47, 50, 53, 56, 59, 62, 65, 68, 
    #                71, 74, 77, 80, 83, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 116, 119, 122, 
    #                125, 128, 131, 134, 137, 140, 143, 146, 149, 152, 155, 158, 160]

    holeNodeIDs = [34, 50, 27058, 27059, 27060, 27061, 27062, 27063, 27064, 27065, 27066, 27067, 27068, 27069, 27070, 27071, 27072, 29150, 
                   29151, 29152, 29153, 29154, 29155, 29156, 29157, 29158, 29159, 29160, 29161, 29162, 29163, 29164]


    # dfOut = rThetaStress(dfNodes, dfStress, hole_center = np.array([3.0, 3.0]), hole_node_ids = holeNodeIDs)
    
    print("stop")