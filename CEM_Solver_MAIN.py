import sys
import os

sys.path.append(r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\CEM_Solver_V2\Lib")

from DatFileParser import DatFileParsing
from DispCalc import DispCalc

def resultsDir(file):
    baseName = os.path.splitext(os.path.basename(file))[0]
    cwd = os.getcwd()
    newFolder = os.path.join("C:\\Documents\\Grad_School\\MastersFinalProject\\CodeBase\\CEM_Solver_V2\\Results", baseName)
    os.makedirs(newFolder, exist_ok=True)
    return newFolder

if __name__ == "__main__":
    # File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\Quad4\TensionPlateFEMAP\TensionPlate-0000_Mod_Mat.dat"
    # File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\2_Quad_Model\2Quad-0000.dat"
    # File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\1_Quad_Model\Single_Quad-0001.dat"
    # File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\Quad4\Non_Square\CQ4_Half_Mixed_Plate-0001.dat"
    # File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\Quad4\Non_Square\CQ4_Half_Mixed_Plate_Sequential-0000.dat"
    
    File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\Mixed_Directions\MixedOrthogonalDir\Single_Ele\Membrane\Single_Square_Membrane-0000.dat"
    # File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\Mixed_Directions\MixedOrthogonalDir\Single_Ele\SingleEle_Plate-0000.dat"
    
    # File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\Mixed_Directions\MixedOrthogonalDir\SquarePlate_Same_Directions-0000.dat"
    # File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\Mixed_Directions\MixedOrthogonalDir\SquarePlate_Diff_Directions_Renumbered-0000.dat"
    # File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\Mixed_Directions\MixedOrthogonalDir\SquarePlate_Diff_Directions-0001.dat"
    # File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\Mixed_Directions\MixedOrthogonalDir\SquarePlate_Same_Directions_VertLoad-0000.dat"
    # File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\Mixed_Directions\MixedOrthogonalDir\NonSquarePlate_Same_Dir-0000.dat"
    # File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\Mixed_Directions\MixedOrthogonalDir\NonSquarePlate_Diff_Dir.dat"
    
    # File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\Quad8\Single\Single_C8-0000.dat"
    # File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\Quad8\Double\Double_C8-0000.dat"
    # File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\Quad8\FullPlate\FullPlate_EdgeLoad.dat"
    # File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\Mixed_Elements\MixedPlates\MixedPlate-0000.dat"
    # File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\Mixed_Elements\SmallMixedPlate\SmallMixedPlate-0000.dat"
    # File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\Mixed_Elements\Q7\Q7_Mixed.dat"
    # File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\Mixed_Elements\MultiDirection\ForwardAndAft.dat"
    # File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\Mixed_Elements\MultiDirection\ForAftTopBottom.dat"

    # File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\Final_Models\Hole_In_Square_Plate\MyFEM\My_Hole_In_Square_Plate.dat"
    # File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\Final_Models\Compounding_Tension_Strip\compoundingtensionsstrip_femap_cq4-0001.dat"
    # File = r"C:\Documents\Grad_School\MastersFinalProject\CodeBase\Nastran\Final_Models\Compounding_Tension_Strip\CompoundingTensionsStrip_MyFEM.dat"
    ''' Creating Results Dir '''
    outDir = resultsDir(File)
    
    ''' Parsing Nastran Dat File '''
    Parsed = DatFileParsing(File)
    dfNodes, dfEle4, dfEle6, dfEle7, dfEle8, dfMatProps, dfForces, dfConstraints = Parsed
    
    ''' Calculating gloabl displacements '''
    dfDisp = DispCalc(dfNodes, dfEle4, dfEle8, dfEle7, dfForces, dfConstraints, dfMatProps)
    dfDisp.to_csv(outDir + "\\01_gloabal_displacements.csv")    
    

    print("stop")