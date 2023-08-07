# @ author: Marco Schaefer
# @ Description:

import traceback
import xml.etree.ElementTree as ET
import numpy as np
from biopandas.pdb import PandasPdb
import pandas as pd
from joblib import Parallel, delayed
from io import TextIOBase
import os
import sys
import time
from pathlib import Path, PurePosixPath
import subprocess
import multiprocessing
import platform
curOS = platform.system()
# Linux
# Windows


# import and install non default modules
#packageNames = ['wget','joblib','biopandas','numpy']
#importCommands =  ['wget','joblib','biopandas','numpy']
installCommad = ''
cnt = 0

#global intercfilesDic
global realProtAtomPositions
global protAtomCnt
global fileContent

realProtAtomPositions = {}
#intercfilesDic = []
protAtomCnt = 0
fileContent = []

# disable/enable file cleaning (delete unecessary/intermediate files))
cleaning = True
verbose = False
dockerOrNative = ""
#dockerOrNative = "docker"

# vina-merger
# from: https://github.com/Rmadeye/vina-merger/tree/master/sample
# improved by: MarcoSchaeferT


class Converter:
    def __init__(self, protein_file: Path, flex_file: Path, output_path: Path):
        self.protein = str(protein_file)
        self.flex = str(flex_file)
        self.output_path = str(output_path)

    def convert(self):

        prot_df = PandasPdb().read_pdb(self.protein)
        flex_df = PandasPdb().read_pdb(self.flex)

        # remove all hydrogens
        prot_df.df['ATOM'] = prot_df.df['ATOM'][prot_df.df['ATOM']
                                                ["element_symbol"] != 'H']

        aalist = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS',
                  'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                  'THR', 'TRP', 'TYR', 'VAL', 'HIE', 'HID', 'ASX', 'GLX', 'HIP']

        flex_df.df['ATOM'].loc[
            ~flex_df.df['ATOM']['residue_name'].isin(aalist), 'record_name'] = 'HETATM'

        for index, row in flex_df.df['ATOM'].iterrows():
            prot_df.df['ATOM']['x_coord'] = np.where((prot_df.df['ATOM']['residue_name'] == row['residue_name']) &
                                                     (prot_df.df['ATOM']['residue_number'] == row['residue_number']) &
                                                     (prot_df.df['ATOM']['atom_name']
                                                      == row['atom_name']),
                                                     row['x_coord'], prot_df.df['ATOM']['x_coord'])
            prot_df.df['ATOM']['y_coord'] = np.where((prot_df.df['ATOM']['residue_name'] == row['residue_name']) &
                                                     (prot_df.df['ATOM']['residue_number'] == row['residue_number']) &
                                                     (prot_df.df['ATOM']['atom_name']
                                                      == row['atom_name']),
                                                     row['y_coord'], prot_df.df['ATOM']['y_coord'])
            prot_df.df['ATOM']['z_coord'] = np.where((prot_df.df['ATOM']['residue_name'] == row['residue_name']) &
                                                     (prot_df.df['ATOM']['residue_number'] == row['residue_number']) &
                                                     (prot_df.df['ATOM']['atom_name'] == row['atom_name']), row['z_coord'], prot_df.df['ATOM']['z_coord'])

        rowCnt = 1
        for index, row in prot_df.df['ATOM'].iterrows():
            prot_df.df['ATOM'].at[index, 'atom_number'] = rowCnt
            rowCnt += 1

        if len(flex_df.df['HETATM']) == 0:
            hetatms = flex_df.df['ATOM'].loc[flex_df.df['ATOM']
                                             ['record_name'] == 'HETATM']
            hetatms = hetatms.drop(['line_idx'], axis=1)
            prot_df.df['ATOM'] = pd.concat(
                [prot_df.df['ATOM'], hetatms], ignore_index=True)

        else:
            hetatm_df = flex_df.df['HETATM']
            hetatm_df['segment_id'] = hetatm_df['segment_id'].iloc[0:0]
            hetatm_df['segment_id'] = hetatm_df['segment_id'].fillna('')
            hetatm_df['blank_4'] = hetatm_df['blank_4'].iloc[0:0]
            hetatm_df['blank_4'] = hetatm_df['blank_4'].fillna('')
            hetatm_df = hetatm_df.drop(['line_idx'], axis=1)

            pd.concat([prot_df.df['ATOM'], hetatm_df], ignore_index=True)

        return prot_df.to_pdb(path=self.output_path,
                              # only ATOMS to remove everythibng else (sugars, other molecules, protein modifications)
                              records=['ATOM'],
                              gz=False,
                              append_newline=True)


#########################################################################
babelPath = ""
multiPDBQTListPath = ""
dataFolder = ""
rawProteinPath = ""
dockerPLIP = ""
proteinPath = ""
multiPDBQTList = ""

PDBQTpathList = multiPDBQTList.split("\n")
if (len(sys.argv) > 1 or False):
    babelPath: Path = Path(r"{}".format(sys.argv[1]))
    multiPDBQTListPath: Path = Path(r"{}".format(sys.argv[2]))
    dataFolder: Path = Path(r"{}".format(sys.argv[3]))
    rawProteinPath: Path = Path(r"{}".format(sys.argv[4]))
    if (curOS == "Linux"):
        babelPath = "obabel"

else:
    babelPath: Path = Path(
        "C:\\Users\\Marco\\Project\\megamol-prolint\\plugins\\prolint\\obabel.exe")
    multiPDBQTListPath: Path = Path(
        "C:/Users/MarcoSchaefer/Projects/results_70/results_list.txt")
    dataFolder: Path = Path(
        "C:\\Users\\MarcoSchaefer\\Projects\\results_70\\InVADo_tmpFiles")
    rawProteinPath: Path = Path(
        "C:\\Users\\MarcoSchaefer\\Projects\\results_70\\7nn9.pdbqt")

dockerPLIP = "plip "
if (dockerOrNative == "docker"):
    dockerPLIP = "docker run --rm --volume ${PWD}:/results --workdir /results --user $(id -u ${USER}):$(id -g ${USER}) pharmai/plip:latest "


proteinPath = rawProteinPath.parent / (rawProteinPath.stem + '.pdb')
print("proteinPath: ", str(proteinPath))

multiPDBQTList = open(multiPDBQTListPath, "r")
multiPDBQTList = multiPDBQTList.read()

PDBQTpathList = multiPDBQTList.split("\n")


def getBasedir(multiPDBQTListPath: Path):
    return multiPDBQTListPath.parent


def toLinuxSubPath(path: Path):
    drive = path.drive
    driveStr = str(drive)
    driveStrlow = driveStr[0].lower()
    out = str(Path(path).as_posix()).replace(driveStr, "/mnt/"+driveStrlow)
    if verbose:
        print(str(out))
    return out


def getZINCID(PDBQTpath, PDBQTpathIndex):
    if verbose:
        print("ID"+str(PDBQTpathIndex))
    third_line = ""
    try:
        with open(PDBQTpath) as f:
            f.readline()
            f.readline()
            third_line = f.readline().strip()
    except:
        print("ERROR: Could not load "+str(PDBQTpath))
        return False

    splitRes = third_line.split("Name = ")
    ZINCID = ""
    if (len(splitRes) > 1):
        ZINCID = third_line.split("Name = ")[-1]
    else:
        ZINCID = "lig"+str(PDBQTpathIndex)
    if (ZINCID != ""):
        return ZINCID
    else:
        print("ERROR: Cloud not determine the ZINCID or LigName!")
        return False


def convertMultiPDBQT_to_singleProtLigComplexes(babelPath, PDBQTpath, index):
    convertMultiPDBQT_to_singlePDBQT(babelPath, PDBQTpath, index)
    ZINCID = getZINCID(PDBQTpath, index+1)

    maxModel = getModelCnt(PDBQTpath)
    for j in range(1, maxModel):
        tmpSinglePDBQT_path = dataFolder / (ZINCID+"_"+str(j)+".pdb")
        if (os.path.exists(tmpSinglePDBQT_path)):
            percentage = str(round(((index+1)/len(PDBQTpathList))*100, 1))
            resultPath = dataFolder / (ZINCID+"_"+str(j)+".inter.txt")
            if (os.path.exists(resultPath) == False):
                if verbose:
                    print(str(tmpSinglePDBQT_path)+" "+percentage+"%")
                tmpComplexPDBQT_targetPath = dataFolder / \
                    ("complex"+ZINCID+"_"+str(j)+".pdb")
                con = Converter(proteinPath, tmpSinglePDBQT_path,
                                tmpComplexPDBQT_targetPath)
                con.convert()
                calcProtLigInteract(babelPath, PDBQTpath, j, index)
            if cleaning:
                try:
                    Path.unlink(tmpSinglePDBQT_path)
                except IOError:
                    print(str(IOError))

    #print("Convert XMLs to .plip")
    # readXMLfile(dataFolder)


def convertMultiPDBQT_to_singlePDBQT(babelPath, PDBQTpath, index):
    ZINCID = getZINCID(PDBQTpath, index+1)
    tmpSinglePDBQT_path = dataFolder / (str(ZINCID)+"_.pdb")
    command = str(babelPath) + " -i pdbqt " + str(PDBQTpath) + \
        " -O " + str(tmpSinglePDBQT_path) + " -m"
    if (curOS != "Linux"):
        silentOutput = ">/dev/null 2>&1"
        if verbose:
            silentOutput = ""
        command = "bash -c 'eval \""+str("obabel")+" -i pdbqt " + toLinuxSubPath(
            PDBQTpath) + " -O " + toLinuxSubPath(tmpSinglePDBQT_path) + " -m "+silentOutput+"\"'"
    if verbose:
        print(command)
    os.system(command)


def processFolder(folderPath, search):
    pngs = []
    if verbose:
        print("searching for: "+str(search))
    for r, d, f in os.walk(str(folderPath)):
        for file in f:
            if file.startswith(search) == True and file.endswith("png"):
                if verbose:
                    print(file)
                pngs.append(file)
    return pngs

# source: https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console


global maxCurIter
global progSign
progSign = "-"
maxCurIter = 0


def printProgressBar(iteration, total, prefix='', suffix='', decimals=2, length=100, fill='█', printEnd="\r"):
    global maxCurIter
    global progSign
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    if maxCurIter < iteration:
        maxCurIter = iteration
    if progSign == "-":
        progSign = "\\"
    elif progSign == "\\":
        progSign = "|"
    elif progSign == "|":
        progSign = "/"
    elif progSign == "/":
        progSign = "-"
    percent = ("{0:." + str(decimals) + "f}").format(100 *
                                                     (maxCurIter / float(total)))
    filledLength = int(length * maxCurIter // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    xOfy = str(maxCurIter)+"/"+str(total)
    empty = " "
    print('\r%s |%s| %s%% %s %s %s %s' %
          (prefix, bar, percent, suffix, progSign, xOfy, empty), end=printEnd)
    # Print New Line on Complete
    if maxCurIter == total:
        print()


def calcProtLigInteract(babelPath, PDBQTpath, i, index):
    #global intercfilesDic
    ZINCIDX = getZINCID(PDBQTpath, index+1)
    #print("ZINCIDX", str(ZINCIDX))
    tmpComplexPDBQT_targetPath = dataFolder
    filename = "complex"+ZINCIDX+"_"+str(i)+".pdb"

    outFilename = ZINCIDX+"_"+str(i)+".inter"
    outFilepath: Path = tmpComplexPDBQT_targetPath / (outFilename+".xml")
    # intercfilesDic.append(outFilepath)
    if (verbose):
        print("--------"+str(outFilepath))
    silentOutput = ">/dev/null 2>&1"
    if verbose:
        print("\n determing interactions...\n")
        silentOutput = ""
    printProgressBar(index, len(PDBQTpathList),
                     prefix='Interaction Force calc... Progress:', suffix='Complete', length=50)

    command = "cd "+str(tmpComplexPDBQT_targetPath)+"; bash -c '"+dockerPLIP + \
        "-f " + filename + " --breakcomposite -txp --name " + \
        outFilename+" "+silentOutput+"'"
    #command = "cd "+str(tmpComplexPDBQT_targetPath)+"; bash -c '"+dockerPLIP +"-f "+ filename +" --maxthreads 1 -txp --name "+outFilename+"'"

    if (verbose):
        print("PLIP: ", command)
    time.sleep(1)
    if (curOS == "Linux"):
        command = "cd "+str(tmpComplexPDBQT_targetPath)+" && "+dockerPLIP + \
            "-f " + filename + " --breakcomposite -txp --name "+outFilename
        os.system(command)
    else:
        completed = subprocess.Popen(["powershell", "-Command", command])
        completed.wait()

    drive = str(tmpComplexPDBQT_targetPath.drive)

    zPath = tmpComplexPDBQT_targetPath
    search = "COMPLEX"+ZINCIDX.upper()+"_"+str(i)
    search = search.replace("-", "_")
    pngs = processFolder(zPath, search)
    if (len(pngs) == 0 and babelPath != "oho"):
        calcProtLigInteract("oho", PDBQTpath, i, index)
        return
    delCommand = str(drive)+" && cd "+str(tmpComplexPDBQT_targetPath) + \
        " && del *"+"COMPLEX"+ZINCIDX+"_"+str(i)+"*.pdb"
    if cleaning:
        os.system(delCommand)

    for png in pngs:
        if (type(png) is str):
            if (("_LIG_" in str(png) or "_PROTEIN_" in str(png)) and not ("_PROTEIN_HIS") in str(png) and not ("_PROTEIN_HIE") in str(png)):
                if verbose:
                    print(str(tmpComplexPDBQT_targetPath))
                path1 = Path(tmpComplexPDBQT_targetPath) / str(png)
                path2 = Path(tmpComplexPDBQT_targetPath) / \
                    ("COMPLEX"+ZINCIDX+"_"+str(i)+".png")
                try:
                    #os.rename(path1, path2)
                    if verbose:
                        print("path1: ", path1)
                        print("path2: ", path2)
                    Path.rename(path1, path2)

                except IOError:
                    print(traceback.print_exc())
            else:
                delPath: Path = tmpComplexPDBQT_targetPath / png
                if verbose:
                    print("del: "+str(delPath))
                if cleaning:
                    try:
                        Path.unlink(delPath)
                    except IOError:
                        print(str(IOError))


def getModelCnt(path):
    # Opening file
    PDBQTfile = open(path, 'r')
    count = 0
    maxModel = 1
    for line in PDBQTfile:
        count += 1
        if line.startswith("MODEL"):
            maxModel = int(line.split(" ")[-1])

    # Closing files
    PDBQTfile.close()
    return maxModel


def getPDBQTpaths(i, PDBQTpathList, dataFolder):
    baseDir = getBasedir(multiPDBQTListPath)
    isCalc = False
    #path = PDBQTpathList[i].split("./")[-1]
    path = baseDir / PDBQTpathList[i]
    if (path != "" and str(path).endswith(".pdbqt")):
        ZINCID = getZINCID(path, i+1)
        maxModel = getModelCnt(path)
        # TODO range should be 0 = maxModel
        for j in range(1, maxModel):
            tmpComplexPDBQT_targetPath = dataFolder / \
                (str(ZINCID)+"_"+str(j)+".inter.txt")
            if ((os.path.exists(tmpComplexPDBQT_targetPath))):
                isCalc
            else:
                if j > 1:
                    if verbose:
                        print("missing: ", str(tmpComplexPDBQT_targetPath))
                    isCalc = True
                    break
        percentage = str(round(((i+1)/len(PDBQTpathList))*100, 1))
        if (not isCalc):
            if verbose:
                print(str(tmpComplexPDBQT_targetPath)+" "+percentage+"%")
        else:
            if verbose:
                print(str(tmpComplexPDBQT_targetPath)+" "+str(i+1) +
                      "/"+str(len(PDBQTpathList))+" "+percentage+"%")
            # if(os.path.exists(tmpComplexPDBQT_targetPath) == False):
            if (isCalc and j > 1):
                if verbose:
                    print("rcalc: ", str(path))
            convertMultiPDBQT_to_singleProtLigComplexes(babelPath, path, i)
    else:
        print("PDBQT path was empty: entry:"+str(i)+" path:"+str(path))


########################################################
########################################################
########################################################

def getRealProtein_atomPositions():
    global realProtAtomPositions
    global protAtomCnt
    prot_withH = PandasPdb().read_pdb(str(proteinPath))
    realProtAtomPositions["empty"] = "empty"

    for index, row in prot_withH.df['ATOM'].iterrows():
       # print(row['element_symbol'])
        if (row['element_symbol'] != "H"):
            realProtAtomPositions[str(len(realProtAtomPositions))] = (index+1)

    protAtomCnt = len(realProtAtomPositions)-1
    if (verbose):
        print("protAtomCnt: ", str(protAtomCnt))


def getLigPos(child: ET.Element):
    ligPos = []
    ligPos.append(child.find("ligcoo").find("x").text)
    ligPos.append(child.find("ligcoo").find("y").text)
    ligPos.append(child.find("ligcoo").find("z").text)
    return ligPos


def getProtPos(child: ET.Element):
    protPos = []
    protPos.append(child.find("protcoo").find("x").text)
    protPos.append(child.find("protcoo").find("y").text)
    protPos.append(child.find("protcoo").find("z").text)
    return protPos


def getProtIDxList(child: ET.Element):
    global realProtAtomPositions
    protIDxlist = []
    for subChild in child.find("prot_idx_list"):
        protIDxlist.append(realProtAtomPositions[subChild.text])
    return protIDxlist


def getLigIDxList(child: ET.Element):
    global protAtomCnt
    ligIDxlist = []
    for subChild in child.find("lig_idx_list"):
        ligIDxlist.append(int(subChild.text)-protAtomCnt)
    return ligIDxlist

# convert list to semicolon seperated string element


def listSemSep(list: list):
    buffer = ""
    for i, ele in enumerate(list):
        buffer += str(ele)
        if (i != len(list)-1):
            buffer += ";"
    return buffer


def parse_hydrophobicInteractions(leaf: ET.Element):
    global realProtAtomPositions
    global protAtomCnt
    global fileContent
    for child in leaf:
        ligPos = getLigPos(child)
        protPos = getProtPos(child)
        protIDx = realProtAtomPositions[child.find("protcarbonidx").text]
        ligIDx = int(child.find("ligcarbonidx").text)-int(protAtomCnt)
        if (verbose):
            print("protIDx: ", protIDx, "\t pos: ", protPos)
            print("ligIDx: ", str(ligIDx), "\t pos: ", ligPos)
        fileLine = "hydrophobicInteraction;protIDx;"+str(protIDx)+";ligIDx;"+str(
            ligIDx)+";protPos;"+listSemSep(protPos)+";ligPos;"+listSemSep(ligPos)
        fileContent.append(fileLine)

# not in use


def parse_waterBridges(leaf: ET.Element):
    global realProtAtomPositions
    global protAtomCnt
    global fileContent
    for child in leaf:
        ligPos = getLigPos(child)
        protPos = getProtPos(child)
        if (verbose):
            print("protIDx: ",
                  realProtAtomPositions[child.find("protcarbonidx").text])
        #fileLine = "halogenBond;protIDx;"+str(protIDx)+";ligIDx;"+str(ligIDx)+";protPos;"+listSemSep(protPos)+";ligPos;"+listSemSep(ligPos)
        # fileContent.append(fileLine)


def parse_saltBridges(leaf: ET.Element):
    global realProtAtomPositions
    global protAtomCnt
    global fileContent
    for child in leaf:
        ligPos = getLigPos(child)
        protPos = getProtPos(child)
        ligIDxlist = getLigIDxList(child)
        protIDxlist = getProtIDxList(child)
        if (verbose):
            print("protIDx: ", protIDxlist, "\t pos: ", protPos)
            print("ligIDx: ", ligIDxlist, "\t pos: ", ligPos)
        fileLine = "saltBridge;protIDxList;"+listSemSep(protIDxlist)+";ligIDxList;"+listSemSep(
            ligIDxlist)+";protPos;"+listSemSep(protPos)+";ligPos;"+listSemSep(ligPos)
        fileContent.append(fileLine)


def parse_piStacks(leaf: ET.Element):
    global realProtAtomPositions
    global protAtomCnt
    global fileContent
    for child in leaf:
        ligPos = getLigPos(child)
        protPos = getProtPos(child)
        ligIDxlist = getLigIDxList(child)
        protIDxlist = getProtIDxList(child)
        if (verbose):
            print("protIDx: ", protIDxlist, "\t pos: ", protPos)
            print("ligIDx: ", ligIDxlist, "\t pos: ", ligPos)
        fileLine = "piStack;protIDxList;"+listSemSep(protIDxlist)+";ligIDxList;"+listSemSep(
            ligIDxlist)+";protPos;"+listSemSep(protPos)+";ligPos;"+listSemSep(ligPos)
        fileContent.append(fileLine)


def parse_piCationInteractions(leaf: ET.Element):
    global realProtAtomPositions
    global protAtomCnt
    global fileContent
    for child in leaf:
        ligPos = getLigPos(child)
        protPos = getProtPos(child)
        ligIDxlist = getLigIDxList(child)
        protIDxlist = getProtIDxList(child)
        if (verbose):
            print("protIDx: ", protIDxlist, "\t pos: ", protPos)
            print("ligIDx: ", ligIDxlist, "\t pos: ", ligPos)
        fileLine = "piCationInteraction;protIDxList;"+listSemSep(protIDxlist)+";ligIDxList;"+listSemSep(
            ligIDxlist)+";protPos;"+listSemSep(protPos)+";ligPos;"+listSemSep(ligPos)
        fileContent.append(fileLine)


def parse_halogenBonds(leaf: ET.Element):
    global realProtAtomPositions
    global protAtomCnt
    global fileContent
    for child in leaf:
        ligPos = getLigPos(child)
        protPos = getProtPos(child)
        ligIDx = -1
        protIDx = -1
        if (verbose):
            print(child)
        idCheck = int(child.find("don_idx").text)-protAtomCnt
        if (idCheck > 0):
            ligIDx = idCheck
        else:
            protIDx = child.find("don_idx").text
        idCheck = int(child.find("acc_idx").text)-protAtomCnt
        if (idCheck > 0):
            ligIDx = idCheck
        else:
            protIDx = child.find("acc_idx").text
        if int(protIDx) != -1:
            protIDx = realProtAtomPositions[protIDx]
        #ligIDx = ligIDx-protAtomCnt
            fileLine = "halogenBond;protIDx;"+str(protIDx)+";ligIDx;"+str(
                ligIDx)+";protPos;"+listSemSep(protPos)+";ligPos;"+listSemSep(ligPos)
            fileContent.append(fileLine)
        if (verbose):
            print("protIDx: ", protIDx, "\t pos: ", protPos)
            print("ligIDx: ", ligIDx, "\t pos: ", ligPos)


def parse_HBonds(leaf: ET.Element):
    global realProtAtomPositions
    global protAtomCnt
    global fileContent
    donorIDx = -1
    ligIDx = -1
    protIDx = -1

    for child in leaf:
        ligPos = getLigPos(child)
        protPos = getProtPos(child)
        idCheck = int(child.find("donoridx").text)-protAtomCnt
        if (idCheck > 0):
            ligIDx = idCheck
        else:
            protIDx = child.find("donoridx").text
        idCheck = int(child.find("acceptoridx").text)-protAtomCnt
        if (idCheck > 0):
            ligIDx = idCheck
        else:
            protIDx = child.find("acceptoridx").text
        if int(protIDx) != -1:
            protisdon = child.find("protisdon").text
            protIDx = realProtAtomPositions[protIDx]
            #ligIDx = ligIDx-protAtomCnt
            fileLine = "HBond;protisdon;"+protisdon+";protIDx;"+str(protIDx)+";ligIDx;"+str(
                ligIDx)+";protPos;"+listSemSep(protPos)+";ligPos;"+listSemSep(ligPos)
            fileContent.append(fileLine)


def parse_metalComplexes(leaf: ET.Element):
    global realProtAtomPositions
    global protAtomCnt
    global fileContent
    for child in leaf:
        metalPos = []
        metalPos.append(child.find("metalcoo").find("x").text)
        metalPos.append(child.find("metalcoo").find("y").text)
        metalPos.append(child.find("metalcoo").find("z").text)
        targetPos = []
        targetPos.append(child.find("targetcoo").find("x").text)
        targetPos.append(child.find("targetcoo").find("y").text)
        targetPos.append(child.find("targetcoo").find("z").text)

        metalIDx = (int(child.find("metal_idx").text)-protAtomCnt, realProtAtomPositions[child.find(
            "metal_idx").text])[int(child.find("metal_idx").text)-protAtomCnt > 0]
        targetIDx = (int(child.find("target_idx").text)-protAtomCnt, realProtAtomPositions[child.find(
            "target_idx").text])[int(child.find("target_idx").text)-protAtomCnt > 0]

        if (verbose):
            print("metalIDx: ", metalIDx, "\t metalPos: ", metalPos)
            print("targetIDx: ", targetIDx, "\t targetPos: ", targetPos)
        fileLine = "metalComplex;metalIDx;"+str(metalIDx)+";targetIDx;"+str(
            targetIDx)+";metalPos;"+listSemSep(metalPos)+";targetPos;"+listSemSep(targetPos)
        fileContent.append(fileLine)


def writeLigPropsJSON(child, path: Path, dataFolder):
    mwt = child.find("lig_properties").find("molweight").text
    logP = child.find("lig_properties").find("logp").text
    fsp3 = -1.0
    ligName = str(path.stem).split("_")[0]
    jsonPath: Path = dataFolder / (ligName+".json")
    data = '{"zinc_id": "'+str(ligName)+'", "logp": '+str(logP) + \
        ', "fractioncsp3": '+str(fsp3)+', "mwt": '+str(mwt)+'}'
    f = open(str(jsonPath), "w")
    f.write(data)
    f.close()


def readXMLfile(dataFolder):
    global fileContent
    for r, d, f in os.walk(dataFolder):
        for file in f:
            if file.endswith(".xml"):
                path = ""
                fileContent = []
                print(file)
                path: Path = dataFolder / file
                interactionType = ['hydrophobic_interactions', 'hydrogen_bonds', 'water_bridges',
                                   'salt_bridges', 'pi_stacks', 'pi_cation_interactions', 'halogen_bonds', 'metal_complexes']
                if (len(sys.argv) <= 1):
                    path = "C:/Users/Marco/Desktop/examplePILPxml/example.xml"
                if (verbose):
                    print(str(path))
                try:
                    root = ET.parse(str(path)).getroot()
                except:
                    print("\n\n ERROR: unable to find "+str(path) + "\n\n")
                    time.sleep(5)
                    return False
                if (verbose):
                    print("\n\n\n\n\n\n\n\n\n\n")
                for type in interactionType:
                    if (verbose):
                        print("\n--------"+type+"---------")
                    intercLeaf = root.findall("bindingsite")
                    for child in intercLeaf:
                        if (child.find("identifiers").find("longname").text != ""):
                            writeLigPropsJSON(child, path, dataFolder)
                            writeLigPropsJSON(child, path, dataFolder)
                            intercLeaf = child.find("interactions").find(type)
                            if (len(intercLeaf) > 0):
                                #print("found interaction of type: ",str(type))
                                if (type == "hydrophobic_interactions"):
                                    parse_hydrophobicInteractions(intercLeaf)
                                elif (type == "hydrogen_bonds"):
                                    parse_HBonds(intercLeaf)
                                elif (type == "water_bridges"):
                                    # should not be present
                                    # parse_waterBridges(intercLeaf)
                                    Text = 0
                                elif (type == "salt_bridges"):
                                    parse_saltBridges(intercLeaf)
                                elif (type == "pi_stacks"):
                                    parse_piStacks(intercLeaf)
                                elif (type == "pi_cation_interactions"):
                                    parse_piCationInteractions(intercLeaf)
                                elif (type == "halogen_bonds"):
                                    parse_halogenBonds(intercLeaf)
                                elif (type == "metal_complexes"):
                                    parse_metalComplexes(intercLeaf)
                            else:
                                if (verbose):
                                    print("no interaction of type: ", str(type))

                newOutFilePath = path.parent / (path.stem+".plip")
                with open(str(newOutFilePath), 'w') as f:
                    for ele in fileContent:
                        if (verbose):
                            print(ele)
                        f.write(ele)
                        f.write("\n")
                if cleaning:
                    try:
                        Path.unlink(path)
                    except IOError:
                        print(str(IOError), str(path))


def writeDoneFile(file):
    f = open(str(file), "w")
    f.write("True")
    f.close()

########################################################
########################################################
########################################################


if (len(sys.argv) > 1):
    if (os.path.exists(dataFolder) == False):
        Path.mkdir(dataFolder)

    # convert protein to pdb

command = str(babelPath) + " " + str(rawProteinPath) + \
    " -O " + str(proteinPath)
if (curOS != "Linux"):
    command = "bash -c '"+str("obabel")+" " + toLinuxSubPath(
        rawProteinPath) + " -O " + toLinuxSubPath(proteinPath)+"'"
os.system(command)

#getPLIPdocker = "docker pull pharmai/plip"
# os.system(getPLIPdocker)


# run the parallel program
print("babelPath = ", str(babelPath))
print("multiPDBQTListPath = ", str(multiPDBQTListPath))
print("dataFolder = ", str(dataFolder))
print("rawProteinPath = ", str(rawProteinPath))
getRealProtein_atomPositions()

#results = Parallel(n_jobs=8, prefer="threads")(delayed(getPDBQTpaths)(i,PDBQTpathList, dataFolder) for i in range(0, len(PDBQTpathList)))
print("Convert XMLs to .plip")
doneFile = dataFolder / Path("interactionDone.txt")

if (os.path.exists(doneFile) == False):
    readXMLfile(dataFolder)
    #results = Parallel(n_jobs=1, prefer="threads")(delayed(getPDBQTpaths)(i,PDBQTpathList, dataFolder) for i in range(0, len(PDBQTpathList)))
    # multiprocessing.cpu_count()
    results = Parallel(n_jobs=int(multiprocessing.cpu_count()/2), prefer="threads")(delayed(getPDBQTpaths)(
        i, PDBQTpathList, dataFolder) for i in range(0, len(PDBQTpathList)))
    results = Parallel(n_jobs=int(multiprocessing.cpu_count()/2), prefer="threads")(delayed(getPDBQTpaths)(
        i, PDBQTpathList, dataFolder) for i in range(0, len(PDBQTpathList)))
    print("Convert XMLs to .plip")
    readXMLfile(dataFolder)
    writeDoneFile(doneFile)
else:
    print("\n interaction forces already calculated!\n")

#results = Parallel(n_jobs=1, prefer="threads")(delayed(readXMLfile)(XMLpaths.get()) for i in range(0, XMLpaths.qsize()))
