# @ author: Marco Schaefer
# @ Description:

import multiprocessing
from tabnanny import verbose
from joblib import Parallel, delayed
import os
from pickle import FALSE
import sys
import time
from pathlib import Path
import platform
curOS = platform.system()
verbose = False
silentOutput = ">/dev/null 2>&1"

# Linux
# Windows


# import and install non default modules
#packageNames = ['wget','joblib']
#importCommands =  ['wget', 'joblib']
installCommad = ''
cnt = 0

# for lib in importCommands:
#    try:
#        exec("import {module}".format(module=lib))
#    except Exception as e:
#        print(e)
#        installCommad = "pip3 install -U "+packageNames[cnt]
#        os.system(installCommad)
#        time.sleep(2)
#    cnt+=1
#########################################################################

if (len(sys.argv) > 1):
    multiPDBQTListPath: Path = Path(sys.argv[1])
    dataFolder: Path = Path(sys.argv[2])
    checkmolPath: Path = Path(sys.argv[3])
else:
    multiPDBQTListPath: Path = Path(
        "C:/Users/MarcoSchaefer/Projects/results_70/results_list.txt")
    dataFolder: Path = Path(
        "C:\\Users\\MarcoSchaefer\\Projects\\results_70\\InVADo_tmpFiles")
    checkmolPath: Path = Path(
        "C:\\Users\\MarcoSchaefer\\Projects\\megamol-prolint_2022Studio\\plugins\\prolint\\checkmol.exe")


if (curOS == "Linux"):
    babelPath = "obabel"
    checkmolPath = checkmolPath.parent / "checkmol_linux"
babelPath = "obabel"

multiPDBQTList = open(str(multiPDBQTListPath), "r")
if verbose:
    print("multiPDBQTList "+str(multiPDBQTList))
multiPDBQTList = multiPDBQTList.read()

#print("PDBQTListPath: "+str(multiPDBQTListPath))
PDBQTpathList = multiPDBQTList.split("\n")

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
    empty = ""
    print('\r%s |%s| %s%% %s %s %s %s' % (prefix, bar, percent, suffix, progSign, xOfy, empty), end=printEnd)
    # Print New Line on Complete
    if maxCurIter == total:
        print()

def getBasedir(multiPDBQTListPath: Path):
    return multiPDBQTListPath.parent


def toLinuxSubPath(path: Path):
    drive = path.drive
    driveStr = str(drive)
    driveStrlow = driveStr[0].lower()
    return str(Path(path).as_posix()).replace(driveStr, "/mnt/"+driveStrlow)


def getZINCID(PDBQTpath, PDBQTpathIndex):
    third_line = ""
    try:
        with open(str(PDBQTpath)) as f:
            f.readline()
            f.readline()
            third_line = f.readline().strip()
    except:
        print("ERROR: Could not load "+str(PDBQTpath))
        return "False"

    splitRes = third_line.split("Name = ")
    ZINCID = ""
    if(len(splitRes) > 1):
        ZINCID = third_line.split("Name = ")[-1]
    else:
        ZINCID = "lig"+str(PDBQTpathIndex)
    if(ZINCID != ""):
        return ZINCID
    else:
        print("ERROR: Cloud not determine the ZINCID or LigName!")
        return False


def convertPDBQT_to_MOL(babelPath, PDBQTpath, MOLtargetPath):
    # print(MOLtargetPath)
    command = str(babelPath) + " " + str(PDBQTpath) + \
        " -O " + str(MOLtargetPath)
    if (curOS != "Linux"):
        command = "bash -c 'eval \""+str("obabel") + " " + toLinuxSubPath(
            PDBQTpath) + " -O " + toLinuxSubPath(MOLtargetPath)+" "+silentOutput+"\"'"
    if verbose:
        print(command)
    # time.sleep(1)
    os.system(command)


def getPDBQTpaths(i, PDBQTpathList, dataFolder):
    baseDir = getBasedir(multiPDBQTListPath)
    printProgressBar(i+1, len(PDBQTpathList),
                     prefix='functional group calc... Progress:', suffix='Complete', length=50)
    path = baseDir / PDBQTpathList[i]
    ZINCID = getZINCID(path, i+1)
    if(ZINCID == FALSE):
        return FALSE
    #print(f"path: {path}")
    if (path != ""):

        MOLtargetPath = dataFolder / (ZINCID+".mol")
        FunctionalGroups_targetPath = dataFolder / (ZINCID+".checkmol")
        percentage = str(round(((i+1)/len(PDBQTpathList))*100, 1))
        if((os.path.exists(str(MOLtargetPath)) and os.path.exists(str(FunctionalGroups_targetPath)))):
            if verbose:
                print(str(MOLtargetPath)+" "+percentage+"%")
        else:
            if verbose:
                print(str(MOLtargetPath)+" "+percentage+"%")
            if(os.path.exists(str(MOLtargetPath)) == False):
                convertPDBQT_to_MOL(babelPath, path, MOLtargetPath)
            if(os.path.exists(str(FunctionalGroups_targetPath)) == False):
                create_functionalGroups_file(
                    checkmolPath, MOLtargetPath, FunctionalGroups_targetPath)
    else:
        print("PDBQT path was empty: entry:"+str(i)+" path:"+str(path))


def create_functionalGroups_file(checkmolPath, MOLPath, FunctionalGroups_targetPath):
    f = open(FunctionalGroups_targetPath, "w")
    command = str(checkmolPath)+" -p " + str(MOLPath)
    # print(command)
    stream = os.popen(command)
    tmp = stream.readline()

    while tmp != "":
        f.write(tmp)
        tmp = stream.readline()
    f.close()


if(os.path.exists(dataFolder) == False):
    Path.mkdir(dataFolder)

results = Parallel(n_jobs=multiprocessing.cpu_count(), prefer="threads")(delayed(getPDBQTpaths)(
    i, PDBQTpathList, dataFolder) for i in range(0, len(PDBQTpathList)))

# downloadSMIFiles(ZINCID)
