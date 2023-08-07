# @ author: Marco Schaefer
# @ Description:

from tabnanny import verbose
import wget
from joblib import Parallel, delayed
import os
import sys
import time
from pathlib import Path
import platform
import multiprocessing
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
#        installCommad = "pip3 install "+packageNames[cnt]
#        os.system(installCommad)
#        time.sleep(2)
#    cnt+=1
#########################################################################
babelPath: Path = Path(sys.argv[1])
baseDir: Path = Path(sys.argv[2])
dataFolder: Path = Path(sys.argv[3])
ligListPath: Path = Path(sys.argv[4])
multiPDBQTListPath: Path = Path(sys.argv[5])

if(curOS == "Linux"):
    babelPath = "obabel"

ligList = open(str(ligListPath), "r")
ligList = ligList.read()

multiPDBQTList = open(str(multiPDBQTListPath), "r")
multiPDBQTList = multiPDBQTList.read()
PDBQTpathList = multiPDBQTList.split("\n")

if verbose:
    print("PDBQTListPath: "+str(ligList))
listOfLigNames = ligList.split(";")
if verbose:
    print("ZINCListNames: "+str(listOfLigNames))

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

def convertSMI2SVG(babelPath, SMIpath, SVGtargetPath):
    command = str(babelPath) + " " + str(SMIpath) + " -O " + str(SVGtargetPath)
    if (curOS != "Linux"):
        command = "bash -c 'eval \""+str("obabel") + " " + toLinuxSubPath(
            SMIpath) + " -O " + toLinuxSubPath(SVGtargetPath)+" "+silentOutput+"\"'"
    os.system(command)


def toLinuxSubPath(path: Path):
    drive = path.drive
    driveStr = str(drive)
    driveStrlow = driveStr[0].lower()
    return str(Path(path).as_posix()).replace(driveStr, "/mnt/"+driveStrlow)


def convertPDBQT2SVG(babelPath, PDBQTpath, SVGtargetPath):
    command = str(babelPath) + " " + str(PDBQTpath) + \
        " --unique -O " + str(SVGtargetPath)
    if (curOS != "Linux"):
        command = "bash -c '"+str("obabel") + " " + toLinuxSubPath(
            PDBQTpath) + " --unique -O " + toLinuxSubPath(SVGtargetPath)+"'"
    if verbose:
        print(str(command))
    os.system(command)


def getBasedir(multiPDBQTListPath):
    return multiPDBQTListPath.parent


def downloadSMIFiles(i, listOfLigNames):
    printProgressBar(i+1, len(listOfLigNames),
                     prefix='get smi files + convert svg... Progress:', suffix='Complete', length=50)
    ligName = listOfLigNames[i]

    # if ligName contains a ZINCID
    if (str(ligName).startswith("ZINC")):
        SMIpath = dataFolder / (ligName+".smi")
        SVGtargetPath = dataFolder / (ligName + ".svg")
        getSMI = "http://zinc15.docking.org/substances/"+str(ligName).split("-")[0]+".smi"
        percentage = str(round(((i+1)/len(listOfLigNames))*100, 1))
        success = 0
        if(os.path.exists(str(SMIpath))):
            if verbose:
                print(str(SMIpath)+" "+percentage+"%")
        else:
            if verbose:
                print(str(SMIpath)+" "+percentage+"%")
            while(success != 1):
                try:
                    wget.download(getSMI, str(SMIpath))
                    success = 1
                except OSError:
                    print("Something went wrong: "+str(getSMI))
                    print("getPDB: "+str(getSMI))
                    print("SMIpath: "+str(SMIpath))
        convertSMI2SVG(babelPath, SMIpath, SVGtargetPath)
    # if the ligName is not a ZINCID
    elif(ligName != ""):
        baseDir = getBasedir(multiPDBQTListPath)
        PDBQTpath = baseDir / PDBQTpathList[i]
        SVGtargetPath = dataFolder / (ligName + ".svg")
        convertPDBQT2SVG(babelPath, PDBQTpath, SVGtargetPath)
    else:
        print("No ligName for element: "+str(i)+" "+str(ligName))


results = Parallel(n_jobs=multiprocessing.cpu_count(), prefer="threads")(delayed(downloadSMIFiles)(
    i, listOfLigNames) for i in range(0, len(listOfLigNames)))

# downloadSMIFiles(ZINCID)
