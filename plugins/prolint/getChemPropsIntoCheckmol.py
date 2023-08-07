# @ author: Marco Schaefer
# @ Description:

from tabnanny import verbose
import wget
from joblib import Parallel, delayed
import os
import sys
import time
import json
from pathlib import Path
import platform
import requests
import re
import multiprocessing
verbose = False

curOS = platform.system()
# Linux
# Windows


# import and install non default modules
# packageNames = ['wget','joblib']
# importCommands =  ['wget', 'joblib']
installCommad = ""
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
ZINCListPath: Path = Path(sys.argv[4])

if curOS == "Linux":
    babelPath = "obabel"

ZINCList = open(ZINCListPath, "r")
ZINCList = ZINCList.read()

if verbose:
    print("PDBQTListPath: " + str(ZINCList))
ZINCListNames = ZINCList.split(";")
if verbose:
    print("ZINCListNames: " + str(ZINCListNames))

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

def add_ChemProps_To_Checkmol(JSON_data, ChemPropTargetPath):
    if JSON_data != "":
        lineToAdd = f"#Properties;mwt:{JSON_data['mwt']};logp:{JSON_data['logp']};fractioncsp3:{JSON_data['fractioncsp3']}"
        if verbose:
            print("ChemPropTargetPath: ", ChemPropTargetPath)
        # print("lineToAdd: ", lineToAdd)
        with open(ChemPropTargetPath, "a") as file:
            file.write(lineToAdd)


def toLinuxSubPath(path: Path):
    drive = path.drive
    driveStr = str(drive)
    driveStrlow = driveStr[0].lower()
    return str(Path(path).as_posix()).replace(driveStr, "/mnt/"+driveStrlow)


def add_ChemProps_To_Checkmol_noZINCID(JSON_data, ChemPropTargetPath):
    lineToAdd = f"#Propteries;mwt:10;logp:1;fractioncsp3:1"
    print("ChemPropTargetPath: ", ChemPropTargetPath)
    try:
        with open(ChemPropTargetPath, "a") as file:
            file.write(lineToAdd)
    except IOError:
        print(str(IOError))


def readJSON(JSONpath):
    data = ""
    try:
        with open(JSONpath) as f:
            data = json.load(f)
    except IOError:
        print(str(IOError))
        return ""
    return data


def downloadJSONfiles(i, ZINCListNames):
    printProgressBar(i+1, len(ZINCListNames),
                     prefix='get chemical properties... Progress:', suffix='Complete', length=50)
    ligName = ZINCListNames[i]
    if verbose:
        print("LigName: " + str(ligName))
    if str(ligName).startswith("ZINC"):
        JSONpath: Path = dataFolder / (ligName + ".json")
        ChemPropTargetPath = dataFolder / (ligName + ".checkmol")
        getJSON: str = (
            "http://zinc15.docking.org/substances/"
            + str(ligName).split("-")[0]
            + ".json:zinc_id+logp+fractioncsp3+mwt"
        )
        percentage = str(round(((i + 1) / len(ZINCListNames)) * 100, 1))
        success = 0
        if os.path.exists(JSONpath):
            # print(str(JSONpath)+" "+percentage+"%")
            Path.unlink(JSONpath)
        if verbose:
            print("wget " + str(JSONpath) + " " + percentage + "%")
        while success != 1:
            JSON_data = ""
            try:
                # wget.download(getJSON, str(JSONpath))
                JSON_data = requests.get(getJSON, stream=True).content
                # print("JSON_data", JSON_data.decode("utf-8"))
                if JSON_data != "":
                    JSON_data = re.findall("{.*}", JSON_data.decode("utf-8"))
                    if(len(JSON_data) > 0):
                        JSON_data = json.loads(JSON_data[0])
                    else:
                        JSON_data = json.loads(
                            '{"zinc_id": "none", "logp": 1.0, "fractioncsp3": 1.0, "mwt": 1.0}')
                success = 1
            except OSError:
                print(OSError)
                print("Something went wrong: " + getJSON)
                print("JSONpath: " + str(JSONpath))

        add_ChemProps_To_Checkmol(JSON_data, ChemPropTargetPath)
    elif ligName != "":
        JSONpath = dataFolder / (ligName + ".json")
        ChemPropTargetPath = dataFolder / (ligName + ".checkmol")

        JSON_data = readJSON(JSONpath)
        add_ChemProps_To_Checkmol(JSON_data, ChemPropTargetPath)
        if verbose:
            print(str(JSON_data), str(ChemPropTargetPath))
    else:
        print("No ligName for element: " + str(i) + " " + str(ligName))


results = Parallel(n_jobs=multiprocessing.cpu_count(), prefer="threads")(
    delayed(downloadJSONfiles)(i, ZINCListNames) for i in range(0, len(ZINCListNames))
)

# downloadSMIFiles(ZINCID)
