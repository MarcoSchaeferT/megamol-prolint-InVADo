#!/usr/bin/env python3
import base64
import io
from PIL import Image
import re
import lxml.etree as ET
import queue
import _thread
from threading import Thread
import struct
import socket
import json
import numpy as np
import pandas as pd
import urllib.request as urllib2
from selenium import webdriver
from flask.cli import FlaskGroup, pass_script_info
import config
import click
import flask
import pathlib
import os
import sys
import time
from screeninfo import get_monitors

numberScreens = 1
numberScreens = len(get_monitors())
print("numberScreens", numberScreens)

# if True just the test data from testData.py will be used with a normal flask server
# if False the data comes from MegaMol using a custom flask server and geckodriver for FireFox
testMode = True
if len(sys.argv) > 1:
    testMode = False
# print("argv:"+str(len(sys.argv))+"testMode: "+str(testMode))

# import and install non default modules
# packageNames = [    'flask==2.0.3', 'click', 'config',  'pandas',       'numpy==1.19.3',     'selenium', 'Pillow',   'lxml']
# importCommands =  [ 'flask',        'click', 'config',  'pandas as pd', 'numpy as np',      'selenium', 'PIL',      'lxml.etree as ET', 'io', 'base64']
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
# requiredPackages = pathlib.Path("../requirements.txt")
# installCommad = "pip install -r "+ str(requiredPackages)
# os.system(installCommad)


DATA_PATH = pathlib.Path().resolve()

verbose = True

######################################
############# TEST DATA ##############
######################################
global data_recv
global clusterData_recv
global funcGroupData_recv
global fgsMapData_recv
global GUIdata_recv
global selectionData_recv
global isIncomingData
global pollTime
global manager
global driver
global app
pollTime = 2.0

threadDataBuffer = queue.Queue()
pauseLifeSignal = queue.Queue()
sendDataBuffer = queue.Queue()
isIncomingData = queue.Queue()
isIncomingData.put("no")
pauseLifeSignal.put("no")

data_dict = {}
if not testMode:
    # configure browser driver
    driver = ""
    driver = webdriver.Firefox()
    if numberScreens >= 2:
        driver.set_window_position(4500, 0)
        driver.fullscreen_window()
    else:
        driver.set_window_position(700, 0)
        driver.set_window_size(2048, 1124)


##################################################
################## SOCKET-SERVER #################
##################################################


def encodeString(data):
    res = data.decode("utf-8")
    # if string is an array of strings seperated with ","
    # otherwise res will be non modified returned
    tmpSplit = res.split(",")
    if len(tmpSplit) > 1:
        res = []
        for element in tmpSplit:
            res.append(str(element))
    return res


def encodeFloat(data):
    res = np.frombuffer(data, np.float32)
    return res


def encodeInt(data):
    res = int.from_bytes(data, byteorder="little", signed=True)
    res = np.frombuffer(data, np.int32)
    return res


def encodeUint(data):
    res = int.from_bytes(data, byteorder="little", signed=False)
    res = np.frombuffer(data, np.int32)
    return res


globalRes = 1

# dataKeys for MegaMol (so that MegaMol knows what to do with that data)


class MegaMolDataKeys:
    clusterLigandModelID = -333
    GUIparam = -444
    gblMdlID = -777


global megaMolDataKeys
megaMolDataKeys = MegaMolDataKeys()


class ClusterDataKeys:
    ligandCount = 0
    modelCount = 1
    energies = 2
    ligandIDs = 3
    zincNames = 4
    clusterAssignment = 5
    clusterSizes = 6
    # (path where svg, smi, checkmol, etc. are located)
    dataPath = 7
    hbonds = 8
    efficiency = 9
    mwt = 10
    logp = 11
    fractioncsp3 = 12
    halogenBonds = 13
    hydrophobicInteractions = 14
    metalComplexes = 15
    piCationInteractions = 16
    piStacks = 17
    saltBridges = 18


global clusterDataKeys
clusterDataKeys = ClusterDataKeys()


class FuncGroupDataKeys:
    funcGroupsPerLigands = 0
    funcGroupsPerClusters = 1


global funcGroupDataKeys
funcGroupDataKeys = FuncGroupDataKeys()


class DataStreamNames:
    clusterData = "clusterData"
    funcGroupData = "functionalGroupData"
    GUIdata = "GUIdata"
    fgsMapData = "selFgsClustDat"
    GUIdata = "GUIdata"
    selectionData = "selectionData"


global streams
streams = DataStreamNames()

dataTypes = {
    "string": {"encode": encodeString, "typeFLAG": "%ss", "megamolCODE": 2},
    "int": {"encode": encodeInt, "typeFLAG": "%si", "megamolCODE": 0},
    "unsigned int": {"encode": encodeUint, "typeFLAG": "%si", "megamolCODE": 3},
    "float": {"encode": encodeFloat, "typeFLAG": "%sf", "megamolCODE": 1},
    "bufferSize": {
        "encode": encodeInt,
    },
    "lifeSignal": {"encode": encodeInt, "megamolCODE": -1},
}

old_dataTypeString = ""


def socketConnectionHandler(old_dataTypeString, dataTypes):
    global pauseLifeSignal
    global isIncomingData
    global pollTime

    result = ""
    old_dataTypeString = ""
    bufferSize = 15
    HOST = "127.0.0.1"  # Standard loopback interface address (localhost)
    PORT = 2040  # Port to listen on (non-privileged ports are > 1023)
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        print("listing at " + HOST + ":" + str(PORT) + "...")
        s.bind((HOST, PORT))
        s.listen(0)

        connection, addr = s.accept()

        with connection:
            print("Connected by", addr)
            # some sample data to perform handshake with socket servers
            data = [5, 36]
            handshake = struct.pack("%si" % len(data), *data)
            connection.sendall(handshake)
            # start a new thread which sends a life signal OR sends data
            _thread.start_new_thread(send_lifeSignalOrData, (connection,))
            while True:

                dataTypeString = ""
                data = connection.recv(bufferSize)
                isIncomingData.put("yes")
                pauseLifeSignal.put("yes")

                if not data:
                    break

                # conn.sendall(data) # just echo all
                # tries to decode incoming data as string (utf-8)
                try:
                    dataTypeString = data.decode("utf-8")
                except:
                    pass

                # checks if dataType is known (decoding defined)
                if dataTypeString in dataTypes.keys():
                    old_dataTypeString = dataTypeString
                    """ old_dataTypeString
                        is a status (if no new dataType is transmitted
                        everything will be decoded with old_dataTypeString)
                    """
                elif old_dataTypeString == "":
                    pass

                # set buffersize for incoming data stream (if old_dataTypeString = 'bufferSize')
                elif str(old_dataTypeString) == str("bufferSize"):
                    bufferSize = 15
                    print("old_dataTypeString:" + str(old_dataTypeString))
                    try:
                        if (
                            dataTypes[old_dataTypeString]["encode"](data)[0]
                            > bufferSize
                        ):
                            bufferSize = dataTypes[old_dataTypeString]["encode"](data)[
                                0
                            ]
                    except:
                        print("wrong endcoding")

                    print("bufferSize: ", str(bufferSize))
                # map the correct decode function to the incoming data stream
                else:
                    result = "eod"
                    try:
                        result = dataTypes[old_dataTypeString]["encode"](data)
                    except:
                        print(
                            "decoding went wrong for: ",
                            dataTypes[old_dataTypeString]["encode"],
                        )
                        isIncomingData.put("no")
                        pauseLifeSignal.put("no")
                    # check if end-of-data
                    # if((not isinstance(result, str)) or (isinstance(result, str) and result != 'eod')):
                    print("res: ", result)
                    if result != "reload" and result != "eod":
                        threadDataBuffer.put(result)
                    elif result == "eod":
                        print("eod")
                        get_data_dict()
                        send_FIN(connection)
                        # time.sleep(1.0)
                        isIncomingData.put("no")
                        pauseLifeSignal.put("no")
                    else:
                        get_data_dict()
                        send_FIN(connection)
                        # time.sleep(1.0)
                        isIncomingData.put("no")
                        pauseLifeSignal.put("no")
                        print("will reload page...")
                        _thread.start_new_thread(
                            reloadPage, (driver, "http://localhost:8050/")
                        )
                        pollTime = 0.2
                        continue

                    print("Received: " + str(result) +
                          "\t :as " + old_dataTypeString)
                # send a signal that the data precssing is done (server is ready for next one!)
                send_FIN(connection)

        return


def send_lifeSignalOrData(connection):
    global pauseLifeSignal
    global isIncomingData
    global pollTime

    maxWaitCount = 2 / pollTime
    pauseCnt = 0

    while True:
        check = "A"
        isIncoming = "A"
        while not pauseLifeSignal.empty():
            check = pauseLifeSignal.get()

        print("is lifeSignal paused: " + str(check))
        if check != "A":
            pauseLifeSignal.put(check)
        if check == "yes":
            pauseCnt += 1
        if pauseCnt >= maxWaitCount:
            pauseCnt = 0
            pauseLifeSignal.put("no")
            isIncomingData.put("no")

        while not isIncomingData.empty():
            isIncoming = isIncomingData.get()
        if isIncoming != "A":
            isIncomingData.put(isIncoming)
        time.sleep(pollTime)

        # sending lifeSignal
        if check == "no":
            print("lifeSignal mode...")
            try:
                send_Header("lifeSignal", 2, connection)
            except:
                # if MegaMol is closed close Firefox window
                global driver
                driver.close()
        # sending real data from 'sendDataBuffer'
        elif isIncoming == "no":
            sendSomething = False
            while not sendDataBuffer.empty():
                print("send data mode...")
                dataType = sendDataBuffer.get()
                bufferElement = sendDataBuffer.get()
                sendSomething = True
                print("type: " + str(dataType) +
                      " data: " + str(bufferElement))
                send_Data(bufferElement, dataType, connection)
            if sendSomething:
                # time.sleep(1.0)
                pauseLifeSignal.put("no")


def reloadPage(driver, url):
    time.sleep(0.1)
    print("calling webpage...")
    driver.get(url)  # tested in combination with scrapy
    # driver.refresh()
    return


def testModeFirefoxWindow(url):
    driver = webdriver.Firefox()
    if numberScreens >= 2:
        driver.set_window_position(4500, 0)
        driver.fullscreen_window()
    else:
        driver.set_window_position(700, 0)
        driver.set_window_size(2048, 1124)
    time.sleep(2.0)
    print("calling webpage...")
    driver.get(url)  # tested in combination with scrapy
    # driver.refresh()
    return


def sendToMegamol():
    global pauseLifeSignal
    # pausing lifeSignal leads to sending of data from 'sendDataBuffer'
    pauseLifeSignal.put("yes")


# *** send_FIN ***
# Info: sends an int value (1) to indicate that the request has been completely processed
# @param connection: the established socket connection
#
def send_FIN(connection):
    value = [-1, -1]
    byteValue = struct.pack(getTypeFLAG("int") % len(value), *value)
    connection.sendall(byteValue)


def getTypeFLAG(t):
    return dataTypes[t]["typeFLAG"]


def getMegamolCODE(t):
    return dataTypes[t]["megamolCODE"]


def send_Data(data, dataType, connection):
    print("try to send HEADER!")
    global threadDataBuffer
    byteValue2 = ""
    if dataType == "string":
        byteValue2 = data.encode()
    else:
        byteValue2 = struct.pack(getTypeFLAG(dataType) % len(data), *data)
    dataLength = len(data)
    send_Header(dataType, dataLength, connection)
    print("HEADER send!")
    # wait for megamol (processing header)
    # data = threadDataBuffer.get()
    # while data != "ready":
    #   data = threadDataBuffer.get()
    # time.sleep(3)
    connection.sendall(byteValue2)


def send_Header(dataType, dataLength, connection):
    header = [getMegamolCODE(dataType), dataLength]
    byteValue = struct.pack(getTypeFLAG("int") % len(header), *header)
    connection.sendall(byteValue)


def runSocket(param):
    cnt = 1
    while True:
        socketConnectionHandler(old_dataTypeString, dataTypes)
        print(str(cnt))
        cnt += 1


#################################################
################## FLASK-SERVER #################
#################################################


def get_data_dict():
    global threadDataBuffer
    global clusterData_recv
    global funcGroupData_recv
    global fgsMapData_recv
    global GUIdata_recv
    global selectionData_recv
    global streams
    global pollTime

    try:
        clusterData_recv
    except:
        clusterData_recv = {}

    try:
        funcGroupData_recv
    except:
        funcGroupData_recv = {}

    try:
        fgsMapData_recv
    except:
        fgsMapData_recv = {}

    try:
        GUIdata_recv
    except:
        GUIdata_recv = {}

    try:
        selectionData_recv
    except:
        selectionData_recv = {}

    clusterData_recvIndex = 0
    funcGroupData_recvIndex = 0
    fgsMapData_recvIndex = 0
    GUIdata_Index = 0
    selectionData_Index = 0

    currentDataStreamName = ""

    while not threadDataBuffer.empty():
        data = threadDataBuffer.get()
        threadDataBuffer.task_done()

        # convert data to list if necessary bcz ndarray is not recognised in JS
        if type(data) is np.ndarray:
            data = data.tolist()

        if data == streams.clusterData:
            pollTime = 2.0
            currentDataStreamName = data
            continue
        elif data == streams.funcGroupData:
            currentDataStreamName = data
            continue
        elif data == streams.fgsMapData:
            currentDataStreamName = data
            continue
        elif data == streams.GUIdata:
            currentDataStreamName = data
            continue
        elif data == streams.selectionData:
            currentDataStreamName = data
            continue

        if currentDataStreamName == streams.clusterData:
            clusterData_recv[clusterData_recvIndex] = data
            clusterData_recvIndex += 1
        if currentDataStreamName == streams.funcGroupData:
            funcGroupData_recv[funcGroupData_recvIndex] = data
            funcGroupData_recvIndex += 1
        if currentDataStreamName == streams.fgsMapData:
            fgsMapData_recv[fgsMapData_recvIndex] = data
            fgsMapData_recvIndex += 1
        if currentDataStreamName == streams.GUIdata:
            GUIdata_recv[GUIdata_Index] = data
            GUIdata_Index += 1
        if currentDataStreamName == streams.selectionData:
            selectionData_recv[selectionData_Index] = data
            selectionData_Index += 1

    return 1


"""
Default app routes for index and favicon
"""


def createAppRoutes():
    if verbose:
        print("createAppRoutes()")
    global app

    @app.route("/")
    def index():
        return flask.render_template("index.html")

    # currently not in use
    @app.route("/favicon.ico")
    def fav():
        print("/favicon.ico")
        b = pathlib.Path("vue_app/public/assets")
        imgPath = DATA_PATH / b
        return flask.send_from_directory(imgPath, "Pocket.svg")

    @app.route("/getClusterData", methods=["GET", "POST"])
    def cluster_stats():
        global clusterData_recv
        if verbose:
            print("clusterData_recv")
            print(clusterData_recv)

        # Datapackage
        # 	[0]	->	ligandCount
        # 	[1]	->	modelCount
        # 	[2]	->	energies
        # 	[3]	->	ligandIDs
        # 	[4]	->	zincNames
        # 	[5]	->	clusterAssignment
        #   [6]	->	clusterSizes
        # 	[7] ->  dataPath (path where svg, smi, checkmol, etc. are located)
        #   [8]	->	hbonds
        #   [9] ->  efficiency
        #   [10] ->  mwt (molecular weigth)
        # 	[11] ->  logp (logP is a quantitative measure of lipophilicit)
        # 	[12] ->  fractioncsp3 (parameter for drug-likeness)

        # print(json.dumps({'clusterData_recv': clusterData_recv}))
        return json.dumps({"clusterData_recv": clusterData_recv, "sort": 0})

    @app.route("/getFGSMapData", methods=["GET", "POST"])
    def getFGSMapData_recv():
        global fgsMapData_recv
        if verbose:
            print("fgsMapData_recv")
            print(fgsMapData_recv)
        # Datapackage (-1 delemiter)
        # 		[0]	 ->	 fgsData [fgsTypeID, x * (LigandID, ModelID), -1, ...]
        #

        # print(json.dumps({'fgsMapData_recv': fgsMapData_recv}))
        dump = json.dumps({"fgsMapData_recv": fgsMapData_recv})
        # fgsMapData_recv = {}
        return dump

    @app.route("/getGUIData", methods=["GET", "POST"])
    def getGUIData_recv():
        global GUIdata_recv
        if verbose:
            print("GUIdata_recv")
            print(GUIdata_recv)
        # Datapackage
        # 	[0]	 ->	 paramName;paramName....
        #   [1]	 ->	 [floatVal, floatVal,...]
        #

        dump = json.dumps({"GUIdata_recv": GUIdata_recv})
        return dump

    @app.route("/getSelectionData", methods=["GET", "POST"])
    def getSelectionData_recv():
        global selectionData_recv
        if verbose:
            print("selectionData_recv")
            print(selectionData_recv)

        # Datapackage
        # 	[0]	 ->  selectionName, selectionName
        #   [1]  ->  [intVal, intVal]
        #

        dump = json.dumps({"selectionData_recv": selectionData_recv})
        return dump

    @app.route("/getFuncGroupData", methods=["GET", "POST"])
    def getFuncGroupData():
        global funcGroupData_recv
        if verbose:
            print("funcGroupData_recv")
            print(funcGroupData_recv)

        # Datapackage
        #       [0]	->	funcGroupsPerLigands
        # 				[LigandID, funcGroupsCnt, X * (funcGroupID, groupCnt), -1, LigandID, ...]
        # 				-1 works as additional delimiter between ligands
        #       [1]	->	funcGroupsPerClusters
        # 				[ClusterID, funcGroupsCnt, X * (funcGroupID, groupCnt, groupCnt * [liganID,modelID], -1), -2, ClusterID, ...]
        #
        # 		*** additional information ***
        #        funcGroupsCnt:	    count of found functional groups
        #        X:				    count of found functional groups
        #        funcGroupID:	    specific ID for a type of functional group
        #        groupCnt:		    count of how many times a specific functional group was found
        # 		[liganID,modelID]:	ligand and model ID of func. group
        # 		-1:				    is used as delimiter between func. groups
        # 		-2:				    is used as delimiter between clusters

        # print(json.dumps({'funcGroupData_recv': funcGroupData_recv}))
        return json.dumps({"funcGroupData_recv": funcGroupData_recv})

    @app.route("/getMarkdownFile/<filename>", methods=["GET", "POST"])
    def getMarkdownFile(filename):
        with open("./vue_app/src/assets/md/" + filename + ".md", "r", encoding="utf-8") as file:
            data = file.read()
        return json.dumps(data)

    @app.route("/sendMegaMolID", methods=["GET", "POST"])
    def setGblMdlID():
        global megaMolDataKeys
        k = megaMolDataKeys
        transferred_data = flask.request.args
        print(transferred_data)

        if not testMode:
            if "clusterLigandModelID" in transferred_data:
                # key, data
                sendDataBuffer.put("int")
                sendDataBuffer.put(
                    [
                        k.clusterLigandModelID,
                        int(transferred_data["clusterID"]),
                        int(transferred_data["ligandID"]),
                        int(transferred_data["modelID"]),
                    ]
                )
                # print(str(int(transferred_data['clusterID'])+" "+ int(transferred_data['ligandID']+" "+ int(transferred_data['modelID']))))
                sendToMegamol()
            elif "gblMdlID" in transferred_data:
                sendDataBuffer.put("int")
                sendDataBuffer.put(
                    [k.gblMdlID, int(transferred_data["gblMdlID"])])
                print(transferred_data["gblMdlID"])
                sendToMegamol()
            elif "clusterLigandModelID" in transferred_data:
                sendDataBuffer.put("int")
                sendDataBuffer.put(
                    [
                        k.clusterLigandModelID,
                        int(transferred_data["clusterID"]),
                        int(
                            transferred_data["ligandID"],
                            int(transferred_data["modelID"]),
                        ),
                    ]
                )
                print(
                    str(
                        int(transferred_data["clusterID"])
                        + " "
                        + int(
                            transferred_data["ligandID"]
                            + " "
                            + int(transferred_data["modelID"])
                        )
                    )
                )
                sendToMegamol()
            else:
                print(
                    "Unkown data recivied! dataKey not registered! Don't know what do!"
                )

            # only a example for sending GUI params to MegaMol
            """
            sendDataBuffer.put("int")
            sendDataBuffer.put([k.GUIparam])
            sendDataBuffer.put("string")
            sendDataBuffer.put("hydrogenBonds")
            sendDataBuffer.put("float")
            sendDataBuffer.put([1.0])
            sendToMegamol()
            """

        return transferred_data

    @app.route("/getZINCsvgs", methods=["GET", "POST"])
    def getZINCsvgs():
        global clusterData_recv
        global clusterDataKeys
        zincSVGs = {}
        namespace = {"W3": "http://www.w3.org/2000/svg"}
        basePath = pathlib.Path("./testData/zinc_svgs")
        # print("Path:"+str(data_recv[clusterDataKeys.dataPath]))
        if not testMode:
            basePath = pathlib.Path(clusterData_recv[clusterDataKeys.dataPath])
        transferred_data = flask.request.json
        # print("transferred_data",transferred_data)
        for element in transferred_data:
            svgPath = basePath / (element + ".svg")

            print(svgPath)
            svg = ET.parse(str(svgPath))
            root = svg.getroot()
            zincText = root.findall(".//W3:text", namespaces=namespace)[-1]
            textParent = zincText.getparent()
            textParent.remove(zincText)

            zincRect = root.findall(".//W3:rect", namespaces=namespace)[0]
            textParent = zincRect.getparent()
            textParent.remove(zincRect)

            zincTitle = root.findall(".//W3:title", namespaces=namespace)[0]
            textParent = zincTitle.getparent()
            textParent.remove(zincTitle)

            xmlstr = ET.tostring(
                svg.getroot(), encoding="unicode", method="xml")
            zincSVGs[element.replace(".svg", "")] = xmlstr
        return json.dumps(zincSVGs)

    @app.route("/getLigandPNGs", methods=["GET", "POST"])
    def getLignadPNGs():
        # partially from https://stackoverflow.com/questions/64065587/how-to-return-multiple-images-with-flask
        basePath = pathlib.Path("./testData/functional_groups_pngs")
        transferred_data = flask.request.json
        byte_pngs = {}
        for imgName in transferred_data:
            img_path = basePath / (imgName.zfill(3) + ".png")
            imported_img = Image.open(img_path, "r")
            byte_arr = io.BytesIO()
            imported_img.save(byte_arr, format="PNG")
            encoded_img = base64.encodebytes(
                byte_arr.getvalue()).decode("ascii")
            byte_pngs[imgName] = "data:image/png;base64," + encoded_img
        return byte_pngs

    @app.route("/getProLigInterPNG", methods=["GET", "POST"])
    def getProLigInterPNG():
        # partially from https://stackoverflow.com/questions/64065587/how-to-return-multiple-images-with-flask
        basePath = pathlib.Path("./testData/plip_pngs")
        transferred_data = flask.request.json
        if not testMode:
            basePath = pathlib.Path(clusterData_recv[clusterDataKeys.dataPath])
        # print("transferred_data",transferred_data)
        byte_pngs = {}
        for element in transferred_data:
            # COMPLEXZINC000039132896_8_PROTEIN_LIG_Z_1 "_PROTEIN_LIG_Z_1"
            pngName = "COMPLEX" + element
            pngPath = basePath / (pngName + ".png")
            print(pngPath)
            imported_img = Image.open(pngPath, "r")
            byte_arr = io.BytesIO()
            imported_img.save(byte_arr, format="PNG")
            encoded_img = base64.encodebytes(
                byte_arr.getvalue()).decode("ascii")
            byte_pngs[element] = "data:image/png;base64," + encoded_img
        return byte_pngs

    @app.route("/navigation", methods=["GET", "POST"])
    def navigation():
        global megaMolDataKeys
        k = megaMolDataKeys
        transferred_data = flask.request.json
        imported_data = transferred_data["exportData"]
        # imported_data is an array ["paramName", value]
        imported_data[1] = float(imported_data[1])
        print(imported_data)
        print(type(imported_data[1]))

        sendDataBuffer.put("int")
        sendDataBuffer.put([k.GUIparam])
        sendDataBuffer.put("string")
        sendDataBuffer.put(imported_data[0])
        sendDataBuffer.put("float")
        sendDataBuffer.put([imported_data[1]])
        sendToMegamol()

        return transferred_data

    @app.route("/getOperationMode", methods=["GET", "POST"])
    def getOperationMode():
        dump = json.dumps({"operatioMode_recv": "megamolMode"})
        if testMode:
            dump = json.dumps({"operatioMode_recv": "testMode"})
        return dump


############# SETUP FLASK APP ##########
# https://stackoverflow.com/questions/27465533/run-code-after-flask-application-has-started

print("here")
if testMode:
    app = flask.Flask(__name__, static_folder="./dist/assets",
                      template_folder="./dist")
    sys.path.append("./testData/")
    from testData import clusterData, funcGroupData, FGSmapData, GUIdata, selectionData

    clusterData_recv = clusterData
    funcGroupData_recv = funcGroupData
    fgsMapData_recv = FGSmapData
    GUIdata_recv = GUIdata
    selectionData_recv = selectionData
    createAppRoutes()
else:
    app = flask.Flask(__name__, static_folder="./dist/assets",
                      template_folder="./dist")

    def create_app(*args, **kwargs):
        global app
        app = flask.Flask(
            __name__, static_folder="./dist/assets", template_folder="./dist"
        )
        ctx = click.get_current_context(silent=True)
        if ctx:
            script_info = ctx.obj
            config_mode = script_info.config_mode
        elif kwargs.get("config_mode"):
            # Production server, e.g., gunincorn
            # We don't have access to the current context, so must
            # read kwargs instead.
            config_mode = kwargs["config_mode"]
        createAppRoutes()
        custom_call()
        return app

    def custom_call():
        # define global_dict here
        thread_ = Thread(
            target=socketConnectionHandler,
            name="Thread1",
            args=[old_dataTypeString, dataTypes],
        )
        thread_.start()
        pass

    @click.group(cls=FlaskGroup, create_app=create_app)
    @click.option("-m", "--config-mode", default="Development")
    @pass_script_info
    def manager(script_info, config_mode):
        script_info.config_mode = config_mode


if testMode:

    if __name__ == "__main__":
        try:
            sys.argv[1]
            app.run(debug=True, host=sys.argv[1])
        except:
            _thread.start_new_thread(
                testModeFirefoxWindow, (("http://localhost:5000/",))
            )
            app.run(debug=False, host="0.0.0.0")

            # app.run(debug=True)
else:
    if __name__ == "__main__":
        manager()
