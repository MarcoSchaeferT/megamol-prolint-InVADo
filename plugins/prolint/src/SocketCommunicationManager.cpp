/*
 * SocketCommunicationManager.cpp
 *
 * Copyright (C) 2020 by University of Tuebingen (BDVA).
 * All rights reserved.
 */


#include "stdafx.h"
#include "SocketCommunicationManager.h"
#include <thread>
#include "MultiPDBQTLoader.h"
#include "mmcore/param/IntParam.h"
#include "mmcore/utility/log/Log.h"
#include "mmcore/utility/sys/Thread.h"
#include "vislib/StringConverter.h"
#include "vislib/net/DNS.h"
#include "vislib/net/IPEndPoint.h"
#include "vislib/net/SocketException.h"
#include <fstream>

#define TIMEOUT 10000000

using namespace megamol;
using namespace megamol::core;
using namespace megamol::prolint;


/***************************************
 ****** SocketCommunicationThread ******
 ***************************************/
/**
 * Constructor
 */
SocketCommunicationThread::SocketCommunicationThread(void)
    : socketValidity(false), sendDataRequested(false), incomingData(false), old_functionalGroup_dataHash(-1) {


    this->header.type = 0;
    this->header.length = 0;
    // Communication with SocketServer
    try {
        // try to start up socket
        vislib::net::Socket::Startup();
        // create socket
        this->socket.Create(
            vislib::net::Socket::FAMILY_INET, vislib::net::Socket::TYPE_STREAM, vislib::net::Socket::PROTOCOL_TCP);
    } catch (vislib::net::SocketException e) {
        megamol::core::utility::log::Log::DefaultLog.WriteMsg(
            megamol::core::utility::log::Log::LEVEL_ERROR, "Socket Exception during startup/create: %s", e.GetMsgA());
    }

    // set sizes for different incoming data types
    this->dataTypes.SetCount(4);
    dataTypes[0] = sizeof(int);
    dataTypes[1] = sizeof(float);
    dataTypes[2] = sizeof(char);
    dataTypes[3] = sizeof(uint);
}

/**
 * Create
 */
bool SocketCommunicationThread::create(void) {
    // TODO allocate variables etc.
    this->host = NULL;
    this->port = NULL;
    this->lmc_CallerSlotPtr = NULL;

    return true;
}

/**
 * Destructor
 */
SocketCommunicationThread::~SocketCommunicationThread(void) { this->release(); }


void SocketCommunicationThread::release(void) {
    socketValidity = false; // flag the socket as being non-functional
    try {
        vislib::net::Socket::Cleanup();
    } catch (vislib::net::SocketException e) {
        megamol::core::utility::log::Log::DefaultLog.WriteMsg(
            megamol::core::utility::log::Log::LEVEL_ERROR, "Socket Exception during cleanup: %s", e.GetMsgA());
    }
}


/*
 * SocketCommunicationThread::OnThreadStarting
 */
void SocketCommunicationThread::OnThreadStarting(void* config) {
   
    // needs to be implemented
}


/*
 * SocketCommunicationThread::Run
 */
DWORD SocketCommunicationThread::Run(void* config) {
    using megamol::core::utility::log::Log;

	startFlaskServer();


    // Initialize or Reinitialize the socket connection
    printf("Initialize Socket --- ip: %s \t port: %d\n", this->host.PeekBuffer(), (int)this->port);
    this->startSocket(this->host, this->port);

    while (this->socketValidity == true) {


        // send data requested
        if (this->sendDataRequested) {
            if (this->lmc_CallerSlotPtr) {
                this->sendClusterData();
                this->sendDataRequested = false;
            }
        } else if (this->incomingData) {

            this->incomingData = false;
            // std::string ready = "ready";
            // send_Type_Size_Data("string", ready.size(), &ready);
            /* this->header.type == megamolCODE send from python server
             * 0: 'int'
             * 1: 'float'
             * 2: 'string'
             * 3: 'unsigned int'
             */

            // this->incomingData = false;
            Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Waiting for incoming data stream...");
            switch (this->header.type) {
            case -1:
                // normal signal of the server that it is still alive
                break;
            case 0:
                // 'int'
                if (this->header.length < 0) break;
                this->intData.SetCount(this->header.length);
                try {
                    if (this->socket.Receive(&this->intData[0], this->header.length * sizeof(int), TIMEOUT, 0, true) !=
                        this->header.length * sizeof(int)) {
                        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "IntData has the worng size!;");
                    }
                } catch (vislib::net::SocketException e) {
                    printf("");
                    return false;
                }
                printf("IntData: ");
                for (int i = 0; i < this->intData.Count(); i++) {
                    printf("%d ", this->intData[i]);
                }
                printf("\n");

                consumeIntData(this->intData);
                break;

            case 1:
                // 'float'
                if (this->header.length < 0) break;
                this->floatData.SetCount(this->header.length);
                try {
                    if (this->socket.Receive(&this->floatData[0], this->header.length * sizeof(float), TIMEOUT, 0,
                            true) != this->header.length * sizeof(float)) {
                        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "FloatData has the worng size!;");
                    }
                } catch (vislib::net::SocketException e) {
                    printf("");
                    return false;
                }
                printf("floatData: ");
                for (int i = 0; i < this->floatData.Count(); i++) {
                    printf("%f ", this->floatData[i]);
                }
                printf("\n");
                if (this->reciving_GUI_data) {
                    this->guiParamValue = this->floatData[0];
                    LigandModelCall* lmc = this->lmc_CallerSlotPtr->CallAs<LigandModelCall>();
                    if (lmc == NULL) return false;

                    if (!(*lmc)(LigandModelCall::CallForGetExtent)) return false;

                    auto a = lmc->getWebData();
                    (*a->GUI)[this->guiParamName.c_str()] = this->guiParamValue;
                    lmc->getWebData()->UpdateWebData();
                    this->reciving_GUI_data = false;
                }
                break;

            case 2:
                // 'string'
                if (this->header.length < 0) break;
                this->stringData.resize(this->header.length);
                try {
                    if (this->socket.Receive(&this->stringData[0], this->header.length * sizeof(char), TIMEOUT, 0,
                            true) != this->header.length * sizeof(char)) {
                        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "StrData has the worng size!;");
                    }
                } catch (vislib::net::SocketException e) {
                    printf("");
                    return false;
                }
                printf("StrData: ");
                printf("%s ", this->stringData.c_str());
                printf("\n");
                if (this->reciving_GUI_data) {
                    this->guiParamName = this->stringData.c_str();
                }
                break;

            case 3:
                // 'uint'
                if (this->header.length < 0) break;
                this->uintData.SetCount(this->header.length);
                try {
                    if (this->socket.Receive(&this->uintData[0], this->header.length * sizeof(uint), TIMEOUT, 0,
                            true) != this->header.length * sizeof(uint)) {
                        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "uIntData has the worng size!;");
                    }
                } catch (vislib::net::SocketException e) {
                    printf("");
                    return false;
                }
                printf("uIntData: ");
                for (int i = 0; i < this->uintData.Count(); i++) {
                    printf("%u ", this->uintData[i]);
                }
                printf("\n");
                break;
            default:
                Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR,
                    "unkown data/code: %d - deserealization for this type is not implemented?", this->header.type);
            }
        } else {
            // waiting for incoming data/header
            if (getHeader()) {
                int signal = this->header.type;
                // printf("Signal: %d", signal);
                if (!(signal >= -1)) {
                    this->socketValidity = false;
                    Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR,
                        "Lifesignal: socket got no life signal from server since %d miliseconds.", TIMEOUT);
                } else if (signal != -1 && signal != 5) {
                    this->incomingData = true;
                } else {
                    // Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Lifesignal recieved");
                }
            }
            // needs to be implemented
        }
        LigandModelCall* lmc = this->lmc_CallerSlotPtr->CallAs<LigandModelCall>();
        if (lmc == NULL) continue;
        if (!(*lmc)(LigandModelCall::CallForGetExtent)) continue;

        auto status = lmc->getLMCStatus();
        if (status.isModelClustRendererBusy == false) {
            sendClusterData();
            sendFuncGroupData();
            sendFgsMap();
            sendGUIdata();
            sendSelectionData();
            if (this->old_WebData_dataHash != lmc->getWebData()->getWebDataHash()) {
                this->old_WebData_dataHash = lmc->getWebData()->getWebDataHash();
            }
        }


        bool modelClustRenderer_busy = 1;
        bool socketCommManager_busy = 0;
        lmc->SetlmcStatus(modelClustRenderer_busy, socketCommManager_busy, status.isFirstRenderCallComplete);
        (*lmc)(LigandModelCall::CallForSetLMCStatus);
    }
    return 0;
}


/*
 * SocketCommunicationThread::consumeIntData
 */
bool SocketCommunicationThread::consumeIntData(vislib::Array<int> intData) {
    DataKeys* k = &dataKeys;
    using megamol::core::utility::log::Log;
    LigandModelCall* lmc = this->lmc_CallerSlotPtr->CallAs<LigandModelCall>();
    auto webData = lmc->getWebData();

    if (intData[0] == k->gblMdlID) {
        printf("Received gblMdlID: %d\n", intData[1]);
        webData->setGlobalMdlID(this->intData[1]);
    } else if (intData[0] == k->ClusterLigandModelID) {
        if (intData.Count() == 4) {
            printf("Received ClusterID, ligandID , ModelID: %d,%d,%d\n", intData[1], intData[2], intData[3]);
            lmc->SetCurrentLigandAndModel(intData[2], intData[3]);
            (*lmc)(LigandModelCall::CallForGetData);
            webData->setGlobalMdlID(lmc->getCurrentGlobalModelID());
            webData->setClusterID(intData[1]);
			webData->setLigandID(intData[2]);
            webData->setModelID(intData[3]);
            
        }
    } else if (intData[0] == k->GUIparam) {
        this->reciving_GUI_data = true;
    } else {
        printf("Received unknown flag!\n");
    }
    webData->UpdateWebData();
    this->old_WebData_dataHash = webData->getWebDataHash();

    return false;
};


/*
 * SocketCommunicationThread::Terminate
 */
bool SocketCommunicationThread::Terminate(void) {
    using megamol::core::utility::log::Log;

    // clear all the requests
    this->sendDataRequested = false;

    return true;
}

/*
 * SocketCommunicationThread::Initialize
 */
void SocketCommunicationThread::Initialize(vislib::TString inHost, int inPort) {
    host = inHost;
    port = inPort;
    this->old_GUI_dataHash = -1;
    this->old_WebData_dataHash = -1;
    this->old_selection_dataHash = -1;
}


/*
 * SocketCommunicationThread::startSocket
 */
bool SocketCommunicationThread::startSocket(const vislib::TString& host, int port) {
    using megamol::core::utility::log::Log;
    // set default values
    header.type = 0;
    header.length = 0;
    bool checkServer = false;
    while (checkServer != true) {
        ;
        try {
            // connect to SocketCommunicationThread to allow real time data transfer
            this->socket.Connect(vislib::net::IPEndPoint::CreateIPv4(T2A(host), port));
            // handshake: receive a header
            if (this->socket.Receive(&this->header, 2 * sizeof(uint), 0, 0, true) != 2 * sizeof(uint)) {
                Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Handshake: did not receive full header to initiate.");

            } else {
                checkServer = true;
            }
        } catch (vislib::net::SocketException e) {
            Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR,
                "Waiting for server to start. Socket Exception during connect/handshake: %s", e.GetMsgA());
            Sleep(1000);
        }
    }

    printf("Header-> type: %d \t length: %d\n", this->header.type, this->header.length);
    this->socketValidity = true; // flag the socket as being functional
    return true;
}

/*
 * SocketCommunicationThread::getData
 */
bool SocketCommunicationThread::getData(void) {
    using megamol::core::utility::log::Log;

    bool retval = true;

    retval &= this->getHeader(); // see what kind of data is arriving

    // has to be implemented

    return retval;
}

/*
 * SocketCommunicationThread::getHeader
 */
bool SocketCommunicationThread::getHeader(void) {
    using megamol::core::utility::log::Log;


    try {

        // handshake: receive a header
        if (this->socket.Receive(&this->header, 2 * sizeof(int), TIMEOUT, 0, true) != 2 * sizeof(int)) {
            Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Handshake: did not receive full header to initiate.");
            printf("Header-> type: %d \t length: %d\n", this->header.type, this->header.length);
        }
    } catch (vislib::net::SocketException e) {
        printf("no header was recieved after %d seconds\n", TIMEOUT / 1000);
        return false;
    }
    // printf("Header-> type: %d \t length: %d\n", this->header.type, this->header.length);
    // printf("Header: %d \t %d\n", test[0], test[1]);

    return true;
}


/*
 * SocketCommunicationThread::send_Type_Size_Data
 */
bool SocketCommunicationThread::send_Type_Size_Data(std::string dataType, SIZE_T size, const void* data) {

    // set up datSize regarding data type
    int dataSize = 1;
    string* tmp = (string*)data;
    if (dataType == "string") {
        dataSize = size * sizeof(const char);
    } else if (dataType == "int") {
        dataSize = size * sizeof(int);
    } else if (dataType == "unsigned int") {
        dataSize = size * sizeof(unsigned int);
    } else if (dataType == "float") {
        dataSize = size * sizeof(float);
    } else if (dataType == "char") {
        dataSize = size * sizeof(char);
    } else {
        printf("ERROR: unkown dataType: %s\n", dataType.c_str());
    }

    // send data size
    if (waitForSignal_FIN()) {
        std::string bufferSize = "bufferSize";
        try {
            this->socket.Send(bufferSize.c_str(), bufferSize.size() * sizeof(const char), 10000, 0, false);
        } catch (vislib::net::SocketException e) {
            megamol::core::utility::log::Log::DefaultLog.WriteMsg(megamol::core::utility::log::Log::LEVEL_ERROR,
                "Socket Exception during %s send: %s", dataType.c_str(), e.GetMsgA());
            return 0;
        }
    }
    if (waitForSignal_FIN()) {
        try {
            this->socket.Send(&dataSize, 1 * sizeof(int), 10000, 0, true);
        } catch (vislib::net::SocketException e) {
            megamol::core::utility::log::Log::DefaultLog.WriteMsg(megamol::core::utility::log::Log::LEVEL_ERROR,
                "Socket Exception during bufferSize send: %s", e.GetMsgA());
            return 0;
        }
    }

    // send data type
    if (waitForSignal_FIN()) {
        try {
            this->socket.Send(dataType.c_str(), dataType.size() * sizeof(const char), 10000, 0, true);
        } catch (vislib::net::SocketException e) {
            megamol::core::utility::log::Log::DefaultLog.WriteMsg(megamol::core::utility::log::Log::LEVEL_ERROR,
                "Socket Exception during %s send: %s", dataType.c_str(), e.GetMsgA());
            return 0;
        }
    }

    // send the real data
    if (waitForSignal_FIN()) {

        if (dataType != "string") {
            try {
                this->socket.Send(data, dataSize, 10000, 0, true);
            } catch (vislib::net::SocketException e) {
                megamol::core::utility::log::Log::DefaultLog.WriteMsg(megamol::core::utility::log::Log::LEVEL_ERROR,
                    "Socket Exception during %s send: %s", dataType.c_str(), e.GetMsgA());
                return 0;
            }
        } else {
            try {
                this->socket.Send(tmp->c_str(), dataSize, 10000, 0, true);
            } catch (vislib::net::SocketException e) {
                megamol::core::utility::log::Log::DefaultLog.WriteMsg(megamol::core::utility::log::Log::LEVEL_ERROR,
                    "Socket Exception during %s send: %s", dataType.c_str(), e.GetMsgA());
                return 0;
            }
        }
    }


    return true;
}

/*
 * SocketCommunicationThread::sendClusterData
 */
bool SocketCommunicationThread::sendClusterData(void) {


    LigandModelCall* lmc = this->lmc_CallerSlotPtr->CallAs<LigandModelCall>();
    if (lmc == NULL) return false;

    if (!(*lmc)(LigandModelCall::CallForGetExtent)) return false;
    if (!(*lmc)(LigandModelCall::CallForGetClusterAndCentroidData)) return false;

    /* Datapackage
        [0]	 ->	 nameOfDataStream
        [1]	 ->	 ligandCount
        [2]	 ->	 modelCount
        [3]	 ->	 energies
        [4]	 ->	 ligandIDs
        [5]	 ->	 zincNames
        [6]	 ->	 clusterAssignment
        [7]	 ->	 clusterSizes
        [8]  ->  dataPath (path where svg, smi, checkmol, etc. are located)
        [9]  ->  hbonds
        [10] ->  efficiency
        [11] ->  mwt (molecular weigth)
        [12] ->  logp (logP is a quantitative measure of lipophilicit)
        [13] ->  fractioncsp3 (parameter for drug-likeness)
        [14] ->  halogenBonds
        [15] ->  hydrophobicInteractions
        [16] ->  metalComplexes
        [17] ->  piCationInteractions
        [18] ->  piStacks
        [19] ->  saltBridges
    */

    if (this->clusterData.clusterDataHash != lmc->getClusterAndCentroidData().clusterDataHash) {
        this->clusterData = lmc->getClusterAndCentroidData();
        if (this->clusterData.clusterSizes->Count() == 0) return false;
        // ###################################################################################
        std::string nameOfDataStream = "clusterData";

        int lgdCount = (int)lmc->getLigandCount();
        auto mdlCount = lmc->getModelCount();

        std::vector<float> allEnergies; // contains all energies of all ligands/models
        std::vector<int> lgdIDs;
        std::string allLgdNames;
        std::vector<float> mdlEfficiency;
        std::vector<float> mwt;
        std::vector<float> logp;
        std::vector<float> fractioncsp3;
        std::vector<int> clustAssigned;
        std::vector<int> clustSizes;
        string dataPath;

		// forces
        std::vector<int> HBonds;
        std::vector<int> halogenBonds;
        std::vector<int> hydrophobicInteractions; 
		std::vector<int> metalComplexes;
        std::vector<int> piCationInteractions;
		std::vector<int> piStacks;
        std::vector<int> saltBridges;

		HBonds.resize(lmc->getTotalModelCount());
        halogenBonds.resize(lmc->getTotalModelCount());
        hydrophobicInteractions.resize(lmc->getTotalModelCount());
        metalComplexes.resize(lmc->getTotalModelCount());
        piCationInteractions.resize(lmc->getTotalModelCount());
        piStacks.resize(lmc->getTotalModelCount());
        saltBridges.resize(lmc->getTotalModelCount());


		allLgdNames = "";
		lgdIDs.resize(lgdCount);
        mwt.resize(lgdCount);
        logp.resize(lgdCount);
        fractioncsp3.resize(lgdCount);
        allEnergies.resize(lmc->getTotalModelCount());
        mdlEfficiency.resize(lmc->getTotalModelCount());

        // get all ligand IDs and energies
        for (int i = 0; i < lgdCount; i++) {
            for (int j = 0; j < mdlCount[i]; j++) {
                lmc->SetCurrentLigandAndModel(i, j);
                if (!(*lmc)(LigandModelCall::CallForGetDataSilent)) return false;
                int gblMdlID = lmc->getCurrentGlobalModelID();

                if (j == 0) {
                    lgdIDs[i] = (int)lmc->getCurrentLigandID();
                    // concatenate all ligand names
                    allLgdNames += std::string((lmc->getCurrentLigandname()).PeekBuffer()) + ",";

                    mwt[i] = lmc->getChemPorpsOfCurLigand()->mwt;
                    logp[i] = lmc->getChemPorpsOfCurLigand()->logp;
                    fractioncsp3[i] = lmc->getChemPorpsOfCurLigand()->fractioncsp3;
                }

                allEnergies[gblMdlID] = lmc->getCurrentBindingenergy();
                mdlEfficiency[gblMdlID] = lmc->getCurrentModelEfficiency();

                // check for h-bonds, push 1 if true, else push 0
                int fCounter = 0;
                fCounter = lmc->getInteractionForces().hProtAcceptors->cnt +
                           lmc->getInteractionForces().hProtDonors->cnt;
                HBonds[gblMdlID] = fCounter > 0 ? 1 : 0;
                halogenBonds[gblMdlID] = lmc->getInteractionForces().halogenBonds->cnt > 0 ? 1 : 0;
                hydrophobicInteractions[gblMdlID] = lmc->getInteractionForces().hydrophobicInteractions->cnt > 0 ? 1 : 0;
                metalComplexes[gblMdlID] = lmc->getInteractionForces().metalComplexes->cnt > 0 ? 1 : 0;
                piCationInteractions[gblMdlID] = lmc->getInteractionForces().piCationInteractions->cnt > 0 ? 1 : 0;
                piStacks[gblMdlID] = lmc->getInteractionForces().piStacks->cnt > 0 ? 1 : 0;
                saltBridges[gblMdlID] = lmc->getInteractionForces().saltBridges->cnt > 0 ? 1 : 0; 
            }
        }

        // transform vislib arrays to vectors for sending
        for (int i = 0; i < clusterData.assignedClusters[0].Count(); i++) {
            clustAssigned.push_back(clusterData.assignedClusters[0][i]);
        }

        for (int i = 0; i < clusterData.clusterSizes[0].Count(); i++) {
            clustSizes.push_back(clusterData.clusterSizes[0][i]);
        }

        // svg paths
        std::vector<string> tmpSplit;
        string tmpString = lmc->getSVGPath();
        split(tmpString, "\\", tmpSplit);

        for (int i = 0; i < tmpSplit.size() - 1; i++) {
            dataPath.append(tmpSplit[i].c_str());
            dataPath.append("\\");
        }

        send_Type_Size_Data("string", nameOfDataStream.size(), &nameOfDataStream); // send data stream name

        send_Type_Size_Data("int", 1, &lgdCount); // send ligand count

        send_Type_Size_Data("int", lgdCount, mdlCount); // send model count

        send_Type_Size_Data("float", allEnergies.size(), allEnergies.data()); // send energies

        send_Type_Size_Data("int", lgdIDs.size(), lgdIDs.data()); // send ligand IDs

        send_Type_Size_Data("string", allLgdNames.size(), &allLgdNames); // send ligand names

        send_Type_Size_Data("int", clustAssigned.size(), clustAssigned.data()); // send assigned clusters

        send_Type_Size_Data("int", clustSizes.size(), clustSizes.data()); // send cluster sizes

        send_Type_Size_Data("string", dataPath.size(), &dataPath); // send svg paths

        send_Type_Size_Data("int", HBonds.size(), HBonds.data()); // send hbonds
    
        send_Type_Size_Data("float", mdlEfficiency.size(), mdlEfficiency.data()); // send efficiency

        send_Type_Size_Data("float", mwt.size(), mwt.data()); // send molecular weight

        send_Type_Size_Data("float", logp.size(), logp.data()); // send logp

        send_Type_Size_Data("float", fractioncsp3.size(), fractioncsp3.data()); // send fractioncsp3

		send_Type_Size_Data("int", halogenBonds.size(), halogenBonds.data()); // send halogenBonds

        send_Type_Size_Data("int", hydrophobicInteractions.size(), hydrophobicInteractions.data()); // send hydrophobicInteractions

        send_Type_Size_Data("int", metalComplexes.size(), metalComplexes.data()); // send

        send_Type_Size_Data("int", piCationInteractions.size(), piCationInteractions.data()); // send piCationInteractions

        send_Type_Size_Data("int", piStacks.size(), piStacks.data()); // send piStacks

        send_Type_Size_Data("int", saltBridges.size(), saltBridges.data()); // send saltBridges


        allEnergies.clear();
        lgdIDs.clear();
        clustAssigned.clear();
        clustSizes.clear();
        HBonds.clear();
        mdlEfficiency.clear();
        mwt.clear();
        logp.clear();
        fractioncsp3.clear();


        // End-of-data-stream
        // std::string eod = "eod";
        // send_Type_Size_Data("string", eod.size(), &eod);


        // ###################################################################################


        // vislib::sys::Thread::Sleep(50);
    }

    return true;
}

/*
 * SocketCommunicationThread::sendfuncGroupData
 */
bool SocketCommunicationThread::sendFuncGroupData(void) {


    LigandModelCall* lmc = this->lmc_CallerSlotPtr->CallAs<LigandModelCall>();
    if (lmc == NULL) return false;


    if (!(*lmc)(LigandModelCall::CallForGetExtent)) return false;
    if (!(*lmc)(LigandModelCall::CallForGetClusterAndCentroidData)) return false;


    /* Datapackage
        [0]	->	nameOfDataStream
        [1]	->	funcGroupsPerLigands
                [LigandID, funcGroupsCnt, X * (funcGroupID, groupCnt), -1, LigandID, ...]
                -1 works as additional delimiter between ligands
        [2]	->	funcGroupsPerClusters
                [ClusterID, funcGroupsCnt, X * (funcGroupID, groupCnt, groupCnt * [liganID,modelID], -1), -2, ClusterID,
       ...]

        *** additional information ***
        funcGroupsCnt:		count of found functional groups
        X:					count of found functional groups
        funcGroupID:		specific ID for a type of functional group
        groupCnt:			count of how many times a specific functional group was found
        [liganID,modelID]:	ligand and model ID of func. group
        -1:					is used as delimiter between func. groups
        -2:					is used as delimiter between clusters
    */

    std::vector<int> funcGroupsPerLigands;
    std::vector<int> funcGroupsPerClusters;

    auto fgs = lmc->getClusteredFunctionalGroups();
    if (fgs->getDataHash() >= 0 && fgs->getDataHash() != this->old_functionalGroup_dataHash) {
        this->old_functionalGroup_dataHash = fgs->getDataHash();


        /*************************************************
         ******** get functional groups per ligand *******
         *************************************************/
        for (int i = 0; i < lmc->getLigandCount(); i++) {
            lmc->SetCurrentLigandAndModel(i, 0);
            (*lmc)(LigandModelCall::CallForGetDataSilent);
            auto fgs_ofCurLigand = lmc->getFunctionalGroupsOfCurLigand();
            // push ligandID
            funcGroupsPerLigands.push_back(i);
            // push funcGroupsCnt
            funcGroupsPerLigands.push_back(fgs_ofCurLigand->fGroups.size());
            // push X * tupel of groupID, speficGroupCnt
            for (int j = 0; j < fgs_ofCurLigand->fGroups.size(); j++) {
                funcGroupsPerLigands.push_back(fgs_ofCurLigand->fGroups[j].groupID);
                funcGroupsPerLigands.push_back(fgs_ofCurLigand->fGroups[j].speficGroupCnt);
            }
            // push -1 as delemiter (fgs_ofCurLigand->fGroups.size() could also be used but with -1 it is more safe)
            funcGroupsPerLigands.push_back(-1);
        }

        /***********************************************************
         ******** get grouped functional groups per clusters *******
         ***********************************************************/
        //

        auto fgsPerClusters = lmc->getClusteredFunctionalGroups();
        auto fgsOneCluster = fgsPerClusters->getPointer_to_functionalGroupClusterSets()[0];
        for (int pocketID = 0; pocketID < fgsPerClusters->getfunctionalGroupSetsCount(); pocketID++) {
            fgsOneCluster = fgsPerClusters->getPointer_to_functionalGroupClusterSets()[pocketID];
            // push clusterID
            funcGroupsPerClusters.push_back(pocketID);

            // push X * combinations of (funcGroupID, groupCnt, Y * [liganID,modelID])
            map<uint, vector<int>> totalFGS_map;
            for (int fgsClust = 0; fgsClust < fgsOneCluster.clusterCnt; fgsClust++) {
                if (fgsOneCluster.clusterSizes.size() == 0) {
                    continue;
                }
                if (fgsOneCluster.clusterSizes[fgsClust] >= 1) {
                    auto fgsMap = fgsOneCluster.fgsMaps_GblMDlIDs[fgsClust];
                    map<uint, std::vector<uint>>::iterator it;
                    for (it = fgsMap.begin(); it != fgsMap.end(); it++) {
                        int fgsTypeID = it->first;
                        int groupCnt = it->second.size();
                        for (int p = 0; p < groupCnt; p++) {
                            lmc->SetCurrentLigandAndModel_byGlobalMdlID(it->second[p]);
                            int ligID = lmc->getCurrentLigandID();
                            int mdlID = lmc->getCurrentModelID();
                            totalFGS_map[fgsTypeID].push_back(ligID);
                            totalFGS_map[fgsTypeID].push_back(mdlID);
                        }
                    }
                }
            }
            // push funcGroupsCnt
            int funcGroupsCnt = 0;
            funcGroupsCnt = totalFGS_map.size();

            funcGroupsPerClusters.push_back(funcGroupsCnt);
            if (funcGroupsCnt == 0) {
                funcGroupsPerClusters.push_back(-1);
            }
            map<uint, vector<int>>::iterator it;
            for (it = totalFGS_map.begin(); it != totalFGS_map.end(); it++) {
                int fgsTypeID = it->first;
                int groupCnt = it->second.size() / 2;
                // push part of combo: funcGroupID
                funcGroupsPerClusters.push_back(fgsTypeID);
                // push part of combo: groupCnt
                funcGroupsPerClusters.push_back(groupCnt);
                // push part of combo: Y * [liganID,modelID]
                for (int z = 0; z < groupCnt; z++) {
                    int ligID = totalFGS_map[fgsTypeID][z * 2 + 0];
                    int mdlID = totalFGS_map[fgsTypeID][z * 2 + 1];
                    funcGroupsPerClusters.push_back(ligID);
                    funcGroupsPerClusters.push_back(mdlID);
                }
                // push -1 as delemiter between functional groups
                funcGroupsPerClusters.push_back(-1);
            }

            // push -2 as delemiter between clusters
            funcGroupsPerClusters.push_back(-2);
        }

        std::string nameOfDataStream = "functionalGroupData";
        send_Type_Size_Data("string", nameOfDataStream.size(), &nameOfDataStream);

        // send functional groups per ligand
        send_Type_Size_Data("int", funcGroupsPerLigands.size(), funcGroupsPerLigands.data());
        // send functional groups per clusters
        send_Type_Size_Data("int", funcGroupsPerClusters.size(), funcGroupsPerClusters.data());

        // End-of-data-stream
        std::string eod = "reload";
        send_Type_Size_Data("string", eod.size(), &eod);
    }

    return true;
}


/*
 * SocketCommunicationThread::sendGUIdata
 */
bool SocketCommunicationThread::sendGUIdata(void) {

    LigandModelCall* lmc = this->lmc_CallerSlotPtr->CallAs<LigandModelCall>();
    if (lmc == NULL) return false;

    if (!(*lmc)(LigandModelCall::CallForGetExtent)) return false;

    auto web = lmc->getWebData()->GUI;
    if (this->old_GUI_dataHash != lmc->getWebData()->getGUIdataHash()) {
        this->old_GUI_dataHash = lmc->getWebData()->getGUIdataHash();

        /* Datapackage
            [0]	 ->	 paramName;paramName....
            [1]	 ->	 [floatVal, floatVal,...]
        */

        std::string paramNames;
        std::vector<float> paramValues;


        // collect data
        paramNames += "GUIdataHash";
        paramValues.push_back(lmc->getWebData()->getGUIdataHash());

        std::map<std::string, float>::iterator it;
        for (it = web->begin(); it != web->end(); it++) {
            paramNames += ";" + it->first;
            paramValues.push_back(it->second);
        }

        // send data
        std::string nameOfDataStream = "GUIdata";
        send_Type_Size_Data("string", nameOfDataStream.size(), &nameOfDataStream); // send name of data stream
        send_Type_Size_Data("string", paramNames.size(), &paramNames);             // send param names
        send_Type_Size_Data("float", paramValues.size(), paramValues.data());      // send param values

        // End-of-data-stream
        std::string eod = "eod";
        send_Type_Size_Data("string", eod.size(), &eod);
    }

    return true;
}

/*
 * SocketCommunicationThread::sendSelectionData
 */
bool SocketCommunicationThread::sendSelectionData(void) {

    LigandModelCall* lmc = this->lmc_CallerSlotPtr->CallAs<LigandModelCall>();
    if (lmc == NULL) return false;

    if (!(*lmc)(LigandModelCall::CallForGetExtent)) return false;

    if (this->old_selection_dataHash != lmc->getWebData()->getSelectionDataHash()) {
        this->old_selection_dataHash = lmc->getWebData()->getSelectionDataHash();

        /* Datapackage
            [0]	 ->  selectionName, selectionName
            [1]  ->  [intVal, intVal]
        */

        std::string selectionNames;
        std::vector<int> selectionValues;

        // collect data
        selectionNames = "selectiodDataHash;ClusterID;LigandID;ModelID";
        selectionValues.push_back(lmc->getWebData()->getSelectionDataHash());
        selectionValues.push_back(lmc->getWebData()->getClusterID());
        selectionValues.push_back(lmc->getWebData()->getLigandID());
        selectionValues.push_back(lmc->getWebData()->getModelID());

        // send data
        std::string nameOfDataStream = "selectionData";
        send_Type_Size_Data("string", nameOfDataStream.size(), &nameOfDataStream);  // send name of data stream
        send_Type_Size_Data("string", selectionNames.size(), &selectionNames);      // send selection names
        send_Type_Size_Data("int", selectionValues.size(), selectionValues.data()); // send selection values


        // End-of-data-stream
        std::string eod = "eod";
        send_Type_Size_Data("string", eod.size(), &eod);
    }

    return true;
}


/*
 * SocketCommunicationThread::sendClusterData
 */
bool SocketCommunicationThread::sendFgsMap(void) {


    LigandModelCall* lmc = this->lmc_CallerSlotPtr->CallAs<LigandModelCall>();
    if (lmc == NULL) return false;

    if (!(*lmc)(LigandModelCall::CallForGetExtent)) return false;

    if (this->old_WebData_dataHash != lmc->getWebData()->getWebDataHash()) {

        auto fgsMap = lmc->getWebData()->curFgsMap;
        if (fgsMap.size() == 0) return true;

        /* Datapackage (-1 delemiter)
            [0]	 ->	 fgsData [fgsTypeID, x * (LigandID, ModelID), -1, ...]
        */

        // collect the data
        std::vector<int> fgsData;
        map<uint, std::vector<uint>>::iterator it;
        for (it = fgsMap.begin(); it != fgsMap.end(); it++) {
            int fgsTypeID = it->first;
            fgsData.push_back(fgsTypeID);

            // push vector of gblMdlIDs
            for (int j = 0; j < it->second.size(); j++) {
                int gblMdlID = it->second[j];
                lmc->SetCurrentLigandAndModel_byGlobalMdlID(gblMdlID);
                int curLig = lmc->getCurrentLigandID();
                int curMdl = lmc->getCurrentModelID();
                fgsData.push_back(curLig);
                fgsData.push_back(curMdl);
            }
            // push delemiter
            fgsData.push_back(-1);
        }


        // send the data
        if (fgsData.size() > 0) {
            std::string nameOfDataStream = "selFgsClustDat";
            send_Type_Size_Data("string", nameOfDataStream.size(), &nameOfDataStream);
            send_Type_Size_Data("int", fgsData.size(), fgsData.data()); // send fgsData
                                                                        // End-of-data-stream
            std::string eod = "eod";
            send_Type_Size_Data("string", eod.size(), &eod);
        }


        fgsData.clear();
    }
    return true;
}


/*
 * SocketCommunicationThread::byteSwap
 */
int SocketCommunicationThread::byteSwap(int input) {
    char output[4];
    output[0] = ((char*)&input)[3];
    output[1] = ((char*)&input)[2];
    output[2] = ((char*)&input)[1];
    output[3] = ((char*)&input)[0];
    return *((int*)output);
}

/*
 * SocketCommunicationThread::waitForSignal_FIN
 */
bool SocketCommunicationThread::waitForSignal_FIN(void) {
    using megamol::core::utility::log::Log;
    int signal[2] = {0, 0};

    try {
        // handshake: receive signal for a completely processed request
        if (this->socket.Receive(&signal, 2 * sizeof(int), 0, 0, true) != 2 * sizeof(int)) {
            Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Signal FIN: did not receive full/correct signal.");
            return false;
        }
    } catch (vislib::net::SocketException e) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Socket Exception during signal (FIN) receive: %s", e.GetMsgA());
        return false;
    }
    return true;
}

void SocketCommunicationThread::setLMC_CallerSlotPtr(megamol::core::CallerSlot* lmc_CallerSlotPtr) {
    this->lmc_CallerSlotPtr = lmc_CallerSlotPtr;
};

void SocketCommunicationThread::setParamPointer(core::param::ParamSlot* serverPathParam) {
    this->serverPathParam = serverPathParam;
}


/**************************************
 ***** SocketCommunicationManager *****
 **************************************/
/**
 * Constructor
 */
SocketCommunicationManager::SocketCommunicationManager(void)
    : Module()
	, ligandModelCall("getData", "Call for a specific model of a specific ligand") 
	, serverPathParam("serverPath", "file path to the flask web server") {
    // set up caller slot
    this->ligandModelCall.SetCompatibleCall<LigandModelCallDescription>();
    this->MakeSlotAvailable(&this->ligandModelCall);

    // set up parameters

    this->scm = NULL;

	 // set up paramter slots
    this->serverPathParam << new param::FilePathParam("../../../plugins/prolint/server/app.py");
    this->MakeSlotAvailable(&this->serverPathParam);
}


/**
 * Destructor
 */
SocketCommunicationManager::~SocketCommunicationManager(void) { this->release(); }

std::string getOsName() {
#ifdef _WIN32
    return "Windows";
#elif _WIN64
    return "Windows";
#elif __APPLE__ || __MACH__
    return "Mac OSX";
#elif __linux__
    return "Linux";
#elif __FreeBSD__
    return "FreeBSD";
#elif __unix || __unix__
    return "Unix";
#else
    return "Other";
#endif
}                     

/**
 * Create
 */
bool SocketCommunicationManager::create(void) {
    // start Websocket on localhost
    vislib::TString ip = "127.0.0.1";
    int port = 2040;
    if (!this->scm) {
        LigandModelCall* lmc = this->ligandModelCall.CallAs<LigandModelCall>();
      

        this->scm = new vislib::sys::RunnableThread<SocketCommunicationThread>;
        this->scm->setLMC_CallerSlotPtr(&this->ligandModelCall);
        this->scm->setParamPointer(&this->serverPathParam);
        this->scm->Initialize(ip, port);
        this->scm->Start();
    }

    return true;
}

bool SocketCommunicationThread::startFlaskServer() {
    // start python server

	// sleep needed because file path has been not loaded from config file yet
	// ToDo: FIXME!
    Sleep(500);    
    
    string tmpFile = this->serverPathParam->Param<param::FilePathParam>()->Value().PeekBuffer();

    std::string baseDir;
    std::string serverName;

	
    std::vector<std::string> tmpSplit;
    split(tmpFile, "\\", tmpSplit);
    if (tmpSplit.size() == 0) {
        split(tmpFile, "/", tmpSplit);
    }
    int tmpSplitCnt = 0;
    for (int i = 0; i < tmpSplit.size() - 1; i++) {
        baseDir.append(tmpSplit[i].c_str());
        baseDir.append("\\");
        tmpSplitCnt++;
    }
    serverName = tmpSplit[tmpSplit.size()-1];


	

	std::string server;
    std:string pythonInstanceCleanUp;
	if (getOsName() == "Linux") {
        server = "cd " + baseDir + " && gnome-terminal -e 'sh -c \"python " + serverName +
                             " "
                             " -m Development run --host=0.0.0.0 --port 8050\"";
        int result = system(server.c_str());
    
	}
	else{
        std::string pathToPython = "..\\..\\..\\build\\install\\bin\\python";
        std::ifstream f(pathToPython.c_str());
        pythonInstanceCleanUp = "taskkill /im python.exe /f";
        server = "cd " + baseDir + " && start cmd.exe /k \"python " + pathToPython + " " + serverName +
                 " "
                 " -m Development run --host=0.0.0.0 --port 8050\"";
        printf("serverPath: %s\n", server.c_str());
        if (!f.good()) {
            server = "cd " + baseDir + " && start cmd.exe /k \"python " + serverName +
                     " "
                     " -m Development run --host=0.0.0.0 --port 8050\"";
        }
        system(pythonInstanceCleanUp.c_str());
        int result = system(server.c_str());
	}
    
    this->serverPathParam->ResetDirty();
	return true;
   
}

/**
 * release
 */
void SocketCommunicationManager::release(void) {}
