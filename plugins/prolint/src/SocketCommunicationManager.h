/*
 * SocketCommunicationManager.h
 *
 * Copyright (C) 2019 by University of Tuebingen (BDVA).
 * All rights reserved.
 */

#ifndef PROLINT_PLUGIN_SocketCommunicationThread_H_INCLUDED
#define PROLINT_PLUGIN_SocketCommunicationThread_H_INCLUDED
typedef unsigned int uint;

#include "LigandModelCall.h"
#include "mmcore/Call.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"
#include "mmcore/param/ParamSlot.h"
#include "mmcore/utility/log/Log.h"

#include "mmcore/utility/sys/RunnableThread.h"
#include "vislib/net/Socket.h"
#include "vislib/sys/CriticalSection.h"
#include "vislib/sys/Runnable.h"
#include <map> 

#define MDD_VERSION 2
using namespace std;

namespace megamol {
namespace prolint {




	
/**************************************
 ***** SocketCommunicationThread ******
 **************************************/
class SocketCommunicationThread : public vislib::sys::Runnable {


public:
    /** Ctor. */
    SocketCommunicationThread();
    /** Dtor. */
    virtual ~SocketCommunicationThread();

    /*
     * The 8 byte header that is transferred to and from socket server (python) with each new set of data.
     */
    struct SCMHeader {
        int type; // this is always associated with a type from the SCMHeaderType enum, but needs to be an int
                           // for byte swapping and to make it easier to send over networks.
        int length; // the size of incoming data
    };

	// maps which type of data is used/is incoming etc.
    // this list sould be the same in the python server (as an dictionary)
	class DataKeys {
    public:
        DataKeys(){
            this->ClusterLigandModelID = -333;
            this->GUIparam = -444;
            this->gblMdlID = -777;
		};
        int ClusterLigandModelID;
        int GUIparam;
        int modelID;
        int ligandModelID;
		int gblMdlID;
        int ligandID;
        int clusterID;
	};
    DataKeys dataKeys;
	 
	string guiParamName;
    float guiParamValue;
    bool reciving_GUI_data = false;


    /**
     * Startup callback of the thread. The Thread class will call that
     * before Run().
     *
     * @param config A pointer to the Configuration, which specifies the
     *               settings of the connector.
     */
    virtual void OnThreadStarting(void* config);

    /**
     * Perform the work of a thread.
     *
     * @param config A pointer to the Configuration, which specifies the
     *               settings of the connector.
     *
     * @return The application dependent return code of the thread. This
     *         must not be STILL_ACTIVE (259).
     */
    virtual DWORD Run(void* config);

    /**
     * Abort the work of the connector by forcefully closing the underlying
     * communication channel.
     *
     * @return true.
     */
    virtual bool Terminate(void);


    /**
     * Checks whether or not the socket is functional.
     *
     * @return 'true' if valid, 'false' otherwise.
     */
    bool IsSocketFunctional(void) { return this->socketValidity; }

    /**
     * Sets the host and the port for the socket connection.
     * This MUST be called before calling Start for the thread, or the socket
     * will never get set up correctly.
     *
     * @param inHost Host name or IP.
     * @param inPort Port number.
     */
    void Initialize(vislib::TString inHost, int inPort);


    /**
     * setCallPointer
     *
     * @ param call		a pointer to the call (should be ligandModelCall)
     */
    
	void setLMC_CallerSlotPtr(megamol::core::CallerSlot* lmc_CallerSlotPtr);

	void setParamPointer(core::param::ParamSlot* serverPathParam);

    /**
     * Sends data.
     *
     * @return True on success.
     */
    bool sendClusterData(void);

	bool sendFuncGroupData(void);

	bool sendGUIdata(void);

	bool sendSelectionData(void);

	/** 
	 *	send fgsMap of current selected functional group cluster of a selected pocket
	 *  this map stores per fgsTypeID a vector of global model IDs
	 */
	bool sendFgsMap(void);

    int  old_functionalGroup_dataHash;

	int old_WebData_dataHash;
    int old_selection_dataHash;
    int old_GUI_dataHash;
	

    /**
     * Requests that the thread sends data to the socket server.
     */
    inline void RequestSendData(void) { this->sendDataRequested = true; }


protected:
    /**
     * Implementation of 'Create'.
     *
     * @return 'true' on success, 'false' otherwise.
     */
    virtual bool create(void);

    /**
     * Implementation of 'Release'.
     */
    virtual void release(void);

private:
    /*
     * This list was taken from the MDDriver code ("imd.h" version 1.2 2008-06-25)
     * The list corresponds to the header value transferred by MDDriver, and specifies
     * what kind of data is to be expected in the following data transfer. It is also
     * used to send commands to MDDriver.
     */
    enum SCMHeaderType {  // currently not in use! (Marco)
        MDD_DISCONNECT, //!< close IMD connection, leaving sim running
        MDD_HANDSHAKE,  //!< endianism and version check message
    };

    vislib::Array<uint> dataTypes;


    /**
     * Starts the socket connection with the given host and port; performs handshaking.
     * Sets the socket validity flag to true if it succeeds.
     *
     * @param host String representing either name or IP of the machine running SocketCommunicationThread.
     * @param port Port number to communicate with on the server machine.
     * @return 'true' on success, 'false' otherwise.
     */
    bool startSocket(const vislib::TString& host, int port);


    /**
     * Receives header data from socket server (python) 
     *
     * @return True on success.
     */
    bool getHeader(void);

    /**
     * Switches order of bytes in an int so it changes endian
     */
    int byteSwap(int input);

    /**
     * Gets data from socket server (python)
     *
     * @return True on success.
     */
    bool getData(void);

    /**
     * Sends  Data Type and Size to the connected socket server
     *
     * @param dataType		tells what type of data is to sen ("int","float","string", etc.)
     * @param size			is the number of elements to send (not the byte size!)
     * @param data			the data that is to send
     *
     * @return				True on success.
     */

    bool send_Type_Size_Data(std::string dataType, SIZE_T size, const void* data);
    bool send_Type_Size_Data(const char* dataType, SIZE_T size, const void* data) {
        return send_Type_Size_Data(std::string(dataType), size, data);
    }

	/**
     * waitForSignal_FIN	waits until the socket server sends a signal that it finished data processing
     *
     * @return				True on success.
     */
    bool waitForSignal_FIN(void);

    bool startFlaskServer();


    // -------------------- variables --------------------

    /** The header */
    SCMHeader header;

    /** The socket for SocketCommunicationThread connection */
    vislib::net::Socket socket;

    /** The socket status */
    bool socketValidity;

    /** The host name to which the thread socket will connect */
    vislib::TString host;

    /** The port for the socket connection */
    int port;

    /** Flag requesting that the thread sends data */
    bool sendDataRequested;

	/** Flag that the thread is waiting for incoming data stream */
    bool incomingData;
    
	// data
    vislib::Array<int> intData;
    vislib::Array<float> floatData;
    std::string stringData;
    vislib::Array<uint> uintData;

	// data functions
    bool consumeIntData(vislib::Array<int> intData);

	/** a pointer to the Call (should be the lignadModelCall)*/
    megamol::core::CallerSlot* lmc_CallerSlotPtr;


	LigandModelCall::ClusterAndCentroidData clusterData;

	core::param::ParamSlot* serverPathParam;

};



/**************************************
 ***** SocketCommunicationManager *****
 **************************************/
class SocketCommunicationManager : public megamol::core::Module {


public:
    /** Ctor. */
    SocketCommunicationManager();
    /** Dtor. */
    virtual ~SocketCommunicationManager();

    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName(void) { return "SocketCommunicationManager"; }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description(void) { return "----------"; }

    /**
     * Answers whether this module is available on the current system.
     *
     * @return 'true' if the module is available, 'false' otherwise.
     */
    static bool IsAvailable(void) { return true; }

	// parmam slots
    core::param::ParamSlot serverPathParam;
    


protected:
    /**
     * Implementation of 'Create'.
     *
     * @return 'true' on success, 'false' otherwise.
     */
    virtual bool create(void);

    /**
     * Implementation of 'Release'.
     */
    virtual void release(void);


private:
    /** LigandModelCall caller slot */
    megamol::core::CallerSlot ligandModelCall;
 
    // -------------------- variables --------------------
    vislib::sys::RunnableThread<SocketCommunicationThread>* scm;

	


};








} // namespace prolint
} /* end namespace megamol */

#endif // PROLINT_PLUGIN_SocketCommunicationThread_H_INCLUDED
