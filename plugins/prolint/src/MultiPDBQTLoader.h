/*
 * MultiPDBQTLoader.h
 *
 * Copyright (C) 2019 by University of Tuebingen (BDVA).
 * All rights reserved.
 */

#ifndef PROLINT_PLUGIN_PDBQT_H_INCLUDED
#define PROLINT_PLUGIN_PDBQT_H_INCLUDED
typedef unsigned int uint;

#include <string>
#include <vector>
#include <filesystem> 
#include "LigandModelCall.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"
#include "mmcore/param/FilePathParam.h"
#include "mmcore/param/FloatParam.h"
#include "mmcore/param/IntParam.h"
#include "mmcore/param/ParamSlot.h"
#include "mmcore/utility/log/Log.h"
#include "mmcore/utility/sys//ASCIIFileBuffer.h"
#include "protein_calls/MolecularDataCall.h"
#include "vislib/Array.h"
#include "vislib/ArrayAllocator.h"
#include "vislib/StringConverter.h"
#include "vislib/math/Vector.h"


using namespace std;
typedef std::filesystem::path Path;

namespace megamol {
namespace prolint {
typedef vislib::Array<vislib::math::Vector<float, 4>> visFloat4;

/****************************
 ***** HELPER FUNCTIONS *****
 ****************************/
// splits a string s into a vector of strings v at the character c
void split(const string& s, const char* c, std::vector<string>& v);


/****************************
 **** PDBQT LOADER CLASS ****
 ****************************/
class MultiPDBQTLoader : public megamol::core::Module {


public:
    /** Ctor. */
    MultiPDBQTLoader();
    /** Dtor. */
    virtual ~MultiPDBQTLoader();

    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName(void) { return "MultiPDBQTLoader"; }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description(void) { return "Loads docking data."; }

    /**
     * Answers whether this module is available on the current system.
     *
     * @return 'true' if the module is available, 'false' otherwise.
     */
    static bool IsAvailable(void) { return true; }

    // TODO move all variables to "private"
    // important sizes			ligandCnt < modelCnt < atomCnt
    int total_ligandCnt;
    int total_modelCnt;
    int total_atomCnt; // CAUTION! counts atoms of all ligands molecules / this value is not multiplied with number
                       // of models per ligand
    int total_atomPosCnt;
    bool PDBQTLoadingStatus;

	bool setLMCStatus(core::Call& call);
	// MPL = MultiPDBQTLoader
	bool MPL_modelClustRenderer_busy;
    bool MPL_SocketCommManager_busy;
    bool MPL_firstRenderCall_complete;


    /*********************
     * storing the data
     *********************/
    /*
     * @ligandsName:	ligandsName[ligandID]
     * @ligandsAtomCnt	ligandsAtomCnt[ligandID]
     * @atomTypes:		atome names (C,H,O,N,...) of a ligand
     *					---> atomTypes[atomIndexSum[lignadID]+i]   (i=0; i < lignadsAtomCnt[ligandID]; i++)
     * @atomTypesAD:	AutoDock atome types (C, OA, NA, HD,...) of a ligand %%% accessing: (see atomTypes)
     * @ligandsModelRange:	stores the index range for accessing values of models of a ligand
     *				    ---> ARRAY[lignadModelRange[ligandID]+1] -to- ARRAY[lignadModelRange[ligandID+1]]
     *                  [0]->8; [1]->20 %%% so if LigandID == 0 you get 8 %%% means Range 0-8!
     *
     ---> accessing arrays below via models of a ligand: (see @lignadsModelRange):
     * @RMSD_lb:		RMSD lower bound of ligand model
     * @RMSD_ub:		RMSD upper bound (see RMSD_lb)
     * @atomPositions:	x,y,z cooridnates of ligand model atoms
     * @atomCharges:	cahrges of of ligand model atoms
     * @bindingEnergy:	binding energy of lignad model
     **/

    vislib::StringA PDBQT_list_Filepath;
    vislib::Array<vislib::StringA> ligandNames;
    vislib::Array<uint> ligandAtomCnts;
    vislib::Array<vislib::StringA> atomTypesAD;
    vislib::Array<megamol::protein_calls::MolecularDataCall::AtomType> atomTypes;
    vislib::Array<uint> atomTypeIdx;
    vislib::Array<uint> ligandsModelRange;
    vislib::Array<uint> modelsPerLigand;
    vislib::Array<float> RMSD_lb;
    vislib::Array<float> RMSD_ub;
    vislib::Array<float> atomPositions;
    vislib::Array<float> atomCharges;
    vislib::Array<float> bindingEnergy;
    vislib::Array<float> modelEfficiency;
    vislib::Array<vislib::math::Vector<uint, 2>> globalMdlID_to_ligID_mdlID;
    vislib::Array<vislib::Array<uint>> ligID_mdlID_to_globalMdlID;

    vislib::Array<vislib::Array<uint>> hBonds_ligMdlDonor;
    vislib::Array<vislib::Array<uint>> hBonds_ligMdlAcceptor;
    vislib::Array<vislib::Array<uint>> hBonds_ligMdlHDonor;
    vislib::Array<vislib::Array<uint>> hBonds_ligMdlHAcceptor;
    vislib::Array<uint> hBondAtomCnt_protein;
	
    std::vector<LigandModelCall::FunctionalGroupsLigand> functionalGroupsStorage;
    std::vector<LigandModelCall::ChemPropsLigand> chemPropsStorage;
    

    SIZE_T frameCapacity;


	void checkWaitPDBQT_data();

    /*********************
     * allocate the data arrays
     *********************/
    void allocateMultiPDBQTLoader(bool needToExtend);
    bool resetMultiPDBQTLoader();


    /*********************
     * MAIN DATA ACCESS: variables and functions
     *********************/

	

    void getLigandData(int lignadID);
    void getLigandData(int lignadID, int modelID);
    /* Returns:
     accessible through 'Ligand' class:
     * -> LigandID
     * -> ligandName
     * -> ligandAtomCnt
     * -> modelCnt
     accessible through 'Model' class:
     * -> atomTypesAD[i]
     * -> atomTypes[i]
     private variables:
     * -> atomIndexSum;
     * -> modelRangeStart
     * -> modelRangeEnd

    /* CALLS:
     * -> CALL: getModelDataOf_currentLigand(modelID) => default first model
     * -> CALL: scanAtomIndex();
     */

    bool getModelDataOf_currentLigand();
    /* Returns:
     accessible through 'Model' class:
     * -> ModelID
     * -> atomTypesAD[i]
     * -> atomTypes[i]
     * -> RMSD_lb
     * -> RMSD_ub
     * -> bindingEnergy
     * -> atomPostitons[i]
     * -> atomCharges[i]
     */


    const vislib::StringA& getFilename();

    void testPrinting(int LigID, int ModID);
    void calculateModelEfficiency(); 

	// data structures
    class PLIPForces {
    public:
        PLIPForces() {};
        // the different forces
        LigandModelCall::HProtAcceptors hProtAcceptors;
        LigandModelCall::HProtDonors hProtDonors;
        LigandModelCall::HydrophobicInteractions hydrophobicInteractions;
        LigandModelCall::SaltBridges saltBridges;
        LigandModelCall::PiStacks piStacks;
        LigandModelCall::PiCationInteractions piCationInteractions;
        LigandModelCall::HalogenBonds halogenBonds;
        LigandModelCall::MetalComplexes metalComplexes;
    };


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

    /**
     * Call callback to get the data
     *
     * @param c The calling call
     *
     * @return True on success
     */
    bool getData(core::Call& call);

    /**
     * Call callback to get the extent of the data
     *
     * @param c The calling call
     *
     * @return True on success
     */
    bool getExtent(core::Call& call);

    bool getLigandAndModelCount(core::Call& call);

    bool loadPDBQTFiles(const vislib::TString& filename);

    bool loadCheckmolFiles();


    /**
     * Call callback to get the data
     *
     * @param call LigandModelCall
     *
     * @return True on success
     */
    bool setRequestedLigandAndModel(core::Call& call);

    bool setRequestedLigandAndModelSilent(core::Call& call);

    bool isSilentDataCall;

    /**
     * Call callback to get the cluster data
     *
     * @param call LigandModelCall
     *
     * @return True on success
     */
    bool getClusterData(core::Call& call);
	
    // stored the currently selected ligand and model
    LigandModelCall::Ligand ligand;

	// store the data set by the web-interface
    LigandModelCall::WebData webData;

	LigandModelCall::ClusteredFunctionalGroups clusteredFunctionalGroups;


    /**
     * Callback function for DBSCAN params
     */
    bool paramsDBSCANchanged(megamol::core::param::ParamSlot& param);

private:
    /** params slots */

    // MultiPDBQTLoader
    core::param::ParamSlot pdbqtListFilenameSlot;

    // DBSCAN
    core::param::ParamSlot epsParam;
    core::param::ParamSlot minPtsParam;
    core::param::ParamSlot minBindEnergyParam;

    /** callee slots */
    core::CalleeSlot dataOutSlot_LigandModelCall;
    core::CalleeSlot dataOutSlot_MolecularDataCall;

    /** caller slots */
    megamol::core::CallerSlot dataInSlot_MolecularDataCall;

    /** The data hashes */
    SIZE_T mdc_dataHash;
    SIZE_T lmc_dataHash;

    /*********************
     * sizes of data arrays
     *********************/
    // frameCapacity == ModelCapacity

    // sizes controlled in allocateMultiPDBQTLoader (needToExtend)
    SIZE_T bindingEnergyCapacity = frameCapacity;
    SIZE_T RMSD_lbCapacity = frameCapacity;
    SIZE_T RMSD_ubCapacity = frameCapacity;
    SIZE_T atomPositionsCapacity = frameCapacity * 100 * 3;
    SIZE_T ligandsNameCapacity = frameCapacity;
    SIZE_T atomTypesCapacity = frameCapacity * 100;

    // size controlled in getAtomPositions_And_Charges
    SIZE_T singleLigand_AtomPositionsCapacity = 100 * 3;


    /*************************************
     * data arrays for MolecluarDataCall *
     *************************************/
    /** The array of residues */
    vislib::Array<megamol::protein_calls::MolecularDataCall::Residue*> residue;
    /** The array of residue type names */
    vislib::Array<vislib::StringA> residueTypeName;
    /** The array of molecules */
    vislib::Array<megamol::protein_calls::MolecularDataCall::Molecule> molecule;
    /** The array of chains */
    vislib::Array<megamol::protein_calls::MolecularDataCall::Chain> chain;
    /**
     * Stores the connectivity information (i.e. subsequent pairs of atom
     * indices)
     */
    vislib::Array<vislib::Array<unsigned int>> connectivity;

    void createHelperDataStructure() ;

    bool getfunctionalGroupDataFiles(const vislib::TString& filename);

    bool getProteinLigand_Interactions(const vislib::TString& filename);

    bool getSVGsAndChemProps(const vislib::TString& filename);

    const char* additionalDataFolderName;
    std::string additionalDataFolderPath;
    std::vector<std::string> SVGPaths;
    /*********************
     * getter functions (structured access to data arrays)
     *********************/
    void scanAtomIndex();
    void getModelCnt();
    uint getGlobalModelID(uint ligandID, uint modelID);
    void getLigandConnections(uint ligandID);
    void getAtomNamesAndTypesAD(uint ligandID);
    uint getLigandAtomCnt(uint ligandID);
    vislib::StringA getLigandsName(uint ligandID);
    float getRMSD_lb(uint globalModelID);
    float getRMSD_ub(uint globalModelID);
    float getBindingEnergy(uint globalModelID);
    float getModelEfficiency(uint globalModelID);
    void getAtomPositions_And_Charges();
    void getSVGPath(uint ligandID);
    void getFunctionalGroups(uint ligandID);
    void getChemProperties(uint ligandID);

	
	/***************
     * XML PARSING *
     ***************/

  
    /*****************************
     * internal helper functions *
     *****************************/

    Path getAdditionalDataFolderPath();

    std::string getFilname_from_filepath(std::string filepath);
    
	Path getDirpath_from_filepath(std::string filepath);

	Path getDirpath_toMegaMol_from_filepath(std::string filepath);


    /**********************************
     * variables for getter functions *
     **********************************/

    vislib::Array<int> atomIndexSum;
    uint modelRangeStart;
    uint modelRangeEnd;

    int ligandIdx;
    int modelIdx;

    float minCharge;
    float maxCharge;
    vislib::math::Cuboid<float> bbox;


    /** Variables */

    /** DBSCAN and Cluster Data
     *
     * eps:					DBSCAB param (neighbour search distance)
     * minPts:				DBSCAB param (minimum of found neighbours to be a clusters)
     * minBindEnergy:		the minimal negative binding energy used for DBSCAN [abs(minBindEnergy) => abs(-5) >
     * abs(-4)] 3 DBSCANParamsChanged: just indicates whether params have changed or not numberOfCluster:		amount
     * of found clusters modelCentroidsVec4:	stores the centroids of the models of the ligands [x,y,z,bindingEnergy]
     * modelCentroids:		stores the centroids of the models of the ligands [x,y,z,radius]
     * clusterCentroids:	stores the centroids of the found clusters [x,y,z,radius]
     * clusterIndexSum:		indexSum over the (the amount of elements) of the different clusters
     * clusterSizes:		stores for every cluster its size
     * assignedClusters:	stores for each model to which cluster it belongs to (IDx = model IDx) [cluster -1 means
     * noise] modelCentroidsCuboids:	growing cuboids used for determination of modelCentroids
     * clusterCentroidsCuboids:	growing cuboids used for determination of clusterCentroids
     * assign_modelID_sphereID: stores the clusterID .GetX() the modelID .GetY() and the sphereID as index [i]
     */
    float eps;
    int minPts;
    float minBindEnergy;
    bool DBSCANParamsChanged;
    int numberOfCluster;

    visFloat4 modelCentroidsVec4;
    vislib::Array<float> modelCentroids;
    vislib::Array<float> clusterCentroids;
    vislib::Array<int> clusterIndexSum;
    vislib::Array<int> clusterSizes;
    int minClusterSize;
    int maxClusterSize;
    int clusterDataHash;
    vislib::Array<int> assignedClusters;
    vislib::Array<vislib::math::Cuboid<float>> modelCentroidsCuboids;
    vislib::Array<vislib::math::Cuboid<float>> clusterCentroidsCuboids;
    vislib::Array<vislib::math::Vector<int, 2>> assign_modelID_sphereID;

	std::vector<PLIPForces> plipForcesStore;
    Path prolintSubPath;
};


/************************
 ***** DBSCAN CLASS *****
 ************************/
class DBSCAN {
public:
    /** Ctor. */
    DBSCAN();
    /** Dtor. */
    virtual ~DBSCAN() { release(); }

    // return the class name
    static const char* ClassName(void) { return "DBSCAN"; }

    // description of the class
    static const char* Description(void) {
        return "clusters data of a float vec3 array using DBSCAN and retruns the cluster assign array.";
    }

    /**
     * Implementation of DBSCAN adapted for spatial clustering (3D)
     * Pseudocode used: https://de.wikipedia.org/wiki/DBSCAN
     */
    vislib::Array<int> DBSCAN_3D(visFloat4& dataVec4, float eps, int minPts, float minValue, int showProgress = 0);

    static float distance2Vectors(vislib::math::Vector<float, 3> v1, vislib::math::Vector<float, 3> v2);

protected:
    /**
     * Implementation of 'Release'.
     */
    void release(void);

    void expandCluster(visFloat4& dataVec4, std::vector<uint>& N_joined, int cluster, int minPts, int showProgress);


    // variables for DBSCAN
    /** Neighbours of a centroid
     * neighbourList	is a 2D array of all neighbours of each point
     * unVistedPoints:	used to mark visited points
     * cluster:			set the cluster
     * minValue:		the abs(minimal value) to be a neighbour candidate
     */

    vislib::Array<int> clusters;
    vislib::Array<int> unVistedPoints;
    int unVistedPointsCnt;
    std::vector<std::vector<uint>> neighbourList_ffrnns;
    int oldProgress;
    float minValue;

};


/**********************************
 **** INTERACTION FORCES CLASS ****
 **********************************/
class InteractionForces {

public:
    /** Ctor. */
    InteractionForces(megamol::core::Call* ligandModelCall, megamol::core::Call* molecularDataCall_receptor,
        megamol::prolint::LigandModelCall::Ligand& ligand);
    /** Dtor. */
    virtual ~InteractionForces() { release(); }

    // return the class name
    static const char* ClassName(void) { return "InteractionForces"; }

    // description of the class
    static const char* Description(void) {
        return "calculates the various interaction forces between protein/receptor and lignads.";
    }

    // functions
    bool calcHBonds(vislib::Array<vislib::Array<uint>>& out_ligMdlTo_prot,
    vislib::Array<vislib::Array<uint>>& out_protTo_LigMdl, vislib::Array<vislib::Array<uint>>& out_ligMdlTo_protH,
    vislib::Array<vislib::Array<uint>>& out_protTo_LigMdlH, vislib::Array<uint>& hBondAtomCnt_protein);

	bool loadPLIPdata(Path filepath, std::vector<MultiPDBQTLoader::PLIPForces>& plipForcesStore);


protected:
    /**
     * Implementation of 'Release'.
     */
    void release(void);


    // variables for H-bond search
    megamol::prolint::LigandModelCall::ClusterAndCentroidData clusterData;
    megamol::prolint::LigandModelCall::Ligand* ligand;

    std::vector<uint> H_donorIndices_prot;
    std::vector<uint> donorIndices_prot;
    std::vector<uint> acceptorIndices_prot;

    visFloat4 proteinData_H_donors;
    visFloat4 proteinData_donors;
    visFloat4 proteinData_acceptors;

    std::vector<uint> H_donorIndices_ligMdl;
    std::vector<uint> donorIndices_ligMdl;
    std::vector<uint> acceptorIndices_ligMdl;

    visFloat4 ligandModelData_H_donors;
    visFloat4 ligandModelData_donors;
    visFloat4 ligandModelData_acceptors;

    std::vector<std::vector<uint>> Hbonds_protTo_LigMdl;
    std::vector<std::vector<uint>> Hbonds_LigMdlTo_prot;

    std::vector<std::vector<uint>> Hbonds_protTo_LigMdl_checked;
    std::vector<std::vector<uint>> Hbonds_LigMdlTo_prot_checked;

	


private:
    megamol::prolint::LigandModelCall* lmc;
    megamol::protein_calls::MolecularDataCall* pmdc;
};


} // namespace prolint
} /* end namespace megamol */

int sortVec2ByX(const vislib::math::Vector<int, 2U>& lhs, const vislib::math::Vector<int, 2U>& rhs);

inline bool fileExists(const std::string& name) {
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

#endif // PROLINT_PLUGIN_PDBQT_H_INCLUDED
