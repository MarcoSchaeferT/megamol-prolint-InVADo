/*
 * LigandModelCall.h
 *
 * Copyright (C) 2010 by University of Stuttgart (VISUS).
 * All rights reserved.
 */

#ifndef MEGAMOL_PROLINT_CALL_LIGANDMODELCALL_H_INCLUDED
#define MEGAMOL_PROLINT_CALL_LIGANDMODELCALL_H_INCLUDED
#if (defined(_MSC_VER) && (_MSC_VER > 1000))
#    pragma once
#endif /* (defined(_MSC_VER) && (_MSC_VER > 1000)) */

#include <glm/glm.hpp>
#include <map>
#include <vector>
#include "mmcore/AbstractGetDataCall.h"
#include "mmcore/Call.h"
#include "mmcore/factories/CallAutoDescription.h"
#include "protein_calls/MolecularDataCall.h"
#include "vislib/Array.h"
#include "vislib/IllegalParamException.h"
#include "vislib/String.h"
#include "vislib/macro_utils.h"
#include "vislib/math/Vector.h"


typedef unsigned int uint;

namespace megamol {
namespace prolint {

/**
 * TODO describe this call
 */

class LigandModelCall : public megamol::core::AbstractGetDataCall {
public:
    /** Index of the 'GetData' function */
    static const unsigned int CallForGetData;

    /** Index of the 'GetExtent' function */
    static const unsigned int CallForGetExtent;

    /** Index of the 'GetDataSilent' function */
    static const unsigned int CallForGetDataSilent;

    /** Index of the 'GetClusterAndCentroidData' function */
    static const unsigned int CallForGetClusterAndCentroidData;

    /** Index of the 'SetLMCStatus' function */
    static const unsigned int CallForSetLMCStatus;

    /** Index of the 'SetLMCStatus' function */
    static const unsigned int CallForSetWebData;

    /** Index of the 'SetFunctionalGroupData' function */
    static const unsigned int CallForSetFunctionalGroupData;

    /**
     * Answer the name of the objects of this description.
     *
     * @return The name of the objects of this description.
     */
    static const char* ClassName(void) { return "LigandModelCall"; }

    /**
     * Gets a human readable description of the module.
     *
     * @return A human readable description of the module.
     */
    static const char* Description(void) { return "Call to get/set ligand and model data/ID"; }

    /**
     * Answer the number of functions used for this call.
     *
     * @return The number of functions used for this call.
     */
    static unsigned int FunctionCount(void) { return 5; }

    /**
     * Answer the name of the function used for this call.
     *
     * @param idx The index of the function to return it's name.
     *
     * @return The name of the requested function.
     */
    static const char* FunctionName(unsigned int idx) {
        switch (idx) {
        case 0: return "GetData";
        case 1: return "GetExtent";
        case 2: return "GetDataSilent";
        case 3: return "GetClusterAndCentroidData";
        case 4: return "SetLMCStatus";
        default:return "";
        }
    }

    /** Ctor. */
    LigandModelCall(void);
    /** Dtor. */
    virtual ~LigandModelCall(void);

    /** DBSCAN and Cluster Data
     *
     * eps:					DBSCAB param (neighbour search distance)
     * minPts:				DBSCAB param (minimum of found neighbours to be a clusters)
     * minBindEnergy:		the minimal negative binding energy used for DBSCAN [abs(minBindEnergy) => abs(-5) >
     * abs(-4)] 3 DBSCANParamsChanged: just indicates whether params have changed or not numberOfCluster:		amount
     * of found
     */
    struct DBSCANParamters {
        float eps;
        int minPts;
        float minBindEnergy;
        bool paramsChanged;
        // set mark when data changed to send new data to python socket
        bool socketSend = false;
    };


    /** Cluster Data (here only pointers to these data arrays)
     * clusters modelCentroidsVec4:	stores the centroids of the models of the ligands [x,y,z,bindingEnergy]
     * modelCentroids:		stores the centroids of the models of the ligands [x,y,z,radius]
     * clusterCentroids:	stores the centroids of the found clusters [x,y,z,radius]
     * clusterIndexSum:		indexSum over the (the amount of elements) of the different clusters
     * assignedClusters:	stores for each model to which cluster it belongs to (IDx = model IDx) [cluster -1 means
     * noise] modelCentroidsCuboids:	growing cuboids used for determination of modelCentroids
     * clusterCentroidsCuboids:	growing cuboids used for determination of clusterCentroids
     */
    struct ClusterAndCentroidData {
        vislib::Array<float>* modelCentroids;
        vislib::Array<float>* clusterCentroids;
        vislib::Array<int>* assignedClusters;
        vislib::Array<int>* clusterIndexSum;
        vislib::Array<int>* clusterSizes;
        int minClusterSize;
        int maxClusterSize;
        vislib::Array<vislib::math::Cuboid<float>>* modelCentroidsCuboids;
        vislib::Array<vislib::math::Cuboid<float>>* clusterCentroidsCuboids;
        // saves for every clusterID the modelID
        vislib::Array<vislib::math::Vector<int, 2>>* assign_modelID_sphereID;
        struct DBSCANParamters DBSCANParams;
        int clusterDataHash = 0;
    };


    /*************************************/
    /********** DATA STRUCTURES **********/
    /*************************************/

    //*** Ligand::Model::forces::HProtAcceptors Class *** //
    class HProtAcceptors {
    public:
        HProtAcceptors() { this->cnt = 0; };

        // variables
        int cnt;
        std::vector<glm::vec3> ligPos;
        std::vector<glm::vec3> protPos;
        std::vector<int> protIDx;
        std::vector<int> ligIDx;
        // functions
        int getNumberInteractions() { return this->cnt; }
        glm::vec3 getLigPos(int i) { return this->ligPos[i]; }
        glm::vec3 getProtPos(int i) { return this->protPos[i]; }
        int getLigID(int i) { return this->ligIDx[i]; }
        int getProtID(int i) { return this->protIDx[i]; }
    };

    //*** Ligand::Model::forces::HProtDonors Class *** //
    class HProtDonors {
    public:
        HProtDonors() = default;

        // variables
        int cnt = 0;
        std::vector<glm::vec3> ligPos;
        std::vector<glm::vec3> protPos;
        std::vector<int> protIDx;
        std::vector<int> ligIDx;
        // functions
        int getNumberInteractions() const { return this->cnt; }
        // check if return of const glm::vec3& is possible
        glm::vec3 getLigPos(int i) const { return this->ligPos[i]; }
        glm::vec3 getProtPos(int i) const { return this->protPos[i]; }
        int getLigID(int i) const { return this->ligIDx[i]; }
        int getProtID(int i) const { return this->protIDx[i]; }
    };


    //*** Ligand::Model::forces::HydrophobicInteractions Class *** //
    class HydrophobicInteractions {
    public:
        HydrophobicInteractions() { this->cnt = 0; };

        // variables
        int cnt;
        std::vector<glm::vec3> ligPos;
        std::vector<glm::vec3> protPos;
        std::vector<int> protIDx;
        std::vector<int> ligIDx;
        // functions
        int getNumberInteractions() { return this->cnt; }
        glm::vec3 getLigPos(int i) { return this->ligPos[i]; }
        glm::vec3 getProtPos(int i) { return this->protPos[i]; }
        int getLigID(int i) { return this->ligIDx[i]; }
        int getProtID(int i) { return this->protIDx[i]; }
    };

    //*** Ligand::Model::forces::SaltBridges Class *** //
    class SaltBridges {
    public:
        SaltBridges() { this->cnt = 0; };

        // variables
        int cnt;
        std::vector<glm::vec3> ligPos;
        std::vector<glm::vec3> protPos;
        std::vector<std::vector<int>> protIDxList;
        std::vector<std::vector<int>> ligIDxList;

        // functions
        int getNumberInteractions() { return this->cnt; }
        glm::vec3 getLigPos(int i) { return this->ligPos[i]; }
        glm::vec3 getProtPos(int i) { return this->protPos[i]; }

        int getLigIDsCnt(int i) { return this->ligIDxList[i].size(); }
        int* getLigIDs(int i) { return this->ligIDxList[i].data(); }

        int getProtIDsCnt(int i) { return this->protIDxList[i].size(); }
        int* getProtIDs(int i) { return this->protIDxList[i].data(); }
    };

    //*** Ligand::Model::forces::PiStacks Class *** //
    class PiStacks {
    public:
        PiStacks() { this->cnt = 0; };

        // variables
        int cnt;
        std::vector<glm::vec3> ligPos;
        std::vector<glm::vec3> protPos;
        std::vector<std::vector<int>> protIDxList;
        std::vector<std::vector<int>> ligIDxList;

        // functions
        int getNumberInteractions() { return this->cnt; }
        glm::vec3 getLigPos(int i) { return this->ligPos[i]; }
        glm::vec3 getProtPos(int i) { return this->protPos[i]; }

        int getLigIDsCnt(int i) { return this->ligIDxList[i].size(); }
        int* getLigIDs(int i) { return this->ligIDxList[i].data(); }

        int getProtIDsCnt(int i) { return this->protIDxList[i].size(); }
        int* getProtIDs(int i) { return this->protIDxList[i].data(); }
    };

    //*** Ligand::Model::forces::PiCationInteractions Class *** //
    class PiCationInteractions {
    public:
        PiCationInteractions() { this->cnt = 0; };

        // variables
        int cnt;
        std::vector<glm::vec3> ligPos;
        std::vector<glm::vec3> protPos;
        std::vector<std::vector<int>> protIDxList;
        std::vector<std::vector<int>> ligIDxList;

        // functions
        int getNumberInteractions() { return this->cnt; }
        glm::vec3 getLigPos(int i) { return this->ligPos[i]; }
        glm::vec3 getProtPos(int i) { return this->protPos[i]; }

        int getLigIDsCnt(int i) { return this->ligIDxList[i].size(); }
        int* getLigIDs(int i) { return this->ligIDxList[i].data(); }

        int getProtIDsCnt(int i) { return this->protIDxList[i].size(); }
        int* getProtIDs(int i) { return this->protIDxList[i].data(); }
    };

    //*** Ligand::Model::forces::HalogenBonds Class *** //
    class HalogenBonds {
    public:
        HalogenBonds() { this->cnt = 0; };

        // variables
        int cnt;
        std::vector<glm::vec3> ligPos;
        std::vector<glm::vec3> protPos;
        std::vector<int> protIDx;
        std::vector<int> ligIDx;
        // functions
        int getNumberInteractions() { return this->cnt; }
        glm::vec3 getLigPos(int i) { return this->ligPos[i]; }
        glm::vec3 getProtPos(int i) { return this->protPos[i]; }
        int getLigID(int i) { return this->ligIDx[i]; }
        int getProtID(int i) { return this->protIDx[i]; }
    };

    //*** Ligand::Model::forces::MetalComplexes Class *** //
    class MetalComplexes {
    public:
        MetalComplexes() { this->cnt = 0; };

        // variables
        int cnt;
        std::vector<glm::vec3> metalPos;
        std::vector<glm::vec3> targetPos;
        std::vector<int> metalIDx;
        std::vector<int> targetIDx;
        // functions
        int getNumberInteractions() { return this->cnt; }
        glm::vec3 getMetalPos(int i) { return this->metalPos[i]; }
        glm::vec3 getTargetPos(int i) { return this->targetPos[i]; }
        int getMetalID(int i) { return this->metalIDx[i]; }
        int getTargetID(int i) { return this->targetIDx[i]; }
    };


public:
    //*** Ligand::Model::forces Class *** //
    class Forces {
    public:
        Forces() = default;
        // the different forces
        HProtAcceptors* hProtAcceptors;
        HProtDonors* hProtDonors;
        HydrophobicInteractions* hydrophobicInteractions;
        SaltBridges* saltBridges;
        PiStacks* piStacks;
        PiCationInteractions* piCationInteractions;
        HalogenBonds* halogenBonds;
        MetalComplexes* metalComplexes;
        bool useMegaMolHBonds = false;
    };

    //*** Ligand::Model::FunctionalGroups::fGroup Class *** //
    class fGroup {
    public:
        fGroup(){};
        // the ID which assigns the specific functional group name
        int groupID;
        // the number of occurrences of this functional group
        int speficGroupCnt;
        // stores the atom IDx (ligand) of the central functional group atom
        std::vector<int> centralAtoms;
    };

    //*** Ligand::Model::FunctionalGroups Class *** //
    class FunctionalGroupsLigand {
    public:
        FunctionalGroupsLigand(){};
        // the total number of functional groups of a molecule (ligand)
        int totalGroupCnt = 0;
        // a vector storing functional group objects (for each group type)
        std::vector<fGroup> fGroups;
    };
    // the vector of strings which assigns the functional groups names to group IDs
    std::vector<std::string> functionalGroupsIDsToWord;
    // store the hierarchy of the functional groups. eg.: amine -> secondary amine -> etc.
    std::vector<int> fgsHierarchy;

    //*** Ligand::Model::ChemPropsLigand Class *** //
    class ChemPropsLigand {
    public:
        ChemPropsLigand() {
            this->mwt = 0.0f;
            this->logp = 0.0f;
            this->fractioncsp3 = 0.0f;
        };
        // molecuilar weigth
        float mwt;
        // logP is a quantitative measure of lipophilicit
        float logp;
        // parameter for drug-likeness
        float fractioncsp3;
    };

    //*** Ligand::Model Class *** //
    class Model {
    public:
        Model(){};
        int modelID;
        uint globalModelID;
        vislib::StringA* atomTypesAD;
        uint* atomTypeIdx;
        vislib::Array<megamol::protein_calls::MolecularDataCall::AtomType>* atomTypes;
        float RMSD_lb;
        float RMSD_ub;
        float bindingEnergy;
        float* atomPostions;
        float* atomCharges;
        float modelEfficiency;

        Forces Forces;
    };

    //*** Ligand Class *** //
    class Ligand {

    public:
        Ligand() = default;
        int ligandID;
        std::string ligandName;
        uint ligandAtomCnt;
        uint modelCnt;
        uint connectionCnt;
        const uint* connections;
        std::string SVGPath = "";

        Model Model;
        FunctionalGroupsLigand* functionalGroupsLigand;
        ChemPropsLigand* chemPropsLigand;
    };


    //*** LmcStatus Class *** //
    class LmcStatus {
    public:
        LmcStatus() = default;

        bool isModelClustRendererBusy;
        bool isSocketCommManagerBusy;
        bool isFirstRenderCallComplete;
    };


    //*** WebData Class *** //

    class WebData {
    public:
        WebData() {
            this->dataHash = 0;
            this->selection_dataHash = 0;
            this->GUIdataHash = 0;

            this->globalMdlID = -1;
            this->ligandID = -1;
            this->modelID = -1;
            this->clusterID = -1;
            this->selectedFgsClusterID = -1;
        };
        ~WebData(){};


    private:
        // data hash
        int dataHash;
        int selection_dataHash;
        int GUIdataHash;

        // selections
        int globalMdlID;
        int ligandID;
        int modelID;
        int clusterID;

    public:
        // points to the stored GUI param values to control render data
        std::map<std::string, float>* GUI = new std::map<std::string, float>();

        // stores per fgsTypeID a vector of global model IDs o the current slected fgsCluster
        std::map<uint, std::vector<uint>> curFgsMap;
        int selectedFgsClusterID;


        // setters
        void setClusterID(int val) {
            this->clusterID = val;
            UpdateSelectionData();
        }

        void setLigandID(int val) {
            this->ligandID = val;
            UpdateSelectionData();
        }

        void setModelID(int val) {
            this->modelID = val;
            UpdateSelectionData();
        }

        void setGlobalMdlID(int val) {
            this->globalMdlID = val;
            // this->GUIandSelection_dataHash++;
        }


        // getters
        int getWebDataHash() { return this->dataHash; }
        int getClusterID() { return this->clusterID; }
        int getLigandID() { return this->ligandID; }
        int getModelID() { return this->modelID; }
        int getGlobalMdlID() { return this->globalMdlID; }
        int getSelectionDataHash() { return this->selection_dataHash; }
        int getGUIdataHash() { return this->GUIdataHash; }


        /* Functions */
        void UpdateWebData();
        void UpdateGUIData();
        void UpdateSelectionData();
    };


    /***************************/
    /***** DATA STRUTURES ******/
    /***************************/

    class FGS_structN {
    public:
        FGS_structN(){};
        ~FGS_structN(){};

        int clusterCnt;
        int clusterID;

        // stores the amount of points for each cluster
        std::vector<int> clusterSizes;

        // the assignment of which dataPoint belongs to which cluster
        std::vector<int> clusterAssignment;

        // the centerPoint of each cluster
        vislib::Array<float> clusterCentroids;

        // the cuboid of each cluster
        std::vector<vislib::math::Cuboid<float>> clusterCuboids;

        // vector storing for each cluster a map<FGStypeID, vector of gblMdlIDs>
        std::vector<std::map<uint, std::vector<uint>>> fgsMaps_GblMDlIDs;
        std::vector<std::map<uint, std::vector<uint>>> fgsMaps_centralAtomArrayIDx;
    };


    // Class FunctionalGroupData
    /* stores for each pocket the grouped/clustered functional groups*/
    class ClusteredFunctionalGroups {
    public:
        ClusteredFunctionalGroups() {
            this->dataHash = -1;

            // number corresponds to cluster count
            this->functionalGroupSetsCount = 0;
        };
        ~ClusteredFunctionalGroups(){};

        // setters
        void setPointer_to_functionalGroupClusterSets(
            FGS_structN* functionalGroupClusterSets, int functionalGroupClusterSetsSize);


        // getters
        FGS_structN* getPointer_to_functionalGroupClusterSets() { return functionalGroupClusterSets; };

        int getfunctionalGroupSetsCount() { return this->functionalGroupSetsCount; };

        int getDataHash() { return this->dataHash; };


    private:
        int dataHash;
        int functionalGroupSetsCount;
        FGS_structN* functionalGroupClusterSets;
    };


    // -------------------- get and set routines --------------------

    /**
     * Set the current ligand and model
     *
     * @param curLig The current ligand
     * @param curMdl The current model
     * @return 'false' if the current ligand or model is out of range, 'true' otherwise.
     */
    bool SetCurrentLigandAndModel(unsigned int curLig, unsigned int curMdl);


    /**
     * Set the global model ID of a model for direct data access
     *
     * info: (Calls in a later step SetCurrentLigandAndModel())
     * @param globalMdlID The current global model ID
     * @return 'false' if the current global ID is out of range (>total_modelCnt), 'true' otherwise.
     */
    bool SetCurrentLigandAndModel_byGlobalMdlID(unsigned int globalMdlID);

    /**
     * Set the total number of ligands and total number of models and number of models of current ligand
     *
     * @param ligCnt The number of ligands
     * @param mdlCnt The number of models per ligand
     * @param totalMdlCnt The number of all models of all ligands
     */
    void SetLigandAndModelCount(unsigned int ligCnt, const unsigned int* mdlCnt, unsigned int totalMdlCnt);
    void SetTotalAtomPosCnt(unsigned int totalAtomPosCnt);
    void SetBindigenergy(float bindingenergy);
    void SetCurrentGlobalModelID(unsigned int globalMdlID) { this->currentGlobalMdlID = globalMdlID; };
    void SetLigandname(vislib::StringA ligandname);
    void SetPDBQT_list_Filepath(vislib::StringA PDBQT_list_Filepath);
    void SetGlobalMdlID_to_ligID_mdlID(vislib::Array<vislib::math::Vector<uint, 2>>* globalMdlID_to_ligID_mdlID);
    void SetSVGPath(std::string SVGPath);
    void SetCurrentModelEfficiency(float modelEfficiency);
    void SetClusterAndCentroidData(vislib::Array<float>* modelCentroids, vislib::Array<float>* clusterCentroids,
        vislib::Array<int>* assignedClusters, vislib::Array<int>* clusterIndexSum,
        vislib::Array<vislib::math::Cuboid<float>>* modelCentroidsCuboids,
        vislib::Array<vislib::math::Cuboid<float>>* clusterCentroidsCuboids,
        vislib::Array<vislib::math::Vector<int, 2>>* assign_modelID_sphereID, vislib::Array<int>* clusterSizes,
        int minClusterSize, int maxClusterSize, int clusterDataHash);

    void SetInteractionForces_hbonds_protAcceptor(HProtAcceptors* hProtAcceptor);
    void SetInteractionForces_hbonds_protDonor(HProtDonors* hProtDonor);
    void SetInteractionForces_hydrophobicInteractions(HydrophobicInteractions* hydrophobicInteractions);
    void SetInteractionForces_saltBridges(SaltBridges* saltBridges);
    void SetInteractionForces_piStacks(PiStacks* piStacks);
    void SetInteractionForces_piCationInteractions(PiCationInteractions* piCationInteractions);
    void SetInteractionForces_halogenBonds(HalogenBonds* halogenBonds);
    void SetInteractionForces_metalComplexes(MetalComplexes* metalComplexes);

    void SetFunctionalGroupsOfLigand(FunctionalGroupsLigand* functionalGroupsLigand);

    void SetChemPropsOfLigand(ChemPropsLigand* chemPropsLigand);

    void SetDBSCANParams(float eps, int minPts, float minBindEnergy);

    void SetlmcStatus(bool modelClustRenderer_busy, bool socketCommManager_busy, bool firstRenderCall_complete);

    void ResetDBSCANChanged();

    void SetWebData(WebData* webData);

    void SetClusteredFunctionalGroups(ClusteredFunctionalGroups* clusteredFunctionalGroups);

    // get ID of current selected ligand
    unsigned int getCurrentLigandID() { return currentLigand; };
    // get ID of current selected model of current ligand
    unsigned int getCurrentModelID() { return currentModel; };
    // get total numnber of ligands
    unsigned int getLigandCount() { return ligandCount; };
    // get numnber of models of current ligand
    const unsigned int* getModelCount() { return modelCount; };
    // get total number of all models of all ligands
    unsigned int getTotalModelCount() { return totalModelCount; };

    // get the unique global model ID of model of a ligand
    int getCurrentGlobalModelID() { return currentGlobalMdlID; };

    // get total numnber of all stored atom positions
    unsigned int getTotalAtomPosCount() { return total_atomPosCnt; };

    float getCurrentBindingenergy() { return currentBindingenrgy; };

    vislib::StringA getCurrentLigandname() { return currentLigandname; };

    const vislib::StringA getPDBQT_list_Filepath() { return PDBQT_list_Filepath; };

    std::string getSVGPath() { return SVGPath; };

    float getCurrentModelEfficiency() { return this->currentModelEfficiency; }

    Forces getInteractionForces() { return this->forces; };

    FunctionalGroupsLigand* getFunctionalGroupsOfCurLigand() { return this->functionalGroupsOfCurLigand; };

    ChemPropsLigand* getChemPorpsOfCurLigand() { return this->chemPropsOfCurLigand; };

    const ClusterAndCentroidData& getClusterAndCentroidData() const { return this->c; };

    LmcStatus getLMCStatus() { return this->lmcStatus; };

    WebData* getWebData() { return this->webData; };

    ClusteredFunctionalGroups* getClusteredFunctionalGroups() { return this->clusteredFunctionalGroups; };


    inline LigandModelCall& operator=(const LigandModelCall& s) {
        AbstractGetDataCall::operator=(s);
        // add variables
        this->currentLigand = s.currentLigand;
        this->currentModel = s.currentModel;
        this->ligandCount = s.ligandCount;
        this->modelCount = s.modelCount;
        this->totalModelCount = s.totalModelCount;
        this->total_atomPosCnt = s.total_atomPosCnt;
        this->currentGlobalMdlID = s.currentGlobalMdlID;
        this->currentBindingenrgy = s.currentBindingenrgy;
        this->currentLigandname = s.currentLigandname;
        this->PDBQT_list_Filepath = s.PDBQT_list_Filepath;
        this->current_hBonds_ligMdlDonorCount = s.current_hBonds_ligMdlDonorCount;
        this->current_hBonds_ligMdlDonor = s.current_hBonds_ligMdlDonor;
        this->current_hBonds_ligMdlAcceptorCount = s.current_hBonds_ligMdlAcceptorCount;
        this->current_hBonds_ligMdlAcceptor = s.current_hBonds_ligMdlAcceptor;
        this->SVGPath = s.SVGPath;
        this->globalMdlID_to_ligID_mdlID = s.globalMdlID_to_ligID_mdlID;
        this->c = s.c;
        this->lmcStatus = s.lmcStatus;
        this->forces = s.forces;
        this->functionalGroupsOfCurLigand = s.functionalGroupsOfCurLigand;
        this->chemPropsOfCurLigand = s.chemPropsOfCurLigand;
        this->globalMdlID_setByWeb = s.globalMdlID_setByWeb;
        this->webData = s.webData;
        this->clusteredFunctionalGroups = s.clusteredFunctionalGroups;
        this->currentModelEfficiency = s.currentModelEfficiency;
        return *this;
    }

private:
    // -------------------- variables --------------------
    unsigned int currentLigand;
    unsigned int currentModel;
    unsigned int ligandCount;
    float currentBindingenrgy;
    vislib::StringA currentLigandname;
    float currentModelEfficiency;
    const unsigned int* modelCount;
    unsigned int totalModelCount;
    unsigned int total_atomPosCnt;
    int currentGlobalMdlID;
    int globalMdlID_setByWeb;
    vislib::TString PDBQT_list_Filepath;
    uint current_hBonds_ligMdlDonorCount;
    uint* current_hBonds_ligMdlDonor;
    uint current_hBonds_ligMdlAcceptorCount;
    uint* current_hBonds_ligMdlAcceptor;
    std::string SVGPath;
    vislib::Array<vislib::math::Vector<uint, 2>>* globalMdlID_to_ligID_mdlID;


    // -------------------- sturctures --------------------
    ClusterAndCentroidData c;
    Forces forces;
    FunctionalGroupsLigand* functionalGroupsOfCurLigand;
    ChemPropsLigand* chemPropsOfCurLigand;
    LmcStatus lmcStatus;
    WebData* webData;
    ClusteredFunctionalGroups* clusteredFunctionalGroups;
};

/** Description class typedef */
typedef megamol::core::factories::CallAutoDescription<LigandModelCall> LigandModelCallDescription;


} // namespace prolint
} /* end namespace megamol */

#endif /* MEGAMOL_PROLINT_CALL_LIGANDMODELCALL_H_INCLUDED */
