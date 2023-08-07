/*
 * MultiPDBQTLoader.cpp
 *
 * Copyright (C) 2019 by University of Tuebingen (BDVA).
 * All rights reserved.
 */


#include "stdafx.h"
#include "MultiPDBQTLoader.h"
#include <fmt/core.h>
#include <fstream>
#include <iostream>
#include "FFRNNS.h"
#include "protein/PDBLoader.h"

#include <Python.h>

#include <filesystem>

bool verbose = false;

using namespace megamol;
using namespace megamol::core;
using namespace megamol::prolint;
using namespace megamol::protein;
using namespace megamol::protein_calls;

    // https://www.oreilly.com/library/view/c-cookbook/0596007612/ch04s07.html // Marco: added more comments
// splits a string -> result stored in a vector
void megamol::prolint::split(const string& s, const char* c, std::vector<string>& v) {

    string::size_type i = 0;
    string::size_type j = s.find(*c);

    while (j != string::npos) {
        v.push_back(s.substr(i, j - i));
        i = ++j;
        j = s.find(*c, j); // if there is no match 'string::npos' will be returned (-1)
        if (j == string::npos) v.push_back(s.substr(i, s.length())); // add the last part after the last split character
    }
}

/**
 * Constructor
 */
MultiPDBQTLoader::MultiPDBQTLoader(void)
    : Module()
    , total_ligandCnt(0)
    , total_modelCnt(0)
    , total_atomCnt(0)
    , total_atomPosCnt(0)
    , PDBQTLoadingStatus(false)
    , MPL_modelClustRenderer_busy(true)
    , MPL_SocketCommManager_busy(false)
    , MPL_firstRenderCall_complete(false)
    , frameCapacity(1000)
    , pdbqtListFilenameSlot("pdbqtListFilename", "The path to the List of PDBQT data files to be loaded")
    , epsParam("DBSCAN::eps", "radius DBSCAN")
    , minPtsParam("DBSCAN::minPts", "minimum of points per cluster DBSCAN")
    , minBindEnergyParam("DBSCAN::minBindEnergy", "minimal binding energy")
    , dataOutSlot_LigandModelCall(
          "LigandModelCall_out", "The slot providing the ligand and model data of docking results (without atom data).")
    , dataOutSlot_MolecularDataCall(
          "MolecularDataCall_out", "The slot is providing the atom related data of one model of one ligand")
   , dataInSlot_MolecularDataCall(
          "MolecularDataCall_in", "The slot is providing molecular data of the receptor protein.")
    , eps(2.f)
    , minPts(10)
    , minBindEnergy(-4.0f)
    , mdc_dataHash(0)
    , lmc_dataHash(0)
    , ligandIdx(0)
    , modelIdx(0)
    , additionalDataFolderName("InVADo_tmpFiles")
    , DBSCANParamsChanged(TRUE)
    , minCharge(FLT_MAX)
    , maxCharge(FLT_MIN)
    , isSilentDataCall(FALSE)
    , minClusterSize(INT_MAX)
    , maxClusterSize(INT_MIN)
    , modelRangeStart(0)
    , modelRangeEnd(0)
    , numberOfCluster(0)
    , clusterDataHash(0) {

    Py_Initialize();

    this->prolintSubPath = Path("plugins\\prolint").make_preferred();
    // set up callee slots
    this->dataOutSlot_MolecularDataCall.SetCallback(MolecularDataCall::ClassName(),
        MolecularDataCall::FunctionName(MolecularDataCall::CallForGetData), &MultiPDBQTLoader::getData);
    this->dataOutSlot_MolecularDataCall.SetCallback(MolecularDataCall::ClassName(),
        MolecularDataCall::FunctionName(MolecularDataCall::CallForGetExtent), &MultiPDBQTLoader::getExtent);
    this->MakeSlotAvailable(&this->dataOutSlot_MolecularDataCall);

    // set up caller slots
    this->dataInSlot_MolecularDataCall.SetCompatibleCall<MolecularDataCallDescription>();
    this->MakeSlotAvailable(&this->dataInSlot_MolecularDataCall);

    // CallForGetData = 0;
    this->dataOutSlot_LigandModelCall.SetCallback(LigandModelCall::ClassName(),
        LigandModelCall::FunctionName(LigandModelCall::CallForGetData), &MultiPDBQTLoader::setRequestedLigandAndModel);

    // CallForGetExtent = 1;
    this->dataOutSlot_LigandModelCall.SetCallback(LigandModelCall::ClassName(),
        LigandModelCall::FunctionName(LigandModelCall::CallForGetExtent), &MultiPDBQTLoader::getLigandAndModelCount);

    // CallForGetDataSilent = 2;
    this->dataOutSlot_LigandModelCall.SetCallback(LigandModelCall::ClassName(),
        LigandModelCall::FunctionName(LigandModelCall::CallForGetDataSilent),
        &MultiPDBQTLoader::setRequestedLigandAndModelSilent);

    // CallForGetClusterAndCentroidData = 3;
    this->dataOutSlot_LigandModelCall.SetCallback(LigandModelCall::ClassName(),
        LigandModelCall::FunctionName(LigandModelCall::CallForGetClusterAndCentroidData),
        &MultiPDBQTLoader::getClusterData);

    // CallForSetLMCStatus = 4;
    this->dataOutSlot_LigandModelCall.SetCallback(LigandModelCall::ClassName(),
        LigandModelCall::FunctionName(LigandModelCall::CallForSetLMCStatus), &MultiPDBQTLoader::setLMCStatus);


    this->MakeSlotAvailable(&this->dataOutSlot_LigandModelCall);

    // set up MultiPDBQTLoader filename slot
    this->pdbqtListFilenameSlot << new param::FilePathParam("");
    this->MakeSlotAvailable(&this->pdbqtListFilenameSlot);

    // set up DBSACAN parameters
    this->epsParam << new param::FloatParam(this->eps);
    this->epsParam.SetUpdateCallback(&MultiPDBQTLoader::paramsDBSCANchanged);
    this->MakeSlotAvailable(&this->epsParam);

    this->minPtsParam << new param::IntParam(this->minPts);
    this->minPtsParam.SetUpdateCallback(&MultiPDBQTLoader::paramsDBSCANchanged);
    this->MakeSlotAvailable(&this->minPtsParam);

    this->minBindEnergyParam << new param::FloatParam(this->minBindEnergy);
    this->minBindEnergyParam.SetUpdateCallback(&MultiPDBQTLoader::paramsDBSCANchanged);
    this->MakeSlotAvailable(&this->minBindEnergyParam);
}

/**
 * Destructor
 */
MultiPDBQTLoader::~MultiPDBQTLoader(void) { this->release(); }

bool MultiPDBQTLoader::create(void) {
    // TODO allocate variables etc.
    return true;
}

void MultiPDBQTLoader::release(void) {
    // destroy variables etc.
    resetMultiPDBQTLoader();
    Py_Finalize();
}

/*
 * MultiPDBQTLoader::getData
 */
bool MultiPDBQTLoader::getData(core::Call& call) {
    using megamol::core::utility::log::Log;

	checkWaitPDBQT_data();
    auto* mdc = dynamic_cast<MolecularDataCall*>(&call);
    if (mdc == nullptr) return false;

    if (this->pdbqtListFilenameSlot.IsDirty() || this->total_ligandCnt < 1) {

        if (this->loadPDBQTFiles(this->pdbqtListFilenameSlot.Param<core::param::FilePathParam>()->Value())) {
            mdc->SetPDBFilename(this->pdbqtListFilenameSlot.Param<core::param::FilePathParam>()->Value());
            this->pdbqtListFilenameSlot.ResetDirty();
        } else {
            return false;
        }
    }

    // this->testPrinting();
    // Sleep(100000);

    // if (ligandIdx != this->ligand.ligandID || modelIdx != this->ligand.Model.modelID) {
    //	testPrinting(ligandIdx, modelIdx);
    //}

    // get current ligand and model
    this->getLigandData(this->ligandIdx, this->modelIdx);

    // write currently selected ligand model to mdc
    mdc->SetDataHash(this->mdc_dataHash);
    mdc->SetFrameCount(1);

    MolecularDataCall::Residue r();

    mdc->SetAtoms(ligand.ligandAtomCnt, this->atomTypes.Count(), ligand.Model.atomTypeIdx, ligand.Model.atomPostions,
        this->atomTypes.PeekElements(), nullptr, nullptr, this->ligand.Model.atomCharges, nullptr);

    mdc->SetBFactorRange(0.0f, 0.0f);
    mdc->SetChargeRange(this->minCharge, this->maxCharge);
    mdc->SetOccupancyRange(0.0f, 0.0f);
    mdc->SetFormerAtomIndices(nullptr);

    mdc->SetConnections(static_cast<unsigned int>(this->connectivity[this->ligandIdx].Count() / 2),
        this->connectivity[this->ligandIdx].PeekElements());
    this->residue[0]->SetPosition(0, this->ligand.ligandAtomCnt);
    this->residue[0]->SetBoundingBox(this->bbox);
    mdc->SetResidues(1, (const MolecularDataCall::Residue**)this->residue.PeekElements());
    mdc->SetSolventResidueIndices(0, 0);
    mdc->SetResidueTypeNames(1, (vislib::StringA*)this->residueTypeName.PeekElements());
    mdc->SetMolecules(1, (MolecularDataCall::Molecule*)this->molecule.PeekElements());
    mdc->SetChains(1, (MolecularDataCall::Chain*)this->chain.PeekElements());

    // Set the filter array for the molecular data call
    mdc->SetFilter(0);
    this->PDBQTLoadingStatus = true;
    return true;
}

/*
 * MultiPDBQTLoader::getExtent
 */
bool MultiPDBQTLoader::getExtent(core::Call& call) {
    MolecularDataCall* dc = dynamic_cast<MolecularDataCall*>(&call);
    if (dc == NULL) return false;

    if (this->pdbqtListFilenameSlot.IsDirty() || this->total_ligandCnt < 1) {
        this->pdbqtListFilenameSlot.ResetDirty();
        if (this->loadPDBQTFiles(this->pdbqtListFilenameSlot.Param<core::param::FilePathParam>()->Value())) {
            dc->SetPDBFilename(this->pdbqtListFilenameSlot.Param<core::param::FilePathParam>()->Value());
        } else {
            return false;
        }
    }

    // grow bounding box by 3.0 Angstrom (for volume rendering / SAS)
    vislib::math::Cuboid<float> bBoxPlus3 = this->bbox;
    bBoxPlus3.Grow(3.0f);

    dc->AccessBoundingBoxes().Clear();
    dc->AccessBoundingBoxes().SetObjectSpaceBBox(bBoxPlus3);
    dc->AccessBoundingBoxes().SetObjectSpaceClipBox(bBoxPlus3);

    dc->SetFrameCount(1);
    dc->SetDataHash(this->mdc_dataHash);

    return true;
}

bool MultiPDBQTLoader::setRequestedLigandAndModelSilent(core::Call& call) {
    LigandModelCall* lmc = dynamic_cast<LigandModelCall*>(&call);
    if (lmc == NULL) return false;
    this->isSilentDataCall = TRUE;
    setRequestedLigandAndModel(*lmc);
    return true;
}

bool MultiPDBQTLoader::setRequestedLigandAndModel(core::Call& call) {
    checkWaitPDBQT_data();
    LigandModelCall* lmc = dynamic_cast<LigandModelCall*>(&call);
    if (lmc == nullptr) return false;

    this->getLigandData(lmc->getCurrentLigandID(), lmc->getCurrentModelID());
    lmc->SetBindigenergy(ligand.Model.bindingEnergy);
    lmc->SetLigandname(this->ligand.ligandName.c_str());
    lmc->SetCurrentGlobalModelID(this->ligand.Model.globalModelID);
    lmc->SetSVGPath(this->ligand.SVGPath);
    lmc->SetFunctionalGroupsOfLigand(this->ligand.functionalGroupsLigand);
    lmc->SetChemPropsOfLigand(this->ligand.chemPropsLigand);
    lmc->SetWebData(&this->webData);

    lmc->SetCurrentModelEfficiency(this->ligand.Model.modelEfficiency);

    if (this->plipForcesStore.size() > this->ligand.Model.globalModelID) {
        lmc->SetInteractionForces_hbonds_protAcceptor(this->ligand.Model.Forces.hProtAcceptors);
        lmc->SetInteractionForces_hbonds_protDonor(this->ligand.Model.Forces.hProtDonors);
        lmc->SetInteractionForces_hydrophobicInteractions(this->ligand.Model.Forces.hydrophobicInteractions);
        lmc->SetInteractionForces_saltBridges(this->ligand.Model.Forces.saltBridges);
        lmc->SetInteractionForces_piStacks(this->ligand.Model.Forces.piStacks);
        lmc->SetInteractionForces_piCationInteractions(this->ligand.Model.Forces.piCationInteractions);
        lmc->SetInteractionForces_halogenBonds(this->ligand.Model.Forces.halogenBonds);
        lmc->SetInteractionForces_metalComplexes(this->ligand.Model.Forces.metalComplexes);
    }


    if (!this->isSilentDataCall) {
        this->ligandIdx = lmc->getCurrentLigandID();
        this->modelIdx = lmc->getCurrentModelID();

    } else {
        this->isSilentDataCall = FALSE;
    }
    lmc->SetDataHash(lmc_dataHash);

    return true;
}


bool MultiPDBQTLoader::getLigandAndModelCount(core::Call& call) {
    checkWaitPDBQT_data();


    LigandModelCall* lmc = dynamic_cast<LigandModelCall*>(&call);
    if (lmc == NULL) return false;
    lmc->SetClusteredFunctionalGroups(&this->clusteredFunctionalGroups);
    lmc->SetlmcStatus(
        this->MPL_modelClustRenderer_busy, this->MPL_SocketCommManager_busy, this->MPL_firstRenderCall_complete);

    lmc->SetLigandAndModelCount(this->total_ligandCnt, this->modelsPerLigand.PeekElements(), this->total_modelCnt);
    lmc->SetTotalAtomPosCnt(atomPositions.Count() / 3);
    lmc->SetDataHash(lmc_dataHash);

    lmc->SetPDBQT_list_Filepath(this->PDBQT_list_Filepath);

    return true;
}

bool MultiPDBQTLoader::getClusterData(core::Call& call) {

    checkWaitPDBQT_data();

    LigandModelCall* lmc = dynamic_cast<LigandModelCall*>(&call);
    if (lmc == NULL) return false;

    getLigandAndModelCount(call);
    lmc->SetGlobalMdlID_to_ligID_mdlID(&this->globalMdlID_to_ligID_mdlID);
    // getRequestedLigandAndModel(call);


    if (this->DBSCANParamsChanged == true) {
        this->DBSCANParamsChanged = false;




        unsigned int ligCnt = this->total_ligandCnt;
        const unsigned int* mdlCnt = this->modelsPerLigand.PeekElements();

        if (this->modelCentroidsVec4.Count() != this->total_modelCnt) {

            this->modelCentroids.Resize(0);
            this->modelCentroidsCuboids.Resize(0);
            this->modelCentroidsVec4.SetCount(this->total_modelCnt);
            this->modelCentroids.SetCount(this->total_modelCnt * 4);
            this->modelCentroidsCuboids.SetCount(this->total_modelCnt);
            vislib::math::Vector<float, 3> currentCentroid;

            unsigned int cnt = 0;

            for (int ligIdx = 0; ligIdx < ligCnt; ligIdx++) {
                for (int mdlIdx = 0; mdlIdx < mdlCnt[ligIdx]; mdlIdx++) {
                    // get the LigandModelData
                    getLigandData(ligIdx, mdlIdx);

                    for (uint z = 0; z < this->ligand.ligandAtomCnt; z++) {
                        if (z != 0) {
                            modelCentroidsCuboids[cnt].GrowToPoint(ligand.Model.atomPostions[z * 3 + 0],
                                ligand.Model.atomPostions[z * 3 + 1], ligand.Model.atomPostions[z * 3 + 2]);
                        } else {
                            modelCentroidsCuboids[cnt].Set(ligand.Model.atomPostions[z * 3 + 0],
                                ligand.Model.atomPostions[z * 3 + 1], ligand.Model.atomPostions[z * 3 + 2],
                                ligand.Model.atomPostions[z * 3 + 0], ligand.Model.atomPostions[z * 3 + 1],
                                ligand.Model.atomPostions[z * 3 + 2]);
                        }
                    }
                    // currentCentroid /= mdc->AtomCount();
                    // MolecularDataCall* pmdc = this->dataInSlot_MolecularDataCall.CallAs<MolecularDataCall>();
                    // auto a = pmdc->GetBoundingBoxes().ObjectSpaceBBox().CalcCenter();
                    // modelCentroidsCuboids[cnt].Move(a.GetX(), a.GetY(), a.GetZ());
                    currentCentroid = modelCentroidsCuboids[cnt].CalcCenter();
                    this->modelCentroids[4 * cnt + 0] = currentCentroid.GetX();
                    this->modelCentroids[4 * cnt + 1] = currentCentroid.GetY();
                    this->modelCentroids[4 * cnt + 2] = currentCentroid.GetZ();
                    this->modelCentroids[4 * cnt + 3] = 1.0f;
                    this->modelCentroidsVec4[cnt].SetX(currentCentroid.GetX());
                    this->modelCentroidsVec4[cnt].SetY(currentCentroid.GetY());
                    this->modelCentroidsVec4[cnt].SetZ(currentCentroid.GetZ());
                    // get binding energy

                    this->modelCentroidsVec4[cnt].SetW(this->ligand.Model.bindingEnergy);
                    cnt++;
                }
            }
        }


        numberOfCluster = 0;
        assignedClusters.SetCount(this->total_modelCnt);
        /**
         * *** Cluster file handling ***
         * @info:	gets PDBQT_list_Filepath + adds "\clustering\" to it to save cluster files under this directory
         * @info:	checks if this file und this path exists
         * @info:	if 'NO' DBSCAN will be executed
         * @info:	if 'YES' clusters will be loaded from file
         */
        Path file;
        file = getDirpath_from_filepath(this->PDBQT_list_Filepath.PeekBuffer());
        file.append(this->additionalDataFolderName);
        file.append("clusterings");
        // create dir
        std::filesystem::create_directory(file);
        file.append(getFilname_from_filepath(this->PDBQT_list_Filepath.PeekBuffer()));

        std::string str_eps = fmt::format("{:.2f}", this->eps).c_str();
        std::string str_minPts = std::to_string(this->minPts).c_str();
        std::string str_minBindEnergy = std::to_string(this->minBindEnergy)
                                            .erase(4, std::to_string(this->minBindEnergy).find_first_of("0") + 1)
                                            .c_str();
        std::string str_underspace = "_";

        file.concat(str_underspace);
        file.concat(str_eps);
        file.concat(str_underspace);
        file.concat(str_minPts);
        file.concat(str_underspace);
        file.concat(str_minBindEnergy);

        ifstream checkFile;

        // only recompute clustering if cached file was not available
        if (!fileExists(file.string())) {
            printf("Cluster file does not exist! calculating...");

            //////////////////////////////////////////////////////////////////

            /** clustering of centroids */
            DBSCAN dbscan;
            assignedClusters =
                dbscan.DBSCAN_3D(this->modelCentroidsVec4, this->eps, this->minPts, this->minBindEnergy, 1);

            ofstream myfile(file.string(), std::ofstream::app);
            if (myfile.is_open()) {
                for (int i = 0; i < assignedClusters.Count(); i++) {
                    myfile << assignedClusters[i];
                    myfile << "\n";
                }
                myfile.close();
            } else
                cout << "Unable to open file";
        } else {
            checkFile.open(file.string());
            int counter = 0;
            std::string line;
            while (std::getline(checkFile, line)) {
                assignedClusters[counter] = std::stoi(line);
                counter++;
            }
            checkFile.close();
        }
        // end of *** Cluster file handling ***

        // TODO: make DBSCAN_3D more memory efficient (directly return a (sorted?) vec2 from DBSCAN)
        /**
         * saves for every clusterID the modelID and sort vec2 by X() clusterID
         * vec2[i] where i is sphereID by the mouse picking
         */

        assign_modelID_sphereID.SetCount(assignedClusters.Count());
#pragma omp parallel for
        for (int i = 0; i < assignedClusters.Count(); i++) {
            assign_modelID_sphereID[i].SetX(assignedClusters[i]);
            assign_modelID_sphereID[i].SetY(i);
        }
        // sort by cluster ID
        assign_modelID_sphereID.Sort(sortVec2ByX);

        int noiseCnt = 0;
        // get number of clusters
        for (int i = 0; i < assignedClusters.Count(); i++) {
            // printf("i: %d \t cluster:%d\n", i, assign_modelID_sphereID[i].GetX());
            if (numberOfCluster < assignedClusters[i]) {
                numberOfCluster = assignedClusters[i];
            }
            if (assignedClusters[i] == -1) {
                noiseCnt++;
            }
        }
        // if -1 there are no clusters (so we need to set to 0)
        // if 0 we have one cluster with ID 0 (so we need to add +1 for allocation things)
        if (numberOfCluster == -1) {
            numberOfCluster = 0;
        } else {
            numberOfCluster++;
        }

        printf("eps:%f \t minPts: %d\n", this->eps, this->minPts);
        (lmc)->SetDBSCANParams(this->eps, this->minPts, this->minBindEnergy);

        /** centroids of clusters*/
        printf("clusterCnt: %d\n", numberOfCluster);
        this->clusterCentroids.Resize(0);
        this->clusterCentroids.SetCount(numberOfCluster * 4);

        this->clusterSizes.SetCount(numberOfCluster);

        this->clusterCentroidsCuboids.Resize(0);
        this->clusterCentroidsCuboids.SetCount(numberOfCluster);

        /** initialize arrays with zero (0) */
        for (int j = 0; j < numberOfCluster; j++) {
            this->clusterSizes[j] = 0;
        }

        /** get number of model centroids for each cluster and sum up x,y,z coordinates */
        vislib::math::Point<float, 3> origin;
        for (int i = 0; i < assignedClusters.Count(); i++) {
            if (assignedClusters[i] >= 0) {
                this->clusterSizes[assignedClusters[i]]++;
                origin = this->clusterCentroidsCuboids[assignedClusters[i]].CalcCenter();
                if (this->clusterCentroidsCuboids[assignedClusters[i]].IsEmpty()) {
                    this->clusterCentroidsCuboids[assignedClusters[i]].Set(this->modelCentroids[4 * i + 0],
                        this->modelCentroids[4 * i + 1], this->modelCentroids[4 * i + 2],
                        this->modelCentroids[4 * i + 0], this->modelCentroids[4 * i + 1],
                        this->modelCentroids[4 * i + 2] + 0.01f);
                } else {
                    this->clusterCentroidsCuboids[assignedClusters[i]].Union(this->modelCentroidsCuboids[i]);
                }
            }
        }

        /** determine cluster centroid*/
        vislib::math::Point<float, 3> center;
        for (int j = 0; j < numberOfCluster; j++) {
            center = clusterCentroidsCuboids[j].CalcCenter();

            this->clusterCentroids[4 * j + 3] = 2;
            this->clusterCentroids[4 * j + 0] = center.GetX();
            this->clusterCentroids[4 * j + 1] = center.GetY();
            this->clusterCentroids[4 * j + 2] = center.GetZ();
        }

        clusterIndexSum.Resize(0);
        clusterIndexSum.SetCount(numberOfCluster + 1);
        clusterIndexSum[0] = noiseCnt;
        minClusterSize = clusterSizes[0];
        maxClusterSize = clusterSizes[0];
        for (int i = 0; i < numberOfCluster; i++) {
            clusterIndexSum[i + 1] = clusterIndexSum[i] + clusterSizes[i];
            // clusterSizes[i] = clusterIndexSum[i+1] - clusterIndexSum[i];
            printf("Cluster:%d Range: %d - %d \t size: %d\n", i, clusterIndexSum[i], clusterIndexSum[i + 1],
                clusterSizes[i]);
            minClusterSize = minClusterSize > clusterSizes[i] ? clusterSizes[i] : minClusterSize;
            maxClusterSize = maxClusterSize < clusterSizes[i] ? clusterSizes[i] : maxClusterSize;
        }
        this->DBSCANParamsChanged = TRUE;
        this->clusterDataHash++;
    }

    // if first call (clusterDataHash == 0) or DBSCANParamsChanged set Cluster and Centroid data
    if ((lmc)->getClusterAndCentroidData().clusterDataHash != this->clusterDataHash || this->DBSCANParamsChanged) {

        (lmc)->SetClusterAndCentroidData(&modelCentroids, &clusterCentroids, &assignedClusters, &clusterIndexSum,
            &modelCentroidsCuboids, &clusterCentroidsCuboids, &assign_modelID_sphereID, &clusterSizes,
            this->minClusterSize, this->maxClusterSize, this->clusterDataHash);


        if (this->DBSCANParamsChanged) {
            this->DBSCANParamsChanged = false;
            Sleep(1000);
            // force calculations
            MolecularDataCall* pmdc = this->dataInSlot_MolecularDataCall.CallAs<MolecularDataCall>();

            InteractionForces intForces(lmc, pmdc, ligand);

            if (lmc->getInteractionForces().useMegaMolHBonds) {
                intForces.calcHBonds(this->hBonds_ligMdlDonor, this->hBonds_ligMdlAcceptor, this->hBonds_ligMdlHDonor,
                    this->hBonds_ligMdlHAcceptor, this->hBondAtomCnt_protein);
            }

            this->plipForcesStore.clear();
            this->plipForcesStore.resize(0);
            intForces.loadPLIPdata(getAdditionalDataFolderPath(), this->plipForcesStore);

            /***************************************
            ************ push HBonds ***************
            ****************************************/
            // TODO: make this more effecient
            // directly save positions of atoms
            // currently only IDs came from calcHBonds()
            // these are then converted into atom positions (in the following code)
            // Reason for conversion: achieve same data and structure like for the PLIP forces
            if (pmdc == NULL) {
                return false;
            }
            if (lmc->getInteractionForces().useMegaMolHBonds) {
                (*pmdc)(megamol::protein_calls::MolecularDataCall::CallForGetExtent);
                (*pmdc)(megamol::protein_calls::MolecularDataCall::CallForGetData);
                (*lmc)(LigandModelCall::CallForGetExtent);
                for (int gblMdlID = 0; gblMdlID < this->total_modelCnt; gblMdlID++) {
                    lmc->SetCurrentLigandAndModel_byGlobalMdlID(gblMdlID);
                    (*lmc)(LigandModelCall::CallForGetData);
                    MultiPDBQTLoader::PLIPForces* fc = &plipForcesStore[gblMdlID];

                    // protein Acceptor
                    auto hAcc = &fc->hProtAcceptors;
                    // protein Donor
                    auto hDon = &fc->hProtDonors;

                    if (this->hBonds_ligMdlDonor.Count() > 0) {
                        for (int i = 0; i < this->hBonds_ligMdlDonor[gblMdlID].Count() / 2; i++) {
                            hAcc->cnt++;
                            int ligID = this->hBonds_ligMdlHDonor[gblMdlID][i * 2 + 1];
                            int protID = this->hBonds_ligMdlHDonor[gblMdlID][i * 2 + 0];
                            glm::vec3 ligPos = glm::vec3(this->ligand.Model.atomPostions[ligID * 3 + 0],
                                this->ligand.Model.atomPostions[ligID * 3 + 1],
                                this->ligand.Model.atomPostions[ligID * 3 + 2]);
                            glm::vec3 protPos = glm::vec3(pmdc->AtomPositions()[protID * 3 + 0],
                                pmdc->AtomPositions()[protID * 3 + 1], pmdc->AtomPositions()[protID * 3 + 2]);
                            hAcc->ligIDx.push_back(ligID);
                            hAcc->protIDx.push_back(protID);
                            hAcc->ligPos.push_back(ligPos);
                            hAcc->protPos.push_back(protPos);
                        }
                    }
                    if (this->hBonds_ligMdlAcceptor.Count() > 0) {

                        for (int i = 0; i < this->hBonds_ligMdlAcceptor[gblMdlID].Count() / 2; i++) {
                            hDon->cnt++;
                            int ligID = this->hBonds_ligMdlHAcceptor[gblMdlID][i * 2 + 1];
                            int protID = this->hBonds_ligMdlHAcceptor[gblMdlID][i * 2 + 0];
                            glm::vec3 ligPos = glm::vec3(this->ligand.Model.atomPostions[ligID * 3 + 0],
                                this->ligand.Model.atomPostions[ligID * 3 + 1],
                                this->ligand.Model.atomPostions[ligID * 3 + 2]);
                            glm::vec3 protPos = glm::vec3(pmdc->AtomPositions()[protID * 3 + 0],
                                pmdc->AtomPositions()[protID * 3 + 1], pmdc->AtomPositions()[protID * 3 + 2]);

                            hDon->ligIDx.push_back(ligID);
                            hDon->protIDx.push_back(protID);
                            hDon->ligPos.push_back(ligPos);
                            hDon->protPos.push_back(protPos);
                        }
                    }
                }
            }
        } // DBSCANParamsChanged


        // lmc->ResetDBSCANChanged();
    }
    return true;
}


/****************************
**** PDBQT LOADER CLASS ****
****************************/

bool MultiPDBQTLoader::loadPDBQTFiles(const vislib::TString& filename) {
    // printf("\n\nFilelocation: %s\n\n", filename.PeekBuffer());
    this->PDBQT_list_Filepath = filename;
    using megamol::core::utility::log::Log;


    // temporary variables
    vislib::StringA listLine = "";
    vislib::StringA PDBQTLine = "";
    unsigned int atomCnt, curList_lineCnt, curPDBQT_lineCnt;

    this->minCharge = FLT_MAX;
    this->maxCharge = -(FLT_MAX - 1.0f);

    bool test = resetMultiPDBQTLoader();
    allocateMultiPDBQTLoader(0);

    vislib::sys::ASCIIFileBuffer file_PDBQT_list; // for results_list.txt
    vislib::sys::ASCIIFileBuffer file_PDBQT;      // for the individual files described in results_list.txt


    // get the root path of the pdbqt files by extracting it from the path of PDBQT list
    string filepath = filename.PeekBuffer();


    if (!filepath.size() > 0) {
        printf("missing PDBQT List Filepath\n");
        return false;
    } else {
        printf("PDBQT List Filepath=%s\n", filepath.c_str());
    }


    Path sysPath;
    sysPath = getDirpath_from_filepath(filepath);

    Log::DefaultLog.WriteMsg(Log::LEVEL_INFO, "Loading list of PDBQT files: %s", T2A(filename.PeekBuffer())); // DEBUG
    // try to load the file
    bool file_loaded = false;
    unsigned int modelCnt = 0;


    /**********************************************************************
    *********************load list of pdbqt files**************************
    **********************************************************************/

    if (file_PDBQT_list.LoadFile(T2A(filename.PeekBuffer()))) {
        // file successfully loaded, read first frame
        file_loaded = true;
        curList_lineCnt = 0;
        curPDBQT_lineCnt = 0;
        Path PathPDBQT;

        this->ligandAtomCnts.SetCount(file_PDBQT_list.Count());
        int modelAtomCnt = 0;
        int oldProgress = 0;

        while (curList_lineCnt < file_PDBQT_list.Count() && !listLine.StartsWith("END")) {

            /**************************************
            *******PDBQT file path preparation*****
            ***************************************/

            // get the current line from the file list
            listLine = file_PDBQT_list.Line(curList_lineCnt);

            if (listLine.StartsWith("./")) {
                // ignore alternate locations
                auto res_string = listLine;
                res_string.TrimSpaces();

                std::string tmp = string(listLine.PeekBuffer()).substr(2, -1);
                Path subPath = Path(tmp).make_preferred();

                PathPDBQT = sysPath / subPath.make_preferred();
                vislib::StringA PathToSinglePDBQT = PathPDBQT.string().c_str();


                /**********************************************************************
                ***********************load one pdbqt file****************************
                **********************************************************************/

                if (file_PDBQT.LoadFile(T2A(PathToSinglePDBQT))) {
                    // store model count of previous ligand
                    if (this->total_ligandCnt > 0) {
                        this->modelsPerLigand.Add(modelCnt);
                    }
                    this->total_ligandCnt++;
                    modelCnt = 0;
                    curPDBQT_lineCnt = 0;

                    // file successfully loaded, read first frame
                    int actualProgress = (int)((float)curList_lineCnt / (float)file_PDBQT_list.Count() * 1000);

                    if (actualProgress != oldProgress) {
                        if (verbose)
                            printf("\r loading PDBQTs %.1f %% done %d/%d: ", (float)actualProgress / 10,
                                curList_lineCnt,
                            file_PDBQT_list.Count());
                        //printf("\r %.1f %% done \t Path PDBQT %d/%d: %s ", (float)actualProgress / 10, curList_lineCnt,
                        //    file_PDBQT_list.Count(), PathPDBQT.string().c_str());
                        fflush(stdout);
                        oldProgress = actualProgress;
                    }
                    // printf("%d \t %d \t %f\n", actualProgress, oldProgress, (float)curList_lineCnt /
                    // (float)file_PDBQT_list.Count());

                    vislib::StringA ligandName = "";
                    while (curPDBQT_lineCnt < file_PDBQT.Count()) {
                        //&& !PDBQTLine.StartsWith("ENDML")) ||  checkNextModel.StartsWith("MODEL")) {

                        // get the current line from the PDBQT file
                        PDBQTLine = file_PDBQT.Line(curPDBQT_lineCnt);
                        if (PDBQTLine.StartsWith("ATOM")) {
                            // this->total_atomPosCnt++;
                            modelAtomCnt++;

                            // the atomTypes are the same for each frame of one ligand
                            // only stores types if the first from a model/ligand

                            if (modelCnt == 1) {
                                // add atomType - atomName
                                auto atomTypeStr = PDBQTLine.Substring(12, 4);
                                atomTypeStr.TrimSpaces();
                                // get the radius of the element
                                float radius = PDBLoader::getElementRadius(atomTypeStr);
                                // get the color of the element
                                vislib::math::Vector<unsigned char, 3> color = PDBLoader::getElementColor(atomTypeStr);
                                // set the new atom type
                                MolecularDataCall::AtomType type(
                                    atomTypeStr, radius, color.X(), color.Y(), color.Z(), atomTypeStr);
                                // search for current atom type in atom type array
                                INT_PTR tmpAtomTypeIdx = this->atomTypes.IndexOf(type);
                                if (tmpAtomTypeIdx == vislib::Array<MolecularDataCall::AtomType>::INVALID_POS) {
                                    this->atomTypeIdx.Add(static_cast<unsigned int>(this->atomTypes.Count()));
                                    this->atomTypes.Add(type);
                                } else {
                                    this->atomTypeIdx.Add(static_cast<unsigned int>(tmpAtomTypeIdx));
                                }

                                // add atomADType
                                auto atomADType = PDBQTLine.Substring(77, 3);
                                atomADType.TrimSpaces();
                                this->atomTypesAD.Add(atomADType);
                            }

                            // add atom postitions
                            auto atom_X = PDBQTLine.Substring(31, 8);
                            auto atom_Y = PDBQTLine.Substring(39, 8);
                            auto atom_Z = PDBQTLine.Substring(47, 8);
                            auto atom_Q = PDBQTLine.Substring(70, 7);

                            atom_X.TrimSpaces();
                            atom_Y.TrimSpaces();
                            atom_Z.TrimSpaces();
                            atom_Q.TrimSpaces();
                            float crg = atof(atom_Q);
                            float x = atof(atom_X);
                            float y = atof(atom_Y);
                            float z = atof(atom_Z);

                            this->atomPositions.Add(x);
                            this->atomPositions.Add(y);
                            this->atomPositions.Add(z);
                            this->atomCharges.Add(crg);

                            // is this the first atom?
                            if (this->minCharge > this->maxCharge) {
                                // if yes, initialize bounding box
                                this->bbox.Set(x, y, z, x, y, z);
                            } else {
                                // if no, update bounding box
                                this->bbox.GrowToPoint(x, y, z);
                            }
                            // store min and max charges
                            this->minCharge = std::fminf(this->minCharge, crg);
                            this->maxCharge = std::fmaxf(this->maxCharge, crg);
                        } else if (PDBQTLine.StartsWith("REMARK")) {
                            auto tmpStr = PDBQTLine.Substring(6);
                            tmpStr.TrimSpaces();
                            if (tmpStr.StartsWith("VINA RESULT:")) {
                                /* gets
                                 * + binding energy
                                 * + RMSD lb/ub (lb = lower bound // ub = upper bound)
                                 * + Atom X, Y, Z, Charge Q
                                 * ---> of the model*/
                                auto bindingEnergy = PDBQTLine.Substring(24, 6);
                                auto RMSD_lb = PDBQTLine.Substring(34, 6);
                                auto RMSD_ub = PDBQTLine.Substring(45, 6);

                                bindingEnergy.TrimSpaces();
                                RMSD_lb.TrimSpaces();
                                RMSD_ub.TrimSpaces();

                                float tmp_bindingEnergy = atof(bindingEnergy);

                                // if energy value is zero or positive something went wrong
                                if (tmp_bindingEnergy < 0) {
                                    this->bindingEnergy.Add(float(atof(bindingEnergy)));

                                } else {
                                    // use energy value from last model as compensation
                                    this->bindingEnergy.Add(this->bindingEnergy[this->bindingEnergy.Count() - 1]);
                                }


                                this->RMSD_lb.Add(float(atof(RMSD_lb)));
                                this->RMSD_ub.Add(float(atof(RMSD_ub)));
                            } else if (tmpStr.StartsWith("Name =") && modelCnt == 1) {
                                // get name of the ligand
                                ligandName = PDBQTLine.Substring(15, 50);
                                ligandName.TrimSpaces();
                                this->ligandNames.Add(ligandName);
                            }
                        } else if (PDBQTLine.StartsWith("MODEL")) {
                            this->total_modelCnt++;
                            modelCnt++;
                            modelAtomCnt = 0;
                            if (this->total_modelCnt >= this->frameCapacity) {
                                // set new Capacity for all arrays
                                bool needToExtend = true;
                                this->allocateMultiPDBQTLoader(needToExtend);
                            }
                        } else if (PDBQTLine.StartsWith("ENDMDL")) {
                            this->ligandAtomCnts[this->total_ligandCnt - 1] = modelAtomCnt;
                            this->total_atomCnt += modelAtomCnt;

                            // fallback if no ZINC name exists
                            if (ligandName == "") {
                                std::string tmpLigName = "lig" + std::to_string(this->total_ligandCnt);
                                this->ligandNames.Add(tmpLigName.c_str());
                                ligandName = tmpLigName.c_str();
                            }
                        }
                        curPDBQT_lineCnt++;
                    }
                } else {
                    Log::DefaultLog.WriteMsg(Log::LEVEL_INFO, "ERROR: Unable to load file: %s", PathToSinglePDBQT);
                }
                this->ligandsModelRange[this->total_ligandCnt - 1] =
                    (this->total_ligandCnt <= 1 ? modelCnt
                                                : modelCnt + this->ligandsModelRange[this->total_ligandCnt - 2]);
            }
            curList_lineCnt++;
        }
        // store model count of previous ligand
        if (this->total_ligandCnt > 0) {
            this->modelsPerLigand.Add(modelCnt);
        }
        Log::DefaultLog.WriteMsg(Log::LEVEL_INFO, "Atom count: %i", this->atomCharges.Count()); // DEBUG
    }
    if (!file_loaded) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Could not load file %s", (const char*)T2A(filename)); // DEBUG
        return false;
    } else {

        //******************************************************************//
        //*********** IF PDBQT LOADING WAS SUCCESSFUL AND PICKING **********//
        //******************************************************************//
        lmc_dataHash++;
        mdc_dataHash++;
        createHelperDataStructure();

        // prepare all
        getfunctionalGroupDataFiles(this->PDBQT_list_Filepath);
        getProteinLigand_Interactions(this->PDBQT_list_Filepath);
        getSVGsAndChemProps(this->PDBQT_list_Filepath);
        loadCheckmolFiles();
    }

    // compute connections (estimate covalent bonds)
    vislib::math::Vector<float, 3> atomPos0, atomPos1;
    this->connectivity.SetCount(this->total_ligandCnt);
    for (unsigned int i = 0; i < this->total_ligandCnt; i++) {
        this->connectivity[i].Clear();
        this->getLigandData(i);
        for (unsigned int j = 0; j < ligand.ligandAtomCnt; j++) {
            // get first atom positions
            atomPos0.Set(this->ligand.Model.atomPostions[3 * j + 0], this->ligand.Model.atomPostions[3 * j + 1],
                this->ligand.Model.atomPostions[3 * j + 2]);
            for (unsigned int k = j + 1; k < ligand.ligandAtomCnt; k++) {
                // get second atom positions
                atomPos1.Set(this->ligand.Model.atomPostions[3 * k + 0], this->ligand.Model.atomPostions[3 * k + 1],
                    this->ligand.Model.atomPostions[3 * k + 2]);
                // check distance
                if ((atomPos0 - atomPos1).Length() <
                    0.58f * (this->atomTypes[this->ligand.Model.atomTypeIdx[j]].Radius() +
                                this->atomTypes[this->ligand.Model.atomTypeIdx[k]].Radius())) {
                    // add connection
                    this->connectivity[i].Add(j);
                    this->connectivity[i].Add(k);
                }
            }
        }
    }
    calculateModelEfficiency();
    this->PDBQTLoadingStatus = true;
    return true;
}

void MultiPDBQTLoader::calculateModelEfficiency() {

    unsigned int ligCnt = this->total_ligandCnt;
    const unsigned int* mdlCnt = this->modelsPerLigand.PeekElements();
    unsigned int totalMdlCnt = this->total_modelCnt;
    uint globalMdlID = 0;

    this->modelEfficiency.SetCount(totalMdlCnt);
    for (int ligIdx = 0; ligIdx < ligCnt; ligIdx++) {
        for (int mdlIdx = 0; mdlIdx < mdlCnt[ligIdx]; mdlIdx++) {
            // get the LigandModelData
            getLigandData(ligIdx, mdlIdx);
            globalMdlID = getGlobalModelID(ligIdx, mdlIdx);
            vislib::Array<vislib::StringA> heavyAtoms;

            for (uint z = 0; z < this->ligand.ligandAtomCnt; z++) {
                auto atom = this->atomTypes[this->ligand.Model.atomTypeIdx[z]].Name();

                // push all heavy atoms into array
                if (atom != vislib::StringA("H")) {
                    heavyAtoms.Append(atom);
                }
            }

            float nrg = this->getBindingEnergy(globalMdlID);
            float eff = nrg / heavyAtoms.Count(); // efficiency of one model

            // efficiencies of all models and ligands into array
            this->modelEfficiency[globalMdlID] = eff;
        }
    }
}


bool MultiPDBQTLoader::loadCheckmolFiles() {

    std::string lingandName;
    for (int i = 0; i < this->ligandNames.Count(); i++) {
        lingandName = this->ligandNames[i].PeekBuffer();


        Path dataFolderPath = getAdditionalDataFolderPath();
        Path checkmolFilePath = dataFolderPath.append(lingandName);
        checkmolFilePath.concat(".checkmol");

        LigandModelCall::FunctionalGroupsLigand functionalGroups;
        LigandModelCall::ChemPropsLigand chemProps;

        string line;
        ifstream file(checkmolFilePath);
        if (checkmolFilePath.filename() == "ZINC000033765216.checkmol") {
            int a = 1;
		}
        if (file.is_open()) {
            try {
            
				while (getline(file, line)) {
					LigandModelCall::fGroup fgroup;


					std::vector<std::string> tmpSplit;
					if (!line._Starts_with("#Properties")) {
						split(line, ":", tmpSplit);
						tmpSplit[0].replace(0, 1, "");
						int groupID = std::stoi(tmpSplit[0]);
						fgroup.groupID = groupID;

						int speficGroupCnt = std::stoi(tmpSplit[1]);
						fgroup.speficGroupCnt = speficGroupCnt;

						if (speficGroupCnt > 1) {
							std::vector<std::string> groupAtoms;
							split(tmpSplit[2], ",", groupAtoms);
							for (int i = 0; i < speficGroupCnt; i++) {
								fgroup.centralAtoms.push_back(std::stoi(groupAtoms[i]));
							}
						} else {
							fgroup.centralAtoms.push_back(std::stoi(tmpSplit[2]));
						}
						functionalGroups.totalGroupCnt++;
						functionalGroups.fGroups.push_back(fgroup);
					} else {
						split(line, ";", tmpSplit);
						std::vector<std::string> oneProp;
						for (int i = 1; i < tmpSplit.size(); i++) {
							split(tmpSplit[i], ":", oneProp);
							if (oneProp[0] == "mwt") {
								chemProps.mwt = std::stof(oneProp[1]);
							} else if (oneProp[0] == "logp") {
								chemProps.logp = std::stof(oneProp[1]);
							} else if (oneProp[0] == "fractioncsp3") {
								chemProps.fractioncsp3 = std::stof(oneProp[1]);
							}
							oneProp.clear();
						}
					}
				}
				file.close();
            }
			catch (const std::exception& e) {
                printf("ERROR while reading line of checkmol file : %s %s\n", e, checkmolFilePath.string().c_str());
			}
		} else {
			printf("Unable to open file: %s\n", checkmolFilePath.string().c_str());
		}
        
        this->functionalGroupsStorage.push_back(functionalGroups);
        this->chemPropsStorage.push_back(chemProps);
    }

    return true;
};

Path MultiPDBQTLoader::getAdditionalDataFolderPath() {

    if (this->additionalDataFolderPath.size() <= 0) {
        std::string PDBQTListFilename =
            this->pdbqtListFilenameSlot.Param<core::param::FilePathParam>()->Value().PeekBuffer();
        Path dataFolderPath = getDirpath_from_filepath(PDBQTListFilename);
        dataFolderPath.append(this->additionalDataFolderName);
        dataFolderPath = dataFolderPath.make_preferred();
        this->additionalDataFolderPath = dataFolderPath.string();
    }
    return Path(this->additionalDataFolderPath);
};

std::string MultiPDBQTLoader::getFilname_from_filepath(std::string filepath) {

    Path filename = Path(filepath).filename();

    return filename.string();
};


std::filesystem::path MultiPDBQTLoader::getDirpath_from_filepath(std::string filepath) {

    std::filesystem::path baseDir;
    baseDir = Path(filepath).make_preferred();

    return baseDir.parent_path();
};

Path MultiPDBQTLoader::getDirpath_toMegaMol_from_filepath(std::string cwd) {

    Path baseDir;
    Path tmp = Path(cwd);

    int i = 0;
    for (const auto& part : tmp) {
        if (part.string() != "build") {
            baseDir.append(part.string());
        } else {
            break;
        }
    }

    return baseDir.make_preferred();
};

void MultiPDBQTLoader::createHelperDataStructure() {
    this->globalMdlID_to_ligID_mdlID.SetCount(this->total_modelCnt);
    this->ligID_mdlID_to_globalMdlID.SetCount(this->total_ligandCnt);

    int cnt = 0;
    for (int i = 0; i < total_ligandCnt; i++) {
        this->ligID_mdlID_to_globalMdlID[i].SetCount(this->modelsPerLigand[i]);
        for (int j = 0; j < this->modelsPerLigand[i]; j++) {
            this->globalMdlID_to_ligID_mdlID[cnt].SetX(i);
            this->globalMdlID_to_ligID_mdlID[cnt].SetY(j);
            this->ligID_mdlID_to_globalMdlID[i][j] = cnt;
            cnt++;
        }
    }
    printf("Count globalMdlID: %d\n", this->globalMdlID_to_ligID_mdlID.Count());
}

const vislib::StringA& MultiPDBQTLoader::getFilename() {
    return this->pdbqtListFilenameSlot.Param<core::param::FilePathParam>()->Value();
}

void MultiPDBQTLoader::checkWaitPDBQT_data() {
    while (this->PDBQTLoadingStatus == false) {
        Sleep(1000);
        //printf("data loading...\n");
        // wait for PDBQT data (loading triggered by MolecularDataCall)
    }
};


/**********
 * allocate the data arrays
 **********/
void MultiPDBQTLoader::allocateMultiPDBQTLoader(bool needToExtend) {

    int incrementStep_frameCapacity = 10000;
    if (needToExtend == true) {
        this->frameCapacity += incrementStep_frameCapacity;
        this->ligandsModelRange.SetCount(this->frameCapacity);
        if (verbose) printf("\nEXTENDING DATA ARRAYS...new frameCapacity: %d\n", frameCapacity);
    } else {
        this->ligandNames.AssertCapacity(this->ligandsNameCapacity);
        this->RMSD_lb.AssertCapacity(this->RMSD_lbCapacity);
        this->RMSD_ub.AssertCapacity(this->RMSD_ubCapacity);
        this->bindingEnergy.AssertCapacity(this->bindingEnergyCapacity);
        this->modelsPerLigand.AssertCapacity(this->ligandsNameCapacity);

        this->ligandsModelRange.SetCount(this->ligandsNameCapacity);

        this->atomTypes.AssertCapacity(this->atomTypesCapacity);
        this->atomTypesAD.AssertCapacity(this->atomTypesCapacity);

        this->atomPositions.AssertCapacity(this->atomPositionsCapacity);
        this->atomCharges.AssertCapacity(this->atomPositionsCapacity);
    }

    // SetCapacityIncrement
    this->ligandNames.SetCapacityIncrement(this->ligandsNameCapacity);
    this->RMSD_lb.SetCapacityIncrement(this->RMSD_lbCapacity);
    this->RMSD_ub.SetCapacityIncrement(this->RMSD_ubCapacity);
    this->bindingEnergy.SetCapacityIncrement(this->bindingEnergyCapacity);
    this->modelsPerLigand.SetCapacityIncrement(this->ligandsNameCapacity);
    this->atomTypes.SetCapacityIncrement(this->atomTypesCapacity);
    this->atomTypesAD.SetCapacityIncrement(this->atomTypesCapacity);
    this->atomPositions.SetCapacityIncrement(this->atomPositionsCapacity);
    this->atomCharges.SetCapacityIncrement(this->atomPositionsCapacity);


    this->bindingEnergyCapacity = incrementStep_frameCapacity;
    this->RMSD_lbCapacity = incrementStep_frameCapacity;
    this->RMSD_ubCapacity = incrementStep_frameCapacity;
    this->atomPositionsCapacity = incrementStep_frameCapacity * 100 * 3;
    this->ligandsNameCapacity = incrementStep_frameCapacity;
    this->atomTypesCapacity = incrementStep_frameCapacity * 100;

    if (this->residue.IsEmpty()) {
        MolecularDataCall::Residue* res = new MolecularDataCall::Residue(0, 0, this->bbox, 0, 0, 0);
        this->residue.Add(res);
    }
    if (this->residueTypeName.IsEmpty()) {
        this->residueTypeName.Add(" ");
    }
    if (this->molecule.IsEmpty()) {
        MolecularDataCall::Molecule m(0, 1, 0);
        this->molecule.Add(m);
    }
    if (this->chain.IsEmpty()) {
        MolecularDataCall::Chain c(0, 1);
        this->chain.Add(c);
    }
}


bool MultiPDBQTLoader::resetMultiPDBQTLoader() {
    if (static_cast<int>(this->ligandNames.Count()) > 0) {
        this->residue.Resize(0);

        this->total_atomCnt = 0;
        this->total_atomPosCnt = 0;
        this->total_ligandCnt = 0;
        this->total_modelCnt = 0;

        this->ligandNames.Clear();
        this->RMSD_lb.Clear();
        this->RMSD_ub.Clear();
        this->bindingEnergy.Clear();
        this->modelsPerLigand.Clear();

        this->ligandsModelRange.Clear();

        this->atomTypes.Clear();
        this->atomTypesAD.Clear();

        this->atomPositions.Clear();
        this->atomCharges.Clear();

        this->globalMdlID_to_ligID_mdlID.Clear();
        return true;
    } else {
        return true;
    }
}


/***********************************
*** MAIN DATA ACCESS functions: ***
***********************************/

void MultiPDBQTLoader::getLigandData(int ligandID) { this->getLigandData(ligandID, 0); }

void MultiPDBQTLoader::getLigandData(int ligandID, int modelID) {
    if (this->ligand.ligandID != ligandID && !this->isSilentDataCall) {
        this->mdc_dataHash++;
    }

    if (atomIndexSum.Count() != total_atomCnt) {
        atomIndexSum.SetCount(total_atomCnt);
        scanAtomIndex();
    }

    // only get data from data arrays if ligandID has changed
    if (ligandID != ligand.ligandID) {
        this->ligand.ligandID = ligandID;
        ligand.ligandAtomCnt = getLigandAtomCnt(ligandID);
        getAtomNamesAndTypesAD(ligandID);
        getModelCnt();
        this->ligand.ligandName = getLigandsName(ligandID);
        getSVGPath(ligandID);
        getLigandConnections(ligandID);
        getFunctionalGroups(ligandID);
        getChemProperties(ligandID);
    }
    this->ligand.Model.modelID = modelID;
    getModelDataOf_currentLigand();
}

/*********************
 * getter functions - NEVER CALL THESE FUNCTIONS, ONLY FOR INTERNAL USE!
 *********************/

bool MultiPDBQTLoader::getModelDataOf_currentLigand() {

    uint globalModelID = getGlobalModelID(this->ligand.ligandID, this->ligand.Model.modelID);

    this->ligand.Model.globalModelID = globalModelID;

    if (this->plipForcesStore.size() > globalModelID) {
        this->ligand.Model.Forces.hProtAcceptors = &this->plipForcesStore[globalModelID].hProtAcceptors;
        this->ligand.Model.Forces.hProtDonors = &this->plipForcesStore[globalModelID].hProtDonors;
        this->ligand.Model.Forces.hydrophobicInteractions =
            &this->plipForcesStore[globalModelID].hydrophobicInteractions;
        this->ligand.Model.Forces.halogenBonds = &this->plipForcesStore[globalModelID].halogenBonds;
        this->ligand.Model.Forces.metalComplexes = &this->plipForcesStore[globalModelID].metalComplexes;
        this->ligand.Model.Forces.piCationInteractions = &this->plipForcesStore[globalModelID].piCationInteractions;
        this->ligand.Model.Forces.piStacks = &this->plipForcesStore[globalModelID].piStacks;
        this->ligand.Model.Forces.saltBridges = &this->plipForcesStore[globalModelID].saltBridges;
    }


    if (!(this->modelRangeStart <= globalModelID) || !(globalModelID <= this->modelRangeEnd)) {
        printf("ERROR: choosen Model %d does not belong to current ligand %d (%d-%d Model)!\n", globalModelID,
            ligand.ligandID, this->modelRangeStart, this->modelRangeEnd);
        printf("lingandID: %d \t modelID: %d\n", this->ligand.ligandID, this->ligand.Model.modelID);
        // exit(EXIT_FAILURE);
    }

    this->ligand.Model.RMSD_lb = getRMSD_lb(globalModelID);
    this->ligand.Model.RMSD_ub = getRMSD_ub(globalModelID);
    this->ligand.Model.bindingEnergy = getBindingEnergy(globalModelID);
    this->ligand.Model.modelEfficiency = getModelEfficiency(globalModelID);

    this->getAtomPositions_And_Charges();
    return true;
}

vislib::StringA MultiPDBQTLoader::getLigandsName(uint ligandID) { return this->ligandNames[ligandID]; }
float MultiPDBQTLoader::getRMSD_lb(uint globalModelID) { return this->RMSD_lb[globalModelID]; }
float MultiPDBQTLoader::getRMSD_ub(uint globalModelID) { return this->RMSD_ub[globalModelID]; }
float MultiPDBQTLoader::getBindingEnergy(uint globalModelID) { return this->bindingEnergy[globalModelID]; }

uint MultiPDBQTLoader::getLigandAtomCnt(uint ligandID) { return this->ligandAtomCnts[ligandID]; }

float MultiPDBQTLoader::getModelEfficiency(uint globalModelID) {
    if (this->modelEfficiency.Count() == 0) {
        return 0.0;
    } else {
        return this->modelEfficiency[globalModelID];
    }
}

void MultiPDBQTLoader::getAtomNamesAndTypesAD(uint ligandID) {
    ligand.Model.atomTypesAD = &this->atomTypesAD[atomIndexSum[ligandID]];
    ligand.Model.atomTypeIdx = &this->atomTypeIdx[atomIndexSum[ligandID]];
    ligand.Model.atomTypes = &this->atomTypes;
}

void MultiPDBQTLoader::scanAtomIndex() {
    atomIndexSum[0] = 0;
    for (int i = 0; i < total_ligandCnt; i++) {
        atomIndexSum[i + 1] = atomIndexSum[i] + ligandAtomCnts[i];
    }
}

void MultiPDBQTLoader::getModelCnt() {
    uint ligandID = ligand.ligandID;
    this->modelRangeStart = ligandID > 0 ? this->ligandsModelRange[ligandID - 1] : 0;
    this->modelRangeEnd = this->ligandsModelRange[ligandID];

    ligand.modelCnt = (this->modelRangeEnd - this->modelRangeStart);
}

uint MultiPDBQTLoader::getGlobalModelID(uint ligandID, uint modelID) {
    return this->ligID_mdlID_to_globalMdlID[ligandID][modelID];
}

void MultiPDBQTLoader::getLigandConnections(uint ligandID) {
    if (connectivity.Count() > 0)
	{
		this->ligand.connectionCnt = this->connectivity[ligandID].Count() / 2;
		this->ligand.connections = this->connectivity[ligandID].PeekElements();
    }
}

void MultiPDBQTLoader::getSVGPath(uint ligandID) {
    // check size of SVGpaths: data is only present after CallForGetClusterData!
    if (this->SVGPaths.size() == this->total_ligandCnt) {
        this->ligand.SVGPath = this->SVGPaths[ligandID];
    } else {
        this->ligand.SVGPath = "";
    }
}
void MultiPDBQTLoader::getFunctionalGroups(uint ligandID) {
    this->ligand.functionalGroupsLigand = &functionalGroupsStorage[ligandID];
}

void MultiPDBQTLoader::getChemProperties(uint ligandID) { this->ligand.chemPropsLigand = &chemPropsStorage[ligandID]; }

void MultiPDBQTLoader::testPrinting(int LigID, int ModID) {
    /********************************************
     *************** TEST PRINTING ***************
     **********************************************/
    printf("total_ligandCnt: %d\ntotal_modelCnt: %d \ntotal_atomCnt: %d \n", this->total_ligandCnt,
        this->total_modelCnt, this->total_atomCnt);
    // Sleep(2000);

    if (LigID > -1 && ModID > -1) {

        this->getLigandData(LigID, ModID);
        printf("Ligand \t\tName \t \t Model \t RMSDlb\tRMSDub \t Frame \t Energy\t Atom \t   x \t   y \t   z "
               "\t q: \t "
               "Type: \t TypeAD: \n");
        for (int z = 0; z < this->ligand.ligandAtomCnt; z++) {
            printf("  %d \t %s \t %d \t %.2f \t %.2f \t   %d \t %.2f \t %d \t %.1f \t %.1f \t \%.1f \t %.2f \t "
                   "  %s "
                   "\t   %s \n",
                this->ligand.ligandID, this->ligand.ligandName, this->ligand.modelCnt, this->ligand.Model.RMSD_lb,
                this->ligand.Model.RMSD_ub, this->ligand.Model.modelID, this->ligand.Model.bindingEnergy, z,
                this->ligand.Model.atomPostions[3 * z + 0], this->ligand.Model.atomPostions[3 * z + 1],
                this->ligand.Model.atomPostions[3 * z + 2], this->ligand.Model.atomCharges[z],
                this->atomTypes[this->ligand.Model.atomTypeIdx[z]].Name(), this->ligand.Model.atomTypesAD[z]);
        }
    } else {
        for (int i = 0; i < this->total_ligandCnt; i++) {

            // sets current ligand - call for ligand name and first model
            // this->getLigandData(i);
            // printf("Ligand: %d   Atoms: %d \n", i, this->ligandsAtomCnt[i]);
            printf("***************************************************************************************************"
                   "***************\n",
                i, this->ligandAtomCnts[i]);

            for (int j = 0; j < this->ligand.modelCnt; j++) {
                // call for the data of a model from the current ligand
                // d.getModelDataOf_currentLigand(j);
                this->getLigandData(i, j);
                printf("Ligand \t\tName \t \t Model \t RMSDlb\tRMSDub \t Frame \t Energy\t Atom \t   x \t   y \t   z "
                       "\t q: \t "
                       "Type: \t TypeAD: \n");
                for (int z = 0; z < this->ligand.ligandAtomCnt; z++) {
                    printf("  %d \t %s \t %d \t %.2f \t %.2f \t   %d \t %.2f \t %d \t %.1f \t %.1f \t \%.1f \t %.2f \t "
                           "  %s "
                           "\t   %s \n",
                        this->ligand.ligandID, this->ligand.ligandName, this->ligand.modelCnt,
                        this->ligand.Model.RMSD_lb, this->ligand.Model.RMSD_ub, this->ligand.Model.modelID,
                        this->ligand.Model.bindingEnergy, z, this->ligand.Model.atomPostions[3 * z + 0],
                        this->ligand.Model.atomPostions[3 * z + 1], this->ligand.Model.atomPostions[3 * z + 2],
                        this->ligand.Model.atomCharges[z], this->atomTypes[this->ligand.Model.atomTypeIdx[z]].Name(),
                        this->ligand.Model.atomTypesAD[z]);
                }
            }
        }
    }
}
// DO NOT TOUCH THIS FUNCTION -> it works!
void MultiPDBQTLoader::getAtomPositions_And_Charges() {

    int offset = 0;
    int numberOfModels = 0;
    int position = 0;

    for (int i = 0; i < ligand.ligandID; i++) {
        numberOfModels = ((this->ligandsModelRange[i]) - (i > 0 ? this->ligandsModelRange[i - 1] : 0));
        offset += this->ligandAtomCnts[i] * numberOfModels;
    }

    position = ligand.ligandAtomCnt * ligand.Model.modelID;

    ligand.Model.atomPostions = &atomPositions[offset * 3 + position * 3];
    ligand.Model.atomCharges = &atomCharges[offset + position];
}


bool MultiPDBQTLoader::paramsDBSCANchanged(param::ParamSlot& param) {
    // set param as clean again
    param.ResetDirty();
    this->DBSCANParamsChanged = TRUE;
    this->eps = this->epsParam.Param<param::FloatParam>()->Value();
    this->minPts = this->minPtsParam.Param<param::IntParam>()->Value();
    this->minBindEnergy = this->minBindEnergyParam.Param<param::FloatParam>()->Value();
    if (this->minBindEnergy > 0.0f) {
        this->minBindEnergy = 0;
        this->minBindEnergyParam.Param<param::FloatParam>()->SetValue(this->minBindEnergy);
    }
    if (this->minPts < 1) {
        this->minPts = 1;
        this->minPtsParam.Param<param::IntParam>()->SetValue(this->minPts);
    }
    if (this->eps < 0.f) {
        this->eps = 0.f;
        this->epsParam.Param<param::FloatParam>()->SetValue(this->eps);
    }


    return true;
};

bool MultiPDBQTLoader::getSVGsAndChemProps(const vislib::TString& filename) {

    auto start = std::chrono::high_resolution_clock::now();

    const char* scriptName = "getSMItoSVG.py";
    const char* scriptName_ChemProps = "getChemPropsIntoCheckmol.py";
    int argvCount = 6;

    Path dataFolderPath;
    Path baseDir;
    Path babelPath;
    Path pythonScriptPath;
    Path pythonScriptPath_ChemProps;

    std::filesystem::path cwd = std::filesystem::current_path();
    baseDir = getDirpath_toMegaMol_from_filepath(cwd.string());


    pythonScriptPath.append(baseDir.string());
    pythonScriptPath.append(this->prolintSubPath.string());
    pythonScriptPath.append(scriptName);
    printf("python: %s\n", pythonScriptPath.string().c_str());

    pythonScriptPath_ChemProps.append(baseDir.string());
    pythonScriptPath_ChemProps.append(this->prolintSubPath.string());
    pythonScriptPath_ChemProps.append(scriptName_ChemProps);
    printf("python: %s\n", pythonScriptPath_ChemProps.string().c_str());

    babelPath.append(baseDir.string());
    babelPath.append(this->prolintSubPath.string());
    babelPath.append("obabel.exe");

    dataFolderPath.append(getAdditionalDataFolderPath().string());
    // create dir
    std::filesystem::create_directory(dataFolderPath);


    wchar_t** argv;

    float percent;


    string ZINCIDListPath;
    ZINCIDListPath.append(getAdditionalDataFolderPath().string());
    ZINCIDListPath.append("ZINCIDList.txt");

    ofstream ZINCNameListFile;
    ZINCNameListFile.open(ZINCIDListPath.c_str());
    for (int i = 0; i < this->ligandNames.Count(); i++) {

        // if (this->int2String.find(lmc->getCurrentLigandID()) == int2String.end() ? 0:1) {

        string ZINCID = this->ligandNames[i].PeekBuffer();
        Path SVGPath;
        SVGPath.append(getAdditionalDataFolderPath().string());
        SVGPath.append(ZINCID + ".svg");

        this->SVGPaths.push_back(SVGPath.string());
        ifstream ifile;
        ifile.open(SVGPath);
        percent = (float)(i + 1) / (float)this->ligandNames.Count() * 100.0f;
        if (ifile) {
            if (verbose) printf("%.2f%% file exists: %s\n", percent, SVGPath.string().c_str());
        } else {
            ZINCNameListFile << ZINCID;
            ZINCNameListFile << ";";
        }
        //}
    }
    ZINCNameListFile.close();


    argv = (wchar_t**)malloc(sizeof(wchar_t*) * argvCount);
    for (int i = 0; i < argvCount; i++) {
        argv[i] = (wchar_t*)malloc(sizeof(wchar_t) * 1);
    }

    string MultiPDBQTList = filename.PeekBuffer();


    size_t size0 = sizeof(scriptName);
    argv[0] = (wchar_t*)malloc(sizeof(wchar_t) * size0);
    argv[0] = Py_DecodeLocale(scriptName, &size0);

    size_t size1 = sizeof(babelPath.c_str());
    argv[1] = (wchar_t*)malloc(sizeof(wchar_t) * size1);
    argv[1] = Py_DecodeLocale(babelPath.string().c_str(), &size1);

    size_t size2 = sizeof(baseDir.c_str());
    argv[2] = (wchar_t*)malloc(sizeof(wchar_t) * size2);
    argv[2] = Py_DecodeLocale(baseDir.string().c_str(), &size2);

    size_t size3 = sizeof(getAdditionalDataFolderPath());
    argv[3] = (wchar_t*)malloc(sizeof(wchar_t) * size3);
    argv[3] = Py_DecodeLocale(getAdditionalDataFolderPath().string().c_str(), &size3);

    size_t size4 = sizeof(ZINCIDListPath.c_str());
    argv[4] = (wchar_t*)malloc(sizeof(wchar_t) * size4);
    argv[4] = Py_DecodeLocale(ZINCIDListPath.c_str(), &size4);

    size_t size5 = sizeof(MultiPDBQTList.c_str());
    argv[5] = (wchar_t*)malloc(sizeof(wchar_t) * size5);
    argv[5] = Py_DecodeLocale(MultiPDBQTList.c_str(), &size5);

    FILE* fp;
    bool success = 0;


    std::string fMode = "r";
    size_t sizeFileString = sizeof(pythonScriptPath.string().c_str());
    size_t sizefileMode = sizeof(fMode.c_str());
    wchar_t* fileString;
    wchar_t* fileMode;
    fileString = (wchar_t*)malloc(sizeof(wchar_t) * sizeFileString);
    fileString = Py_DecodeLocale(pythonScriptPath.string().c_str(), &sizeFileString);
    fileMode = (wchar_t*)malloc(sizeof(wchar_t) * sizefileMode);
    fileMode = Py_DecodeLocale(fMode.c_str(), &sizefileMode);

   
    while (success != 1) {
        try {
            fp = _Py_wfopen(fileString, fileMode);

            PySys_SetArgv(argvCount, argv);

            PyRun_SimpleFile(fp, pythonScriptPath.string().c_str());
            success = 1;
        } catch (const std::exception& e) {
            printf("Python excecution failed: %s\n", e);
        }
    }

    std::string fMode2 = "r";
    size_t sizeFileString2 = sizeof(pythonScriptPath_ChemProps.string().c_str());
    size_t sizefileMode2 = sizeof(fMode2.c_str());
    wchar_t* fileString2;
    wchar_t* fileMode2;
    fileString2 = (wchar_t*)malloc(sizeof(wchar_t) * sizeFileString2);
    fileString2 = Py_DecodeLocale(pythonScriptPath_ChemProps.string().c_str(), &sizeFileString2);
    fileMode2 = (wchar_t*)malloc(sizeof(wchar_t) * sizefileMode2);
    fileMode2 = Py_DecodeLocale(fMode2.c_str(), &sizefileMode2);

    success = 0;
    while (success != 1) {
        try {
            fp = _Py_wfopen(fileString2, fileMode2);

            PySys_SetArgv(argvCount, argv);

            PyRun_SimpleFile(fp, pythonScriptPath_ChemProps.string().c_str());
            success = 1;
        } catch (const std::exception& e) {
            printf("Python excecution failed: %s\n", e);
        }
    }

    // Py_Finalize();

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

    return true;
}


bool MultiPDBQTLoader::getfunctionalGroupDataFiles(const vislib::TString& filename) {


    auto start = std::chrono::high_resolution_clock::now();

    const char* scriptName = "getFunctionalGroupData.py";

    int argvCount = 4;

    Path dataFolderPath;
    Path baseDir;
    Path babelPath;
    Path checkmolPath;
    Path pythonScriptPath;

    std::filesystem::path cwd = std::filesystem::current_path();
    baseDir = getDirpath_toMegaMol_from_filepath(cwd.string());

    pythonScriptPath.append(baseDir.string());
    pythonScriptPath.append(this->prolintSubPath.string());
    pythonScriptPath.append(scriptName);
    printf("python: %s\n", pythonScriptPath.string().c_str());


    checkmolPath.append(baseDir.string());
    checkmolPath.append(this->prolintSubPath.string());
    checkmolPath.append("checkmol.exe");


    wchar_t** argv;

    float percent;


    string MultiPDBQTList = filename.PeekBuffer();


    argv = (wchar_t**)malloc(sizeof(wchar_t*) * argvCount);
    for (int i = 0; i < argvCount; i++) {
        argv[i] = (wchar_t*)malloc(sizeof(wchar_t) * 1);
    }

    size_t size0 = sizeof(scriptName);
    argv[0] = (wchar_t*)malloc(sizeof(wchar_t) * size0);
    argv[0] = Py_DecodeLocale(scriptName, &size0);

    size_t size1 = sizeof(MultiPDBQTList.c_str());
    argv[1] = (wchar_t*)malloc(sizeof(wchar_t) * size1);
    argv[1] = Py_DecodeLocale(MultiPDBQTList.c_str(), &size1);

    size_t size2 = sizeof(getAdditionalDataFolderPath());
    argv[2] = (wchar_t*)malloc(sizeof(wchar_t) * size2);
    argv[2] = Py_DecodeLocale(getAdditionalDataFolderPath().string().c_str(), &size2);

    size_t size3 = sizeof(checkmolPath.c_str());
    argv[3] = (wchar_t*)malloc(sizeof(wchar_t) * size3);
    argv[3] = Py_DecodeLocale(checkmolPath.string().c_str(), &size3);

    FILE* fp;
    bool success = 0;

    std::string fMode = "r";
    size_t sizeFileString = sizeof(pythonScriptPath.string().c_str());
    size_t sizefileMode = sizeof(fMode.c_str());
    wchar_t* fileString;
    wchar_t* fileMode;
    fileString = (wchar_t*)malloc(sizeof(wchar_t) * sizeFileString);
    fileString = Py_DecodeLocale(pythonScriptPath.string().c_str(), &sizeFileString);
    fileMode = (wchar_t*)malloc(sizeof(wchar_t) * sizefileMode);
    fileMode = Py_DecodeLocale(fMode.c_str(), &sizefileMode);

    while (success != 1) {
        try {
            fp = _Py_wfopen(fileString, fileMode);

            PySys_SetArgv(argvCount, argv);

            PyRun_SimpleFile(fp, pythonScriptPath.string().c_str());
            success = 1;
        } catch (const std::exception& e) {
            printf("Python excecution failed: %s\n", e);
        }
    }

    // Py_Finalize();

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

    return true;
};

bool MultiPDBQTLoader::getProteinLigand_Interactions(const vislib::TString& filename) {
    MolecularDataCall* pmdc = this->dataInSlot_MolecularDataCall.CallAs<MolecularDataCall>();
    if (pmdc == NULL) {
        return false;
    }
    (*pmdc)(megamol::protein_calls::MolecularDataCall::CallForGetExtent);
    (*pmdc)(megamol::protein_calls::MolecularDataCall::CallForGetData);


    auto start = std::chrono::high_resolution_clock::now();

    const char* scriptName = "getProteinLigand_Interactions.py";

    int argvCount = 5;

    Path dataFolderPath;
    Path baseDir;
    Path babelPath;
    Path proteinFilePath = Path(pmdc->GetPDBFilename().PeekBuffer());
    Path pythonScriptPath;


    std::filesystem::path cwd = std::filesystem::current_path();
    baseDir = getDirpath_toMegaMol_from_filepath(cwd.string());

    pythonScriptPath.append(baseDir.string());
    pythonScriptPath.append(this->prolintSubPath.string());
    pythonScriptPath.append(scriptName);
    printf("python: %s\n", pythonScriptPath.string().c_str());

    babelPath.append(baseDir.string());
    babelPath.append(this->prolintSubPath.string());
    babelPath.append("obabel.exe");

    wchar_t** argv;

    float percent;


    string MultiPDBQTList = filename.PeekBuffer();


    argv = (wchar_t**)malloc(sizeof(wchar_t*) * argvCount);
    for (int i = 0; i < argvCount; i++) {
        argv[i] = (wchar_t*)malloc(sizeof(wchar_t) * 1);
    }

    size_t size0 = sizeof(scriptName);
    argv[0] = (wchar_t*)malloc(sizeof(wchar_t) * size0);
    argv[0] = Py_DecodeLocale(scriptName, &size0);

    size_t size1 = sizeof(babelPath.c_str());
    argv[1] = (wchar_t*)malloc(sizeof(wchar_t) * size1);
    argv[1] = Py_DecodeLocale(babelPath.string().c_str(), &size1);

    size_t size2 = sizeof(MultiPDBQTList.c_str());
    argv[2] = (wchar_t*)malloc(sizeof(wchar_t) * size2);
    argv[2] = Py_DecodeLocale(MultiPDBQTList.c_str(), &size2);

    size_t size3 = sizeof(getAdditionalDataFolderPath());
    argv[3] = (wchar_t*)malloc(sizeof(wchar_t) * size3);
    argv[3] = Py_DecodeLocale(getAdditionalDataFolderPath().string().c_str(), &size3);

    size_t size4 = sizeof(proteinFilePath.c_str());
    argv[4] = (wchar_t*)malloc(sizeof(wchar_t) * size4);
    argv[4] = Py_DecodeLocale(proteinFilePath.string().c_str(), &size4);

    FILE* fp;
    bool success = 0;

    std::string fMode = "r";
    size_t sizeFileString = sizeof(pythonScriptPath.string().c_str());
    size_t sizefileMode = sizeof(fMode.c_str());
    wchar_t* fileString;
    wchar_t* fileMode;
    fileString = (wchar_t*)malloc(sizeof(wchar_t) * sizeFileString);
    fileString = Py_DecodeLocale(pythonScriptPath.string().c_str(), &sizeFileString);
    fileMode = (wchar_t*)malloc(sizeof(wchar_t) * sizefileMode);
    fileMode = Py_DecodeLocale(fMode.c_str(), &sizefileMode);
    while (success != 1) {
        try {
            fp = _Py_wfopen(fileString, fileMode);

            PySys_SetArgv(argvCount, argv);

            PyRun_SimpleFile(fp, pythonScriptPath.string().c_str());
            success = 1;
        } catch (const std::exception& e) {
            printf("Python excecution failed: %s\n", e);
        }
    }

    // Py_Finalize();

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

    return true;
};


bool MultiPDBQTLoader::setLMCStatus(core::Call& call) {

    LigandModelCall* lmc = dynamic_cast<LigandModelCall*>(&call);
    if (lmc == NULL) return false;
    auto status = lmc->getLMCStatus();

    this->MPL_modelClustRenderer_busy = status.isModelClustRendererBusy;
    this->MPL_SocketCommManager_busy = status.isSocketCommManagerBusy;
    this->MPL_firstRenderCall_complete = status.isFirstRenderCallComplete;
}


/************************
 ***** DBSCAN CLASS *****
 ************************/
DBSCAN::DBSCAN(void) { this->minValue = 0; };

void DBSCAN::release() {
    if (this->neighbourList_ffrnns.size() > 0) {
        for (int i = 0; i < this->neighbourList_ffrnns.size(); i++) {
            this->neighbourList_ffrnns[i].resize(0);
        }
    }
    this->neighbourList_ffrnns.resize(0);
    this->unVistedPoints.Resize(0);
    this->clusters.Resize(0);
}

vislib::Array<int> DBSCAN::DBSCAN_3D(visFloat4& dataVec4, float eps, int minPts, float minValue, int showProgress) {
    if (minValue != NULL) {
        this->minValue = minValue;
    }
    release();
    this->unVistedPointsCnt = 0;
    this->unVistedPoints.SetCount(dataVec4.Count());
    this->clusters.SetCount(dataVec4.Count());

    // initilaize
    for (int i = 0; i < dataVec4.Count(); i++) {
        this->unVistedPoints[i] = 0;
        this->clusters[i] = -1;
    }


    /** Cluster: index of cluster
     *
     */
    int cluster = 0;

    ////////////////////////////
    // getAllNeighbours

    auto start = std::chrono::high_resolution_clock::now();
    FFRNNS rnss;
    rnss.setData(dataVec4, eps, 0.0f);
    this->neighbourList_ffrnns = rnss.getNeigbours();

    if (showProgress) {
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "Elapsed time neighbour search - FFRNNS: " << elapsed.count() << " s\n";
    }
    // filter for energy
#pragma omp parallel for
    for (int i = 0; i < dataVec4.Count(); i++) {

        if ((abs(dataVec4[i].GetW()) < abs(this->minValue))) {
            this->neighbourList_ffrnns[i].resize(0);
            this->unVistedPoints[i] = 1;
            this->clusters[i] = -1;
        } else {
            std::vector<uint> tmp;
            for (int j = 0; j < neighbourList_ffrnns[i].size(); j++) {
                int index = neighbourList_ffrnns[i][j];
                if ((abs(dataVec4[index].GetW()) >= abs(this->minValue))) {
                    tmp.push_back(index);
                }
            }
            neighbourList_ffrnns[i].resize(tmp.size());
            neighbourList_ffrnns[i] = tmp;
        }
        if (showProgress) {
            if (i % 100 == 0) {
                printf("filtering.../      \r");
                printf("filtering...--      \r");
                printf("filtering...\      \r");
            }
        }
    }

    ////////////////////////////


    // start of DBSCAN => iterating over all points
    for (int i = 0; i < dataVec4.Count(); i++) {
        if (this->unVistedPoints[i] != 1) {
            this->unVistedPoints[i] = 1;
            this->unVistedPointsCnt++;
            // sort(this->neighbourList_ffrnns[i].begin(), this->neighbourList_ffrnns[i].end());
            if (this->neighbourList_ffrnns[i].size() < minPts || this->neighbourList_ffrnns[i].size() == 0) {
                this->clusters[i] = -1; // mark as noise
            } else {
                if (this->clusters[i] <= -1) {
                    this->clusters[i] = cluster;
                }
                expandCluster(dataVec4, this->neighbourList_ffrnns[i], cluster, minPts, showProgress);
                cluster++; // next cluster
            }
        }
    }
    return this->clusters;
};


void DBSCAN::expandCluster(
    visFloat4& dataVec4, std::vector<uint>& N_joined, int cluster, int minPts, int showProgress) {
    /**
     * N_star / neighbourList[NJ_Idx]:		Neighbours of a point from the centroid neighbours in the cluster
     * expansion
     * step N_joined:							Neigbhbours of the centroid and all N_star for current centroid in
     * expansion step
     */

    // if the cluster expands (N_joined grows) the loop range grows automatically!
    for (int a = 0; a < N_joined.size(); a++) {
        int NJ_Idx = N_joined[a]; // N_joinedIdx

        if (showProgress) {
            int actualProgress = (int)((float)this->unVistedPointsCnt / (float)dataVec4.Count() * 10000);
            if (this->oldProgress != actualProgress) {
                this->oldProgress = actualProgress;
                printf(
                    "searching clusters...  found: %d \t %.2f %% done      \r", cluster, (float)actualProgress / 100);
            }
            if (actualProgress == 10000) {
                printf(
                    "searching clusters...  found: %d \t %.2f %% done      \n", cluster, (float)actualProgress / 100);
            }
        }


        //  if (this->unVistedPoints[NJ_Idx] != 1) {
        // mark as visited
        this->unVistedPointsCnt++;
        this->unVistedPoints[NJ_Idx] = 1;
        // this->neighbourList[NJ_Idx] is equal to N_star
        if (this->clusters[NJ_Idx] <= -1) {
            this->clusters[NJ_Idx] = cluster;
        }

        if (this->neighbourList_ffrnns[NJ_Idx].size() >= minPts) {
            // Sort the vector
            std::vector<uint> t1;
            t1.resize(this->neighbourList_ffrnns[NJ_Idx].size());
            t1 = this->neighbourList_ffrnns[NJ_Idx];
            sort(t1.begin(), t1.end());

            std::vector<uint> t0;
            t0.resize(N_joined.size());
            t0 = N_joined;
            sort(t0.begin(), t0.end());

            // execution policy festlegen
            // set_difference
            std::vector<uint> dest1;
            std::set_difference(t1.begin(), t1.end(), t0.begin(), t0.end(), std::back_inserter(N_joined));
        }
    }
}

// returns the distance between to vec3 vectors
float DBSCAN::distance2Vectors(vislib::math::Vector<float, 3> v1, vislib::math::Vector<float, 3> v2) {
    float dis = sqrtf(powf(v2.X() - v1.X(), 2) + powf(v2.Y() - v1.Y(), 2) + powf(v2.Z() - v1.Z(), 2));
    return dis;
}


/**********************************
 **** INTERACTION FORCES CLASS ****
 **********************************/

// *Ctor* of InteractionForces
InteractionForces::InteractionForces(megamol::core::Call* ligandModelCall,
    megamol::core::Call* molecularDataCall_receptor, megamol::prolint::LigandModelCall::Ligand& ligand) {

    this->lmc = dynamic_cast<megamol::prolint::LigandModelCall*>(ligandModelCall);
    if (this->lmc == NULL) {
        printf("Missing or wrong Call for constructor of InteractionForces class: needs LigandModelCall!");
        assert(this->lmc != NULL);
    }
    (*lmc)(megamol::prolint::LigandModelCall::CallForGetExtent);
    (*lmc)(megamol::prolint::LigandModelCall::CallForGetExtent);

    this->pmdc = dynamic_cast<megamol::protein_calls::MolecularDataCall*>(molecularDataCall_receptor);
    if (this->pmdc == NULL) {
        printf("Missing or wrong Call for constructor of InteractionForces class: needs MolecularDataCall for receptor "
               "protein!");
        assert(this->pmdc != NULL);
    }
    (*pmdc)(megamol::protein_calls::MolecularDataCall::CallForGetExtent);
    (*pmdc)(megamol::protein_calls::MolecularDataCall::CallForGetData);

    this->ligand = &ligand;
};

// *Dtor* of InteractionForces
void InteractionForces::release(void){
    // needs to implemented
};


bool InteractionForces::calcHBonds(vislib::Array<vislib::Array<uint>>& out_ligMdlTo_prot,
    vislib::Array<vislib::Array<uint>>& out_protTo_LigMdl, vislib::Array<vislib::Array<uint>>& out_ligMdlTo_protH,
    vislib::Array<vislib::Array<uint>>& out_protTo_LigMdlH, vislib::Array<uint>& hBondAtomCnt_protein) {

    /*
        // index % 2 != 0 --> proteinAtomID; index % 2 == 0 --> ligandModelAtomID
        uint hBonds_ligMdlHDonorCount;
        uint* hBonds_ligMdlHDonor;

        // index % 2 != 0 --> proteinAtomID; index % 2 == 0 --> ligandModelAtomID
        uint hBonds_ligMdlHAcceptorCount;
        uint* hBonds_ligMdlHAcceptor;

        uint hBonds_ligMdlDonorCount;
        // index % 2 != 0 --> proteinAtomID; index % 2 == 0 --> ligandModelAtomID
        uint* hBonds_ligMdlDonor;
        uint hBonds_ligMdlAcceptorCount;
        // index % 2 != 0 --> proteinAtomID; index % 2 == 0 --> ligandModelAtomID
        uint* hBonds_ligMdlAcceptor;
        // number of H-bonds counted for every protein atom
        uint* hBondAtomCnt_protein;
        // number of protein atoms
        uint hBondAtomCnt_proteinCount;
        // max number of H-bonds counted of every protein atom
        uint max_hBondAtomCnt_protein;
    */


    if (this->pmdc == NULL) return false;
    if (this->lmc == NULL) return false;

    float distanceToFromHBond = 3.0f;
    float angleToFromHBond = 30.0f;

    this->H_donorIndices_prot.resize(0);
    this->donorIndices_prot.resize(0);
    this->acceptorIndices_prot.resize(0);

    this->proteinData_H_donors.Resize(0);
    this->proteinData_donors.Resize(0);
    this->proteinData_acceptors.Resize(0);

    this->H_donorIndices_ligMdl.resize(0);
    this->donorIndices_ligMdl.resize(0);
    this->acceptorIndices_ligMdl.resize(0);

    this->ligandModelData_H_donors.Resize(0);
    this->ligandModelData_donors.Resize(0);
    this->ligandModelData_acceptors.Resize(0);

    hBondAtomCnt_protein.Resize(0);

    if (Hbonds_protTo_LigMdl_checked.size() > 0) {
        for (int i = 0; i < Hbonds_protTo_LigMdl_checked.size(); i++) {
            Hbonds_protTo_LigMdl_checked[i].resize(0);
        }
    }
    this->Hbonds_protTo_LigMdl_checked.resize(0);

    if (Hbonds_LigMdlTo_prot_checked.size() > 0) {
        for (int i = 0; i < Hbonds_LigMdlTo_prot_checked.size(); i++) {
            Hbonds_LigMdlTo_prot_checked[i].resize(0);
        }
    }
    this->Hbonds_LigMdlTo_prot_checked.resize(0);

    if (out_ligMdlTo_prot.Count() > 0) {
        for (int i = 0; i < out_ligMdlTo_prot.Count(); i++) {
            out_ligMdlTo_prot[i].Resize(0);
        }
    }
    out_ligMdlTo_prot.SetCount(lmc->getTotalModelCount());

    if (out_protTo_LigMdl.Count() > 0) {
        for (int i = 0; i < out_protTo_LigMdl.Count(); i++) {
            out_protTo_LigMdl[i].Resize(0);
        }
    }
    out_protTo_LigMdl.SetCount(lmc->getTotalModelCount());

    if (out_ligMdlTo_protH.Count() > 0) {
        for (int i = 0; i < out_ligMdlTo_protH.Count(); i++) {
            out_ligMdlTo_protH[i].Resize(0);
        }
    }
    out_ligMdlTo_protH.SetCount(lmc->getTotalModelCount());

    if (out_protTo_LigMdlH.Count() > 0) {
        for (int i = 0; i < out_protTo_LigMdlH.Count(); i++) {
            out_protTo_LigMdlH[i].Resize(0);
        }
    }
    out_protTo_LigMdlH.SetCount(lmc->getTotalModelCount());


    // electro-negativities of non metals
    std::map<std::string, float> electroNegativities{{"H", 2.20f}, {"C", 2.50f}, {"N", 3.04f}, {"P", 2.19f}, {"O", 3.44f},
        {"S", 2.58f}, {"Se", 2.55f}, {"F", 3.98f}, {"Cl", 3.16f}, {"Br", 2.96f}, {"I", 2.10f}, {"A", 2.50f}};


    // get protein data
    if (!(*pmdc)(MolecularDataCall::CallForGetData)) {
        return 0;
    }

    // stores the globalModelID and the corresponding intern atom index of Donor/Acceptor atom
    vislib::Array<vislib::math::Vector<int, 2>> invGblMdlIDs_p_to_m;
    vislib::Array<vislib::math::Vector<int, 2>> invGblMdlIDs_m_to_p;
    // stores the globalModelID and the corresponding intern atom index of Hydrogen Donor atom
    vislib::Array<vislib::math::Vector<int, 2>> invGblMdlIDs_mH_to_p;
    vislib::Array<uint> invProtIDs_donor;
    vislib::Array<uint> invProtIDs_acceptor;


    vislib::Array<uint> checkAtomDonor;
    vislib::Array<uint> checkAtomAcceptor;
    checkAtomDonor.SetCount(pmdc->AtomCount());
    checkAtomAcceptor.SetCount(pmdc->AtomCount());
    for (int i = 0; i < pmdc->AtomCount(); i++) {
        checkAtomDonor[i] = 0;
        checkAtomAcceptor[i] = 0;
    }

    uint exc = 0;
    //#pragma omp parallel for
    for (int i = 0; i < pmdc->ConnectionCount(); i++) {
        const uint* cons = pmdc->Connection();
        int atom0Index = cons[i * 2 + 0];
        int atom1Index = cons[i * 2 + 1];

        if (atom1Index > atom0Index) {
            const megamol::protein_calls::MolecularDataCall::AtomType tmpAtomType_0 =
                pmdc->AtomTypes()[pmdc->AtomTypeIndices()[atom0Index]];
            const megamol::protein_calls::MolecularDataCall::AtomType tmpAtomType_1 =
                pmdc->AtomTypes()[pmdc->AtomTypeIndices()[atom1Index]];

            string check = tmpAtomType_0.Element().PeekBuffer();
            string check2 = tmpAtomType_1.Element().PeekBuffer();

            // if hydrogen is present it will be always atom0
            if (check == "H") {

            } else if (check2 == "H") {
                atom0Index = cons[i * 2 + 1];
                atom1Index = cons[i * 2 + 0];
            }


            const megamol::protein_calls::MolecularDataCall::AtomType atom0 =
                pmdc->AtomTypes()[pmdc->AtomTypeIndices()[atom0Index]];
            const megamol::protein_calls::MolecularDataCall::AtomType atom1 =
                pmdc->AtomTypes()[pmdc->AtomTypeIndices()[atom1Index]];

            vislib::StringA test0 = atom0.Element();
            vislib::StringA test1 = atom1.Element();
            string t0 = test0.PeekBuffer();
            string t1 = test1.PeekBuffer();
            float delatEN = std::abs(electroNegativities[t0] - electroNegativities[t1]);

            check = atom0.Element().PeekBuffer();
            // check if hydrogen is present
            if (check == "H") {
                if (delatEN > 0.5) {

                    if (!checkAtomDonor[atom0Index]) {
                        checkAtomDonor[atom0Index] = exc;
                        H_donorIndices_prot.push_back(atom0Index);
                        proteinData_H_donors.Append({
                            pmdc->AtomPositions()[atom0Index * 3 + 0], pmdc->AtomPositions()[atom0Index * 3 + 1],
                            pmdc->AtomPositions()[atom0Index * 3 + 2],
                            0.0f // atom0.Radius()
                        });
                    }

                    if (!checkAtomDonor[atom1Index]) {
                        checkAtomDonor[atom1Index] = exc;
                        donorIndices_prot.push_back(atom1Index);
                        proteinData_donors.Append({
                            pmdc->AtomPositions()[atom1Index * 3 + 0], pmdc->AtomPositions()[atom1Index * 3 + 1],
                            pmdc->AtomPositions()[atom1Index * 3 + 2],
                            0.0f // atom1.Radius()
                        });
                    }

                    if (!checkAtomAcceptor[atom1Index]) {
                        acceptorIndices_prot.push_back(atom1Index);
                        checkAtomAcceptor[atom1Index] = exc;
                        proteinData_acceptors.Append({
                            pmdc->AtomPositions()[atom1Index * 3 + 0], pmdc->AtomPositions()[atom1Index * 3 + 1],
                            pmdc->AtomPositions()[atom1Index * 3 + 2],
                            0.0f // atom1.Radius()
                        });
                    }
                }
                // if hydrogen is not present it only can be an accpetor
            } else {
                if (delatEN > 0.5) {
                    // printf("%s - %s\n", t0.c_str(), t1.c_str());
                    // who has the charge? push that one with higher electroNegativity
                    if (electroNegativities[t0] > electroNegativities[t1]) {
                        if (!checkAtomAcceptor[atom0Index]) {
                            checkAtomAcceptor[atom0Index] = exc;
                            acceptorIndices_prot.push_back(atom0Index);
                            proteinData_acceptors.Append({
                                pmdc->AtomPositions()[atom0Index * 3 + 0], pmdc->AtomPositions()[atom0Index * 3 + 1],
                                pmdc->AtomPositions()[atom0Index * 3 + 2],
                                0.0f // atom0.Radius()
                            });
                        }
                    } else {
                        if (!checkAtomAcceptor[atom1Index]) {
                            checkAtomAcceptor[atom1Index] = exc;
                            acceptorIndices_prot.push_back(atom1Index);
                            proteinData_acceptors.Append({
                                pmdc->AtomPositions()[atom1Index * 3 + 0], pmdc->AtomPositions()[atom1Index * 3 + 1],
                                pmdc->AtomPositions()[atom1Index * 3 + 2],
                                0.0f // atom1.Radius()
                            });
                        }
                    }
                }
            }
        }
    }

    ///////////////////////////////////////////////////////////////

    // get ligand model data

    this->clusterData = lmc->getClusterAndCentroidData();
    auto clusterIndexSum = this->clusterData.clusterIndexSum->PeekElements();
    for (int i = 0; i < this->clusterData.clusterCentroids->Count() / 4; i++) {
        int clusterIDx_start = i > 0 ? clusterIndexSum[i] : clusterIndexSum[0];
        int clusterIDx_end = clusterIndexSum[i + 1];

        for (int j = clusterIDx_start; j < clusterIDx_end; j++) {
            lmc->SetCurrentLigandAndModel_byGlobalMdlID(
                this->clusterData.assign_modelID_sphereID->PeekElements()[j].GetY());
            (*lmc)(LigandModelCall::CallForGetData);
            uint gblMdlID = lmc->getCurrentGlobalModelID();

            //(*mdc)(MolecularDataCall::CallForGetData);
            vislib::Array<uint> checkAtomDonor;
            vislib::Array<uint> checkAtomAcceptor;
            checkAtomDonor.SetCount(this->ligand->ligandAtomCnt);
            checkAtomAcceptor.SetCount(this->ligand->ligandAtomCnt);
            for (int i = 0; i < this->ligand->ligandAtomCnt; i++) {
                checkAtomDonor[i] = 0;
                checkAtomAcceptor[i] = 0;
            }
            //********************************* old code from above for protein data
            for (int z = 0; z < this->ligand->connectionCnt; z++) {
                const uint* cons = this->ligand->connections;
                int atom0Index = cons[z * 2 + 0];
                int atom1Index = cons[z * 2 + 1];


                if (atom1Index > atom0Index) {
                    const megamol::protein_calls::MolecularDataCall::AtomType tmpAtomType_0 =
                        this->ligand->Model.atomTypes->PeekElements()[this->ligand->Model.atomTypeIdx[atom0Index]];
                    const megamol::protein_calls::MolecularDataCall::AtomType tmpAtomType_1 =
                        this->ligand->Model.atomTypes->PeekElements()[this->ligand->Model.atomTypeIdx[atom1Index]];
                    /*
                    printf("i:%d=%d \t %d->%d \t %s->%s\n",
                        this->clusterData.assign_modelID_sphereID->PeekElements()[j].GetY(), gblMdlID, atom0Index,
                        atom1Index, tmpAtomType_0.Element().PeekBuffer(), tmpAtomType_1.Element().PeekBuffer());
                    */

                    string check = tmpAtomType_0.Element().PeekBuffer();
                    string check2 = tmpAtomType_1.Element().PeekBuffer();

                    // if hydrogen is present it will be always atom0
                    if (check == "H") {

                    } else if (check2 == "H") {
                        atom0Index = cons[z * 2 + 1];
                        atom1Index = cons[z * 2 + 0];
                    }


                    const megamol::protein_calls::MolecularDataCall::AtomType atom0 =
                        this->ligand->Model.atomTypes->PeekElements()[this->ligand->Model.atomTypeIdx[atom0Index]];
                    const megamol::protein_calls::MolecularDataCall::AtomType atom1 =
                        this->ligand->Model.atomTypes->PeekElements()[this->ligand->Model.atomTypeIdx[atom1Index]];

                    vislib::StringA test0 = atom0.Element();
                    vislib::StringA test1 = atom1.Element();
                    string t0 = test0.PeekBuffer();
                    string t1 = test1.PeekBuffer();
                    float delatEN = std::abs(electroNegativities[t0] - electroNegativities[t1]);

                    check = atom0.Element().PeekBuffer();
                    vislib::math::Vector<uint, 2> checkVec;
                    // check if hydrogen is present
                    if (check == "H") {
                        if (delatEN > 0.5) {
                            // check if already included not neccessary because H can only have one connection
                            H_donorIndices_ligMdl.push_back(atom0Index);
                            if (!checkAtomDonor[atom0Index]) {

                                invGblMdlIDs_mH_to_p.Append({(int)gblMdlID, atom0Index});

                                ligandModelData_H_donors.Append({
                                    this->ligand->Model.atomPostions[atom0Index * 3 + 0],
                                    this->ligand->Model.atomPostions[atom0Index * 3 + 1],
                                    this->ligand->Model.atomPostions[atom0Index * 3 + 2],
                                    0.0f // atom0.Radius()
                                });
                                checkAtomDonor[atom0Index] = exc;
                            }

                            // check if already included
                            if (!checkAtomDonor[atom1Index]) {

                                invGblMdlIDs_m_to_p.Append({(int)gblMdlID, atom1Index});
                                donorIndices_ligMdl.push_back(atom1Index);

                                ligandModelData_donors.Append({
                                    this->ligand->Model.atomPostions[atom1Index * 3 + 0],
                                    this->ligand->Model.atomPostions[atom1Index * 3 + 1],
                                    this->ligand->Model.atomPostions[atom1Index * 3 + 2],
                                    0.0f // atom1.Radius()
                                });
                                checkAtomDonor[atom1Index] = exc;
                            }

                            // check if already included
                            if (!checkAtomAcceptor[atom1Index]) {
                                invGblMdlIDs_p_to_m.Append({(int)gblMdlID, atom1Index});
                                acceptorIndices_ligMdl.push_back(atom1Index);
                                ligandModelData_acceptors.Append({
                                    this->ligand->Model.atomPostions[atom1Index * 3 + 0],
                                    this->ligand->Model.atomPostions[atom1Index * 3 + 1],
                                    this->ligand->Model.atomPostions[atom1Index * 3 + 2],
                                    0.0f // atom1.Radius()
                                });
                                checkAtomAcceptor[atom1Index] = exc;
                            }
                        }
                        // if hydrogen is not present it only can be an accpetor
                    } else {
                        if (delatEN > 0.5) {
                            // who has the charge? push that one with higher electroNegativity
                            if (electroNegativities[t0] > electroNegativities[t1]) {
                                // check if already included
                                if (!checkAtomAcceptor[atom0Index]) {

                                    acceptorIndices_ligMdl.push_back(atom0Index);

                                    invGblMdlIDs_p_to_m.Append({(int)gblMdlID, atom0Index});

                                    ligandModelData_acceptors.Append({
                                        this->ligand->Model.atomPostions[atom0Index * 3 + 0],
                                        this->ligand->Model.atomPostions[atom0Index * 3 + 1],
                                        this->ligand->Model.atomPostions[atom0Index * 3 + 2],
                                        0.0f // atom0.Radius()
                                    });
                                    checkAtomAcceptor[atom0Index] = exc;
                                }
                            } else {
                                // check if already included
                                if (!checkAtomAcceptor[atom1Index]) {
                                    acceptorIndices_ligMdl.push_back(atom1Index);
                                    invGblMdlIDs_p_to_m.Append({(int)gblMdlID, atom1Index});

                                    ligandModelData_acceptors.Append({
                                        this->ligand->Model.atomPostions[atom1Index * 3 + 0],
                                        this->ligand->Model.atomPostions[atom1Index * 3 + 1],
                                        this->ligand->Model.atomPostions[atom1Index * 3 + 2],
                                        0.0f // atom1.Radius()
                                    });
                                    checkAtomAcceptor[atom1Index] = exc;
                                }
                            }
                        }
                    }
                }
            }
            //*********************************
        }
    }
    /*
    for (int i = 0; i < invGblMdlIDs_p_to_m.Count(); i++) {
        printf("ModelEntry:%d \t mdl:%d id:%d\n", i, invGblMdlIDs_p_to_m[i].GetX(), invGblMdlIDs_p_to_m[i].GetY());

    }
    */

    //********* neighbour serach - distance to form a hydrogen bond ********* //
    FFRNNS ffrnns_protTo_LigMdl;
    FFRNNS ffrnns_LigMdlTo_prot;

    ffrnns_protTo_LigMdl.setData(proteinData_donors, distanceToFromHBond);
    if (ligandModelData_acceptors.Count() == 0) return false;
    this->Hbonds_protTo_LigMdl = ffrnns_protTo_LigMdl.getNeigbours(ligandModelData_acceptors);

    ffrnns_protTo_LigMdl.setData(proteinData_acceptors, distanceToFromHBond);
    if (ligandModelData_donors.Count() == 0) return false;
    Hbonds_LigMdlTo_prot = ffrnns_protTo_LigMdl.getNeigbours(ligandModelData_donors);

    // only for Debugging
    /*
    for (int i = 0; i < Hbonds_protTo_LigMdl.size(); i++) {
        printf("ModelEntry:%d\n",i);
        glm::vec3 mdlAcc(ligandModelData_acceptors[i].GetX(), ligandModelData_acceptors[i].GetY(),
            ligandModelData_acceptors[i].GetZ());

        for (int j = 0; j < Hbonds_protTo_LigMdl[i].size(); j++) {
            int index = Hbonds_protTo_LigMdl[i][j];
            glm::vec3 protDon(
                proteinData_donors[index].GetX(), proteinData_donors[index].GetY(), proteinData_donors[index].GetZ());
            float distance = glm::length(mdlAcc-protDon);
            printf("mdlAcc: %d; n:%d protDon: %d \t distance:%f \n",i, j,Hbonds_protTo_LigMdl[i][j],distance);
        }

        printf("\n");
    }
    */

    /*********************************************************************
     ******************** CHECK CORRECT ANGLE CRITERIA *******************
     *********************************************************************/
    float angle = angleToFromHBond * static_cast<float>(vislib::math::PI_DOUBLE / 180.0);


    Hbonds_protTo_LigMdl_checked.resize(ligandModelData_acceptors.Count());
    int cntDeletedHBonds_protTo_LigMdl = 0;
    int cntHBonds_protTo_LigMdl = 0;
    int old_internMdlAtom, old_protDon, old_HprotDon = -1;
    // i = model atomID in Buffer (all lignads all models)
    // value = protein donor atomIDs
    for (int i = 0; i < ligandModelData_acceptors.Count(); i++) {
        auto acceptorPos = ligandModelData_acceptors[i];
        if (Hbonds_protTo_LigMdl[i].size() != 0) {
            // sort(Hbonds_protTo_LigMdl[i].begin(), Hbonds_protTo_LigMdl[i].end());
            for (int j = 0; j < Hbonds_protTo_LigMdl[i].size(); j++) {
                auto hydrogenPos = proteinData_H_donors[Hbonds_protTo_LigMdl[i][j]];
                auto donorPos = proteinData_donors[Hbonds_protTo_LigMdl[i][j]];
                // bool checkAngle = (hydrogenPos - donorPos).Angle(acceptorPos - donorPos) <= angle;
                bool checkAngle = (acceptorPos - donorPos).Angle(hydrogenPos - donorPos) <= angle;
                if (checkAngle) {
                    Hbonds_protTo_LigMdl_checked[i].push_back(Hbonds_protTo_LigMdl[i][j]);
                    // back assignment to models
                    int mdlID = invGblMdlIDs_p_to_m[i].GetX();
                    int internModelAtomIndex = invGblMdlIDs_p_to_m[i].GetY();
                    // avoid doubled-multiple entries; TODO: check reason for mutliple identical entries
                    if (old_internMdlAtom == internModelAtomIndex &&
                        old_protDon == donorIndices_prot[Hbonds_protTo_LigMdl[i][j]] &&
                        old_HprotDon == H_donorIndices_prot[Hbonds_protTo_LigMdl[i][j]]) {
                        continue;
                    }
                    cntHBonds_protTo_LigMdl++;
                    // printf("protDon-HprotDon-mdlAtom:%d - %d - %d\n", donorIndices_prot[Hbonds_protTo_LigMdl[i][j]],
                    //     H_donorIndices_prot[Hbonds_protTo_LigMdl[i][j]], internModelAtomIndex);
                    // push protein atomID
                    out_protTo_LigMdl[mdlID].Append(donorIndices_prot[Hbonds_protTo_LigMdl[i][j]]);
                    out_protTo_LigMdlH[mdlID].Append(H_donorIndices_prot[Hbonds_protTo_LigMdl[i][j]]);
                    // push model atomID
                    out_protTo_LigMdl[mdlID].Append(internModelAtomIndex);
                    out_protTo_LigMdlH[mdlID].Append(internModelAtomIndex);
                    old_internMdlAtom = internModelAtomIndex;
                    old_protDon = donorIndices_prot[Hbonds_protTo_LigMdl[i][j]];
                    old_HprotDon = H_donorIndices_prot[Hbonds_protTo_LigMdl[i][j]];

                } else {
                    cntDeletedHBonds_protTo_LigMdl++;
                }
            }
        }
    }


    Hbonds_LigMdlTo_prot_checked.resize(ligandModelData_donors.Count());
    int cntDeletedHBonds_LigMdlTo_prot = 0;
    int cntHBonds_LigMdlTo_prot = 0;
    int old_protAcc, old_mdlAtom, old_HmdlAtom = -1;
    // i = model atomID in Buffer (all lignads all models)
    // value = protein acceptor atomIDs
    for (int i = 0; i < ligandModelData_donors.Count(); i++) {
        auto hydrogenPos = ligandModelData_H_donors[i];
        auto donorPos = ligandModelData_donors[i];
        // sort(Hbonds_LigMdlTo_prot[i].begin(), Hbonds_LigMdlTo_prot[i].end());
        if (Hbonds_LigMdlTo_prot[i].size() != 0) {
            for (int j = 0; j < Hbonds_LigMdlTo_prot[i].size(); j++) {
                auto acceptorPos = proteinData_acceptors[Hbonds_LigMdlTo_prot[i][j]];
                bool checkAngle = (hydrogenPos - donorPos).Angle(acceptorPos - donorPos) <= angle;
                if (checkAngle) {
                    Hbonds_LigMdlTo_prot_checked[i].push_back(Hbonds_LigMdlTo_prot[i][j]);
                    // back assignment to models
                    int mdlID = invGblMdlIDs_m_to_p[i].GetX();
                    int internModelAtomIndex = invGblMdlIDs_m_to_p[i].GetY();
                    int mdlID_H = invGblMdlIDs_mH_to_p[i].GetX();
                    int internModelAtomIndex_H = invGblMdlIDs_mH_to_p[i].GetY();
                    // avoid doubled-multiple entries; TODO: check reason for mutliple identical entries
                    if (old_protAcc == acceptorIndices_prot[Hbonds_LigMdlTo_prot[i][j]] &&
                        old_mdlAtom == internModelAtomIndex && old_HmdlAtom == internModelAtomIndex_H) {
                        continue;
                    }
                    cntHBonds_LigMdlTo_prot++;
                    // printf("protAcc-mdlDon-HmdlDon:%d - %d - %d\n", acceptorIndices_prot[Hbonds_LigMdlTo_prot[i][j]],
                    //     internModelAtomIndex, internModelAtomIndex_H);
                    // push protein atomID
                    out_ligMdlTo_prot[mdlID].Append(acceptorIndices_prot[Hbonds_LigMdlTo_prot[i][j]]);
                    out_ligMdlTo_protH[mdlID_H].Append(acceptorIndices_prot[Hbonds_LigMdlTo_prot[i][j]]);
                    // push model atomID
                    out_ligMdlTo_prot[mdlID].Append(internModelAtomIndex);
                    out_ligMdlTo_protH[mdlID_H].Append(internModelAtomIndex_H);
                    old_protAcc = acceptorIndices_prot[Hbonds_LigMdlTo_prot[i][j]];
                    old_mdlAtom = internModelAtomIndex;
                    old_HmdlAtom = internModelAtomIndex_H;
                } else {
                    cntDeletedHBonds_LigMdlTo_prot++;
                }
            }
        }
    }
    printf("not valid H-Bonds \t protTo_LigMdl: %.2f%% \t LigMdlTo_prot: %.2f%% \n",
        (float)((float)cntDeletedHBonds_protTo_LigMdl /
                (float)(cntHBonds_protTo_LigMdl + cntDeletedHBonds_protTo_LigMdl) * 100.0f),
        (float)((float)cntDeletedHBonds_LigMdlTo_prot /
                (float)(cntHBonds_LigMdlTo_prot + cntDeletedHBonds_LigMdlTo_prot) * 100.0f));
    printf(
        "valid H-Bonds \t protTo_LigMdl: %d \t LigMdlTo_prot: %d \n", cntHBonds_protTo_LigMdl, cntHBonds_LigMdlTo_prot);
    printf("not valid H-Bonds \t protTo_LigMdl: %d \t LigMdlTo_prot: %d \n", cntDeletedHBonds_protTo_LigMdl,
        cntDeletedHBonds_LigMdlTo_prot);

    hBondAtomCnt_protein.SetCount(pmdc->AtomCount());
#pragma omp parallel for
    for (int i = 0; i < hBondAtomCnt_protein.Count(); i++) {
        hBondAtomCnt_protein[i] = 0;
    }

    // determine the H-Bonds per protein atom
    for (int o = 0; o < ligandModelData_acceptors.Count(); o++) {
        for (int i = 0; i < Hbonds_protTo_LigMdl_checked[o].size(); i++) {
            int z = donorIndices_prot[Hbonds_protTo_LigMdl_checked[o][i]];
            hBondAtomCnt_protein[z]++;
        }
    }

    for (int o = 0; o < ligandModelData_acceptors.Count(); o++) {
        for (int i = 0; i < Hbonds_protTo_LigMdl_checked[o].size(); i++) {
            int z = H_donorIndices_prot[Hbonds_protTo_LigMdl_checked[o][i]];
            hBondAtomCnt_protein[z]++;
        }
    }
    for (int o = 0; o < ligandModelData_donors.Count(); o++) {
        for (int i = 0; i < Hbonds_LigMdlTo_prot_checked[o].size(); i++) {
            int z = acceptorIndices_prot[Hbonds_LigMdlTo_prot_checked[o][i]];
            hBondAtomCnt_protein[z]++;
        }
    }

    return true;
}

bool InteractionForces::loadPLIPdata(Path dataPath, std::vector<MultiPDBQTLoader::PLIPForces>& plipForcesStore) {

    if (this->lmc == NULL) return false;
    Path dataFolderPath = Path(dataPath);
    std::string lingandName;
    plipForcesStore.resize(lmc->getTotalModelCount());
    for (int i = 0; i < lmc->getLigandCount(); i++) {

        // get ligandName
        lmc->SetCurrentLigandAndModel(i, 0);
        (*lmc)(LigandModelCall::CallForGetExtent);
        (*lmc)(LigandModelCall::CallForGetDataSilent);

        for (int j = 0; j < lmc->getModelCount()[lmc->getCurrentLigandID()]; j++) {
            lmc->SetCurrentLigandAndModel(i, j);
            (*lmc)(LigandModelCall::CallForGetData);
            int glbMdlID = lmc->getCurrentGlobalModelID();
            MultiPDBQTLoader::PLIPForces* fc = &plipForcesStore[glbMdlID];
            lingandName = lmc->getCurrentLigandname();


            if (lmc->getClusterAndCentroidData().assignedClusters[0][lmc->getCurrentGlobalModelID()] != -1) {


                Path plipXMLFilePathBase = dataFolderPath;
                plipXMLFilePathBase = plipXMLFilePathBase.append(lingandName);

                Path plipXMLFilePath = plipXMLFilePathBase;
                //****************************
                plipXMLFilePath = plipXMLFilePathBase;
                plipXMLFilePath = plipXMLFilePath.concat("_");
                plipXMLFilePath.concat(std::to_string(j + 1));
                plipXMLFilePath.concat(".inter.plip");

                string line;

                // serach for different models of a ligand

                if (verbose) printf("File: %s\n", plipXMLFilePath.string().c_str());
                //****************************


                ifstream file(plipXMLFilePath.string());
                if (file.is_open()) {

                    while (getline(file, line)) {
                        // lineSplit = ls
                        std::vector<std::string> ls;
                        split(line, ";", ls);
                        if (ls[0] == "hydrophobicInteraction" || ls[0] == "halogenBond" || ls[0] == "HBond") {
                            int protIDx, ligIDx;
                            bool protisdon;
                            glm::vec3 protPos, ligPos;
                            for (int i = 1; i < ls.size(); i++) {
                                // parse data lines
                                // TODO: may divide into functions
                                if (ls[i] == "protIDx") {
                                    protIDx = std::stoi(ls[i + 1]);
                                    i++;
                                } else if (ls[i] == "ligIDx") {
                                    ligIDx = std::stoi(ls[i + 1]);
                                    i++;
                                } else if (ls[i] == "protPos") {
                                    protPos =
                                        glm::vec3(std::stof(ls[i + 1]), std::stof(ls[i + 2]), std::stof(ls[i + 3]));
                                    i += 3;
                                } else if (ls[i] == "ligPos") {
                                    ligPos =
                                        glm::vec3(std::stof(ls[i + 1]), std::stof(ls[i + 2]), std::stof(ls[i + 3]));
                                    i += 3;
                                } else if (ls[i] == "protisdon") {
                                    std::string tmp = ls[i + 1];
                                    protisdon = tmp == "True" ? true : false;
                                    i++;
                                }
                            }
                            // push parsed data to storage
                            if (ls[0] == "hydrophobicInteraction") {
                                fc->hydrophobicInteractions.cnt++;
                                fc->hydrophobicInteractions.protIDx.push_back(protIDx);
                                fc->hydrophobicInteractions.ligIDx.push_back(ligIDx);
                                fc->hydrophobicInteractions.protPos.push_back(protPos);
                                fc->hydrophobicInteractions.ligPos.push_back(ligPos);
                            } else if (ls[0] == "halogenBond") {
                                fc->halogenBonds.cnt++;
                                fc->halogenBonds.protIDx.push_back(protIDx);
                                fc->halogenBonds.ligIDx.push_back(ligIDx);
                                fc->halogenBonds.protPos.push_back(protPos);
                                fc->halogenBonds.ligPos.push_back(ligPos);
                            } else if (ls[0] == "HBond" && !lmc->getInteractionForces().useMegaMolHBonds) {
                                if (protisdon) {
                                    fc->hProtDonors.cnt++;
                                    fc->hProtDonors.protIDx.push_back(protIDx);
                                    fc->hProtDonors.ligIDx.push_back(ligIDx);
                                    fc->hProtDonors.protPos.push_back(protPos);
                                    fc->hProtDonors.ligPos.push_back(ligPos);
                                } else {
                                    fc->hProtAcceptors.cnt++;
                                    fc->hProtAcceptors.protIDx.push_back(protIDx);
                                    fc->hProtAcceptors.ligIDx.push_back(ligIDx);
                                    fc->hProtAcceptors.protPos.push_back(protPos);
                                    fc->hProtAcceptors.ligPos.push_back(ligPos);
                                }
                            }

                        } else if (ls[0] == "saltBridge" || ls[0] == "piStack" || ls[0] == "piCationInteraction") {
                            std::vector<int> protIDxList, ligIDxList;
                            glm::vec3 protPos, ligPos;
                            for (int i = 1; i < ls.size(); i++) {
                                // parse data lines
                                if (ls[i] == "protIDxList") {
                                    while (ls[i + 1] != "ligIDxList") {
                                        protIDxList.push_back(std::stoi(ls[i + 1]));
                                        i++;
                                    }
                                } else if (ls[i] == "ligIDxList") {
                                    while (ls[i + 1] != "protPos") {
                                        ligIDxList.push_back(std::stoi(ls[i + 1]));
                                        i++;
                                    }
                                } else if (ls[i] == "protPos") {
                                    protPos =
                                        glm::vec3(std::stof(ls[i + 1]), std::stof(ls[i + 2]), std::stof(ls[i + 3]));
                                    i += 3;
                                } else if (ls[i] == "ligPos") {
                                    ligPos =
                                        glm::vec3(std::stof(ls[i + 1]), std::stof(ls[i + 2]), std::stof(ls[i + 3]));
                                    i += 3;
                                }
                            }
                            // push the data to storage
                            if (ls[0] == "saltBridge") {
                                fc->saltBridges.cnt++;
                                fc->saltBridges.protIDxList.push_back(protIDxList);
                                fc->saltBridges.ligIDxList.push_back(ligIDxList);
                                fc->saltBridges.protPos.push_back(protPos);
                                fc->saltBridges.ligPos.push_back(ligPos);
                            } else if (ls[0] == "piStack") {
                                fc->piStacks.cnt++;
                                fc->piStacks.protIDxList.push_back(protIDxList);
                                fc->piStacks.ligIDxList.push_back(ligIDxList);
                                fc->piStacks.protPos.push_back(protPos);
                                fc->piStacks.ligPos.push_back(ligPos);
                            } else if (ls[0] == "piCationInteraction") {
                                fc->piCationInteractions.cnt++;
                                fc->piCationInteractions.protIDxList.push_back(protIDxList);
                                fc->piCationInteractions.ligIDxList.push_back(ligIDxList);
                                fc->piCationInteractions.protPos.push_back(protPos);
                                fc->piCationInteractions.ligPos.push_back(ligPos);
                            }

                        } else if (ls[0] == "metalComplex") {
                            int metalIDx = -1;
                            int targetIDx = -1;
                            glm::vec3 metalPos, targetPos;
                            for (int i = 1; i < ls.size(); i++) {
                                // parse data lines
                                if (ls[i] == "metalIDx") {
                                    metalIDx = std::stoi(ls[i + 1]);
                                    i++;
                                } else if (ls[i] == "targetIDx") {
                                    targetIDx = std::stoi(ls[i + 1]);
                                    i++;
                                } else if (ls[i] == "metalPos") {
                                    metalPos =
                                        glm::vec3(std::stof(ls[i + 1]), std::stof(ls[i + 2]), std::stof(ls[i + 3]));
                                    i += 3;
                                } else if (ls[i] == "targetPos") {
                                    targetPos =
                                        glm::vec3(std::stof(ls[i + 1]), std::stof(ls[i + 2]), std::stof(ls[i + 3]));
                                    i += 3;
                                }
                            }
                            if (metalIDx == -1 || targetIDx == -1) {
                                __debugbreak();
                            }
                                // push parsed data to storage if (ls[0] == "hydrophobicInteraction") {
                                fc->metalComplexes.cnt++;
                            fc->metalComplexes.metalIDx.push_back(metalIDx);
                            fc->metalComplexes.targetIDx.push_back(targetIDx);
                            fc->metalComplexes.metalPos.push_back(metalPos);
                            fc->metalComplexes.targetPos.push_back(targetPos);
                        }
                    }
                } else {
                    if (verbose) printf("Unable to open file:%s\n", plipXMLFilePath.string().c_str());
                }
            }
        }
    }

    return true;
};


//////////////////////////////////


/**
 * sort function for a int vec2 by X()
 */
int sortVec2ByX(const vislib::math::Vector<int, 2U>& lhs, const vislib::math::Vector<int, 2U>& rhs) {
    if (lhs.GetX() < rhs.GetX()) {
        return -1;
    } else if (lhs.GetX() > rhs.GetX()) {
        return 1;
    } else {
        return 0;
    }
}