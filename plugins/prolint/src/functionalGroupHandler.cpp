
#include "functionalGroupHandler.h"


using namespace megamol;
using namespace megamol::core;
using namespace megamol::prolint;
using namespace megamol::protein_calls;
using namespace std;



/**
 * collects,filters and groups functional groups of a cluster/pocket
 *
 * @param	clusterID:		the cluster/pocket ID for which the func. groups are to be collected and grouped
 * @param	minPTS:			the min number of a certain func. group to build cluster/hotSpot
 * @param	searchRadius:	(distance for searching neighbouring func. groups of same type around a func. group)
 * @param	FGS_struct:		data structure to store all the results of collecting and grouping func. groups
 * @param	Call->ligand_molecularDataCall:		serves the positional information of ligand atoms
 * @param	Call->ligandModelCall:	serves the clusterData of a pocket and the func. groups of each ligand
 *									(without positional inforamtions)
 * @src		prolint::FunctionalGroupHandler::handle
 */
bool FGS_Handler::handle(int clusterID, int minPTS, float searchRadius, LigandModelCall::FGS_structN& fgsStruct,
    core::Call* ligMDC_, core::Call* lmc_) {

	this->clusterParamsChanged = 0;
    if (this->searchRadius != searchRadius || this->minPTS != minPTS) {
        this->searchRadius = searchRadius;
        this->minPTS = minPTS;
        this->clusterParamsChanged = true;
        this->dataChanged = true;
    } else {
        this->dataChanged = false;
    }

    MolecularDataCall* ligMDC0 = dynamic_cast<MolecularDataCall*>(ligMDC_);
    if (ligMDC0 == NULL) return false;

    LigandModelCall* lmc0 = dynamic_cast<LigandModelCall*>(lmc_);
    if (lmc0 == NULL) return false; 

    this->ligMDC = ligMDC0;
    this->lmc = lmc0;

    if ((fgsStruct.clusterID != clusterID || this->clusterParamsChanged) && clusterID >= -1) {

        // filter and collect(sum up) functional groups to clusters
        locate_functionalGroup_clusters(fgsStruct, clusterID);
        fgsStruct.clusterID = clusterID;
	}

    return true;
};

/**
 *  dertmines min, max count of all groupings of func. groups from FuncGroupData of a given cluster/pocket
 *
 * @param	FGS_struct:		data structure to store all the results of collecting and grouping func. groups
 *			(without positional inforamtions)
 * @return	returns min and max count of func.groups as a vector of two floats
 * @src		prolint::FunctionalGroupHandler::getMinMaxCntFGS
 */
vislib::math::Vector<float, 2> FGS_Handler::getMinMaxCntFGS(LigandModelCall::FGS_structN& fgsClust) {
    vislib::math::Vector<float, 2> minMax;
    minMax.SetX(MAXINT);
    minMax.SetY(MININT);

    if (fgsClust.clusterCentroids.Count() == 0) {
        printf("No FunctionalGroupData for this cluster/pocket to dertmine min max count of func. groups!\n");
        return {0.f, 0.f};
    } else {
        for (int i = 0; i < fgsClust.clusterSizes.size(); i++) {
            minMax.SetY(minMax.GetY() < fgsClust.clusterSizes[i] ? fgsClust.clusterSizes[i] : minMax.GetY());
            minMax.SetX(minMax.GetX() > fgsClust.clusterSizes[i] ? fgsClust.clusterSizes[i] : minMax.GetX());
        }
        return minMax;
    }
};


/**
 *  collects,clusters functional groups of a cluster/pocket
 *
 * @param	clusteredFGS:	data structure of the clustered,filtered functional groups
 * @param	clusterID:		the cluster/pocket ID for which the func. groups are to be collected, clustered
 *
 * @return	bool
 * @src		prolint::FunctionalGroupHandler::locate_functionalGroup_clusters
 */
bool FGS_Handler::locate_functionalGroup_clusters(LigandModelCall::FGS_structN& fgsStruct, int clusterID) {

	int ID = clusterID;
    this->fgsClustCentroids.Resize(0);
    this->fgsTypeID_per_ClustCentroid.resize(0);
    this->clustRes_per_type.clear();

	this->FGSpositions_per_type.clear();
    this->gblMdlIDs_per_type.clear();
    this->centralAtomArrayIDx_per_type.clear();

	MolecularDataCall* ligMDC = dynamic_cast<MolecularDataCall*>(this->ligMDC);
    if (ligMDC == NULL) return false;

    LigandModelCall* lmc = dynamic_cast<LigandModelCall*>(this->lmc);
    if (lmc == NULL) return false; 

    /******************************
     ***** collect FGS data *******
     *****************************/

    // call the cluster atom data
    vector<visFloat4> funcGroupPositions;
    funcGroupPositions.resize(lmc->functionalGroupsIDsToWord.size());

    auto clusterDat = lmc->getClusterAndCentroidData();
    auto clusterIndexSum = clusterDat.clusterIndexSum->PeekElements();

    // get start and end index for iterating over ligands of the cluster (clusterID)
    int clusterIDx_start =
        ID > 0 ? clusterDat.clusterIndexSum->PeekElements()[ID] : clusterDat.clusterIndexSum->PeekElements()[0];
    int clusterIDx_end = clusterDat.clusterIndexSum->PeekElements()[ID + 1];
    int clusterSize = clusterDat.clusterSizes->PeekElements()[ID];


	

    // collect func groups of the cluster/pocket
    // iterating over ligands
    for (int j = clusterIDx_start; j < clusterIDx_end; j++) {

        // set current ligand
        lmc->SetCurrentLigandAndModel_byGlobalMdlID(clusterDat.assign_modelID_sphereID->PeekElements()[j].GetY());
        // get data of current ligand
        (*lmc)(LigandModelCall::CallForGetData);
        (*ligMDC)(MolecularDataCall::CallForGetData);
        auto FGS_ofCurLig = lmc->getFunctionalGroupsOfCurLigand();

		// check for min binding energy
        if (lmc->getCurrentBindingenergy() > lmc->getClusterAndCentroidData().DBSCANParams.minBindEnergy) {
            continue;
        }

        collectHierarchyHighest(FGS_ofCurLig, ligMDC, lmc);
		
	}


	// filtering
    map<uint, visFloat4>::iterator it;
    DBSCAN dbs;
    for (it = FGSpositions_per_type.begin(); it != FGSpositions_per_type.end(); it++) {
        uint typeID = it->first;
        auto res = dbs.DBSCAN_3D(it->second, (float)this->searchRadius, this->minPTS, 0.0f);
        LigandModelCall::FGS_structN conRes = consumeClusterResults(it->second, res, false);

		if (conRes.clusterCnt <= 0) {
            continue;
        }
        if (conRes.clusterCnt > 1 && typeID == 201) {
            int a = 0;
        }
		for (int i = 0; i < conRes.clusterCnt; i++) {
            vislib::math::Vector<float, 4> a(conRes.clusterCentroids[i * 4 + 0], conRes.clusterCentroids[i * 4 + 1],
                conRes.clusterCentroids[i * 4 + 2], 0.0f);
            fgsClustCentroids.Append(a);
            // typeID + internal cluster ID
            vislib::math::Vector<uint, 2> tmp;
            tmp[0] = typeID;
            tmp[1] = i;
            fgsTypeID_per_ClustCentroid.push_back(tmp);
		}
		int cnter = 0;
       
        for (int j = 0; j < conRes.clusterAssignment.size(); j++) {
            if (conRes.clusterAssignment[j] != -1) {
                conRes.fgsMaps_GblMDlIDs[conRes.clusterAssignment[j]][typeID].push_back(gblMdlIDs_per_type[typeID][j]);
                conRes.fgsMaps_centralAtomArrayIDx[conRes.clusterAssignment[j]][typeID].push_back(
                    centralAtomArrayIDx_per_type[typeID][j]); 
            }
        }
        
        clustRes_per_type[typeID] = conRes;
    }
	
	int dsa = 5;
	
    DBSCAN dbscan;
    float minValue = 0.0f;

    vislib::Array<int> clustRes;
    // cluster it
    this->searchRadius;

    if (fgsClustCentroids.Count() == 0) {
        fgsStruct.clusterCnt = -1;
        return true;
	}
    clustRes = dbscan.DBSCAN_3D(fgsClustCentroids, 0.5f, 1, minValue);

    // consume the cluster results
    LigandModelCall::FGS_structN consumedClusterRes = consumeClusterResults(fgsClustCentroids, clustRes, true);
    for (int i = 0; i < consumedClusterRes.clusterCnt; i++) {
        consumedClusterRes.clusterSizes[i] = 0;
	}

	// set fgs type IDs and regarding gblMdlIDs
    for (int i = 0; i < consumedClusterRes.clusterAssignment.size(); i++) {
        int curClust = consumedClusterRes.clusterAssignment[i];
        if (curClust != -1) {
            uint curTypeID = fgsTypeID_per_ClustCentroid[i].GetX();
            uint curInternClustID = fgsTypeID_per_ClustCentroid[i].GetY();
            auto curFGSmap_gblMdl = clustRes_per_type[curTypeID].fgsMaps_GblMDlIDs[curInternClustID];
            auto curFGSmap_centAtomIDx = clustRes_per_type[curTypeID].fgsMaps_centralAtomArrayIDx[curInternClustID];

            consumedClusterRes.fgsMaps_GblMDlIDs[curClust][curTypeID] = curFGSmap_gblMdl[curTypeID];
            consumedClusterRes.fgsMaps_centralAtomArrayIDx[curClust][curTypeID] = curFGSmap_centAtomIDx[curTypeID];
            consumedClusterRes.clusterSizes[curClust] += curFGSmap_centAtomIDx[curTypeID].size();
		}
	}


    // if a cluster was found set the clusterResults
    if (consumedClusterRes.clusterCnt > -1) {
        fgsStruct = consumedClusterRes;
    } else {
        fgsStruct.clusterCnt = -1;
    }
       
    printf("%d funcGroups clusters for pocketID: %d \n", fgsStruct.clusterCnt, ID);
    return true;
};


bool FGS_Handler::collectHierarchyHighest(
    LigandModelCall::FunctionalGroupsLigand* FGS_ofCurLig, MolecularDataCall* ligMDC, LigandModelCall* lmc) {

    // parse fgsHierarchy
    // reverse order (start at basic leaves)
    for (int i = FGS_ofCurLig->totalGroupCnt - 1; i >= 0; i--) {
        auto f = FGS_ofCurLig->fGroups[i];

        for (int j = FGS_ofCurLig->fGroups[i].speficGroupCnt - 1; j >= 0; j--) {
            const float* fAtom = &ligMDC->AtomPositions()[(f.centralAtoms[j] - 1) * 3 + 0];
            // check if it is root in hierarchy
            if (lmc->fgsHierarchy[f.groupID] == -1) {

                // push data
                vislib::math::Vector<float, 4> atomPos(fAtom[0], fAtom[1], fAtom[2], 0.0f);
                this->FGSpositions_per_type[f.groupID].Append(atomPos);
                this->gblMdlIDs_per_type[f.groupID].push_back(lmc->getCurrentGlobalModelID());
                this->centralAtomArrayIDx_per_type[f.groupID].push_back(j);
            }
        }
    }


    return true;
}

/**
 *  generic function to process the results of the DBSCAN based clustering
 *
 * @param	points:		the dataPoints which were used for clustering
 * @param	clusterAssignment:		the assignment of which dataPoint belongs to which cluster
 *
 * @return	returns ClusteringStats_struct (clusterSizes, clusterAssignment, clusterCentroids, clusterCuboids)
 * @src		prolint::FunctionalGroupHandler::locate_functionalGroup_clusters
 */
LigandModelCall::FGS_structN FGS_Handler::consumeClusterResults(
    visFloat4& points, vislib::Array<int>& clusterAssignment, bool setFGSmap) {
    LigandModelCall::FGS_structN fgs;
    if (points.Count() == 0) {
        return fgs;
	}

    // copy cluster assignment;
    fgs.clusterAssignment.resize(clusterAssignment.Count());
    std::memcpy(&fgs.clusterAssignment[0], clusterAssignment.PeekElements(), clusterAssignment.Count() * sizeof(int));

    // get number of clusters
    int clusterCnt = -1;
    for (int i = 0; i < clusterAssignment.Count(); i++) {
        if (clusterAssignment[i] >= 0) {
            clusterCnt = clusterCnt < clusterAssignment[i] ? clusterAssignment[i] : clusterCnt;
        }
    }

    // set number of clusters
    fgs.clusterCnt = clusterCnt;

    if (clusterCnt != -1) {
        clusterCnt++;

        fgs.clusterCnt = clusterCnt;
        fgs.clusterCuboids.resize(clusterCnt);
        fgs.clusterCentroids.Resize(clusterCnt*4);
        fgs.clusterCentroids.SetCount(clusterCnt * 4);
        fgs.clusterSizes.resize(clusterCnt);
        fgs.fgsMaps_GblMDlIDs.resize(clusterCnt);
        fgs.fgsMaps_centralAtomArrayIDx.resize(clusterCnt);

        // initialize clusterSizes with zero
        for (int i = 0; i < clusterCnt; i++) {
            fgs.clusterSizes[i] = 0;
        }

        // create the cluster cuboids
        for (int i = 0; i < clusterAssignment.Count(); i++) {
            if (clusterAssignment[i] > -1) {
                // initilaize cuboid
                if (fgs.clusterCuboids[clusterAssignment[i]].IsEmpty()) {
                    fgs.clusterCuboids[clusterAssignment[i]].Set(points[i].GetX(), points[i].GetY(), points[i].GetZ(),
                        points[i].GetX(), points[i].GetY(), points[i].GetZ());
                    // extend cuboid
                } else {
                    fgs.clusterCuboids[clusterAssignment[i]].GrowToPoint(
                        points[i].GetX(), points[i].GetY(), points[i].GetZ());
                }
                fgs.clusterSizes[clusterAssignment[i]]++;
            }
        }

        // calc the cluster centroids
        for (int i = 0; i < fgs.clusterCuboids.size(); i++) {
            vislib::math::Vector<float, 3> center = fgs.clusterCuboids[i].CalcCenter();
            fgs.clusterCentroids[i * 4 + 0] = center.GetX();
            fgs.clusterCentroids[i * 4 + 1] = center.GetY();
            fgs.clusterCentroids[i * 4 + 2] = center.GetZ();
            fgs.clusterCentroids[i * 4 + 3] = 0.4f;
        } 
    }
    return fgs;
};

