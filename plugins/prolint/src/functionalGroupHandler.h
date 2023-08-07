#include "protein_calls/MolecularDataCall.h"
#include "MultiPDBQTLoader.h"
#include "vislib/math/Vector.h"
#include "vislib/Array.h"
#include <map>
#include "LigandModelCall.h"

// typedefs
typedef unsigned int uint;

// namespaces
namespace megamol {
namespace prolint {




	



/*******************************************/
/***** FUNCTINAL GROOUP HANLDER CLASS ******/
/*******************************************/

class FGS_Handler {

public:
    FGS_Handler(){};
    ~FGS_Handler(){};


    // public varibales:

    // shows if new func. group data was calculated
    bool dataChanged = false;

	visFloat4 fgsClustCentroids;
	// fgsType and internal ClusterID per cluster centroid
    std::vector<vislib::math::Vector<uint,2>> fgsTypeID_per_ClustCentroid;

    std::map<uint, LigandModelCall::FGS_structN> clustRes_per_type;


    // public functions:

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
    bool handle(int clusterID, int minPTS, float searchRadius, LigandModelCall::FGS_structN& fgsClust,
        core::Call* ligMDC_, core::Call* lmc_);


    /**
     *  dertmines min, max count of all groupings of func. groups from FuncGroupData of a given cluster/pocket
     *
     * @param	FGS_struct:		data structure to store all the results of collecting and grouping func. groups
     *			(without positional inforamtions)
     * @return	returns min and max count of func.groups as a vector of two floats
     * @src		prolint::FunctionalGroupHandler::getMinMaxCntFGS
     */
    vislib::math::Vector<float, 2> getMinMaxCntFGS(LigandModelCall::FGS_structN& fgsClust);


private:
    // private variables:

    // call pointers
    core::Call *ligMDC;
    core::Call *lmc;


    // DBSCAN params (for clustering)
    bool clusterParamsChanged;
    int minPTS;
    float searchRadius;

    


    // private functions:

    /**
     *  collects,clusters functional groups of a cluster/pocket
     *
     * @param	clusteredFGS:	data structure the clustered,filtered functional groups
     * @param	clusterID:		the cluster/pocket ID for which the func. groups are to be collected, clustered
     *
     * @return	returns min and max count of func.groups as a vector of two floats
     * @src		prolint::FunctionalGroupHandler::locate_functionalGroup_clusters
     */
    bool locate_functionalGroup_clusters(LigandModelCall::FGS_structN& fgsStruct, int clusterID);

	map<uint, visFloat4> FGSpositions_per_type;
    /* map:
	*	key: typeID <vector>
			[0]: gblMdlIDs <int vector>
			[1]: centralAtomIDs (which of may multiple central atoms of this fgs type) <int vector>
	*/
    map<uint,std::vector<uint>> gblMdlIDs_per_type;
    map<uint,std::vector<uint>> centralAtomArrayIDx_per_type;

	bool collectHierarchyHighest(LigandModelCall::FunctionalGroupsLigand* FGS_ofCurLig,
        megamol::protein_calls::MolecularDataCall* ligMDC, megamol::prolint::LigandModelCall* lmc);


    /**
     *  generic function to process the results of the DBSCAN based clustering
     *
     * @param	points:		the dataPoints which were used for clustering
     * @param	clusterAssignment:		the assignment of which dataPoint belongs to which cluster
     *
     * @return	returns ClusteringStats_struct (clusterSizes, clusterAssignment, clusterCentroids, clusterCuboids)
     * @src		prolint::FunctionalGroupHandler::locate_functionalGroup_clusters
     */
    LigandModelCall::FGS_structN consumeClusterResults(
        visFloat4& points, vislib::Array<int>& clusterAssignment, bool setFGSmap);

};







} // namespace prolint
} // namespace megamol
