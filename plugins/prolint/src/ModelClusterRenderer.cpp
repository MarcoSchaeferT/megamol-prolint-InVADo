/*
 * ModelClusterRenderer.cpp
 *
 * Copyright (C) 2019 by University of Tuebingen (BDVA).
 * All rights reserved.
 */


#include "stdafx.h"
#include "ModelClusterRenderer.h"
#include <algorithm>
#include <cuda_gl_interop.h>
#include <direct.h>
#include <omp.h>
#include <regex>
#include "FFRNNS.h"
#include "MultiPDBQTLoader.h"
#include "mmcore/CoreInstance.h"
#include "mmcore/utility/ShaderSourceFactory.h"
#include "vislib/graphics/gl/CameraOpenGL.h"
#include "vislib/math/Matrix.h"
#include "vislib/math/mathfunctions.h"
#define NUM_THREADS 256


#define CHECK_FOR_OGL_ERROR()                                                                                          \
    do {                                                                                                               \
        GLenum err;                                                                                                    \
        err = glGetError();                                                                                            \
        if (err != GL_NO_ERROR) {                                                                                      \
            fprintf(stderr, "%s(%d) glError: %s\n", __FILE__, __LINE__, gluErrorString(err));                          \
        }                                                                                                              \
    } while (0)

#ifdef _DEBUG
#    define SVGs_onDemand true
#else
#    define SVGs_preLoad true
#endif

// using namespace std;
using namespace megamol;
using namespace megamol::core;
using namespace megamol::prolint;
using namespace megamol::protein_calls;
using namespace megamol::geocalls;

using namespace glm;

// Enable High Performance Graphics while using Integrated Graphics. 
__declspec(dllexport) DWORD NvOptimusEnablement = 0x00000001; // Nvidia 

/**
 * Constructor
 */
ModelClusterRenderer::ModelClusterRenderer(void)
    : Renderer3DModule()
    , ligandModel_callerSlot("getData", "Call for a specific model of a specific ligand")
    , ligandMolecularData_callerSlot(
          "getMolecularData", "in this case the atom data of a specific model of a specific ligand")
    , triMeshData_callerSlot("getTriMeshData ", "calling for triangle mesh data of the protein")
    , proteinMolecularData_callerSlot("getProteinMolecularData ", "in this case the atom data of the protein")
    , volumetricData_callerSlot("getAOvolumeData", "call for density volume of protein data")
    , transferFunction_interactionForce_callerSlot("TransferFunction_interactionForce",
          "(for interactionForce intensity on protein surface) get the color "
          "transfer function from a serialized JSON string.")
    , transferFunction_clusterSpheres_callerSlot("TransferFunction_clusterSpheres",
          "(to encode the amount of models in a cluster sphere) get the color "
          "transfer function from a serialized JSON string.")
    , transferFunction_functGroups_callerSlot("TransferFunction_functGroups",
          "(to encode the number of functional groups of centroid spheres)get the "
          "color transfer function from a serialized JSON string.")
    , dataSlot_MousePos("getMouseData", "gets mouse position and mouse event flags")
    , dataOutSlot_MDC_protein("MolecularDataCall_protein",
          "The slot provides protein resiudes which are involved in protein<->ligand interactions")
    , dataOutSlot_MDC_ligand("MolecularDataCall_ligand",
          "The slot provides ligand pose(s)/model(s) which are request from MultiPDBQTLoader")
    , getClipPlane_callerSlot("getClipPlaneSlot", "Connects to a clipping plane module")
    , starMenuTextSizeParam("star Menu::MenuTextSize", "Textsize of star menu")
    , interactionForceParam("InterForces->  Surface Coloring::interactionForce (on/off)",
          "Enable/Disable protein surface coloring based on chosen interaction force")
    , interactionForceAmountColoringParam("InterForces->  Surface Coloring::interactionForceAmountCol (on/off)",
          "Enable/Disable protein surface coloring based on chosen interaction force")
    , surfaceColoring_dropdownParam("InterForces->  Surface Coloring::Force",
          "drop-down for interaction force selection (to color the protein surface) ")
    , pocketAreaParam("Protein Coloring::pocketArea", "show protein pockets surface")
    , opaquePocketParam("Protein Coloring::opaquePocket", "disables transparency for protein pockets surface")
    , transparencyParam(
          "Protein Coloring::transparencyRenderning", "Enable/Disable transparent renderning of protein surface")
    , chemicalPropertyParam("Protein Coloring::chemicalProperty",
          "Enable/Disable protein surface coloring based on a chemical property choosen in MSMS-meshLoader")
    , starMenuCntParam("star Menu::starMenuCircleCnt", "number of elements in the star menu")
    , starMenuSizeParam("star Menu::starMenuSize", "size of the star menu")
    , starMenuSorting_dropdownParam("star Menu::Sorting", "interaction force used for sorting")
    //, starMenuInnerStartParam("star Menu::starMenuInnerStart", "controls where the inner star lines begin")
    //, starMenuInnerStart(7)
    , subSphereSizefactorParam("star Menu::minCircleSizeParam", "controls the size of sub circles")
    , subCircleSizefactor(0.14f)
    , curSelectedStarMenuModelIDx(-1)
    , subCircleCnt(3)
    , clusterNoiseParam("Render Data:: show noise (not clustered models)",
          "enable/disable showing noise centroids - models which not fullfil the DBSCAN params")
    , clusterSpheresParam(
          "Render Data:: show cluster spheres", "enable/disable showing cluster centroids - clusters found with DBSCAN")
    , residuesParam("Render Data:: protein residues (on/off)", "show protein residues")
    , residuesBackboneParam("Render Data:: protein backbone (on/off)", "show protein backbone")
    , residuesLabelsParam("Render Data:: protein residue labels", "show residues labels")
    , interactionResiduesParam("Render Data:: protein residues (interaction only)",
          "show residues which involved in chosen type of protein-ligand interaction")
    , clusterAtomsParam(
          "Render Data:: show all clusterAtoms", "show all atoms that belongs to the current selected cluster")
    , funcGroupCentroidsParam("Render Data:: show functional Groups",
          "show all cluster centroids of clustered functional groups (see params of \"functional group clustering\") ")
    , improveTextReadabilityParam("Render Data:: improve text readability", " draws an outline around the text")
    , barchartGylpheParam("Render Data:: bar chart glyphs (on/off)",
          " draws bar chart glyphs for each pocket (top: contacted residues counts; bot: surface areas)")
    , aoSampFactParam(
          "ambient Occlusion::aoSampleFactor", "shadowing amount factor multiplied with the ambient occlusion factor")
    , aoSampFact(1.0f)
    , aoSampDistParam("ambient Occlusion::aoSampleDistance", "the access step length in voxels (Ambient Occlusion")
    , aoSampDist(0.9f)
    , aoIntensityMinParam("ambient Occlusion::aoIntensityMin", "lowest value for interpolation of AO-Factor")
    , aoIntensityMin(0.1f)
    , aoIntensityMaxParam("ambient Occlusion::aoIntensityMax", "highest value for interpolation of AO-Factor")
    , aoIntensityMax(1.3f)
    , searchRadius_fgcParam("Clustering->  functional groups::searchRadius",
          "sets the searchRadius in angstrom used for clustering (DBSCAN)")
    , searchRadius_fgc(3.f)
    , minPTS_fgcParam("Clustering->  functional groups::minPTS",
          "sets the minimum of points/groups that have to be present of the "
          "same functional group to form a cluster (DBSCAN)")
    , withoutSocketParam("WebServer::WithoutWebServer", "selects if the connection to the web server is started or not")
    , hydrogenBondsParam(
          "InterForces-> Sticks ::hydrogenBonds", "show hydrogen bonds - red protein donor - green ligand donor")
    , hydrogenConesParam("InterForces-> Sticks ::hydrogenCones",
          "show hydrogen bonds as cone - red protein donor - green ligand donor")
    , halogenBondsParam("InterForces-> Sticks ::halogenBonds", "halogenBonds")
    , hydrophobicInteractionsParam("InterForces-> Sticks ::hydrophobicInteractions", "hydrophobicInteractions")
    , metalComplexesParam("InterForces-> Sticks ::metalComplexes", "metalComplexes")
    , piCationInteractionsParam("InterForces-> Sticks ::piCationInteractions", "piCationInteractions")
    , piStacksParam("InterForces-> Sticks ::piStacks", "piStacks")
    , saltBridgesParam("InterForces-> Sticks ::saltBridges", "saltBridges")
    , minPTS_fgc(5)
    , mouseX(0)
    , mouseY(0)
    , noiseModel_CentroidRadius(0.5f)
    , cluster_CentroidRadius(2.f)
    , higlight_cluster_CentroidRadius(2.f)
    , higlight_noiseModel_CentroidRadius(0.5)
    , searchRadius_InteractionResidues(3.0f)
    , searchRadius_PocketPatch(1.0f)
    , old_tmc_dataHash(-1)
    , old_pmdc_DataHash(-1)
    , old_vdc_DataHash(-1)
    , old_clipPlaneOrientation(-1)
    , old_clipPlaneDist(-1)
    , SVGcnt(0)
    , arrowUpDown(-1)
    , lmdcOut_newData(false)
    , sphereShader()
    , firstRenderCall(true)
    , scaleWorld(1.0f)
    , font(utility::SDFFont::ROBOTO_SANS, utility::SDFFont::RENDERTYPE_FILL)
    , withoutSocket(true) {


    this->starMenu.textSize = 0.05f;
    this->starMenu.circleCnt = 8;
    this->starMenu.size = 13;
    this->starMenu.maxCircleCnt = 15;
    this->starMenu.hoverIDx = -1;
    this->starMenu.selectedIDx = -1;
    this->starMenu.old_selectedIDx = -1;
    this->starMenu.fontColor = {253.f, 180.f, 98.f, 255.f};
    this->starMenu.curPageStartIDx = 0;
    this->starMenu.curPage = 0;
    this->starMenu.curSortingFroceType = this->interForceTypes.bestScores;


    /**************************************************
     ************ SET DEFAULT PARAM VALUES ************
     **************************************************/

    // set interaction force colors
    this->forceColorsMap[this->interForceTypes.HBonds] =
        glm::vec4(0.407f, 0.196f, 0.658f, 1.0f); // rgb( 103,785	49,98	167,79
    this->forceColorsMap[this->interForceTypes.halogenBonds] =
        glm::vec4(0.216f, 0.494f, 0.722f, 1.0f); // rgb( 55,08	125,97	184,11
    this->forceColorsMap[this->interForceTypes.hydrophobicInteractions] =
        glm::vec4(0.302f, 0.686f, 0.29f, 1.0f); // rgb( 77,01	174,93	73,95
    this->forceColorsMap[this->interForceTypes.metalComplexes] =
        glm::vec4(0.650f, 0.337f, 0.639f, 0.157f); // rgb(166, 86, 40)
    this->forceColorsMap[this->interForceTypes.piCationInteractions] =
        glm::vec4(1.0f, 0.498f, 0.0f, 1.0f);                                                       // rgb( 255	126,99	0
    this->forceColorsMap[this->interForceTypes.piStacks] = glm::vec4(0.96f, 0.867f, 0.098f, 1.0f); // rgb( 245	221	25
    this->forceColorsMap[this->interForceTypes.saltBridges] =
        glm::vec4(0.89f, 0.102f, 0.11f, 1.0f); // rgb( 226,95	26,01	28,05

    this->forceColorsMap[this->interForceTypes.Hbonds_protDonor] =
        glm::vec4(118.f / 255.f, 42.f / 255.f, 131.f / 255.f, 1.0f); // rgb(118, 42, 131)
    this->forceColorsMap[this->interForceTypes.Hbonds_protAcceptor] =
        glm::vec4(194.f / 255.f, 165.f / 255.f, 207.f / 255.f, 1.0f); // rgb(194, 165, 207)
                                                                      // 244,165,130
                                                                      // 146,197,222

    // set color for rendering of all cluster atoms
    auto col = &this->forceColorsMap[this->interForceTypes.piCationInteractions];
    this->renVarsACAS.baseColor.Set(col->x, col->y, col->z, col->w);

    // control protein mesh coloring
    this->colorControl.interactionForce = 1;
    this->colorControl.interactionForceAmountColoring = 1;
    this->colorControl.pocketArea = 1;
    this->colorControl.opaquePocket = 0;
    this->colorControl.transparency = 0;
    this->colorControl.chemicalProperty = 0;
    this->colorControl.surfaceColoring_dropdown = 0;

    // control render data
    this->renderDataControl.clusterNoise = 0;
    this->renderDataControl.clusterSpheres = 1;
    this->renderDataControl.clusterAtoms = 0;
    this->renderDataControl.funcGroupCentroids = 1;
    this->renderDataControl.interactionResidues = 1;
    this->renderDataControl.residues = 1;
    this->renderDataControl.residuesBackbone = 0;
    this->renderDataControl.residueLabels = 1;
    this->renderDataControl.improveTextReadability = 1;
    this->renderDataControl.barchartGylphe = 0;

    // control render data => interactionForces
    this->renderDataControl.intForces.hydrogenBonds = 0;
    this->renderDataControl.intForces.hydrogenCones = 1;
    this->renderDataControl.intForces.halogenBonds = 0;
    this->renderDataControl.intForces.hydrophobicInteractions = 0;
    this->renderDataControl.intForces.metalComplexes = 0;
    this->renderDataControl.intForces.piCationInteractions = 0;
    this->renderDataControl.intForces.piStacks = 0;
    this->renderDataControl.intForces.saltBridges = 1;

    // indices for interaction forces
    this->subCircleContentIDxes["orderNum"] = 0;
    this->subCircleContentIDxes["score"] = 1;
    this->subCircleContentIDxes[this->interForceTypes.HBonds] = 2;
    this->subCircleContentIDxes[this->interForceTypes.halogenBonds] = 3;
    this->subCircleContentIDxes[this->interForceTypes.hydrophobicInteractions] = 4;
    this->subCircleContentIDxes[this->interForceTypes.metalComplexes] = 5;
    this->subCircleContentIDxes[this->interForceTypes.piCationInteractions] = 6;
    this->subCircleContentIDxes[this->interForceTypes.piStacks] = 7;
    this->subCircleContentIDxes[this->interForceTypes.saltBridges] = 8;
    this->subCircleContentIDxes[this->interForceTypes.bestScores] = 9;

    this->interForceTypes.forceTypeVec.resize(this->subCircleContentIDxes.size(), "");
    this->interForceTypes.forceTypeVec[this->subCircleContentIDxes[this->interForceTypes.HBonds]] =
        this->interForceTypes.HBonds;
    this->interForceTypes.forceTypeVec[this->subCircleContentIDxes[this->interForceTypes.halogenBonds]] =
        this->interForceTypes.halogenBonds;
    this->interForceTypes.forceTypeVec[this->subCircleContentIDxes[this->interForceTypes.hydrophobicInteractions]] =
        this->interForceTypes.hydrophobicInteractions;
    this->interForceTypes.forceTypeVec[this->subCircleContentIDxes[this->interForceTypes.metalComplexes]] =
        this->interForceTypes.metalComplexes;
    this->interForceTypes.forceTypeVec[this->subCircleContentIDxes[this->interForceTypes.piCationInteractions]] =
        this->interForceTypes.piCationInteractions;
    this->interForceTypes.forceTypeVec[this->subCircleContentIDxes[this->interForceTypes.piStacks]] =
        this->interForceTypes.piStacks;
    this->interForceTypes.forceTypeVec[this->subCircleContentIDxes[this->interForceTypes.saltBridges]] =
        this->interForceTypes.saltBridges;
    this->interForceTypes.forceTypeVec[this->subCircleContentIDxes[this->interForceTypes.bestScores]] =
        this->interForceTypes.bestScores;

    // pointer to subCircle text colors
    this->subCircleColorPointers.resize(this->subCircleContentIDxes.size());
    this->subCircleColorPointers[this->subCircleContentIDxes["orderNum"]] = &this->starMenu.fontColor.x;
    this->subCircleColorPointers[this->subCircleContentIDxes["score"]] = &this->starMenu.fontColor.x;
    this->subCircleColorPointers[this->subCircleContentIDxes[this->interForceTypes.HBonds]] =
        &this->forceColorsMap[this->interForceTypes.HBonds].x;
    this->subCircleColorPointers[this->subCircleContentIDxes[this->interForceTypes.halogenBonds]] =
        &this->forceColorsMap[this->interForceTypes.halogenBonds].x;
    this->subCircleColorPointers[this->subCircleContentIDxes[this->interForceTypes.hydrophobicInteractions]] =
        &this->forceColorsMap[this->interForceTypes.hydrophobicInteractions].x;
    this->subCircleColorPointers[this->subCircleContentIDxes[this->interForceTypes.metalComplexes]] =
        &this->forceColorsMap[this->interForceTypes.metalComplexes].x;
    this->subCircleColorPointers[this->subCircleContentIDxes[this->interForceTypes.piCationInteractions]] =
        &this->forceColorsMap[this->interForceTypes.piCationInteractions].x;
    this->subCircleColorPointers[this->subCircleContentIDxes[this->interForceTypes.piCationInteractions]] =
        &this->forceColorsMap[this->interForceTypes.piCationInteractions].x;
    this->subCircleColorPointers[this->subCircleContentIDxes[this->interForceTypes.piStacks]] =
        &this->forceColorsMap[this->interForceTypes.piStacks].x;
    this->subCircleColorPointers[this->subCircleContentIDxes[this->interForceTypes.saltBridges]] =
        &this->forceColorsMap[this->interForceTypes.saltBridges].x;

    // control hoverings and selections
    this->s.selectedGlobalMdlID = -1;
    this->s.old_selectedGlobalMdlID = -1;
    this->s.selectedIDx = -1;
    this->s.hoveringIDx = -1;
    this->s.old_SelectedIDx = -1;
    this->s.old_HoveringIDx = -1;
    this->s.curClusterID = -1;
    this->s.old_curClusterID = -2;
    this->s.changedcurClusterID = 0;
    this->s.selectedFgsClusterID = -1;
    this->s.changedCurForceType = 1;
    this->s.curForceType_surfaceColoring = this->interForceTypes.HBonds;


    // set extern param names
    this->paramName_proteinColoring0 = "::inst::MSMSMeshLoader1::color::coloringMode0";
    this->paramName_proteinColoring1 = "::inst::MSMSMeshLoader1::color::coloringMode1";
    this->paramName_modelClustering_eps = "::inst::MultiPDBQTLoader1::DBSCAN::eps";
    this->paramName_modelClustering_minBindEnergy = "::inst::MultiPDBQTLoader1::DBSCAN::minBindEnergy";
    this->paramName_modelClustering_minPts = "::inst::MultiPDBQTLoader1::DBSCAN::minPts";
    this->paramName_clipPlaneEnable = "::inst::ClipPlane1::enable";
    this->paramName_clipPlaneNormal = "::inst::ClipPlane1::normal";
    this->paramName_clipPlaneDist = "::inst::ClipPlane1::dist";
    this->paramName_cartoonRnederer = "::inst::Mux4Renderer3D1::renderer2active";
    this->planeOrientations["XY"] = glm::vec3(0.0, 0.0, 1.0);
    this->planeOrientations["XZ"] = glm::vec3(0.0, 1.0, 0.0);
    this->planeOrientations["YZ"] = glm::vec3(1.0, 0.0, 0.0);


    // set up caller slots
    this->ligandMolecularData_callerSlot.SetCompatibleCall<MolecularDataCallDescription>();
    this->MakeSlotAvailable(&this->ligandMolecularData_callerSlot);

    this->ligandModel_callerSlot.SetCompatibleCall<LigandModelCallDescription>();
    this->MakeSlotAvailable(&this->ligandModel_callerSlot);

    this->triMeshData_callerSlot.SetCompatibleCall<CallTriMeshDataDescription>();
    this->MakeSlotAvailable(&this->triMeshData_callerSlot);

    this->proteinMolecularData_callerSlot.SetCompatibleCall<MolecularDataCallDescription>();
    this->MakeSlotAvailable(&this->proteinMolecularData_callerSlot);

    this->volumetricData_callerSlot.SetCompatibleCall<misc::VolumetricDataCallDescription>();
    this->MakeSlotAvailable(&this->volumetricData_callerSlot);

    // set up caller slots TRANSFER FUNCTIONS
    this->transferFunction_interactionForce_callerSlot
        .SetCompatibleCall<core::view::CallGetTransferFunctionDescription>();
    this->MakeSlotAvailable(&this->transferFunction_interactionForce_callerSlot);

    this->transferFunction_clusterSpheres_callerSlot
        .SetCompatibleCall<core::view::CallGetTransferFunctionDescription>();
    this->MakeSlotAvailable(&this->transferFunction_clusterSpheres_callerSlot);

    this->transferFunction_functGroups_callerSlot.SetCompatibleCall<core::view::CallGetTransferFunctionDescription>();
    this->MakeSlotAvailable(&this->transferFunction_functGroups_callerSlot);

    // set up caller slot CLIPPING PLANE
    this->getClipPlane_callerSlot.SetCompatibleCall<core::view::CallClipPlaneDescription>();
    this->MakeSlotAvailable(&this->getClipPlane_callerSlot);

    // set up callee slots
    // Mouse Pos
    this->dataSlot_MousePos.SetCallback(prolint::MouseInputCall::ClassName(), prolint::MouseInputCall::FunctionName(0),
        &ModelClusterRenderer::getMousePos);
    this->MakeSlotAvailable(&this->dataSlot_MousePos);

    // MDC_protein
    this->dataOutSlot_MDC_protein.SetCallback(MolecularDataCall::ClassName(),
        MolecularDataCall::FunctionName(MolecularDataCall::CallForGetData), &ModelClusterRenderer::getDataPMDC);

    this->dataOutSlot_MDC_protein.SetCallback(MolecularDataCall::ClassName(),
        MolecularDataCall::FunctionName(MolecularDataCall::CallForGetExtent), &ModelClusterRenderer::getExtentPMDC);
    this->MakeSlotAvailable(&this->dataOutSlot_MDC_protein);

    // MDC_ligand
    this->dataOutSlot_MDC_ligand.SetCallback(MolecularDataCall::ClassName(),
        MolecularDataCall::FunctionName(MolecularDataCall::CallForGetData), &ModelClusterRenderer::getDataLMDC);

    this->dataOutSlot_MDC_ligand.SetCallback(MolecularDataCall::ClassName(),
        MolecularDataCall::FunctionName(MolecularDataCall::CallForGetExtent), &ModelClusterRenderer::getExtentLMDC);
    this->MakeSlotAvailable(&this->dataOutSlot_MDC_ligand);


    // set up FONT RENDERING parameters
    this->starMenuTextSizeParam << new param::IntParam(this->starMenu.textSize * 100);
    this->starMenuTextSizeParam.SetUpdateCallback(&ModelClusterRenderer::paramsStarMenuChanged);
    this->MakeSlotAvailable(&this->starMenuTextSizeParam);

    this->starMenuCntParam << new param::IntParam(this->starMenu.circleCnt);
    this->starMenuCntParam.SetUpdateCallback(&ModelClusterRenderer::paramsStarMenuChanged);
    this->MakeSlotAvailable(&this->starMenuCntParam);

    this->starMenuSizeParam << new param::IntParam(this->starMenu.size);
    this->starMenuSizeParam.SetUpdateCallback(&ModelClusterRenderer::paramsStarMenuChanged);
    this->MakeSlotAvailable(&this->starMenuSizeParam);

    /*
    this->starMenuInnerStartParam << new param::IntParam(this->starMenuInnerStart);
    this->starMenuInnerStartParam.SetUpdateCallback(&ModelClusterRenderer::paramsStarMenuChanged);
    this->MakeSlotAvailable(&this->starMenuInnerStartParam);
    */

    this->subSphereSizefactorParam << new param::IntParam(this->subCircleSizefactor * 100);
    this->subSphereSizefactorParam.SetUpdateCallback(&ModelClusterRenderer::paramsStarMenuChanged);
    this->MakeSlotAvailable(&this->subSphereSizefactorParam);

    // set up RENDER DATA CONTROL paramters
    param::BoolParam* clusterNoisePara = new param::BoolParam(this->renderDataControl.clusterNoise);
    this->clusterNoiseParam << clusterNoisePara;
    this->MakeSlotAvailable(&this->clusterNoiseParam);

    param::BoolParam* clusterSpheresPara = new param::BoolParam(this->renderDataControl.clusterSpheres);
    this->clusterSpheresParam << clusterSpheresPara;
    this->MakeSlotAvailable(&this->clusterSpheresParam);

    param::BoolParam* hydrogenBondsPara = new param::BoolParam(this->renderDataControl.intForces.hydrogenBonds);
    this->hydrogenBondsParam << hydrogenBondsPara;
    this->MakeSlotAvailable(&this->hydrogenBondsParam);

    param::BoolParam* hydrogenConesPara = new param::BoolParam(this->renderDataControl.intForces.hydrogenCones);
    this->hydrogenConesParam << hydrogenConesPara;
    this->MakeSlotAvailable(&this->hydrogenConesParam);

    param::BoolParam* residuesPara = new param::BoolParam(this->renderDataControl.residues);
    this->residuesParam << residuesPara;
    this->MakeSlotAvailable(&this->residuesParam);

    param::BoolParam* residuesBackbonePara = new param::BoolParam(this->renderDataControl.residuesBackbone);
    this->residuesBackboneParam << residuesBackbonePara;
    this->MakeSlotAvailable(&this->residuesBackboneParam);

    param::BoolParam* residuesLabelsPara = new param::BoolParam(this->renderDataControl.residueLabels);
    this->residuesLabelsParam << residuesLabelsPara;
    this->MakeSlotAvailable(&this->residuesLabelsParam);

    param::BoolParam* interactionResiduesPara = new param::BoolParam(this->renderDataControl.interactionResidues);
    this->interactionResiduesParam << interactionResiduesPara;
    this->MakeSlotAvailable(&this->interactionResiduesParam);

    param::BoolParam* clusterAtomsPara = new param::BoolParam(this->renderDataControl.clusterAtoms);
    this->clusterAtomsParam << clusterAtomsPara;
    this->MakeSlotAvailable(&this->clusterAtomsParam);

    param::BoolParam* funcGroupCentroidsPara = new param::BoolParam(this->renderDataControl.funcGroupCentroids);
    this->funcGroupCentroidsParam << funcGroupCentroidsPara;
    this->MakeSlotAvailable(&this->funcGroupCentroidsParam);

    param::BoolParam* improveTextReadabilityPara = new param::BoolParam(this->renderDataControl.improveTextReadability);
    this->improveTextReadabilityParam << improveTextReadabilityPara;
    this->MakeSlotAvailable(&this->improveTextReadabilityParam);

    param::BoolParam* barchartPara = new param::BoolParam(this->renderDataControl.barchartGylphe);
    this->barchartGylpheParam << barchartPara;
    this->MakeSlotAvailable(&this->barchartGylpheParam);

    // set up RENDER DATA CONTROL paramters => PLIP FORCES
    param::BoolParam* halogenBondsPara = new param::BoolParam(this->renderDataControl.intForces.halogenBonds);
    this->halogenBondsParam << halogenBondsPara;
    this->MakeSlotAvailable(&this->halogenBondsParam);

    param::BoolParam* hydrophobicInteractionsPara =
        new param::BoolParam(this->renderDataControl.intForces.hydrophobicInteractions);
    this->hydrophobicInteractionsParam << hydrophobicInteractionsPara;
    this->MakeSlotAvailable(&this->hydrophobicInteractionsParam);

    param::BoolParam* metalComplexesPara = new param::BoolParam(this->renderDataControl.intForces.metalComplexes);
    this->metalComplexesParam << metalComplexesPara;
    this->MakeSlotAvailable(&this->metalComplexesParam);

    param::BoolParam* piCationInteractionsPara =
        new param::BoolParam(this->renderDataControl.intForces.piCationInteractions);
    this->piCationInteractionsParam << piCationInteractionsPara;
    this->MakeSlotAvailable(&this->piCationInteractionsParam);

    param::BoolParam* piStacksPara = new param::BoolParam(this->renderDataControl.intForces.piStacks);
    this->piStacksParam << piStacksPara;
    this->MakeSlotAvailable(&this->piStacksParam);

    param::BoolParam* saltBridgesPara = new param::BoolParam(this->renderDataControl.intForces.saltBridges);
    this->saltBridgesParam << saltBridgesPara;
    this->MakeSlotAvailable(&this->saltBridgesParam);


    // set up PROTEIN MESH COLORING CONTROL paramter
    param::BoolParam* interactionForcePara = new param::BoolParam(this->colorControl.interactionForce);
    this->interactionForceParam << interactionForcePara;
    this->MakeSlotAvailable(&this->interactionForceParam);

    param::BoolParam* interactionForceAmountColringPara =
        new param::BoolParam(this->colorControl.interactionForceAmountColoring);
    this->interactionForceAmountColoringParam << interactionForceAmountColringPara;
    this->MakeSlotAvailable(&this->interactionForceAmountColoringParam);

    // enum param: interaction force surface coloring
    param::EnumParam* cm0 =
        new param::EnumParam(int(this->subCircleContentIDxes[this->s.curForceType_surfaceColoring]));
    cm0->SetTypePair(this->subCircleContentIDxes[this->interForceTypes.HBonds], this->interForceTypes.HBonds.c_str());
    cm0->SetTypePair(
        this->subCircleContentIDxes[this->interForceTypes.halogenBonds], this->interForceTypes.halogenBonds.c_str());
    cm0->SetTypePair(this->subCircleContentIDxes[this->interForceTypes.hydrophobicInteractions],
        this->interForceTypes.hydrophobicInteractions.c_str());
    cm0->SetTypePair(this->subCircleContentIDxes[this->interForceTypes.metalComplexes],
        this->interForceTypes.metalComplexes.c_str());
    cm0->SetTypePair(this->subCircleContentIDxes[this->interForceTypes.piCationInteractions],
        this->interForceTypes.piCationInteractions.c_str());
    cm0->SetTypePair(
        this->subCircleContentIDxes[this->interForceTypes.piStacks], this->interForceTypes.piStacks.c_str());
    cm0->SetTypePair(
        this->subCircleContentIDxes[this->interForceTypes.saltBridges], this->interForceTypes.saltBridges.c_str());
    this->surfaceColoring_dropdownParam << cm0;
    this->MakeSlotAvailable(&this->surfaceColoring_dropdownParam);

    // enum param: interaction force starMenu sorting
    param::EnumParam* cm1 = new param::EnumParam(int(this->subCircleContentIDxes[this->starMenu.curSortingFroceType]));
    cm1->SetTypePair(this->subCircleContentIDxes[this->interForceTypes.HBonds], this->interForceTypes.HBonds.c_str());
    cm1->SetTypePair(
        this->subCircleContentIDxes[this->interForceTypes.halogenBonds], this->interForceTypes.halogenBonds.c_str());
    cm1->SetTypePair(this->subCircleContentIDxes[this->interForceTypes.hydrophobicInteractions],
        this->interForceTypes.hydrophobicInteractions.c_str());
    cm1->SetTypePair(this->subCircleContentIDxes[this->interForceTypes.metalComplexes],
        this->interForceTypes.metalComplexes.c_str());
    cm1->SetTypePair(this->subCircleContentIDxes[this->interForceTypes.piCationInteractions],
        this->interForceTypes.piCationInteractions.c_str());
    cm1->SetTypePair(
        this->subCircleContentIDxes[this->interForceTypes.piStacks], this->interForceTypes.piStacks.c_str());
    cm1->SetTypePair(
        this->subCircleContentIDxes[this->interForceTypes.saltBridges], this->interForceTypes.saltBridges.c_str());
    cm1->SetTypePair(
        this->subCircleContentIDxes[this->interForceTypes.bestScores], this->interForceTypes.bestScores.c_str());
    this->starMenuSorting_dropdownParam << cm1;
    this->MakeSlotAvailable(&this->starMenuSorting_dropdownParam);

    param::BoolParam* pocketAreaPara = new param::BoolParam(this->colorControl.pocketArea);
    this->pocketAreaParam << pocketAreaPara;
    this->MakeSlotAvailable(&this->pocketAreaParam);

    param::BoolParam* opaquePocketPara = new param::BoolParam(this->colorControl.opaquePocket);
    this->opaquePocketParam << opaquePocketPara;
    this->MakeSlotAvailable(&this->opaquePocketParam);

    param::BoolParam* transparencyPara = new param::BoolParam(this->colorControl.transparency);
    this->transparencyParam << transparencyPara;
    this->MakeSlotAvailable(&this->transparencyParam);

    param::BoolParam* chemicalPropertyPara = new param::BoolParam(this->colorControl.chemicalProperty);
    this->chemicalPropertyParam << chemicalPropertyPara;
    this->MakeSlotAvailable(&this->chemicalPropertyParam);


    // set up MESH RENDERING AMBIENT OCCLUSION parameters
    this->aoSampFactParam << new param::FloatParam(this->aoSampFact);
    this->MakeSlotAvailable(&this->aoSampFactParam);

    this->aoSampDistParam << new param::FloatParam(this->aoSampDist);
    this->MakeSlotAvailable(&this->aoSampDistParam);

    this->aoIntensityMinParam << new param::FloatParam(this->aoIntensityMin);
    this->MakeSlotAvailable(&this->aoIntensityMinParam);

    this->aoIntensityMaxParam << new param::FloatParam(this->aoIntensityMax);
    this->MakeSlotAvailable(&this->aoIntensityMaxParam);

    // set up CONTROL FUNCTIONAL GROUP CLUSTERING
    this->searchRadius_fgcParam << new param::FloatParam(this->searchRadius_fgc);
    this->MakeSlotAvailable(&this->searchRadius_fgcParam);

    this->minPTS_fgcParam << new param::IntParam(this->minPTS_fgc);
    this->MakeSlotAvailable(&this->minPTS_fgcParam);

    this->withoutSocketParam << new param::BoolParam(this->withoutSocket);
    this->MakeSlotAvailable(&this->withoutSocketParam);


    /************************************************
     *************** SET WEB-GUI DATA ***************
     ************************************************/

    // bool
    GUI["clusterNoise"] = (float)this->renderDataControl.clusterNoise;
    GUI["clusterSpheres"] = (float)this->renderDataControl.clusterSpheres;
    GUI["hydrogenBonds"] = (float)this->renderDataControl.intForces.hydrogenBonds;
    GUI["hydrogenCones"] = (float)this->renderDataControl.intForces.hydrogenCones;
    GUI["halogenBonds"] = (float)this->renderDataControl.intForces.halogenBonds;
    GUI["hydrophobicInteractions"] = (float)this->renderDataControl.intForces.hydrophobicInteractions;
    GUI["metalComplexes"] = (float)this->renderDataControl.intForces.metalComplexes;
    GUI["piCationInteractions"] = (float)this->renderDataControl.intForces.piCationInteractions;
    GUI["piStacks"] = (float)this->renderDataControl.intForces.piStacks;
    GUI["saltBridges"] = (float)this->renderDataControl.intForces.saltBridges;
    GUI["residues"] = (float)this->renderDataControl.residues;
    GUI["residuesBackbone"] = (float)this->renderDataControl.residuesBackbone;
    GUI["residuesLabels"] = (float)this->renderDataControl.residueLabels;
    GUI["interactionResidues"] = (float)this->renderDataControl.interactionResidues;
    GUI["clusterAtoms"] = (float)this->renderDataControl.clusterAtoms;
    GUI["funcGroupCentroids"] = (float)this->renderDataControl.funcGroupCentroids;
    GUI["improveTextReadability"] = (float)this->renderDataControl.improveTextReadability;
    GUI["barchartGylphe"] = (float)this->renderDataControl.barchartGylphe;
    GUI["pocketArea"] = (float)this->colorControl.pocketArea;
    GUI["opaquePocket"] = (float)this->colorControl.opaquePocket;
    GUI["transparency"] = (float)this->colorControl.transparency;
    GUI["chemicalProperty"] = (float)this->colorControl.chemicalProperty;
    GUI["interactionForce"] = (float)this->colorControl.interactionForce;
    GUI["interactionForceAmountColoring"] = (float)this->colorControl.interactionForceAmountColoring;
    GUI["minPTS_fgc"] = (float)this->minPTS_fgc;
    GUI["clipPlaneEnable"] = false;


    // int
    GUI["modelClustering_minPts"] = 10;
    GUI["minPTS_fgc"] = (float)this->minPTS_fgc;
    GUI["starMenuCircleCnt"] = (float)this->starMenu.circleCnt;
    GUI["starMenuSize"] = (float)this->starMenu.size;
    GUI["proteinColoringMode0"] = 1;
    GUI["proteinColoringMode1"] = 1;
    GUI["clipPlaneOrientation"] = 0;
    GUI["surfaceColoring_dropdown"] = this->surfaceColoring_dropdownParam.Param<param::EnumParam>()->Value();
    GUI["starMenuSorting_dropdown"] = this->starMenuSorting_dropdownParam.Param<param::EnumParam>()->Value();

    // float
    GUI["modelClustering_eps"] = 3.5f;
    GUI["modelClustering_minBindEnergy"] = -4.0f;
    GUI["searchRadius_fgc"] = (float)this->searchRadius_fgc;
    GUI["starMenuTextSize"] = (float)this->starMenu.textSize * 100;
    GUI["subSphereSizefactor"] = (float)this->subCircleSizefactor * 100;
    GUI["clipPlaneDist"] = 50.0f; // precent


    // text colors
    this->fontOutline_color = {0.0f, 0.0f, 0.0f, 255.0f};

    // this->fontResidueLables_color = {231.f, 138.f, 195.f, 255.f};
    this->fontResidueLables_color = {0.f, 0.f, 0.f, 255.f};
    this->fontFGSclustLables_color = {166.f, 216.f, 84.f, 255.0f};

    this->fontOutline_color /= 255.0f;
    this->starMenu.fontColor /= 255.0f;
    this->fontResidueLables_color /= 255.0f;
    this->fontFGSclustLables_color /= 255.0f;
}

/**
 * Destructor
 */
ModelClusterRenderer::~ModelClusterRenderer(void) { this->Release(); }

bool ModelClusterRenderer::create(void) {
    if (!vislib::graphics::gl::GLSLShader::InitialiseExtensions()) return false;

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glEnable(GL_VERTEX_PROGRAM_TWO_SIDE);
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_ARB);
    using megamol::core::utility::log::Log;

    using namespace vislib::sys;
    using namespace vislib::graphics::gl;

    ShaderSource triangleMesh_vertSrc;
    ShaderSource triangleMesh_fragSrc;

    ShaderSource dodecahedron_vertSrc;
    ShaderSource dodecahedron_fragSrc;

    ShaderSource simpleAmbient_vertSrc;
    ShaderSource simpleAmbient_fragSrc;

    ShaderSource circle_vertSrc;
    ShaderSource circle_geomSrc;
    ShaderSource circle_fragSrc;

    ShaderSource circle_Billboard_vertSrc;
    ShaderSource circle_Billboard_geomSrc;
    ShaderSource circle_Billboard_fragSrc;

    ShaderSource sphere_vertSrc;
    ShaderSource sphere_geomSrc;
    ShaderSource sphere_fragSrc;

    ShaderSource cone_vertSrc;
    ShaderSource cone_fragSrc;

    ShaderSource cylinder_vertSrc;
    ShaderSource cylinder_fragSrc;

    using namespace vislib::sys;
    using namespace vislib::graphics::gl;

    /************************************************************
     ******************* TRIANGLE MESH SHADER  ******************
     ************************************************************/
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource(
            "prolint::docking::withTransferFunction::triangleMesh::vertex", triangleMesh_vertSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to load vertex shader source for triangleMesh shader");
        return false;
    }
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource(
            "prolint::docking::withTransferFunction::triangleMesh::fragment", triangleMesh_fragSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to load fragment shader source for triangleMesh shader");
        return false;
    }
    try {
        if (!this->triangleMeshShader.Create(triangleMesh_vertSrc.Code(), triangleMesh_vertSrc.Count(),
                triangleMesh_fragSrc.Code(), triangleMesh_fragSrc.Count())) {
            throw vislib::Exception("Generic creation failure", __FILE__, __LINE__);
        }
    } catch (vislib::Exception e) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to create striangleMesh shader: %s\n", e.GetMsgA());
        return false;
    }

    /************************************************************
     ************* DODECAHEDRON GEOMTERY SHADER  ****************
     ************************************************************/

    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource(
            "dodecahedronShader::dodecahedronVertex", dodecahedron_vertSrc)) {
        Log::DefaultLog.WriteMsg(
            Log::LEVEL_ERROR, "Unable to load vertex shader source for geometry dodecahedron shader");
        return false;
    }
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource(
            "dodecahedronShader::dodecahedronFragment", dodecahedron_fragSrc)) {
        Log::DefaultLog.WriteMsg(
            Log::LEVEL_ERROR, "Unable to load fragment shader source for geometry dodecahedron shader");
        return false;
    }
    try {
        if (!this->dodecahedron_Shader.Create(dodecahedron_vertSrc.Code(), dodecahedron_vertSrc.Count(),
                dodecahedron_fragSrc.Code(), dodecahedron_fragSrc.Count())) {
            throw vislib::Exception("Generic creation failure", __FILE__, __LINE__);
        }
    } catch (vislib::Exception e) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to create dodecahedron_Shader: %s\n", e.GetMsgA());
        return false;
    }


    /*****************************************************
     ************* SIMPLE AMBINET SHADER  ****************
     *****************************************************/

    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource(
            "simpleAmbientShader::simpleAmbientShader_vertex", simpleAmbient_vertSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to load vertex shader source for simpleAmbientShader");
        return false;
    }
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource(
            "simpleAmbientShader::simpleAmbientShader_fragment", simpleAmbient_fragSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to load fragment shader source for  simpleAmbientShader");
        return false;
    }
    try {
        if (!this->simpleAmbientShader.Create(simpleAmbient_vertSrc.Code(), simpleAmbient_vertSrc.Count(),
                simpleAmbient_fragSrc.Code(), simpleAmbient_fragSrc.Count())) {
            throw vislib::Exception("Generic creation failure", __FILE__, __LINE__);
        }
    } catch (vislib::Exception e) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to create simpleAmbientShader: %s\n", e.GetMsgA());
        return false;
    }


    // setup geometry shader
    // set GL_TRIANGLES_ADJACENCY_EXT primitives as INPUT
    // this->dodecahedron_ShaderGeom.SetProgramParameter(GL_GEOMETRY_INPUT_TYPE_EXT, GL_POINTS);
    // set TRIANGLE_STRIP as OUTPUT
    // this->dodecahedron_ShaderGeom.SetProgramParameter(GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_TRIANGLE_STRIP);
    // set maximum number of vertices to be generated by geometry shader to
    // this->dodecahedron_ShaderGeom.SetProgramParameter(GL_GEOMETRY_VERTICES_OUT_EXT, 9*12);

    try {
        this->dodecahedron_Shader.Link();
    } catch (vislib::Exception e) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to link dodecahedron shader: %s\n", e.GetMsgA());
        return false;
    }


    /************************************************************
     ********** SPHERE GEOMTERY SHADER (for Circles) ************
     ************************************************************/
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource(
            "prolint::docking::sphereVertexGeom", circle_vertSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to load vertex shader source for geometry sphere shader");
        return false;
    }
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource(
            "prolint::docking::sphereGeom", circle_geomSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to load geometry shader source for geometry sphere shader");
        return false;
    }
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource(
            "prolint::docking::sphereFragmentGeom", circle_fragSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to load fragment shader source for geometry sphere shader");
        return false;
    }
    try {
        if (!this->circleShader.Compile(circle_vertSrc.Code(), circle_vertSrc.Count(), circle_geomSrc.Code(),
                circle_geomSrc.Count(), circle_fragSrc.Code(), circle_fragSrc.Count())) {
            throw vislib::Exception("Generic creation failure", __FILE__, __LINE__);
        }
    } catch (vislib::Exception e) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to create circleShader: %s\n", e.GetMsgA());
        return false;
    }

    // setup geometry shader
    // set GL_TRIANGLES_ADJACENCY_EXT primitives as INPUT
    this->circleShader.SetProgramParameter(GL_GEOMETRY_INPUT_TYPE_EXT, GL_POINTS);
    // set TRIANGLE_STRIP as OUTPUT
    this->circleShader.SetProgramParameter(GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_TRIANGLE_STRIP);
    // set maximum number of vertices to be generated by geometry shader to
    this->circleShader.SetProgramParameter(GL_GEOMETRY_VERTICES_OUT_EXT, 4);

    try {
        this->circleShader.Link();
    } catch (vislib::Exception e) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to link circleShader: %s\n", e.GetMsgA());
        return false;
    }


    /**********************************************
     ********** SPHERE GEOMTERY SHADER ************
     **********************************************/

    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource(
            "sphereShaderGeo::sphereVertexGeom", sphere_vertSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to load vertex shader source for geometry sphere shader");
        return false;
    }

    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource(
            "sphereShaderGeo::sphereGeom", sphere_geomSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to load geometry shader source for geometry sphere shader");
        return false;
    }

    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource(
            "sphereShaderGeo::sphereFragmentGeom_withLight", sphere_fragSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to load fragment shader source for geometry sphere shader");
        return false;
    }

    try {
        if (!this->sphereShader.Compile(sphere_vertSrc.Code(), sphere_vertSrc.Count(), sphere_geomSrc.Code(),
                sphere_geomSrc.Count(), sphere_fragSrc.Code(), sphere_fragSrc.Count())) {
            throw vislib::Exception("Generic creation failure", __FILE__, __LINE__);
        }
    } catch (vislib::Exception e) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to create sphereShader: %s\n", e.GetMsgA());
        return false;
    }

    // setup geometry shader
    // set GL_TRIANGLES_ADJACENCY_EXT primitives as INPUT
    this->sphereShader.SetProgramParameter(GL_GEOMETRY_INPUT_TYPE_EXT, GL_POINTS);
    // set TRIANGLE_STRIP as OUTPUT
    this->sphereShader.SetProgramParameter(GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_TRIANGLE_STRIP);
    // set maximum number of vertices to be generated by geometry shader to
    this->sphereShader.SetProgramParameter(GL_GEOMETRY_VERTICES_OUT_EXT, 4);

    try {
        this->sphereShader.Link();
    } catch (vislib::Exception e) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to link sphereShader: %s\n", e.GetMsgA());
        return false;
    }


    /******************************************************
     ********** CIRCLE TEXTURE GEOMTERY SHADER ************
     ******************************************************/

    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource(
            "circleBillboardShader::circleVertex", circle_Billboard_vertSrc)) {
        Log::DefaultLog.WriteMsg(
            Log::LEVEL_ERROR, "Unable to load vertex shader source for circleBillboardShader::circleVertex");
        return false;
    }
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource(
            "circleBillboardShader::circleGeom", circle_Billboard_geomSrc)) {
        Log::DefaultLog.WriteMsg(
            Log::LEVEL_ERROR, "Unable to load geometry shader source for circleBillboardShader::circleGeom");
        return false;
    }
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource(
            "circleBillboardShader::circleFragment", circle_Billboard_fragSrc)) {
        Log::DefaultLog.WriteMsg(
            Log::LEVEL_ERROR, "Unable to load fragment shader source for circleBillboardShader::circleFragment");
        return false;
    }
    try {
        if (!this->circleBillTextureShader.Compile(circle_Billboard_vertSrc.Code(), circle_Billboard_vertSrc.Count(),
                circle_Billboard_geomSrc.Code(), circle_Billboard_geomSrc.Count(), circle_Billboard_fragSrc.Code(),
                circle_Billboard_fragSrc.Count())) {
            throw vislib::Exception("Generic creation failure", __FILE__, __LINE__);
        }
    } catch (vislib::Exception e) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to create circleBillTextureShader: %s\n", e.GetMsgA());
        return false;
    }


    // setup geometry shader
    // set GL_TRIANGLES_ADJACENCY_EXT primitives as INPUT
    this->circleBillTextureShader.SetProgramParameter(GL_GEOMETRY_INPUT_TYPE_EXT, GL_POINTS);
    // set TRIANGLE_STRIP as OUTPUT
    this->circleBillTextureShader.SetProgramParameter(GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_TRIANGLE_STRIP);
    // set maximum number of vertices to be generated by geometry shader to
    this->circleBillTextureShader.SetProgramParameter(GL_GEOMETRY_VERTICES_OUT_EXT, 4);

    try {
        this->circleBillTextureShader.Link();
    } catch (vislib::Exception e) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to link circleBillTextureShader: %s\n", e.GetMsgA());
        return false;
    }


    /************************************************************
     ************************ CONE SHADER  **********************
     ************************************************************/
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource("cone::vertex", cone_vertSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to load vertex shader source for cone::vertex");
        return false;
    }
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource("cone::fragment", cone_fragSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to load fragment shader source for cone::fragment");
        return false;
    }
    try {
        if (!this->coneShader.Create(
                cone_vertSrc.Code(), cone_vertSrc.Count(), cone_fragSrc.Code(), cone_fragSrc.Count())) {
            throw vislib::Exception("Generic creation failure", __FILE__, __LINE__);
        }
    } catch (vislib::Exception e) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to create coneShader: %s\n", e.GetMsgA());
        return false;
    }


    /*
    // setup geometry shader
    // set GL_TRIANGLES_ADJACENCY_EXT primitives as INPUT
    this->coneShader.SetProgramParameter(GL_GEOMETRY_INPUT_TYPE_EXT, GL_POINTS);
    // set TRIANGLE_STRIP as OUTPUT
    this->coneShader.SetProgramParameter(GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_TRIANGLE_STRIP);
    // set maximum number of vertices to be generated by geometry shader to
    this->segments = 20;
    this->coneShader.SetProgramParameter(GL_GEOMETRY_VERTICES_OUT_EXT, (int)this->segments * 2 + 1);
    */

    try {
        this->coneShader.Link();
    } catch (vislib::Exception e) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to link coneShader: %s\n", e.GetMsgA());
        return false;
    }


    /************************************************************
     ********************** CYLINDER SHADER  ********************
     ************************************************************/
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource("cylinder::vertex", cylinder_vertSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to load vertex shader source for cone::vertex");
        return false;
    }
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource("cylinder::fragment", cylinder_fragSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to load fragment shader source for cone::fragment");
        return false;
    }
    try {
        if (!this->cylinderShader.Create(
                cylinder_vertSrc.Code(), cylinder_vertSrc.Count(), cylinder_fragSrc.Code(), cylinder_fragSrc.Count())) {
            throw vislib::Exception("Generic creation failure", __FILE__, __LINE__);
        }
    } catch (vislib::Exception e) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to create coneShader: %s\n", e.GetMsgA());
        return false;
    }


    /*
    // setup geometry shader
    // set GL_TRIANGLES_ADJACENCY_EXT primitives as INPUT
    this->coneShader.SetProgramParameter(GL_GEOMETRY_INPUT_TYPE_EXT, GL_POINTS);
    // set TRIANGLE_STRIP as OUTPUT
    this->coneShader.SetProgramParameter(GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_TRIANGLE_STRIP);
    // set maximum number of vertices to be generated by geometry shader to
    this->segments = 20;
    this->coneShader.SetProgramParameter(GL_GEOMETRY_VERTICES_OUT_EXT, (int)this->segments * 2 + 1);
    */

    try {
        this->cylinderShader.Link();
    } catch (vislib::Exception e) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to link coneShader: %s\n", e.GetMsgA());
        return false;
    }


    // initialize font object
    if (!this->font.IsInitialised()) {
        this->font.Initialise(this->GetCoreInstance());
        this->font.SetBillboardMode(true);
    }

    // TODO allocate variables etc.

    return true;
}

void ModelClusterRenderer::release(void) {}


bool ModelClusterRenderer::MouseEvent(int x, int y, core::view::MouseFlags flags) { return true; }


bool ModelClusterRenderer::getMousePos(core::Call& call) {
    MouseInputCall* cmi = dynamic_cast<MouseInputCall*>(&call);
    if (cmi == NULL) {
        printf("\nWrong Call!\n");
        return false;
    }
    this->mouseX = cmi->GetMouseX();
    this->mouseY = cmi->GetMouseY();
    this->cam = cmi->GetCam();
    this->mouseFlag = cmi->GetMouseFlags();
    if (this->mouseFlag == 0) {
        // printf("Clicked\n");
        this->isClicked = 1;
        // printf("--mouseX: %d \t mouseY: %d \t flag:%d \n", cmi->GetMouseX(), cmi->GetMouseY(), cmi->GetMouseFlags());
    } else {
        isClicked = 0;
    }


    return true;
}


//###################################################################
//###################################################################
//###################################################################

bool ModelClusterRenderer::GetExtents(megamol::core::Call& call) {

    view::AbstractCallRender3D* cr3d = dynamic_cast<view::AbstractCallRender3D*>(&call);
    if (cr3d == NULL) return false;

    MolecularDataCall* mdc = this->ligandMolecularData_callerSlot.CallAs<MolecularDataCall>();
    if (mdc == NULL) return false;
    if (!(*mdc)(MolecularDataCall::CallForGetExtent)) return false;

    MolecularDataCall* pmdc = this->proteinMolecularData_callerSlot.CallAs<MolecularDataCall>();
    if (pmdc == NULL) return false;
    if (!(*pmdc)(MolecularDataCall::CallForGetExtent)) return false;

    misc::VolumetricDataCall* vdc = this->volumetricData_callerSlot.CallAs<misc::VolumetricDataCall>();
    if (vdc == NULL) return false;
    // if (!(*vdc)(1)) return false;
    // if (!(*vdc)(0)) return false;

    CallTriMeshData* tmc = this->triMeshData_callerSlot.CallAs<CallTriMeshData>();
    float scale = 1.0f;
    if (tmc != NULL) {
        if (!(*tmc)(1)) return false;
        if (!(*tmc)(0)) return false;
        std::vector<float> edgeOfCalls;
        std::vector<float> volOfCalls;
        edgeOfCalls.resize(4);
        volOfCalls.resize(4);
        edgeOfCalls[0] = mdc->AccessBoundingBoxes().ObjectSpaceBBox().LongestEdge();
        edgeOfCalls[1] = pmdc->AccessBoundingBoxes().ObjectSpaceBBox().LongestEdge();
        edgeOfCalls[2] = tmc->AccessBoundingBoxes().ObjectSpaceBBox().LongestEdge();
        edgeOfCalls[3] = vdc->AccessBoundingBoxes().ObjectSpaceBBox().LongestEdge();

        volOfCalls[0] = mdc->AccessBoundingBoxes().ObjectSpaceBBox().Volume();
        volOfCalls[1] = pmdc->AccessBoundingBoxes().ObjectSpaceBBox().Volume();
        volOfCalls[2] = tmc->AccessBoundingBoxes().ObjectSpaceBBox().Volume();
        volOfCalls[3] = vdc->AccessBoundingBoxes().ObjectSpaceBBox().Volume();

        float maxVol = volOfCalls[0];
        scale = 2.0f / edgeOfCalls[0];
        if (vislib::math::IsEqual(edgeOfCalls[0], 0.0f)) return false;
        for (size_t i = 1; i < edgeOfCalls.size(); i++) {
            if (!vislib::math::IsEqual(edgeOfCalls[i], 0.0f)) {
                if (volOfCalls[i] > maxVol) {
                    scale = 2.0f / edgeOfCalls[i];
                    maxVol = volOfCalls[i];
                }
            }
        }
        cr3d->AccessBoundingBoxes().MakeScaledWorld(scale);
        cr3d->SetTimeFramesCount(mdc->FrameCount());
        this->scaleWorld = scale;
    }

    return true;
}

bool ModelClusterRenderer::Render(megamol::core::Call& call) {

    int changedRenderData = 0;
    // cast the call to Render3D
    view::AbstractCallRender3D* cr3d = dynamic_cast<view::AbstractCallRender3D*>(&call);
    if (cr3d == nullptr) {
        return false;
    }

    // get camera information
    this->cameraInfo = cr3d->GetCameraParameters();

    // get pointer to MolecularDataCall
    MolecularDataCall* mdc = this->ligandMolecularData_callerSlot.CallAs<MolecularDataCall>();
    if (mdc == nullptr) 
	{
        return false;
    }

    // get pointer to LigandModelCall
    LigandModelCall* lmc = this->ligandModel_callerSlot.CallAs<LigandModelCall>();
    if (lmc == nullptr) {
        return false;
    }

    if (!(*lmc)(LigandModelCall::CallForGetExtent))
	{
        return false;
    }

    auto status = lmc->getLMCStatus();

    while (status.isSocketCommManagerBusy) {
        // wait for SocketCommunicationManager to finish
        if (!(*lmc)(LigandModelCall::CallForGetExtent)) {
            return false;
        }
        status = lmc->getLMCStatus();
    }

    this->withoutSocket = this->withoutSocketParam.Param<param::BoolParam>()->Value();

    if (status.isFirstRenderCallComplete == false) {

        /**************************************************
         ************** GET EXTERN PARAM SLOTs ************
         **************************************************/

        // get protein coloringMode from  "MSMSMeshLoader"
        this->proteinColoringMode0_Param =
            dynamic_cast<param::ParamSlot*>(this->FindNamedObject(this->paramName_proteinColoring0, true).get());
        this->proteinColoringMode1_Param =
            dynamic_cast<param::ParamSlot*>(this->FindNamedObject(this->paramName_proteinColoring1, true).get());
        int paramColoringMode0_value = this->proteinColoringMode0_Param->Param<param::EnumParam>()->Value();
        int paramColoringMode1_value = this->proteinColoringMode1_Param->Param<param::EnumParam>()->Value();
        GUI["proteinColoringMode0"] = static_cast<float>(paramColoringMode0_value);
        GUI["proteinColoringMode1"] = static_cast<float>(paramColoringMode1_value);

        // get DBSCAN clsutering parms from  "MultiPDBQTLoader"
        this->modelClustering_eps_Param =
            dynamic_cast<param::ParamSlot*>(this->FindNamedObject(this->paramName_modelClustering_eps, true).get());
        this->modelClustering_minBindEnergy_Param = dynamic_cast<param::ParamSlot*>(
            this->FindNamedObject(this->paramName_modelClustering_minBindEnergy, true).get());
        this->modelClustering_minPts_Param =
            dynamic_cast<param::ParamSlot*>(this->FindNamedObject(this->paramName_modelClustering_minPts, true).get());

        float eps = this->modelClustering_eps_Param->Param<param::FloatParam>()->Value();
        int minBindEnergy = this->modelClustering_minBindEnergy_Param->Param<param::FloatParam>()->Value();
        int minPts = this->modelClustering_minPts_Param->Param<param::IntParam>()->Value();

        GUI["modelClustering_eps"] = static_cast<float>(eps);
        GUI["modelClustering_minBindEnergy"] = static_cast<float>(minBindEnergy);
        GUI["modelClustering_minPts"] = static_cast<float>(minPts);

        // get clip lane params from "ClipPlane1"

        this->clipPlaneEnable_Param =
            dynamic_cast<param::ParamSlot*>(this->FindNamedObject(this->paramName_clipPlaneEnable, true).get());
        this->clipPlaneNormal_Param =
            dynamic_cast<param::ParamSlot*>(this->FindNamedObject(this->paramName_clipPlaneNormal, true).get());
        this->clipPlaneDist_Param =
            dynamic_cast<param::ParamSlot*>(this->FindNamedObject(this->paramName_clipPlaneDist, true).get());

        bool clipPlaneEnable = this->clipPlaneEnable_Param->Param<param::BoolParam>()->Value();
        float clipPlaneDist = this->clipPlaneDist_Param->Param<param::FloatParam>()->Value();

        GUI["clipPlaneEnable"] = (float)clipPlaneEnable;
        GUI["clipPlaneDist"] = (float)clipPlaneDist;


        this->cartoonRenderer_Param =
            dynamic_cast<param::ParamSlot*>(this->FindNamedObject(this->paramName_cartoonRnederer, true).get());

        bool isCartoonRenderer = this->cartoonRenderer_Param->Param<param::BoolParam>()->Value();
        GUI["isCartoonRenderer"] = (float)isCartoonRenderer;

        LigandModelCall* lmc = this->ligandModel_callerSlot.CallAs<LigandModelCall>();
        (*lmc)(LigandModelCall::CallForGetExtent);
        auto tmpWebdata = lmc->getWebData();
        tmpWebdata->GUI = &this->GUI;


// nasty way: preload all SVGs to have quicker star menu building
// : there is an out commented code fragment like this later in the code to
// only load currently needed SVGs and store them for later use
// but for stability reasons for the moment it is put on this place here
#ifdef SVGs_preLoad
        for (int ID = 0; ID < lmc->getLigandCount(); ID++) {
            lmc->SetCurrentLigandAndModel(ID, 0);
            if (!(*lmc)(LigandModelCall::CallForGetDataSilent)) return false;
            if (this->int2String.find(ID) == int2String.end()) {
                this->int2String.insert({ID, this->SVGcnt});
                int mapID = int2String[ID];
                SVG1 a;
                this->svg1.push_back(a);
            std:
                string SVGPath = lmc->getSVGPath();
                if (this->verbose) printf("SVGPath: %s\n", SVGPath.c_str());
                parseSVG(lmc->getSVGPath(), &svg1[this->SVGcnt]);
                this->SVGcnt++;
            }
        }
#endif // SVGs_preLoad
    }

    /** trigger clustering once if not yet done for this data set
     *  lmc->DataHash() is by default 0
     *  after a sucessfull file loading of PDBQT files the variable lmc_dataHash is >0
     *  but lmc->SetdataHash(lmc_dataHash) will be only executed for the first time if there was an lmc CALL
     *  for exmaple through this renderer
     *  so the "if" statment would be (0 != 0 ) for the first
     *  and without any change of the data (1 != 0)
     */

    if (!(*lmc)(LigandModelCall::CallForGetClusterAndCentroidData)) return false; //# FIXME
    this->clusterData = lmc->getClusterAndCentroidData();
    if (this->clusterData.DBSCANParams.paramsChanged) {
        forcesLineStorage.resetAll();
        this->markedVerticesIndicesMap.clear();
        this->maxForceCntPerProtAtomMap.clear();
        while (this->clusterData.clusterDataHash == this->old_clusterDataHash) {
            this->clusterData = lmc->getClusterAndCentroidData();
            Sleep(100);
        }
    }

    if (clusterData.DBSCANParams.paramsChanged == TRUE) {
        changedRenderData = 1;
        this->pmdcOut_newData = true;

        // update starMenu.circleCnt if cluster size (minPts) is below it
        if (this->starMenu.circleCnt > this->clusterData.DBSCANParams.minPts) {
            if (this->clusterData.DBSCANParams.minPts < 1) {
                this->starMenuCntParam.Param<param::IntParam>()->SetValue(1);
            } else {
                this->starMenuCntParam.Param<param::IntParam>()->SetValue(this->clusterData.DBSCANParams.minPts);
            }
        }
    }

    /**************************************************
     ************** UPDATE BOOL PARAMETERs ************
     **************************************************/

    LigandModelCall::WebData* webdata = lmc->getWebData();
    auto webGUI = (*webdata->GUI);
    bool changedParam = false;

    this->webDataChanged = false;
    if (this->old_getWebData_dataHash != lmc->getWebData()->getWebDataHash() &&
        lmc->getWebData()->getWebDataHash() != NULL) {
        this->webDataChanged = true;
        this->old_getWebData_dataHash = lmc->getWebData()->getWebDataHash();
    }
    // auto webGUI = this->GUI;

    // control render data
    if (updateBoolParams(&this->clusterNoiseParam, "clusterNoise", webdata, &this->renderDataControl.clusterNoise)) {
        changedParam = true;
        changedRenderData = 1;
    };
    if (updateBoolParams(
            &this->clusterSpheresParam, "clusterSpheres", webdata, &this->renderDataControl.clusterSpheres)) {
        changedParam = true;
        changedRenderData = 1;
    };
    if (updateBoolParams(&this->interactionResiduesParam, "interactionResidues", webdata,
            &this->renderDataControl.interactionResidues)) {
        changedParam = true;
        this->pmdcOut_newData = true;
    };
    if (updateBoolParams(&this->residuesParam, "residues", webdata, &this->renderDataControl.residues)) {
        changedParam = true;
        this->pmdcOut_newData = true;
    };
    if (updateBoolParams(
            &this->residuesBackboneParam, "residuesBackbone", webdata, &this->renderDataControl.residuesBackbone)) {
        changedParam = true;
        this->pmdcOut_newData = true;
    };
    if (updateBoolParams(
            &this->residuesLabelsParam, "residuesLabels", webdata, &this->renderDataControl.residueLabels)) {
        changedParam = true;
    };
    if (updateBoolParams(&this->clusterAtomsParam, "clusterAtoms", webdata, &this->renderDataControl.clusterAtoms)) {
        changedParam = true;
    };
    if (updateBoolParams(&this->funcGroupCentroidsParam, "funcGroupCentroids", webdata,
            &this->renderDataControl.funcGroupCentroids)) {
        changedParam = true;
    };
    if (updateBoolParams(&this->improveTextReadabilityParam, "improveTextReadability", webdata,
            &this->renderDataControl.improveTextReadability)) {
        changedParam = true;
    };
    if (updateBoolParams(
            &this->barchartGylpheParam, "barchartGylphe", webdata, &this->renderDataControl.barchartGylphe)) {
        changedParam = true;
    };
    if (updateBoolParams(
            &this->hydrogenBondsParam, "hydrogenBonds", webdata, &this->renderDataControl.intForces.hydrogenBonds)) {
        changedParam = true;
    };
    if (updateBoolParams(
            &this->hydrogenConesParam, "hydrogenCones", webdata, &this->renderDataControl.intForces.hydrogenCones)) {
        changedParam = true;
    };


    // control render data => plipForces
    if (updateBoolParams(
            &this->halogenBondsParam, "halogenBonds", webdata, &this->renderDataControl.intForces.halogenBonds)) {
        changedParam = true;
    };
    if (updateBoolParams(&this->hydrophobicInteractionsParam, "hydrophobicInteractions", webdata,
            &this->renderDataControl.intForces.hydrophobicInteractions)) {
        changedParam = true;
    };
    if (updateBoolParams(
            &this->metalComplexesParam, "metalComplexes", webdata, &this->renderDataControl.intForces.metalComplexes)) {
        changedParam = true;
    };
    if (updateBoolParams(&this->piCationInteractionsParam, "piCationInteractions", webdata,
            &this->renderDataControl.intForces.piCationInteractions)) {
        changedParam = true;
    };
    if (updateBoolParams(&this->piStacksParam, "piStacks", webdata, &this->renderDataControl.intForces.piStacks)) {
        changedParam = true;
    };
    if (updateBoolParams(
            &this->saltBridgesParam, "saltBridges", webdata, &this->renderDataControl.intForces.saltBridges)) {
        changedParam = true;
    };

    // control color
    if (updateBoolParams(&this->pocketAreaParam, "pocketArea", webdata, &this->colorControl.pocketArea)) {
        changedParam = true;
    };
    if (updateBoolParams(&this->opaquePocketParam, "opaquePocket", webdata, &this->colorControl.opaquePocket)) {
        changedParam = true;
    };
    if (updateBoolParams(&this->transparencyParam, "transparency", webdata, &this->colorControl.transparency)) {
        changedParam = true;
    };
    if (updateBoolParams(
            &this->chemicalPropertyParam, "chemicalProperty", webdata, &this->colorControl.chemicalProperty)) {
        changedParam = true;
    };
    if (updateBoolParams(
            &this->interactionForceParam, "interactionForce", webdata, &this->colorControl.interactionForce)) {
        changedParam = true;
    };
    if (updateBoolParams(&this->interactionForceAmountColoringParam, "interactionForceAmountColoring", webdata,
            &this->colorControl.interactionForceAmountColoring)) {
        changedParam = true;
    };


    /**************************************************
     ************** UPDATE INT PARAMETERs *************
     **************************************************/


    // star menu params
    int placeHolder;
    if (updateIntParams(&this->starMenuCntParam, "starMenuCircleCnt", webdata, &placeHolder)) {
        changedParam = true;
    };
    if (updateIntParams(&this->starMenuSizeParam, "starMenuSize", webdata, &placeHolder)) {
        changedParam = true;
    };
    if (updateIntParams(&this->starMenuTextSizeParam, "starMenuTextSize", webdata, &placeHolder)) {
        changedParam = true;
    };
    if (updateIntParams(&this->subSphereSizefactorParam, "subSphereSizefactor", webdata, &placeHolder)) {
        changedParam = true;
    };

    /**************************************************
     ************** UPDATE ENUM PARAMETERs ************
     **************************************************/
    if (updateEnumParams(&this->surfaceColoring_dropdownParam, "surfaceColoring_dropdown", webdata,
            &this->colorControl.surfaceColoring_dropdown)) {
        changedParam = true;
        s.curForceType_surfaceColoring =
            this->interForceTypes.forceTypeVec[this->colorControl.surfaceColoring_dropdown];
        this->s.changedCurForceType = true;
    };
    int curForceID = 0;
    if (updateEnumParams(&this->starMenuSorting_dropdownParam, "starMenuSorting_dropdown", webdata, &curForceID)) {
        changedParam = true;
        this->starMenu.curSortingFroceType = this->interForceTypes.forceTypeVec[curForceID];
        this->starMenu.changed_curSortingFroceType = true;
        this->starMenu.curPage = 0;
        s.selectedGlobalMdlID = -1;
    };


    /*********************************************************************
    ************** UPDATE Functional Group Clustering Params ************
    *********************************************************************/
    if (updateIntParams(&this->minPTS_fgcParam, "minPTS_fgc", webdata, &this->minPTS_fgc)) {
        changedParam = true;
        // only valid values
        if (this->minPTS_fgc < 1) {
            this->minPTS_fgc = 1;
            this->minPTS_fgcParam.Param<param::IntParam>()->SetValue(this->minPTS_fgc);
        }
        this->fgcClusterParamsChanged = true;
    };

    /**************************************************
     ************* UPDATE FLOAT PARAMETERs ************
     **************************************************/

    if (this->searchRadius_fgcParam.IsDirty() || this->webDataChanged) {
        auto p = this->searchRadius_fgcParam.Param<param::FloatParam>();
        if (webGUI["searchRadius_fgc"] != (float)p->Value() && this->webDataChanged) {
            p->SetValue((float)webGUI["searchRadius_fgc"]);
            this->searchRadius_fgcParam.ForceSetDirty();
        }
        if (this->searchRadius_fgcParam.IsDirty()) {
            this->searchRadius_fgc = p->Value();
            GUI["searchRadius_fgc"] = (float)p->Value();
            changedParam = true;

            // only valid values
            if (this->searchRadius_fgc < 0.1f) {
                this->searchRadius_fgc = 0.1f;
                this->searchRadius_fgcParam.Param<param::FloatParam>()->SetValue(this->searchRadius_fgc);
            }
            this->fgcClusterParamsChanged = true;
            this->searchRadius_fgcParam.ResetDirty();
        }
    }


    /****************************************************
     ************** UPDATE EXTREN PARAMETERs ************
     ****************************************************/
    // check for dirty is not possible for extern params!
    if (this->proteinColoringMode0_Param->Param<param::EnumParam>()->Value() != (int)webGUI["proteinColoringMode0"] ||
        this->webDataChanged) {
        auto p = this->proteinColoringMode0_Param->Param<param::EnumParam>();
        if (webGUI["proteinColoringMode0"] != (float)p->Value() && this->webDataChanged) {
            p->SetValue((int)webGUI["proteinColoringMode0"]);
            this->proteinColoringMode0_Param->ForceSetDirty();
        } else if (webGUI["proteinColoringMode0"] != (float)p->Value()) {
            GUI["proteinColoringMode0"] = (float)p->Value();
            changedParam = true;
            // this->proteinColoringMode0_Param->ResetDirty();
        }
    }

    if (this->proteinColoringMode1_Param->Param<param::EnumParam>()->Value() != (int)webGUI["proteinColoringMode1"] ||
        this->webDataChanged) {
        auto p = this->proteinColoringMode1_Param->Param<param::EnumParam>();
        if (webGUI["proteinColoringMode1"] != (float)p->Value() && this->webDataChanged) {
            p->SetValue((int)webGUI["proteinColoringMode1"]);
            this->proteinColoringMode1_Param->ForceSetDirty();
        } else if (webGUI["proteinColoringMode1"] != (float)p->Value()) {
            GUI["proteinColoringMode1"] = (float)p->Value();
            changedParam = true;
            // this->proteinColoringMode1_Param->ResetDirty();
        }
    }

    if (this->modelClustering_eps_Param->Param<param::FloatParam>()->Value() != (float)webGUI["modelClustering_eps"] ||
        this->webDataChanged) {
        auto p = this->modelClustering_eps_Param->Param<param::FloatParam>();
        if (webGUI["modelClustering_eps"] != (float)p->Value() && this->webDataChanged) {
            float tmpVal = (float)webGUI["modelClustering_eps"];
            p->SetValue(tmpVal);
            changedParam = true;
            this->modelClustering_eps_Param->ForceSetDirty();
        } else if (webGUI["modelClustering_eps"] != (float)p->Value()) {
            GUI["modelClustering_eps"] = (float)p->Value();
            changedParam = true;
            // this->modelClustering_eps_Param->ResetDirty();
        }
    }

    if (this->modelClustering_minBindEnergy_Param->Param<param::FloatParam>()->Value() !=
            (int)webGUI["modelClustering_minBindEnergy"] ||
        this->webDataChanged) {
        auto p = this->modelClustering_minBindEnergy_Param->Param<param::FloatParam>();
        if (webGUI["modelClustering_minBindEnergy"] != (float)p->Value() && this->webDataChanged) {
            p->SetValue((float)webGUI["modelClustering_minBindEnergy"]);
            changedParam = true;
            this->modelClustering_minBindEnergy_Param->ForceSetDirty();
        } else if (webGUI["modelClustering_minBindEnergy"] != (float)p->Value()) {
            GUI["modelClustering_minBindEnergy"] = (float)p->Value();
            changedParam = true;
            // this->modelClustering_minBindEnergy_Param->ResetDirty();
        }
    }

    if (this->modelClustering_minPts_Param->Param<param::IntParam>()->Value() !=
            (int)webGUI["modelClustering_minPts"] ||
        this->webDataChanged) {
        auto p = this->modelClustering_minPts_Param->Param<param::IntParam>();
        if (webGUI["modelClustering_minPts"] != (float)p->Value() && this->webDataChanged) {
            p->SetValue((int)webGUI["modelClustering_minPts"]);
            changedParam = true;
            this->modelClustering_minPts_Param->ForceSetDirty();
        } else if (webGUI["modelClustering_minPts"] != (float)p->Value()) {
            GUI["modelClustering_minPts"] = (float)p->Value();
            changedParam = true;
            // this->modelClustering_minPts_Param->ResetDirty();
        }
    }

    if (this->clipPlaneEnable_Param->Param<param::BoolParam>()->Value() != (bool)webGUI["clipPlaneEnable"] ||
        this->webDataChanged) {
        auto p = this->clipPlaneEnable_Param->Param<param::BoolParam>();
        if (webGUI["clipPlaneEnable"] != (float)p->Value() && this->webDataChanged) {
            p->SetValue((int)webGUI["clipPlaneEnable"]);
            changedParam = true;
            this->clipPlaneEnable_Param->ForceSetDirty();
        } else if (webGUI["clipPlaneEnable"] != (float)p->Value()) {
            GUI["clipPlaneEnable"] = (float)p->Value();
            changedParam = true;
            // this->clipPlaneEnable_Param->ResetDirty();
        }
    }


    if (this->webDataChanged || (float)webGUI["clipPlaneDist"] != (float)this->old_clipPlaneDist ||
        this->clipPlaneDist_Param->Param<param::FloatParam>()->Value() != this->old_clipPlane_distVal) {
        auto p = this->clipPlaneDist_Param->Param<param::FloatParam>();
        if ((float)webGUI["clipPlaneDist"] != (float)this->old_clipPlaneDist ||
            this->clipPlaneDist_Param->Param<param::FloatParam>()->Value() != this->old_clipPlane_distVal ||
            (int)webGUI["clipPlaneOrientation"] != (int)this->old_clipPlaneOrientation) {
            int check = (int)webGUI["clipPlaneOrientation"];
            auto a = cr3d->AccessBoundingBoxes().ObjectSpaceBBox().GetOrigin();
            float oriVal = a.GetZ();
            auto dim = cr3d->AccessBoundingBoxes().ObjectSpaceBBox().Depth();
            if (check == 1) {
                dim = cr3d->AccessBoundingBoxes().ObjectSpaceBBox().Height();
                oriVal = a.GetY();
            } else if (check == 2) {
                oriVal = a.GetX();
                dim = cr3d->AccessBoundingBoxes().ObjectSpaceBBox().Width();
            }

            float distVal = oriVal + (float)webGUI["clipPlaneDist"] / 100.0f * dim;
            p->SetValue(distVal);
            this->old_clipPlane_distVal = distVal;
            this->old_clipPlaneDist = (float)webGUI["clipPlaneDist"];
            changedParam = true;
            this->clipPlaneDist_Param->ForceSetDirty();
        }
    }


    if (this->webDataChanged || (int)webGUI["clipPlaneOrientation"] != (int)this->old_clipPlaneOrientation) {
        if ((int)webGUI["clipPlaneOrientation"] != (int)this->old_clipPlaneOrientation) {
            int check = (int)webGUI["clipPlaneOrientation"];
            this->old_clipPlaneOrientation = check;
            changedParam = true;
            glm::vec3 ori;
            if (check == 0) {
                ori = this->planeOrientations["XY"];
            } else if (check == 1) {
                ori = this->planeOrientations["XZ"];
            } else if (check == 2) {
                ori = this->planeOrientations["YZ"];
            }

            this->clipPlaneNormal_Param->Param<param::Vector3fParam>()->SetValue(
                vislib::math::Point<float, 3>(ori.x, ori.y, ori.z));
            this->clipPlaneNormal_Param->ForceSetDirty();
        }
    }

     if (this->cartoonRenderer_Param->Param<param::BoolParam>()->Value() != (bool)webGUI["isCartoonRenderer"] ||
        this->webDataChanged) {
        auto p = this->cartoonRenderer_Param->Param<param::BoolParam>();
        if (webGUI["isCartoonRenderer"] != (float)p->Value() && this->webDataChanged) {
            p->SetValue((int)webGUI["isCartoonRenderer"]);
            changedParam = true;
            this->cartoonRenderer_Param->ForceSetDirty();
        } else if (webGUI["isCartoonRenderer"] != (float)p->Value()) {
            GUI["isCartoonRenderer"] = (float)p->Value();
            changedParam = true;
            // this->clipPlaneEnable_Param->ResetDirty();
        }
    }


    if (changedParam) {
        // lmc->UpdateWebData();
        lmc->getWebData()->UpdateGUIData();
    }

    /************************************************************
    ************** GET and UPDATE trnasfer functions ************
    *************************************************************/

    // tf_interactionForce
    this->tf.tf_interactionForce =
        this->transferFunction_interactionForce_callerSlot.CallAs<core::view::CallGetTransferFunction>();
    if (this->tf.tf_interactionForce == nullptr || !(*this->tf.tf_interactionForce)()) {
        megamol::core::utility::log::Log::DefaultLog.WriteMsg(megamol::core::utility::log::Log::LEVEL_ERROR,
            "ModelClusterRenderer requires a transfer function for tf_interactionForce!");
    }

    // tf_clusterSpheres
    this->tf.tf_clusterSpheres =
        this->transferFunction_clusterSpheres_callerSlot.CallAs<core::view::CallGetTransferFunction>();
    if (this->tf.tf_clusterSpheres == nullptr || !(*this->tf.tf_clusterSpheres)()) {
        megamol::core::utility::log::Log::DefaultLog.WriteMsg(megamol::core::utility::log::Log::LEVEL_ERROR,
            "ModelClusterRenderer requires a transfer function for tf_clusterSpheres!");
    }

    // tf_functGroups
    this->tf.tf_functGroups =
        this->transferFunction_functGroups_callerSlot.CallAs<core::view::CallGetTransferFunction>();
    if (this->tf.tf_functGroups == nullptr || !(*this->tf.tf_functGroups)()) {
        megamol::core::utility::log::Log::DefaultLog.WriteMsg(megamol::core::utility::log::Log::LEVEL_ERROR,
            "ModelClusterRenderer requires a transfer function for tf_functGroups!");
    }


    // collect data for rendering (sphere data and color)
    if (changedRenderData) {

        this->s.selectedIDx = -1;
        this->s.hoveringIDx = -1;
        this->s.old_SelectedIDx = -1;
        this->s.old_HoveringIDx = -1;
        this->noiseModelCentroidsIDXs.Resize(0);
        this->noiseModelCentroidsIDXs.AssertCapacity(1);
        printf("changed render data...\n");
        for (int z = 0; z < this->clusterData.assignedClusters->Count(); z++) {
            if (this->clusterData.assignedClusters->PeekElements()[z] == -1 && this->renderDataControl.clusterNoise) {
                noiseModelCentroidsIDXs.Append(z);
            }
        }


        this->clusterDodeca_verts.Resize(0);
        this->clusterDodeca_normals.Resize(0);


        int clusterCentroidCnt = 0;
        if (this->renderDataControl.clusterSpheres) {
            clusterCentroidCnt = this->clusterData.clusterCentroids->Count();
        }
        int combineSize = (noiseModelCentroidsIDXs.Count() * 4 + clusterCentroidCnt);
        combined_noiseModelCentroids_clusterCentroids.Resize(0);
        combined_noiseModelCentroids_clusterCentroids.SetCount(combineSize);
        sphereColorsScales.Resize(0);
        sphereColorsScales.SetCount(combineSize / 4);

        int noiseCnter = noiseModelCentroidsIDXs.Count();


// set sphere data and color
#pragma omp parallel for
        for (int i = 0; i < combineSize / 4; i++) {
            if (i < noiseCnter) {
                int IDx = this->clusterData.assign_modelID_sphereID->PeekElements()[i].GetY() * 4;
                auto modelCentroids = this->clusterData.modelCentroids->PeekElements();
                combined_noiseModelCentroids_clusterCentroids[i * 4 + 0] = modelCentroids[IDx + 0];
                combined_noiseModelCentroids_clusterCentroids[i * 4 + 1] = modelCentroids[IDx + 1];
                combined_noiseModelCentroids_clusterCentroids[i * 4 + 2] = modelCentroids[IDx + 2];
                combined_noiseModelCentroids_clusterCentroids[i * 4 + 3] = this->noiseModel_CentroidRadius;
                sphereColorsScales[i] = 0.0f;
            } else {
                int IDx = (i - noiseCnter) * 4;

                auto clusterCentroids = this->clusterData.clusterCentroids->PeekElements();
                combined_noiseModelCentroids_clusterCentroids[i * 4 + 0] = (clusterCentroids[IDx + 0]);
                combined_noiseModelCentroids_clusterCentroids[i * 4 + 1] = (clusterCentroids[IDx + 1]);
                combined_noiseModelCentroids_clusterCentroids[i * 4 + 2] = (clusterCentroids[IDx + 2]);
                combined_noiseModelCentroids_clusterCentroids[i * 4 + 3] = this->cluster_CentroidRadius;
                float scale = 1.0f - (float)this->clusterData.clusterSizes[0][i - noiseCnter] /
                                         (float)this->clusterData.maxClusterSize;
                sphereColorsScales[i] = (float)this->clusterData.clusterSizes[0][i - noiseCnter];
            }
        }
        printf("Spheres: %d\n", combined_noiseModelCentroids_clusterCentroids.Count() / 4);
    } // endif changed render data

    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glEnable(GL_VERTEX_PROGRAM_TWO_SIDE);
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_ARB);
    glPushMatrix();


    // compute scale factor and scale world
    glScalef(this->scaleWorld, this->scaleWorld, this->scaleWorld);


    this->transformMatrices = this->mouseRayIntersection.getTransformMatrices(this->cam);
    glm::vec4 camPos4_cam = glm::vec4(0, 0, 0, 1);
    this->camPos_local = glm::vec3(camPos4_cam * this->transformMatrices.camera_to_local);


    // get viewpoint parameters for raycasting
    float viewportStuff[4] = {this->cameraInfo->TileRect().Left(), this->cameraInfo->TileRect().Bottom(),
        this->cameraInfo->TileRect().Width(), this->cameraInfo->TileRect().Height()};
    if (viewportStuff[2] < 1.0f) viewportStuff[2] = 1.0f;
    if (viewportStuff[3] < 1.0f) viewportStuff[3] = 1.0f;
    viewportStuff[2] = 2.0f / viewportStuff[2];
    viewportStuff[3] = 2.0f / viewportStuff[3];

    // Get GL_MODELVIEW matrix
    glGetFloatv(GL_MODELVIEW_MATRIX, this->modelMatrix_column);
    vislib::math::Matrix<GLfloat, 4, vislib::math::COLUMN_MAJOR> modelMatrix_intern(&this->modelMatrix_column[0]);
    this->modelMatrix = &modelMatrix_intern;
    // Get GL_PROJECTION matrix
    glGetFloatv(GL_PROJECTION_MATRIX, this->projMatrix_column);
    vislib::math::Matrix<GLfloat, 4, vislib::math::COLUMN_MAJOR> projMatrix_intern(&this->projMatrix_column[0]);
    this->projMatrix = &projMatrix_intern;

    //***********************************************************************


    //************************************************************************************//
    //******************************* HOVERING AND PICKING *******************************//
    //************************************************************************************//
    (*lmc)(LigandModelCall::CallForGetDataSilent);


    // draw selected Model
    int selectedStarMenuModelIDx = -1;
    if (starCircles.Count() > 0 && s.curClusterID != -1) {
        selectedStarMenuModelIDx = mouseRayIntersection.getFirstMouseRayIntersection(
            starCircles, this->cam, this->mouseX, this->mouseY, this->isClicked);
        // printf("selectedStarMenuModelIDx: %d\n", selectedStarMenuModelIDx);
    }

    // gives priority to star menu circle selections and disables hovering/picking for underneath spheres
    if (selectedStarMenuModelIDx >= 0) {
        if (this->isClicked == 1) {
            this->curSelectedStarMenuModelIDx = selectedStarMenuModelIDx;
            // check if 'close menu' is selected
            if (this->curSelectedStarMenuModelIDx == this->indexOfLastStarCircle + 1) {
                this->s.selectedIDx = -1;
                this->isClicked = -1;
                this->s.selectedGlobalMdlID = -1;
                lmc->getWebData()->setClusterID(-1);
            } else {
                this->s.selectedGlobalMdlID =
                    this->oneClusterData_sortedWithEnergy[this->starMenu.curPageStartIDx + selectedStarMenuModelIDx]
                        .GetX();


                lmc->SetCurrentLigandAndModel_byGlobalMdlID(this->s.selectedGlobalMdlID);
                (*lmc)(LigandModelCall::CallForGetData);
                auto web = lmc->getWebData();
                web->setLigandID(lmc->getCurrentLigandID());
                web->setModelID(lmc->getCurrentModelID());
                this->s.hoveringIDx = -1;
            }
        }
        // gives priority for up/down arrows if star menu and disables hovering/picking for underneath spheres
    } else if (this->arrowUpDown >= 0) {
        this->s.hoveringIDx = -1;
    } else {
        this->s.hoveringIDx = mouseRayIntersection.getFirstMouseRayIntersection(
            combined_noiseModelCentroids_clusterCentroids, this->cam, this->mouseX, this->mouseY, this->isClicked);
        if (this->s.hoveringIDx != -1 && this->isClicked == 1) {
            this->s.selectedIDx = this->s.hoveringIDx;

            lmc->SetCurrentLigandAndModel_byGlobalMdlID(this->s.hoveringIDx);
            lmc->getWebData()->setLigandID(lmc->getCurrentLigandID());
            lmc->getWebData()->setModelID(lmc->getCurrentModelID());
        }
    }


    // load from web-interface requested model
    if (this->webDataChanged) {
        auto webData = lmc->getWebData();

        this->s.selectedGlobalMdlID =
            webData->getGlobalMdlID() != -1 ? webData->getGlobalMdlID() : this->s.selectedGlobalMdlID;
        this->s.selectedIDx = webData->getClusterID() + this->noiseModelCentroidsIDXs.Count();
        if (webData->getClusterID() == -1) {
            this->s.selectedGlobalMdlID = -1;
        }
    }


    // detect change of clusterID
    this->s.curClusterID = this->s.selectedIDx - this->noiseModelCentroidsIDXs.Count();
    if (this->s.curClusterID >= 0 && this->s.curClusterID < this->clusterData.assignedClusters->Count()) {
        this->s.curClusterID = this->s.selectedIDx - this->noiseModelCentroidsIDXs.Count();
        // lmc->getWebData()->setClusterID(s.curClusterID);
    } else {
        this->s.curClusterID = -1;
    }


    if (this->s.old_curClusterID != this->s.curClusterID) {
        this->s.changedcurClusterID = 1;
        this->s.old_curClusterID = this->s.curClusterID;
        lmc->getWebData()->setClusterID(this->s.curClusterID);
    }


    /*
    bool deselect = false;
    if (this->isClicked == 1 && this->s.selectedIDx == this->s.old_SelectedIDx && this->s.hoveringIDx != -1 &&
        this->s.selectedIDx != -1) {
        deselect = true;
        printf("deselect %d\n", this->s.selectedIDx);
        this->s.selectedIDx = -1;
    }
    */

    // printf("hover: %d", this->hoveringIDx);
    // set back radius -> when selection has changed
    if ((this->s.selectedIDx != this->s.old_SelectedIDx && this->s.old_SelectedIDx > -1)) {
        if (this->s.old_SelectedIDx < noiseModelCentroidsIDXs.Count()) {
            this->curSelectedStarMenuModelIDx = -1;
            combined_noiseModelCentroids_clusterCentroids[this->s.old_SelectedIDx * 4 + 3] =
                this->noiseModel_CentroidRadius;
        } else {
            int clusterIDx = this->s.old_SelectedIDx - this->noiseModelCentroidsIDXs.Count();
            this->curSelectedStarMenuModelIDx = -1;
            combined_noiseModelCentroids_clusterCentroids[this->s.old_SelectedIDx * 4 + 3] =
                this->cluster_CentroidRadius;
        }
    }

    // set back radius -> when hovering has changed
    if (this->s.hoveringIDx != this->s.old_HoveringIDx && this->s.old_HoveringIDx > -1) {
        if (this->s.old_HoveringIDx < noiseModelCentroidsIDXs.Count()) {
            combined_noiseModelCentroids_clusterCentroids[this->s.old_HoveringIDx * 4 + 3] =
                this->noiseModel_CentroidRadius;
        } else {
            int clusterIDx = this->s.old_HoveringIDx - this->noiseModelCentroidsIDXs.Count();
            combined_noiseModelCentroids_clusterCentroids[this->s.old_HoveringIDx * 4 + 3] =
                this->cluster_CentroidRadius;
        }
    }


    // set color and radius for current hovering
    if (this->s.hoveringIDx >= 0) {
        if (this->s.hoveringIDx < noiseModelCentroidsIDXs.Count()) {
            combined_noiseModelCentroids_clusterCentroids[this->s.hoveringIDx * 4 + 3] =
                this->higlight_noiseModel_CentroidRadius;
        } else {
            combined_noiseModelCentroids_clusterCentroids[this->s.hoveringIDx * 4 + 3] =
                this->higlight_cluster_CentroidRadius;
        }
    }

    // set color and radius for current selection
    if (this->s.selectedIDx >= 0) {
        if (this->s.selectedIDx < noiseModelCentroidsIDXs.Count()) {
            combined_noiseModelCentroids_clusterCentroids[this->s.selectedIDx * 4 + 3] = 0;
            // this->higlight_noiseModel_CentroidRadius;
        } else {
            combined_noiseModelCentroids_clusterCentroids[this->s.selectedIDx * 4 + 3] = 0;
            // this->higlight_cluster_CentroidRadius;
        }
    }


    /******************************************************
     ******** Dertermine functional groups clusters********
     ******************************************************/

    // calc and store the functional groups for each cluster
    if ((this->clusterData.DBSCANParams.paramsChanged && this->clusterData.clusterDataHash > 1) ||
        this->fgcClusterParamsChanged) {
        // handle sending data to python server
        auto fgs = lmc->getClusteredFunctionalGroups();
        bool modelClustRenderer_busy = 1;
        bool socketCommManager_busy = 0;
        bool firstRenderCall_complete = 1;
        if (fgs->getDataHash() >= 0 && fgs->getDataHash() != this->old_functionalGroup_dataHash && !withoutSocket) {
            this->old_functionalGroup_dataHash = fgs->getDataHash();
            lmc->SetlmcStatus(modelClustRenderer_busy, socketCommManager_busy, firstRenderCall_complete);
            (*lmc)(LigandModelCall::CallForSetLMCStatus);
        }
        // allocate storage to store a seperate functional-group-class-object of each cluster
        if (this->fgsPerClusterComplete.size() > 0) {
            for (int i = 0; i < fgsPerClusterComplete.size(); i++) {
                this->fgsPerClusterComplete[i].clusterID = -1;
                this->fgsPerClusterFiltered[i].clusterID = -1;
            }
            this->funcGroupDodeca_verts.Resize(0);
            this->funcGroupDodeca_normals.Resize(0);
        }

        this->fgsPerClusterComplete.resize(0);
        this->fgsPerClusterComplete.resize(this->clusterData.clusterSizes->Count());
        this->fgsPerClusterFiltered.resize(0);
        this->fgsPerClusterFiltered.resize(this->clusterData.clusterSizes->Count());

        printf("determining functional groups...\t radius: %f minPTS: %d\n", this->searchRadius_fgc, this->minPTS_fgc);
        for (int i = 0; i < this->clusterData.clusterSizes->Count(); i++) {
            // set back dodecahedron selection
            this->renVarsFGS.selectedSphereID = -1;
            this->fgsHandler.handle(
                i, this->minPTS_fgc, this->searchRadius_fgc, this->fgsPerClusterFiltered[i], mdc, lmc);
            // int clustCnt = this->fgsPerClusterFiltered[i].clusterCnt;
            this->fgsPerClusterComplete[i].clusterCnt = 1;
            this->fgsPerClusterComplete[i].fgsMaps_GblMDlIDs.resize(1);
            this->fgsPerClusterComplete[i].clusterSizes.resize(1);
            this->fgsPerClusterComplete[i].clusterSizes[0] = 0;

            // iterate over all fgs clusters of current pocket
            for (int j = 0; j < this->fgsPerClusterFiltered[i].clusterCnt; j++) {
                // set FGSmap of current fgs cluster
                // iterate of typedID and their gblMdlIDs vectors
                map<uint, std::vector<uint>>::iterator it;
                auto curFGSmap_gblMdl = this->fgsPerClusterFiltered[i].fgsMaps_GblMDlIDs[j];
                auto curFGSmap_centalAtomIDx = this->fgsPerClusterFiltered[i].fgsMaps_centralAtomArrayIDx[j];
                for (it = curFGSmap_gblMdl.begin(); it != curFGSmap_gblMdl.end(); it++) {
                    uint typeID = it->first;
                    printf("pocket:%d subClust:%d tyepID:%d \n", i, j, typeID);
                    // iterate over  gblMdlIDs vector
                    for (int y = 0; y < it->second.size(); y++) {
                        // set current ligand
                        lmc->SetCurrentLigandAndModel_byGlobalMdlID(it->second[y]);
                        // get data of current ligand
                        (*lmc)(LigandModelCall::CallForGetData);
                        (*mdc)(MolecularDataCall::CallForGetData);
                        auto FGS_ofCurLig = lmc->getFunctionalGroupsOfCurLigand();
                        // iterate of fgs of current ligand model
                        for (int z = FGS_ofCurLig->totalGroupCnt - 1; z >= 0; z--) {
                            auto f = FGS_ofCurLig->fGroups[z];
                            int curFgsTypeID = f.groupID;
                            for (int u = FGS_ofCurLig->fGroups[z].speficGroupCnt - 1; u >= 0; u--) {
                                bool check = false;
                                if (f.groupID == typeID) {
                                    check = true;
                                } else {
                                    while (curFgsTypeID != -1) {
                                        int higher_lvl = lmc->fgsHierarchy[curFgsTypeID];
                                        if (higher_lvl == typeID) {
                                            check = true;
                                            break;
                                        }
                                        curFgsTypeID = higher_lvl;
                                    }
                                }
                                // printf("curFGSmap_centalAtomIDx[typeID][y] == u :%d == %d\n",
                                // curFGSmap_centalAtomIDx[typeID][y], u);
                                if (check && curFGSmap_centalAtomIDx[typeID][y] == u) {
                                    this->fgsPerClusterComplete[i].fgsMaps_GblMDlIDs[0][f.groupID].push_back(
                                        it->second[y]);
                                    this->fgsPerClusterComplete[i].clusterSizes[0]++;
                                }
                            }
                        }
                    }
                }
            }
        }
        this->clusteredFunctionalGroups = lmc->getClusteredFunctionalGroups();
        this->clusteredFunctionalGroups->setPointer_to_functionalGroupClusterSets(
            this->fgsPerClusterComplete.data(), this->fgsPerClusterComplete.size());
    }

    if (webdata->selectedFgsClusterID != this->s.selectedFgsClusterID) {
        webdata->selectedFgsClusterID = this->s.selectedFgsClusterID;

        if (this->s.selectedFgsClusterID != -1 && this->renVarsFGS.selectedSphereID != -1) {

            auto fgsMap =
                this->fgsPerClusterFiltered[this->s.curClusterID].fgsMaps_GblMDlIDs[this->renVarsFGS.selectedSphereID];
            // set fgsMap for webdata
            webdata->curFgsMap = fgsMap;
            webdata->UpdateWebData();
        }
    }


    // CHECK IF NEEDED
    // if (this->functionalGroupHandler.dataChanged) {
    //    this->renVarsFGS.selectedSphereID = -1;
    //}

    // set arrays back if selected cluster changed
    if (this->s.changedcurClusterID && this->s.curClusterID > -1 || this->fgcClusterParamsChanged) {
        this->allClusterAtoms.Resize(0);
        this->funcGroupDodeca_verts.Resize(0);
        this->funcGroupDodeca_normals.Resize(0);
        this->starSubSpheres_onSelect.Resize(0);
        this->starMenu.curPageStartIDx = 0;
        this->starMenu.curPage = 0;

        /*
         if (this->s.curClusterID > -1) {
            vislib::graphics::SceneSpacePoint3D lookAt;
            lookAt.SetX(lmc->getClusterAndCentroidData().clusterCentroids->PeekElements()[this->s.curClusterID * 4 +
         0]); lookAt.SetY(lmc->getClusterAndCentroidData().clusterCentroids->PeekElements()[this->s.curClusterID * 4 +
         1]); lookAt.SetZ(lmc->getClusterAndCentroidData().clusterCentroids->PeekElements()[this->s.curClusterID * 4 +
         2]);

            this->cameraInfo->SetPosition(lookAt);
            cr3d->SetCameraParameters(this->cameraInfo);
            this->cameraInfo->SetView(this->cameraInfo->Position(), this->cameraInfo->LookAt(), this->cameraInfo->Up());
         }
         */
    }


    if (this->tf.tf_clusterSpheres->IsDirty() || changedRenderData) {
        this->tf.tf_clusterSpheres->SetRange(
            {(float)this->clusterData.minClusterSize, (float)this->clusterData.maxClusterSize});
        this->tf_clusterSizes.Resize(0);
        this->tf_clusterSizes.SetCount(this->sphereColorsScales.Count());
        for (int i = 0; i < this->sphereColorsScales.Count(); i++) {
            this->tf_clusterSizes[i] = (float)this->sphereColorsScales[i];
        }

        glGenBuffers(1, &this->gl_sphereColorsScales);
        glDeleteBuffers(1, &this->gl_sphereColorsScales);
        glGenBuffers(1, &this->gl_sphereColorsScales);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, this->gl_sphereColorsScales);
        glBufferStorage(GL_SHADER_STORAGE_BUFFER, this->tf_clusterSizes.Count() * 1 * sizeof(float),
            this->tf_clusterSizes.PeekElements(), GL_DYNAMIC_STORAGE_BIT);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
        this->tf.tf_clusterSpheres->ResetDirty();
        this->renVarsCNS.gl_sphereValuesSSBO = this->gl_sphereColorsScales;
        this->renVarsCNS.hideOnHover = false;
    }


    // render all atoms of selected cluster
    if (this->renderDataControl.clusterAtoms) {
        if (this->allClusterAtoms.Count() <= 0) {
            auto clusterDat = lmc->getClusterAndCentroidData();
            auto clusterIndexSum = clusterDat.clusterIndexSum->PeekElements();

            int clusterIDx_start = this->s.curClusterID > 0
                                       ? this->clusterData.clusterIndexSum->PeekElements()[this->s.curClusterID]
                                       : this->clusterData.clusterIndexSum->PeekElements()[0];
            int clusterIDx_end = this->clusterData.clusterIndexSum->PeekElements()[this->s.curClusterID + 1];
            int clusterSize = this->clusterData.clusterSizes->PeekElements()[this->s.curClusterID];
            this->allClusterAtoms.Resize(0);
            for (int j = clusterIDx_start; j < clusterIDx_end; j++) {

                lmc->SetCurrentLigandAndModel_byGlobalMdlID(
                    clusterDat.assign_modelID_sphereID->PeekElements()[j].GetY());
                (*lmc)(LigandModelCall::CallForGetData);
                (*mdc)(MolecularDataCall::CallForGetData);


                for (int i = 0; i < mdc->AtomCount(); i++) {
                    this->allClusterAtoms.Append(mdc->AtomPositions()[i * 3 + 0]);
                    this->allClusterAtoms.Append(mdc->AtomPositions()[i * 3 + 1]);
                    this->allClusterAtoms.Append(mdc->AtomPositions()[i * 3 + 2]);
                    this->allClusterAtoms.Append(0.3f);
                }
            }
        }
        /*************************************
         ******** Render Cluster Atoms ********
         *************************************/
        if (this->allClusterAtoms.Count() > 0) {
            sphereRenderer((float*)this->allClusterAtoms.PeekElements(), this->allClusterAtoms.Count() / 4,
                this->tf.tf_functGroups, this->renVarsACAS, viewportStuff);
        }
    }

    if (this->s.curClusterID > -1) {
        if (this->fgsPerClusterFiltered[this->s.curClusterID].clusterCnt > 0 &&
            this->renderDataControl.funcGroupCentroids) {
            int hover_IDx = this->mouseRayIntersection.getFirstMouseRayIntersection(
                this->fgsPerClusterFiltered[this->s.curClusterID].clusterCentroids, this->cam, this->mouseX,
                this->mouseY, 1);
            this->renVarsFGS.hoverSphereID = hover_IDx;
            if (hover_IDx != -1 && this->isClicked == 1) {
                this->renVarsFGS.selectedSphereID = hover_IDx;
            }

            if (this->s.changedcurClusterID || this->tf.tf_functGroups->IsDirty() || this->fgcClusterParamsChanged) {


                // workaround: TODO check why the data from the vector can't be directly set to SSBO?
                this->fgsClustSizes.Resize(0);
                this->fgsClustSizes.SetCount(this->fgsPerClusterFiltered[this->s.curClusterID].clusterSizes.size());

                for (int i = 0; i < this->fgsClustSizes.Count(); i++) {
                    this->fgsClustSizes[i] = this->fgsPerClusterFiltered[this->s.curClusterID].clusterSizes[i];
                }


                glGenBuffers(1, &this->gl_funczGroupColorsScales);
                glDeleteBuffers(1, &this->gl_funczGroupColorsScales);
                glGenBuffers(1, &this->gl_funczGroupColorsScales);
                glBindBuffer(GL_SHADER_STORAGE_BUFFER, this->gl_funczGroupColorsScales);
                glBufferStorage(GL_SHADER_STORAGE_BUFFER, this->fgsClustSizes.Count() * sizeof(float),
                    this->fgsClustSizes.PeekElements(), GL_DYNAMIC_STORAGE_BIT);
                glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
                this->renVarsFGS.gl_sphereValuesSSBO = this->gl_funczGroupColorsScales;
                this->renVarsFGS.hoverColor = {1, 0.647f, 0.0784f, 1};

                vislib::math::Vector<float, 2> te;
                te = this->fgsHandler.getMinMaxCntFGS(this->fgsPerClusterFiltered[this->s.curClusterID]);
                this->tf.tf_functGroups->SetRange({(float)te.GetX(), (float)te.GetY()});
                this->tf.tf_functGroups->ResetDirty();
            }
            /******************************************
             ******** Render FUNCTIONAL GROUPs ********
             ******************************************/

            if (this->fgsPerClusterFiltered[this->s.curClusterID].clusterCentroids.Count() > 0) {

                if (this->s.changedcurClusterID) {
                    this->funcGroupDodeca_verts.Resize(0);
                    this->funcGroupDodeca_normals.Resize(0);
                }

                // sphereRenderer(this->functionalGroupsPerCluster[this->curClusterID].grpClusterFGSCentroid,
                //    this->tf.tf_functGroups, this->renVarsFGS, viewportStuff);
                this->dodecahedronRenderer(this->fgsPerClusterFiltered[this->s.curClusterID].clusterCentroids,
                    &this->funcGroupDodeca_verts, &this->funcGroupDodeca_normals, this->tf.tf_functGroups,
                    this->renVarsFGS, viewportStuff);
            }
        }
    }
    this->s.selectedFgsClusterID = this->renVarsFGS.selectedSphereID;


    /********************************************************
     ******** Render Cluster and NoiseModelCentroids ********
     ********************************************************/

    this->renVarsCNS.hoverSphereID = this->s.hoveringIDx;


    // if (this->combined_noiseModelCentroids_clusterCentroids.Count() / 4 >
    // this->clusterData.clusterSizes->Count()) {
    //   this->renVarsCNS.hideOnSelection = false;
    //      sphereRenderer(this->combined_noiseModelCentroids_clusterCentroids, this->tf.tf_clusterSpheres,
    //      this->renVarsCNS,   viewportStuff);
    // }


    this->renVarsCNS.selectedSphereID = this->s.selectedIDx;
    this->renVarsCNS.hideOnSelection = true;
    this->dodecahedronRenderer(this->combined_noiseModelCentroids_clusterCentroids, &this->clusterDodeca_verts,
        &this->clusterDodeca_normals, this->tf.tf_clusterSpheres, this->renVarsCNS, viewportStuff);

    ::glDisable(GL_LIGHTING);
    ::glEnable(GL_DEPTH_TEST);
    /**********************************************
     ************ PLIP FORCE RENDERING ************
     **********************************************/

    lmc->SetCurrentLigandAndModel_byGlobalMdlID(s.selectedGlobalMdlID);
    if (lmc->getLMCStatus().isFirstRenderCallComplete && this->renderDataControl.intForces.saltBridges) {
        auto list = forcesLineStorage.getCurPosArray(this->interForceTypes.saltBridges, lmc);
        auto size = forcesLineStorage.getForceCnt(this->interForceTypes.saltBridges, lmc);
        auto startPositions = forcesLineStorage.getCurStartPosArray(this->interForceTypes.saltBridges, lmc);
        auto directions = forcesLineStorage.getCurDirArray(this->interForceTypes.saltBridges, lmc);
        auto lenghts = forcesLineStorage.getCurLengthArray(this->interForceTypes.saltBridges, lmc);

        auto col = &this->forceColorsMap[this->interForceTypes.saltBridges];
        this->renVarsIFC.baseColor.Set(col->x, col->y, col->z, col->w);
        // StickRenderer(list, size, this->tf.tf_functGroups, this->renVarsACAS, viewportStuff);
        cylinderRenderer(startPositions, directions, lenghts, size, this->renVarsIFC, viewportStuff);
    }
    if (this->renderDataControl.intForces.halogenBonds) {
        auto list = forcesLineStorage.getCurPosArray(this->interForceTypes.halogenBonds, lmc);
        auto size = forcesLineStorage.getForceCnt(this->interForceTypes.halogenBonds, lmc);
        auto startPositions = forcesLineStorage.getCurStartPosArray(this->interForceTypes.halogenBonds, lmc);
        auto directions = forcesLineStorage.getCurDirArray(this->interForceTypes.halogenBonds, lmc);
        auto lenghts = forcesLineStorage.getCurLengthArray(this->interForceTypes.halogenBonds, lmc);

        auto col = &this->forceColorsMap[this->interForceTypes.halogenBonds];
        this->renVarsIFC.baseColor.Set(col->x, col->y, col->z, col->w);
        cylinderRenderer(startPositions, directions, lenghts, size, this->renVarsIFC, viewportStuff);
    }
    if (this->renderDataControl.intForces.hydrophobicInteractions) {
        auto list = forcesLineStorage.getCurPosArray(this->interForceTypes.hydrophobicInteractions, lmc);
        auto size = forcesLineStorage.getForceCnt(this->interForceTypes.hydrophobicInteractions, lmc);
        auto startPositions = forcesLineStorage.getCurStartPosArray(this->interForceTypes.hydrophobicInteractions, lmc);
        auto directions = forcesLineStorage.getCurDirArray(this->interForceTypes.hydrophobicInteractions, lmc);
        auto lenghts = forcesLineStorage.getCurLengthArray(this->interForceTypes.hydrophobicInteractions, lmc);

        auto col = &this->forceColorsMap[this->interForceTypes.hydrophobicInteractions];
        this->renVarsIFC.baseColor.Set(col->x, col->y, col->z, col->w);
        cylinderRenderer(startPositions, directions, lenghts, size, this->renVarsIFC, viewportStuff);
    }
    if (this->renderDataControl.intForces.metalComplexes) {
        auto list = forcesLineStorage.getCurPosArray(this->interForceTypes.metalComplexes, lmc);
        auto size = forcesLineStorage.getForceCnt(this->interForceTypes.metalComplexes, lmc);
        auto startPositions = forcesLineStorage.getCurStartPosArray(this->interForceTypes.metalComplexes, lmc);
        auto directions = forcesLineStorage.getCurDirArray(this->interForceTypes.metalComplexes, lmc);
        auto lenghts = forcesLineStorage.getCurLengthArray(this->interForceTypes.metalComplexes, lmc);

        auto col = &this->forceColorsMap[this->interForceTypes.metalComplexes];
        this->renVarsIFC.baseColor.Set(col->x, col->y, col->z, col->w);
        cylinderRenderer(startPositions, directions, lenghts, size, this->renVarsIFC, viewportStuff);
    }
    if (this->renderDataControl.intForces.piCationInteractions) {
        auto list = forcesLineStorage.getCurPosArray(this->interForceTypes.piCationInteractions, lmc);
        auto size = forcesLineStorage.getForceCnt(this->interForceTypes.piCationInteractions, lmc);
        auto startPositions = forcesLineStorage.getCurStartPosArray(this->interForceTypes.piCationInteractions, lmc);
        auto directions = forcesLineStorage.getCurDirArray(this->interForceTypes.piCationInteractions, lmc);
        auto lenghts = forcesLineStorage.getCurLengthArray(this->interForceTypes.piCationInteractions, lmc);

        auto col = &this->forceColorsMap[this->interForceTypes.piCationInteractions];
        this->renVarsIFC.baseColor.Set(col->x, col->y, col->z, col->w);
        cylinderRenderer(startPositions, directions, lenghts, size, this->renVarsIFC, viewportStuff);
    }
    if (this->renderDataControl.intForces.piStacks) {
        auto list = forcesLineStorage.getCurPosArray(this->interForceTypes.piStacks, lmc);
        auto size = forcesLineStorage.getForceCnt(this->interForceTypes.piStacks, lmc);
        auto startPositions = forcesLineStorage.getCurStartPosArray(this->interForceTypes.piStacks, lmc);
        auto directions = forcesLineStorage.getCurDirArray(this->interForceTypes.piStacks, lmc);
        auto lenghts = forcesLineStorage.getCurLengthArray(this->interForceTypes.piStacks, lmc);

        auto col = &this->forceColorsMap[this->interForceTypes.piStacks];
        this->renVarsIFC.baseColor.Set(col->x, col->y, col->z, col->w);
        cylinderRenderer(startPositions, directions, lenghts, size, this->renVarsIFC, viewportStuff);
        /*
        for (int j = 0; j < size; j++) {
            glBegin(GL_LINES);
            glColor4f(1.0f, 1.0f, 0.2f, 1.0f);
            glVertex3f(list[j * 4 + 0], list[j * 4 + 1], list[j * 4 + 2]);
            j++;
            glVertex3f(list[j * 4 + 0], list[j * 4 + 1], list[j * 4 + 2]);
            glEnd();
        }
        */
    }
    if (this->renderDataControl.intForces.hydrogenBonds) {
        auto list = forcesLineStorage.getCurPosArray(this->interForceTypes.Hbonds_protAcceptor, lmc);
        auto size = forcesLineStorage.getForceCnt(this->interForceTypes.Hbonds_protAcceptor, lmc);
        auto startPositions = forcesLineStorage.getCurStartPosArray(this->interForceTypes.Hbonds_protAcceptor, lmc);
        auto directions = forcesLineStorage.getCurDirArray(this->interForceTypes.Hbonds_protAcceptor, lmc);
        auto lenghts = forcesLineStorage.getCurLengthArray(this->interForceTypes.Hbonds_protAcceptor, lmc);

        auto col = &this->forceColorsMap[this->interForceTypes.Hbonds_protAcceptor];
        this->renVarsIFC.baseColor.Set(col->x, col->y, col->z, col->w);
        cylinderRenderer(startPositions, directions, lenghts, size, this->renVarsIFC, viewportStuff);
    }
    if (this->renderDataControl.intForces.hydrogenBonds) {
        auto list = forcesLineStorage.getCurPosArray(this->interForceTypes.Hbonds_protDonor, lmc);
        auto size = forcesLineStorage.getForceCnt(this->interForceTypes.Hbonds_protDonor, lmc);
        auto startPositions = forcesLineStorage.getCurStartPosArray(this->interForceTypes.Hbonds_protDonor, lmc);
        auto directions = forcesLineStorage.getCurDirArray(this->interForceTypes.Hbonds_protDonor, lmc);
        auto lenghts = forcesLineStorage.getCurLengthArray(this->interForceTypes.Hbonds_protDonor, lmc);

        auto col = &this->forceColorsMap[this->interForceTypes.Hbonds_protDonor];
        this->renVarsIFC.baseColor.Set(col->x, col->y, col->z, col->w);
        cylinderRenderer(startPositions, directions, lenghts, size, this->renVarsIFC, viewportStuff);
    }
    //::glEnable(GL_LIGHTING);


    /**********************************************
     *********** H-BOND CONE RENDERING ************
     **********************************************/
    prepareHBond_Cones();

    if (this->renderDataControl.intForces.hydrogenCones) {
        vislib::Array<float> color;
        vislib::Array<float> color2;

        color.Append(255.0f / 255.0f);
        color.Append(114.0f / 255.0f);
        color.Append(114.0f / 255.0f);
        coneRenderer(this->protHBonds_HPos, this->protHBonds_HDirection, this->protHBonds_coneLengths,
            this->protHBonds_coneRadii, this->forceColorsMap[this->interForceTypes.Hbonds_protDonor], this->renVarsSMS,
            viewportStuff);
        color2.Append(167.0f / 255.0f);
        color2.Append(255.0f / 255.0f);
        color2.Append(114.0f / 255.0f);
        coneRenderer(this->modelHBonds_HPos, this->modelHBonds_HDirection, this->modelHBonds_coneLengths,
            this->modelHBonds_coneRadii, this->forceColorsMap[this->interForceTypes.Hbonds_protAcceptor],
            this->renVarsSMS, viewportStuff);
    }
    // circleTextureRenderer(testCone, this->tf.tf_functGroups, this->renVarsSMS, viewportStuff);


    /********************************************
     ********** PROTEIN MESH RENDERING **********
     ********************************************/

    // triangle mesh call
    CallTriMeshData* tmc = this->triMeshData_callerSlot.CallAs<CallTriMeshData>();
    if (tmc != NULL) {
        proteinMeshRenderer(cr3d, tmc, call, this->tf.tf_interactionForce);
    }
    if (this->s.changedcurClusterID) {
        this->renVarsCNS.selectedSphereID = -1;
        this->renVarsFGS.selectedSphereID = -1;
        this->pmdcOut_newData = true;
    }


    /********************************************
     *********** BAR-CHART RENDERING ************
     ********************************************/
    if (this->renderDataControl.barchartGylphe) {
        float liftSize = 10.0f;
        auto c = this->cameraInfo;
        const auto& z = c->Right();
        const auto& clusterCentroids = this->clusterData.clusterCentroids->PeekElements();
        glm::vec3 rightVec = glm::vec3(z.GetX(), z.GetY(), z.GetZ());
        glm::vec3 eyeDirection =
            glm::vec3(c->EyeDirection().GetX(), c->EyeDirection().GetY(), c->EyeDirection().GetZ());
        glm::vec3 CameraUp = vec3(c->Up().GetX(), c->Up().GetY(), c->Up().GetZ());
        glm::vec3 upVec = glm::cross(rightVec, eyeDirection);
        upVec = normalize(upVec);
        vislib::Array<vislib::math::Vector<float, 3>> drawOrder;
        const int ClusterCentroidsCnt = lmc->getClusterAndCentroidData().clusterCentroids->Count() / 4;
        drawOrder.SetCount(ClusterCentroidsCnt);

        for (size_t i = 0; i < ClusterCentroidsCnt; i++) {
            glm::vec3 clusterPos = glm::vec3(clusterCentroids[i * 4 + 0] + upVec.x * liftSize,
                clusterCentroids[i * 4 + 1] + upVec.y * liftSize, clusterCentroids[i * 4 + 2] + upVec.z * liftSize);
            vislib::math::Vector<float, 2> drawOderEntry;
            glm::vec3 dir = this->camPos_local - clusterPos;
            float distance = glm::length(dir) * (-1);
            drawOderEntry.SetX((float)i);
            drawOderEntry.SetY(distance);
            drawOrder[i] = drawOderEntry;
        }
        drawOrder.Sort(sortVecf3ByY);

        vislib::Array<float> sColors;
        sColors.SetCount(6);

        // outline color
        sColors[0] = 0.302f;
        sColors[1] = 0.686f;
        sColors[2] = 0.29f;
        // circle color
        sColors[3] = 1.0f;
        sColors[4] = 1.0f;
        sColors[5] = 1.0f;

        vislib::Array<float> sSpheres;
        sSpheres.SetCount(8);

        for (size_t i = 0; i < ClusterCentroidsCnt; i++) {
            const int ID = (int)drawOrder[i].GetX();
            // outline circle
            sSpheres[0] = clusterCentroids[ID * 4 + 0] + upVec.x * liftSize;
            sSpheres[1] = clusterCentroids[ID * 4 + 1] + upVec.y * liftSize;
            sSpheres[2] = clusterCentroids[ID * 4 + 2] + upVec.z * liftSize;
            sSpheres[3] = 6.5f;
            // circle
            sSpheres[4] = clusterCentroids[ID * 4 + 0] + upVec.x * liftSize;
            sSpheres[5] = clusterCentroids[ID * 4 + 1] + upVec.y * liftSize;
            sSpheres[6] = clusterCentroids[ID * 4 + 2] + upVec.z * liftSize;
            sSpheres[7] = 6.3f;

            glDisable(GL_DEPTH_TEST);
            CircleRenderer(sSpheres, sColors, viewportStuff);
            barChartRenderer(
                this->clusterForceAreas, this->clusterForceResidueCnts, this->clusterAreas, lmc, ID, liftSize);
            glEnable(GL_DEPTH_TEST);
        }
    }


    /********************************************
     ************** TEXT RENDERING **************
     ********************************************/
    // initialize the font object
    if (!this->font.IsInitialised()) {
        // initialize the font object
        this->font.Initialise(this->GetCoreInstance());
        this->font.SetBillboardMode(true);
    }


    /****** render star menu things ******/
    if (this->s.selectedIDx >= 0) {

        // Make sure that the text is the last thing that is rendered, otherwise, it might be covered!
        glDisable(GL_DEPTH_TEST);
        glm::vec3 center = vec3(combined_noiseModelCentroids_clusterCentroids[this->s.selectedIDx * 4 + 0],
            combined_noiseModelCentroids_clusterCentroids[this->s.selectedIDx * 4 + 1],
            combined_noiseModelCentroids_clusterCentroids[this->s.selectedIDx * 4 + 2]);

        glm::vec3 upVec = glm::vec3(cameraInfo->Up().GetX(), cameraInfo->Up().GetY(), cameraInfo->Up().GetZ());
        glm::vec3 rightVec =
            glm::vec3(cameraInfo->Right().GetX(), cameraInfo->Right().GetY(), cameraInfo->Right().GetZ());
        glm::vec3 tmpDir = glm::cross(upVec, rightVec);

        // TODO: fix this! the text should be moved towards the user, but should stay within the bounding box
        //       Clip tmpDir with Bounding Box?
        // tmpDir *= mdc->AccessBoundingBoxes().ObjectSpaceBBox().LongestEdge() / 2.0f;

        glm::vec3 eyeDirection = normalize(tmpDir);
        glm::vec3 CameraUp = vec3(cameraInfo->Up().GetX(), cameraInfo->Up().GetY(), cameraInfo->Up().GetZ());
        glm::vec3 UpVecStarMenu = glm::cross(rightVec, eyeDirection);
        UpVecStarMenu = normalize(UpVecStarMenu);

        // little correction to get the up vector of the star menu
        // UpVecStarMenu = glm::rotate(CameraRight, (float)(-90 * vislib::math::PI_DOUBLE / 180), eyeDirection);

        // if selected sphere is a noise sphere (!= cluster sphere)
        float correction = this->scaleWorld / 0.0117368829;
        float postScale = 1.0f;
        float postsubTextScale = 0.47f * correction;
        float postSVGtextScale = 0.46f * correction;
        if (this->s.selectedIDx < noiseModelCentroidsIDXs.Count()) {
            lmc->SetCurrentLigandAndModel_byGlobalMdlID(
                this->clusterData.assign_modelID_sphereID->PeekElements()[this->s.selectedIDx].GetY());
            (*lmc)(LigandModelCall::CallForGetData);
            float bindergy = lmc->getCurrentBindingenergy();

            vislib::StringA tmpStr;
            tmpStr.Format("%.1f kcal/mol", bindergy);

            float textScaleFunctionValue1 = 0.1937 * powf(sqrt(this->starMenu.textSize * 100 * 6.7), 1.8923);
            float textScaleFunctionValue2 = 0.1937 * powf(sqrt(this->starMenu.textSize * 100 * 4.0), 1.8923);

            utility::SDFFont::Alignment fontAlignment = utility::SDFFont::Alignment::ALIGN_CENTER_BOTTOM;
            if (this->renderDataControl.improveTextReadability) {
                font.SetRenderType(font.RENDERTYPE_OUTLINE);
                font.DrawString(&fontOutline_color.x,
                    center.x - eyeDirection.x + UpVecStarMenu.x * textScaleFunctionValue1,
                    center.y - eyeDirection.y + UpVecStarMenu.y * textScaleFunctionValue1,
                    center.z - eyeDirection.z + UpVecStarMenu.z * textScaleFunctionValue1, this->starMenu.textSize,
                    false, tmpStr.PeekBuffer(), fontAlignment);
            }

            font.SetRenderType(font.RENDERTYPE_FILL);
            font.DrawString(&starMenu.fontColor.x,
                center.x - eyeDirection.x + UpVecStarMenu.x * textScaleFunctionValue1,
                center.y - eyeDirection.y + UpVecStarMenu.y * textScaleFunctionValue1,
                center.z - eyeDirection.z + UpVecStarMenu.z * textScaleFunctionValue1, this->starMenu.textSize, false,
                tmpStr.PeekBuffer(), fontAlignment);

            // single info line
            tmpStr.Format(lmc->getCurrentLigandname());
            if (this->renderDataControl.improveTextReadability) {
                font.SetRenderType(font.RENDERTYPE_OUTLINE);
                font.DrawString(&fontOutline_color.x,
                    center.x - eyeDirection.x + UpVecStarMenu.x * textScaleFunctionValue2,
                    center.y - eyeDirection.y + UpVecStarMenu.y * textScaleFunctionValue2,
                    center.z - eyeDirection.z + UpVecStarMenu.z * textScaleFunctionValue2, this->starMenu.textSize,
                    false, tmpStr.PeekBuffer(), fontAlignment);
            }

            font.SetRenderType(font.RENDERTYPE_FILL);
            font.DrawString(&starMenu.fontColor.x,
                center.x - eyeDirection.x + UpVecStarMenu.x * textScaleFunctionValue2,
                center.y - eyeDirection.y + UpVecStarMenu.y * textScaleFunctionValue2,
                center.z - eyeDirection.z + UpVecStarMenu.z * textScaleFunctionValue2, this->starMenu.textSize, false,
                tmpStr.PeekBuffer(), fontAlignment);

            // if selected sphere is a cluster sphere (!=noise sphere)
            /** TODO: improve performance -> recalculation of star menu only
             * - if selectedIDx is changed
             * - if star menu params are changed
             */
        } else {

            int tmpCluster = this->s.selectedIDx - noiseModelCentroidsIDXs.Count();
            int clusterIDx_start = tmpCluster > 0 ? this->clusterData.clusterIndexSum->PeekElements()[tmpCluster]
                                                  : this->clusterData.clusterIndexSum->PeekElements()[0];
            int clusterIDx_end = this->clusterData.clusterIndexSum->PeekElements()[tmpCluster + 1];
            int clusterSize =
                this->clusterData.clusterSizes->PeekElements()[tmpCluster]; //(clusterIDx_end - clusterIDx_start + 1);

            // sort models of cluster by binding energy
            if (this->s.selectedIDx != this->s.old_SelectedIDx || this->starMenu.changed_curSortingFroceType) {
                // set starMenu.selectedIDx if clusterSphere is changed
                this->starMenu.selectedIDx = -1;
                this->oneClusterData_sortedWithEnergy.SetCount(clusterSize);
                int cnt = 0;
                for (int i = clusterIDx_start; i < clusterIDx_end; i++) {
                    // set query
                    lmc->SetCurrentLigandAndModel_byGlobalMdlID(
                        this->clusterData.assign_modelID_sphereID->PeekElements()[i].GetY());
                    // call for data
                    (*lmc)(LigandModelCall::CallForGetDataSilent);
                    float bindergy = lmc->getCurrentBindingenergy();
                    this->oneClusterData_sortedWithEnergy[cnt].SetX(
                        this->clusterData.assign_modelID_sphereID->PeekElements()[i].GetY());
                    this->oneClusterData_sortedWithEnergy[cnt].SetY(bindergy);
                    int counter = 0;
                    if (this->starMenu.curSortingFroceType == this->interForceTypes.HBonds) {
                        counter = lmc->getInteractionForces().hProtAcceptors->cnt +
                                  lmc->getInteractionForces().hProtDonors->cnt;
                    } else if (this->starMenu.curSortingFroceType == this->interForceTypes.halogenBonds) {
                        counter = lmc->getInteractionForces().halogenBonds->cnt;
                    } else if (this->starMenu.curSortingFroceType == this->interForceTypes.hydrophobicInteractions) {
                        counter = lmc->getInteractionForces().hydrophobicInteractions->cnt;
                    } else if (this->starMenu.curSortingFroceType == this->interForceTypes.metalComplexes) {
                        counter = lmc->getInteractionForces().metalComplexes->cnt;
                    } else if (this->starMenu.curSortingFroceType == this->interForceTypes.piCationInteractions) {
                        counter = lmc->getInteractionForces().piCationInteractions->cnt;
                    } else if (this->starMenu.curSortingFroceType == this->interForceTypes.piStacks) {
                        counter = lmc->getInteractionForces().piStacks->cnt;
                    } else if (this->starMenu.curSortingFroceType == this->interForceTypes.saltBridges) {
                        counter = lmc->getInteractionForces().saltBridges->cnt;
                    }

                    this->oneClusterData_sortedWithEnergy[cnt].SetZ((float)counter);
                    cnt++;
                }
                // TODO: may use THRUST sort for large clusters?
                this->oneClusterData_sortedWithEnergy.Sort(sortVecf3ByYifZisNoneZero);
            }

            utility::SDFFont::Alignment fontAlignment[20];
            vislib::Array<vislib::math::Vector<float, 3>> starDirections;
            starDirections.SetCount(this->starMenu.maxCircleCnt);

            double starMenuRotationAngle = 360 / (double)this->starMenu.circleCnt;
            // iterate over star menu points
            this->lines.Resize(0);
            for (int i = 0; i < this->starMenu.circleCnt; i++) {
                double rotAngle = (starMenuRotationAngle * i);
                /**********************************************************
                 * Handle right text alignment depending on rotation angle
                 **********************************************************/
                /**
                 * TOP [AROUND 0°/360°]
                 */
                if ((rotAngle >= 0 && rotAngle < 45) || (rotAngle > 315 && rotAngle <= 360)) { // +/- 45° von 0°/360°;
                    fontAlignment[i] = utility::SDFFont::Alignment::ALIGN_CENTER_BOTTOM;
                    if (rotAngle >= 22.5 && rotAngle < 45) {
                        fontAlignment[i] = utility::SDFFont::Alignment::ALIGN_LEFT_MIDDLE;
                        // printf("i: %d \t ALIGN_LEFT_MIDDLE, \n", i);
                    } else if (rotAngle >= 315 && rotAngle < 337.5) {
                        fontAlignment[i] = utility::SDFFont::Alignment::ALIGN_RIGHT_MIDDLE;
                        // printf("i: %d \t ALIGN_RIGHT_MIDDLE, \n", i);
                    } else if (rotAngle >= 10 && rotAngle < 22.5) {
                        fontAlignment[i] = utility::SDFFont::Alignment::ALIGN_RIGHT_BOTTOM;
                        // printf("i: %d \t ALIGN_RIGHT_BOTTOM, \n", i);
                    } else if (rotAngle >= 337.5 && rotAngle < 350) {
                        fontAlignment[i] = utility::SDFFont::Alignment::ALIGN_LEFT_BOTTOM;
                        // printf("i: %d \t ALIGN_LEFT_BOTTOM, \n", i);
                    }
                }

                /**
                 * RIGHT [AROUND 90°]
                 */
                if ((rotAngle >= 45 && rotAngle <= 135)) { // +/- 45° of 90°;
                    fontAlignment[i] = utility::SDFFont::Alignment::ALIGN_LEFT_MIDDLE;
                    // printf("i: %d \t ALIGN_LEFT_MIDDLE, \n", i);
                }

                /**
                 * BOTTOM [AROUND 180°]
                 */
                if ((rotAngle > 135 && rotAngle < 225)) { // +/- 45° von 180°;
                    fontAlignment[i] = utility::SDFFont::Alignment::ALIGN_CENTER_TOP;
                    if (rotAngle >= 135 && rotAngle < 157.5) {
                        fontAlignment[i] = utility::SDFFont::Alignment::ALIGN_LEFT_MIDDLE;
                        // printf("i: %d \t ALIGN_RIGHT_MIDDLE, \n", i);
                    } else if (rotAngle >= 202.5 && rotAngle < 225) {
                        fontAlignment[i] = utility::SDFFont::Alignment::ALIGN_RIGHT_MIDDLE;
                        // printf("i: %d \t ALIGN_LEFT_MIDDLE, \n", i);
                    } else if (rotAngle >= 190 && rotAngle < 202.5) {
                        fontAlignment[i] = utility::SDFFont::Alignment::ALIGN_RIGHT_TOP;
                        // printf("i: %d \t ALIGN_RIGHT_TOP, \n", i);
                    } else if (rotAngle >= 157.5 && rotAngle < 170) {
                        fontAlignment[i] = utility::SDFFont::Alignment::ALIGN_LEFT_TOP;
                        // printf("i: %d \t ALIGN_LEFT_TOP, \n", i);
                    }
                }

                /**
                 * LEFT [AROUND 270°]
                 */
                if ((rotAngle >= 225 && rotAngle <= 315)) { // +/- 45° of 270°;
                    fontAlignment[i] = utility::SDFFont::Alignment::ALIGN_RIGHT_MIDDLE;
                    // printf("i: %d \t ALIGN_RIGHT_MIDDLE, \n", i);
                }
                // fontAlignment[i - 1] = utility::SDFFont::Alignment::ALIGN_CENTER_MIDDLE;
                // rotate up vector of star menu to get next text orientation/position
                rotAngle = rotAngle * vislib::math::PI_DOUBLE / 180;
                glm::vec3 rotatedUpVec_StarMenu = glm::rotate(UpVecStarMenu, (float)rotAngle, eyeDirection);
                starDirections[i] = {rotatedUpVec_StarMenu.x, rotatedUpVec_StarMenu.y, rotatedUpVec_StarMenu.z};

                // draw a line of the star menu
                /*
                 float lineWidth = 0.13f;
                 getLineCoords(center, this->starMenuSize, lineWidth, rotatedUpVec_StarMenu, eyeDirection,
                 this->lines, this->starMenuInnerStart);
                 */
            }


            // scaling factor
            float starMenuOuterEnd = this->starMenu.size + 1;


            /************************************************/
            /********* START: STAR MENU PAGE HANDLE *********/
            /************************************************/
            int starMenuPageCnt = clusterSize % this->starMenu.circleCnt == 0
                                      ? (int)(clusterSize / this->starMenu.circleCnt)
                                      : (int)(clusterSize / this->starMenu.circleCnt) + 1;
            // do not run out of bounds
            if (this->starMenu.curPage < 0) {
                this->starMenu.curPage = starMenuPageCnt - 1;
            } else if (this->starMenu.curPage >= starMenuPageCnt) {
                this->starMenu.curPage = 0;
            }

            // get page for cur gblMdlID
            if (this->webDataChanged) {
                (*lmc)(LigandModelCall::CallForGetDataSilent);
                auto webData = lmc->getWebData();
                for (int i = 0; i < this->oneClusterData_sortedWithEnergy.Count(); i++) {
                    if (oneClusterData_sortedWithEnergy[i].GetX() == webData->getGlobalMdlID()) {
                        lmc->SetCurrentLigandAndModel_byGlobalMdlID(webData->getGlobalMdlID());
                        (*lmc)(LigandModelCall::CallForGetDataSilent);
                        if (oneClusterData_sortedWithEnergy[i].GetY() == lmc->getCurrentBindingenergy()) {
                            this->starMenu.curPage = (int)(float(i) / (float)this->starMenu.circleCnt);
                            printf("starMenu.curPage: %d \t i:%d cnt:%d \n", this->starMenu.curPage, i,
                                (float)this->starMenu.circleCnt);
                            this->starMenu.selectedIDx = i - (this->starMenu.curPage * this->starMenu.circleCnt);
                            this->curSelectedStarMenuModelIDx = this->starMenu.selectedIDx;
                        }
                    }
                }
            }


            this->starMenu.curPageStartIDx = this->starMenu.circleCnt * starMenu.curPage;
            int curStarMenuEndClusterIDx = 0;

            if ((this->starMenu.curPageStartIDx + this->starMenu.circleCnt) >= clusterSize) {
                curStarMenuEndClusterIDx = clusterSize;
            } else {
                curStarMenuEndClusterIDx = this->starMenu.curPageStartIDx + this->starMenu.circleCnt;
            }
            // * 4: x,y,z and radius;  (sphere == circle)
            int circleDif = this->starMenu.curPageStartIDx + 1 - curStarMenuEndClusterIDx;
            if (this->starCircles.Count() == 0 || this->starMenu.circleCnt != this->starMenu.old_circleCnt ||
                this->starMenu.curPage != this->starMenu.old_curPage || this->s.changedcurClusterID) {

                int allocSize = (this->starMenu.curPageStartIDx + this->starMenu.circleCnt) > clusterSize
                                    ? (clusterSize - this->starMenu.curPageStartIDx)
                                    : this->starMenu.circleCnt;

                this->starCircles.Resize(0);
                this->starSubSpheres_1_orderNumber.Resize(0);
                this->starSubSpheres_2_dockScore.Resize(0);
                this->starSubSpheres_3_Hbond.Resize(0);
                this->starSubSpheres_4_halogenBond.Resize(0);
                this->starSubSpheres_5_hydrophobicInters.Resize(0);
                this->starSubSpheres_6_metalComplexes.Resize(0);
                this->starSubSpheres_7_piCationInters.Resize(0);
                this->starSubSpheres_8_piSticks.Resize(0);
                this->starSubSpheres_9_saltBridges.Resize(0);
                this->starSubSpheres_onSelect.Resize(0);

                this->starCircles.SetCount((allocSize + 1) * 4);
                this->starSubSpheres_1_orderNumber.SetCount(allocSize * 4);
                this->starSubSpheres_2_dockScore.SetCount(allocSize * 4);
                this->starSubSpheres_3_Hbond.SetCount(allocSize * 4);
                this->starSubSpheres_4_halogenBond.SetCount(allocSize * 4);
                this->starSubSpheres_5_hydrophobicInters.SetCount(allocSize * 4);
                this->starSubSpheres_6_metalComplexes.SetCount(allocSize * 4);
                this->starSubSpheres_7_piCationInters.SetCount(allocSize * 4);
                this->starSubSpheres_8_piSticks.SetCount(allocSize * 4);
                this->starSubSpheres_9_saltBridges.SetCount(allocSize * 4);

                subCirclePointers.clear();
                this->subCirclePointers.push_back(&starSubSpheres_1_orderNumber[0]);
                this->subCirclePointers.push_back(&starSubSpheres_2_dockScore[0]);
                this->subCirclePointers.push_back(&starSubSpheres_3_Hbond[0]);
                this->subCirclePointers.push_back(&starSubSpheres_4_halogenBond[0]);
                this->subCirclePointers.push_back(&starSubSpheres_5_hydrophobicInters[0]);
                this->subCirclePointers.push_back(&starSubSpheres_6_metalComplexes[0]);
                this->subCirclePointers.push_back(&starSubSpheres_7_piCationInters[0]);
                this->subCirclePointers.push_back(&starSubSpheres_8_piSticks[0]);
                this->subCirclePointers.push_back(&starSubSpheres_9_saltBridges[0]);
            }

            // change content for subCircles
            if (this->starMenu.curSubCircleContent.size() == 0 ||
                this->starMenu.circleCnt != this->starMenu.old_circleCnt ||
                this->starMenu.curPage != this->starMenu.old_curPage || this->starMenu.changed_curSortingFroceType ||
                this->s.changedcurClusterID) {

                this->starMenu.old_circleCnt = this->starMenu.circleCnt;
                this->starMenu.old_curPage = this->starMenu.curPage;

                auto cont = &this->starMenu.curSubCircleContent;

                cont->clear();
                cont->resize(this->starMenu.circleCnt);
                vislib::StringA stringCont;
                for (size_t x = 0; x < cont->size(); x++) {
                    cont[0][x].resize(this->subCircleContentIDxes.size() + 1);
                    int starIDx = this->starMenu.curPageStartIDx;
                    int gblMdlID = oneClusterData_sortedWithEnergy[starIDx + x].GetX();

                    lmc->SetCurrentLigandAndModel_byGlobalMdlID(gblMdlID);
                    (*lmc)(LigandModelCall::CallForGetDataSilent);
                    auto lmcCont = lmc->getInteractionForces();
                    cont[0][x][this->subCircleContentIDxes["orderNum"]] = std::to_string(starIDx + x + 1).c_str();
                    vislib::StringA tmpString;
                    tmpString.Format("%.1f", lmc->getCurrentBindingenergy());
                    cont[0][x][this->subCircleContentIDxes["score"]] = tmpString;
                    cont[0][x][this->subCircleContentIDxes[this->interForceTypes.HBonds]] =
                        (lmcCont.hProtAcceptors->cnt + lmcCont.hProtDonors->cnt) > 0 ? "H" : "";
                    cont[0][x][this->subCircleContentIDxes[this->interForceTypes.halogenBonds]] =
                        lmcCont.halogenBonds->cnt > 0 ? "HA" : "";
                    cont[0][x][this->subCircleContentIDxes[this->interForceTypes.hydrophobicInteractions]] =
                        lmcCont.hydrophobicInteractions->cnt > 0 ? "HydI" : "";
                    cont[0][x][this->subCircleContentIDxes[this->interForceTypes.metalComplexes]] =
                        lmcCont.metalComplexes->cnt > 0 ? "M" : "";
                    cont[0][x][this->subCircleContentIDxes[this->interForceTypes.piCationInteractions]] =
                        lmcCont.piCationInteractions->cnt > 0 ? "pi+" : "";
                    cont[0][x][this->subCircleContentIDxes[this->interForceTypes.piStacks]] =
                        lmcCont.piStacks->cnt > 0 ? "pi" : "";
                    cont[0][x][this->subCircleContentIDxes[this->interForceTypes.saltBridges]] =
                        lmcCont.saltBridges->cnt > 0 ? "S" : "";
                }
            }


            //****************DONT TOUCH ANYTHING HERE (UNDER BLOOD AND SWEAT FOUND MAGIC)******************//


            /**
             * page infos (info line)
             */
            float pageStrScale = this->starMenu.size + 6.5;

            vislib::StringA pageStr;
            pageStr.Format("(%d-%d)/%d page %d/%d", this->starMenu.curPageStartIDx + 1, curStarMenuEndClusterIDx,
                clusterSize, this->starMenu.curPage + 1, starMenuPageCnt);
            float textScaleFunctionValue1 = 0.1937 * pow(sqrt(this->starMenu.textSize * 100 * 8.7), 1.8923) * postScale;


            if (this->renderDataControl.improveTextReadability) {
                font.SetRenderType(font.RENDERTYPE_OUTLINE);
                this->font.DrawString(&fontOutline_color.x,
                    center.x + starDirections[0].GetX() * pageStrScale +
                        (starDirections[0].GetX() * textScaleFunctionValue1),
                    center.y + starDirections[0].GetY() * pageStrScale +
                        (starDirections[0].GetY() * textScaleFunctionValue1),
                    center.z + starDirections[0].GetZ() * pageStrScale +
                        (starDirections[0].GetZ() * textScaleFunctionValue1),
                    this->starMenu.textSize * postsubTextScale, false, pageStr.PeekBuffer(), fontAlignment[0]);
            }

            font.SetRenderType(font.RENDERTYPE_FILL);
            this->font.DrawString(&starMenu.fontColor.x,
                center.x + starDirections[0].GetX() * pageStrScale +
                    (starDirections[0].GetX() * textScaleFunctionValue1),
                center.y + starDirections[0].GetY() * pageStrScale +
                    (starDirections[0].GetY() * textScaleFunctionValue1),
                center.z + starDirections[0].GetZ() * pageStrScale +
                    (starDirections[0].GetZ() * textScaleFunctionValue1),
                this->starMenu.textSize * postsubTextScale, false, pageStr.PeekBuffer(), fontAlignment[0]);

            // ligandName
            if (this->curSelectedStarMenuModelIDx != -1) {
                pageStr.Format(lmc->getCurrentLigandname());
                if (this->renderDataControl.improveTextReadability) {
                    font.SetRenderType(font.RENDERTYPE_OUTLINE);
                    this->font.DrawString(&fontOutline_color.x,
                        center.x - starDirections[0].GetX() * pageStrScale -
                            (starDirections[0].GetX() * textScaleFunctionValue1 * 1.1),
                        center.y - starDirections[0].GetY() * pageStrScale -
                            (starDirections[0].GetY() * textScaleFunctionValue1 * 1.1),
                        center.z - starDirections[0].GetZ() * pageStrScale -
                            (starDirections[0].GetZ() * textScaleFunctionValue1 * 1.1),
                        this->starMenu.textSize * postsubTextScale, false, pageStr.PeekBuffer(), fontAlignment[0]);
                }
                font.SetRenderType(font.RENDERTYPE_FILL);
                this->font.DrawString(&starMenu.fontColor.x,
                    center.x - starDirections[0].GetX() * pageStrScale -
                        (starDirections[0].GetX() * textScaleFunctionValue1 * 1.1),
                    center.y - starDirections[0].GetY() * pageStrScale -
                        (starDirections[0].GetY() * textScaleFunctionValue1 * 1.1),
                    center.z - starDirections[0].GetZ() * pageStrScale -
                        (starDirections[0].GetZ() * textScaleFunctionValue1 * 1.1),
                    this->starMenu.textSize * postsubTextScale, false, pageStr.PeekBuffer(), fontAlignment[0]);
            }

            /**
             * arrow elements (clicking through starMenu pages)
             */

            // get orthogonal vector for up direction to  +/- move elements horizontal
            glm::vec3 tmpOrtoVec =
                glm::rotate(glm::vec3(starDirections[0].GetX(), starDirections[0].GetY(), starDirections[0].GetZ()),
                    (float)(-90 * vislib::math::PI_DOUBLE / 180), eyeDirection);

            // set with textSize scaling distance between the arrow elements
            tmpOrtoVec *= (0.1937 * powf(sqrt(this->starMenu.textSize * 10 * 5), 1.8923)) * 2.2;

            float arrowDis = DBSCAN::distance2Vectors(
                {tmpOrtoVec.x, tmpOrtoVec.y, tmpOrtoVec.z}, {-tmpOrtoVec.x, -tmpOrtoVec.y, -tmpOrtoVec.z});
            arrowSpheres.SetCount(2 * 4);

            float textScaleFunctionValue2 = 0.1937 * pow(sqrt(this->starMenu.textSize * 100 * 8), 1.8923) * postScale;
            //
            arrowSpheres[0 * 4 + 0] = center.x + starDirections[0].GetX() * pageStrScale +
                                      (starDirections[0].GetX() * textScaleFunctionValue2) + tmpOrtoVec.x;
            arrowSpheres[0 * 4 + 1] = center.y + starDirections[0].GetY() * pageStrScale +
                                      (starDirections[0].GetY() * textScaleFunctionValue2) + tmpOrtoVec.y;
            arrowSpheres[0 * 4 + 2] = center.z + starDirections[0].GetZ() * pageStrScale +
                                      (starDirections[0].GetZ() * textScaleFunctionValue2) + tmpOrtoVec.z;
            arrowSpheres[0 * 4 + 3] = arrowDis / 2;

            arrowSpheres[1 * 4 + 0] = center.x + starDirections[0].GetX() * pageStrScale +
                                      (starDirections[0].GetX() * textScaleFunctionValue2) - tmpOrtoVec.x;
            arrowSpheres[1 * 4 + 1] = center.y + starDirections[0].GetY() * pageStrScale +
                                      (starDirections[0].GetY() * textScaleFunctionValue2) - tmpOrtoVec.y;
            arrowSpheres[1 * 4 + 2] = center.z + starDirections[0].GetZ() * pageStrScale +
                                      (starDirections[0].GetZ() * textScaleFunctionValue2) - tmpOrtoVec.z;
            arrowSpheres[1 * 4 + 3] = arrowDis / 2;
            //
            float arrowCol[] = {starMenu.fontColor.x, starMenu.fontColor.y, starMenu.fontColor.z, starMenu.fontColor.w};
            float arrowCol2[] = {
                starMenu.fontColor.x, starMenu.fontColor.y, starMenu.fontColor.z, starMenu.fontColor.w};

            // check for click on arrow elements
            this->arrowUpDown = mouseRayIntersection.getFirstMouseRayIntersection(
                this->arrowSpheres, this->cam, this->mouseX, this->mouseY, this->isClicked);


            // handle color for highlighting the hovered arrow element
            if (arrowUpDown > -1) {
                if (arrowUpDown == 0) {
                    arrowCol2[0] = starMenu.fontColor.x;
                    arrowCol2[1] = starMenu.fontColor.y;
                    arrowCol2[2] = starMenu.fontColor.z;
                    arrowCol2[3] = starMenu.fontColor.w;

                    arrowCol[0] = 0.0f;
                    arrowCol[1] = 1.0f;
                    arrowCol[2] = 0.0f;
                    arrowCol[3] = 1.0f;
                } else {
                    arrowCol[0] = starMenu.fontColor.x;
                    arrowCol[1] = starMenu.fontColor.y;
                    arrowCol[2] = starMenu.fontColor.z;
                    arrowCol[3] = starMenu.fontColor.w;

                    arrowCol2[0] = 0.0f;
                    arrowCol2[1] = 1.0f;
                    arrowCol2[2] = 0.0f;
                    arrowCol2[3] = 1.0f;
                }
                if (this->isClicked == 1) {
                    arrowUpDown == 0 ? this->starMenu.curPage-- : this->starMenu.curPage++;
                    // set selection back if page is changed
                    this->starMenu.selectedIDx = -1;
                }
            }

            // do not run out of bounds
            if (this->starMenu.curPage < 0) {
                this->starMenu.curPage = starMenuPageCnt - 1;
            } else if (this->starMenu.curPage >= starMenuPageCnt) {
                this->starMenu.curPage = 0;
            }

            // draw the arrow elements
            vislib::StringA pageStrArrowDown = "<";
            vislib::StringA pageStrArrowUp = ">";

            font.SetRenderType(font.RENDERTYPE_OUTLINE);
            this->font.DrawString(&fontOutline_color.x, arrowSpheres[0 * 4 + 0], arrowSpheres[0 * 4 + 1],
                arrowSpheres[0 * 4 + 2], this->starMenu.textSize * postsubTextScale, false,
                pageStrArrowDown.PeekBuffer(), utility::SDFFont::Alignment::ALIGN_CENTER_MIDDLE);
            font.SetRenderType(font.RENDERTYPE_FILL);
            this->font.DrawString(arrowCol, arrowSpheres[0 * 4 + 0], arrowSpheres[0 * 4 + 1], arrowSpheres[0 * 4 + 2],
                this->starMenu.textSize * postsubTextScale, false, pageStrArrowDown.PeekBuffer(),
                utility::SDFFont::Alignment::ALIGN_CENTER_MIDDLE);

            font.SetRenderType(font.RENDERTYPE_OUTLINE);
            this->font.DrawString(&fontOutline_color.x, arrowSpheres[1 * 4 + 0], arrowSpheres[1 * 4 + 1],
                arrowSpheres[1 * 4 + 2], this->starMenu.textSize * postsubTextScale, false, pageStrArrowUp.PeekBuffer(),
                utility::SDFFont::Alignment::ALIGN_CENTER_MIDDLE);
            font.SetRenderType(font.RENDERTYPE_FILL);
            this->font.DrawString(arrowCol2, arrowSpheres[1 * 4 + 0], arrowSpheres[1 * 4 + 1], arrowSpheres[1 * 4 + 2],
                this->starMenu.textSize * postsubTextScale, false, pageStrArrowUp.PeekBuffer(),
                utility::SDFFont::Alignment::ALIGN_CENTER_MIDDLE);


            /*********************************************/
            /********* END: STAR MENU PAGE HANDLE ********/
            /*********************************************/
            /**
             * draw star lines and text
             */
            vislib::math::Vector<float, 3> vec1;
            vislib::math::Vector<float, 3> vec2;
            /**
             *  Draw the star menu
             */
            vec1 = {center.x + starDirections[0].GetX() * (starMenuOuterEnd + 2),
                center.y + starDirections[0].GetY() * (starMenuOuterEnd + 2),
                center.z + starDirections[0].GetZ() * (starMenuOuterEnd + 2)};
            vec2 = {center.x + starDirections[1].GetX() * (starMenuOuterEnd + 2),
                center.y + starDirections[1].GetY() * (starMenuOuterEnd + 2),
                center.z + starDirections[1].GetZ() * (starMenuOuterEnd + 2)};

            // starSpheres scaling
            float tmp = 0.1937 * pow(sqrt(this->starMenu.textSize * 100 * 6.5), 1.8923) * postScale;
            float gDistance;
            (float)tmp * 0.9 > (float)(DBSCAN::distance2Vectors(vec1, vec2) / 2)* 0.95
                ? gDistance = DBSCAN::distance2Vectors(vec1, vec2) / 2 * 0.95
                : gDistance = tmp * 0.9;

            for (int i = 0; i < this->starMenu.circleCnt; i++) {
                if (this->starMenu.curPageStartIDx + i < clusterSize) {
                    this->indexOfLastStarCircle = i;

                    if (this->starMenu.old_selectedIDx != this->starMenu.selectedIDx) {
                        this->starMenu.old_selectedIDx = this->starMenu.selectedIDx;
                        this->lmdcOut_newData = true;
                    }

                    int subCircleCounter = 1;
                    // main sphere
                    starCircles[i * 4 + 0] = center.x + starDirections[i].GetX() * (starMenuOuterEnd + 2);
                    starCircles[i * 4 + 1] = center.y + starDirections[i].GetY() * (starMenuOuterEnd + 2);
                    starCircles[i * 4 + 2] = center.z + starDirections[i].GetZ() * (starMenuOuterEnd + 2);
                    starCircles[i * 4 + 3] = gDistance;


                    // sub-spheres right
                    float r, b, alpha;
                    glm::vec3 tmpVec;
                    for (int j = 0; j < this->subCirclePointers.size(); j++) {
                        auto subS = this->subCirclePointers[j];

                        // bogenStrecke b, winkel alpha, radius r
                        // alpha = b/r
                        r = gDistance;
                        b = gDistance * subCircleSizefactor * 1.8f * j;
                        alpha = b / r;
                        tmpVec = glm::rotate(UpVecStarMenu, alpha, eyeDirection);
                        tmpVec = normalize(tmpVec);
                        subS[i * 4 + 0] = starCircles[i * 4 + 0] + tmpVec.x * (1.02 + subCircleSizefactor) * gDistance;
                        subS[i * 4 + 1] = starCircles[i * 4 + 1] + tmpVec.y * (1.02 + subCircleSizefactor) * gDistance;
                        subS[i * 4 + 2] = starCircles[i * 4 + 2] + tmpVec.z * (1.02 + subCircleSizefactor) * gDistance;
                        subS[i * 4 + 3] = gDistance * subCircleSizefactor;
                    }


                    if (this->starMenu.selectedIDx == i) {
                        starSubSpheres_onSelect.SetCount(3 * 4);
                        for (int j = 1; j <= 3; j++) {
                            r = gDistance;
                            b = gDistance * subCircleSizefactor * -1.8f * j;
                            alpha = b / r;
                            tmpVec = glm::rotate(UpVecStarMenu, alpha, eyeDirection);
                            tmpVec = normalize(tmpVec);
                            int jj = j - 1;
                            starSubSpheres_onSelect[jj * 4 + 0] =
                                starCircles[i * 4 + 0] + tmpVec.x * (1.02 + subCircleSizefactor) * gDistance;
                            starSubSpheres_onSelect[jj * 4 + 1] =
                                starCircles[i * 4 + 1] + tmpVec.y * (1.02 + subCircleSizefactor) * gDistance;
                            starSubSpheres_onSelect[jj * 4 + 2] =
                                starCircles[i * 4 + 2] + tmpVec.z * (1.02 + subCircleSizefactor) * gDistance;
                            starSubSpheres_onSelect[jj * 4 + 3] = gDistance * subCircleSizefactor;
                        }
                    }
                } // end if (< clusterSize)
            }     // end for (<starMenu.circleCnt)


            // add one more sphere/circle as click-close menu element
            glm::vec3 tmpVec = glm::rotate(UpVecStarMenu, (float)(90 * vislib::math::PI_DOUBLE / 180), eyeDirection);
            starCircles[(indexOfLastStarCircle + 1) * 4 + 0] =
                center.x + tmpVec.x * (starMenuOuterEnd + 2 + gDistance * 0.5) +
                (starDirections[0].GetX() * (starMenuOuterEnd + 2 + gDistance * 0.65));
            starCircles[(indexOfLastStarCircle + 1) * 4 + 1] =
                center.y + tmpVec.y * (starMenuOuterEnd + 2 + gDistance * 0.5) +
                (starDirections[0].GetY() * (starMenuOuterEnd + 2 + gDistance * 0.65));
            starCircles[(indexOfLastStarCircle + 1) * 4 + 2] =
                center.z + tmpVec.z * (starMenuOuterEnd + 2 + gDistance * 0.5) +
                (starDirections[0].GetZ() * (starMenuOuterEnd + 2 + gDistance * 0.65));
            starCircles[(indexOfLastStarCircle + 1) * 4 + 3] = gDistance * 0.3;

            /*************************************
             ************ Render SVGs ************
             *************************************/

            this->SVGLineCoords.Resize(0);
            this->SVGLineColors.Resize(0);
            this->SVGPolygons.Resize(0);
            this->SVGTextContent.Resize(0);
            this->SVGTextSizes.Resize(0);
            this->SVGTextCoords.Resize(0);
            this->SVGTextColors.Resize(0);
            this->SVGScales.Resize(0);
            float svgScale;
            // this->SVGcnt = 0;
            for (int i = 0; i < this->starMenu.circleCnt; i++) {
                if (this->starMenu.curPageStartIDx + i < clusterSize) {
                    vec3 a, b, direction;
                    float length;

                    vec3 centerX = vec3(starCircles[i * 4 + 0], starCircles[i * 4 + 1], starCircles[i * 4 + 2]);
                    vec3 direct = centerX - center;
                    direct = normalize(direct);
                    // centerX = centerX + direct *25.0f;

                    auto oldLig = lmc->getCurrentLigandID();
                    auto oldMdl = lmc->getCurrentModelID();

                    lmc->SetCurrentLigandAndModel_byGlobalMdlID(
                        oneClusterData_sortedWithEnergy[starMenu.curPageStartIDx + i].GetX());
                    if (!(*lmc)(LigandModelCall::CallForGetDataSilent)) return false;
                    int ID = lmc->getCurrentLigandID();


// check if SVG already parsed (if not load it)
#ifdef SVGs_onDemand
                    if (this->int2String.find(ID) == int2String.end()) {
                        this->int2String.insert({ID, this->SVGcnt});
                        int mapID = int2String[ID];
                        SVG1 a;
                        this->svg1.push_back(a);
                        std::string SVGPath = lmc->getSVGPath();
                        printf("SVGPath: %s\n", SVGPath.c_str());
                        parseSVG(lmc->getSVGPath(), &svg1[this->SVGcnt]);
                        this->SVGcnt++;
                    }
#endif // SVGs_onDemand

                    int mapID = int2String[ID];


                    // determine the svg scale factor to fit svg into the circles
                    float xRange = this->svg1[mapID].svgMAX_xy.GetX() - this->svg1[mapID].svgMIN_xy.GetX();
                    float yRange = this->svg1[mapID].svgMAX_xy.GetY() - this->svg1[mapID].svgMIN_xy.GetY();
                    svgScale = (2 * gDistance) / xRange;
                    if (svgScale > (2 * gDistance) / yRange) {
                        svgScale = (2 * gDistance) / yRange;
                    }
                    svgScale = svgScale * 0.80f;
                    auto middle = this->svg1[mapID].svgCentroid;
                    // svgScale = float(((int)svgScale)) / 100.f;
                    // svgScale = 1 / 10.f;


                    /******************
                     **** SVG TEXT ****
                     ******************/
                    for (int i = 0; i < this->svg1[mapID].textCnt; i++) {

                        auto textElements = this->svg1[mapID].svg_textElements;

                        a = vec3((textElements[i].x - middle.GetX()) * svgScale,
                            (textElements[i].y - middle.GetY()) * svgScale, centerX.z);

                        float factorB = 1.03f;
                        float factor = 1.00f; // move to left with lower as 1.0

                        a = centerX + rightVec * a.x * factor - upVec * a.y * factorB;

                        vislib::StringA tmpStr;
                        tmpStr.Format(this->svg1[mapID].svg_textElements[i].text.c_str());
                        this->SVGTextContent.Append(tmpStr);
                        this->SVGTextSizes.Append(this->svg1[mapID].svg_textElements[i].fontSize * postSVGtextScale);

                        vislib::math::Vector<float, 4> fontColor(
                            this->svg1[mapID].svg_textElements[i].color.GetX() / 255,
                            this->svg1[mapID].svg_textElements[i].color.GetY() / 255,
                            this->svg1[mapID].svg_textElements[i].color.GetZ() / 255, 1.0f);
                        this->SVGTextColors.Append(fontColor);
                        this->SVGTextCoords.Append({a.x, a.y, a.z});
                        this->SVGScales.Append(svgScale);
                    }

                    /******************
                     **** SVG LINE ****
                     ******************/
                    for (int i = 0; i < this->svg1[mapID].lineCnt; i++) {
                        auto lineElements = this->svg1[mapID].svg_lineElements;


                        a = vec3(lineElements[i].x1 - middle.GetX(), lineElements[i].y1 - middle.GetY(), centerX.z);

                        vec3 lineStart = vec3((lineElements[i].x1 - middle.GetX()) * svgScale,
                            (lineElements[i].y1 - middle.GetY()) * svgScale, centerX.z);

                        b = vec3(lineElements[i].x2 - middle.GetX(), lineElements[i].y2 - middle.GetY(), centerX.z);

                        a = centerX + rightVec * a.x - upVec * a.y;
                        b = centerX + rightVec * b.x - upVec * b.y,
                        lineStart = centerX + rightVec * lineStart.x - upVec * lineStart.y,

                        direction = b - a;
                        length = glm::distance(b, a) * svgScale;
                        length = length * 1.01f;

                        getLineCoords(
                            lineStart, length, svgScale * 0.7f, direction, eyeDirection, this->SVGLineCoords, 0);
                        vislib::math::Vector<float, 4> color(this->svg1[mapID].svg_lineElements[i].color.GetX() / 255,
                            this->svg1[mapID].svg_lineElements[i].color.GetY() / 255,
                            this->svg1[mapID].svg_lineElements[i].color.GetZ() / 255,
                            this->svg1[mapID].svg_lineElements[i].opacity);
                        this->SVGLineColors.Append(color);
                    }

                    /*********************
                     **** SVG POLYGON ****
                     *********************/
                    for (int i = 0; i < this->svg1[mapID].polygonCnt; i++) {
                        // for (int i = 0; i < 0; i++) {

                        auto polygon = this->svg1[mapID].svg_polygons[i].points;

                        a = vec3(polygon[0].GetX() - middle.GetX(), polygon[0].GetY() - middle.GetY(), centerX.z);
                        vec3 a_real = vec3((polygon[0].GetX() - middle.GetX()) * svgScale,
                            (polygon[0].GetY() - middle.GetY()) * svgScale, centerX.z);
                        // move in x,y --> rightVec, upVec direction
                        a = centerX + rightVec * a.x - upVec * a.y;
                        a_real = centerX + rightVec * a_real.x - upVec * a_real.y;

                        this->SVGPolygons.Append({a_real.x, a_real.y, a_real.z});


                        for (int j = 1; j < this->svg1[mapID].svg_polygons[i].points.Count(); j++) {

                            b = vec3(polygon[j].GetX() - middle.GetX(), polygon[j].GetY() - middle.GetY(), centerX.z);

                            // move in x,y --> rightVec, upVec direction
                            b = centerX + rightVec * b.x - upVec * b.y,

                            // get direction vector
                                direction = b - a;

                            // determine length between a,b
                            length = glm::distance(b, a) * svgScale;
                            length = length * 1.01f;

                            //
                            vec3 lineDirection = normalize(direction);
                            eyeDirection = normalize(eyeDirection);
                            vec3 rotVec = glm::cross(lineDirection, eyeDirection);
                            float angle = glm::acos(glm::dot(lineDirection, eyeDirection)) -
                                          (float)(90 * vislib::math::PI_DOUBLE / 180);
                            lineDirection = glm::rotate(lineDirection, angle, rotVec);

                            vec3 b_real = a_real + (lineDirection * length);
                            this->SVGPolygons.Append({b_real.x, b_real.y, b_real.z});
                        }
                    }
                    lmc->SetCurrentLigandAndModel(oldLig, oldMdl);
                    if (!(*lmc)(LigandModelCall::CallForGetDataSilent)) return false;
                }
            }

            // do not check star spheres if no cluster is selected
            if (s.curClusterID > -1) {
                this->starMenu.hoverIDx = mouseRayIntersection.getFirstMouseRayIntersection(
                    this->starCircles, this->cam, this->mouseX, this->mouseY, this->isClicked);
                if (this->isClicked == 1 && this->starMenu.hoverIDx > -1) {
                    this->starMenu.selectedIDx = this->starMenu.hoverIDx;
                }
            }


            // glLightfv(GL_LIGHT0, GL_POSITION, value_ptr(camPos_local));
            if (starCircles.Count() > 1) {

                vislib::Array<float> starSphereColors;
                vislib::Array<float> starSubSphereColors;

                int sphereCnt = (this->starCircles.Count() / 4);
                starSphereColors.SetCount(sphereCnt * 3);
                starSubSphereColors.SetCount(sphereCnt * 3);

                for (int i = 0; i < sphereCnt - 1; i++) {

                    // TODO: make this more efficient and dynamic!
                    if (this->starMenu.hoverIDx == i && this->starMenu.hoverIDx != this->starMenu.selectedIDx) {
                        starSphereColors[i * 3 + 0] = 173.f / 255.f;
                        starSphereColors[i * 3 + 1] = 190.f / 255.f;
                        starSphereColors[i * 3 + 2] = 247.f / 255.f;
                    } else if (this->starMenu.selectedIDx == i) {
                        starSphereColors[i * 3 + 0] = 197.f / 255.f;
                        starSphereColors[i * 3 + 1] = 227.f / 255.f;
                        starSphereColors[i * 3 + 2] = 217.f / 255.f;
                    } else {
                        starSphereColors[i * 3 + 0] = 1.f;
                        starSphereColors[i * 3 + 1] = 1.f;
                        starSphereColors[i * 3 + 2] = 1.f;
                    }
                    // first sub-sphere
                    starSubSphereColors[i * 3 + 0] = 0.85f;
                    starSubSphereColors[i * 3 + 1] = 0.85f;
                    starSubSphereColors[i * 3 + 2] = 0.85f;
                }
                // set color of exit button/circle
                starSphereColors[(sphereCnt - 1) * 3 + 0] = this->starMenu.fontColor.x;
                starSphereColors[(sphereCnt - 1) * 3 + 1] = this->starMenu.fontColor.y;
                starSphereColors[(sphereCnt - 1) * 3 + 2] = this->starMenu.fontColor.z;

                // *** Render Circles ***
                CircleRenderer(starCircles, starSphereColors, &viewportStuff[0]);
                CircleRenderer(starSubSpheres_1_orderNumber, starSubSphereColors, &viewportStuff[0]);
                CircleRenderer(starSubSpheres_2_dockScore, starSubSphereColors, &viewportStuff[0]);
                CircleRenderer(starSubSpheres_3_Hbond, starSubSphereColors, &viewportStuff[0]);
                CircleRenderer(starSubSpheres_4_halogenBond, starSubSphereColors, &viewportStuff[0]);
                CircleRenderer(starSubSpheres_5_hydrophobicInters, starSubSphereColors, &viewportStuff[0]);
                CircleRenderer(starSubSpheres_6_metalComplexes, starSubSphereColors, &viewportStuff[0]);
                CircleRenderer(starSubSpheres_7_piCationInters, starSubSphereColors, &viewportStuff[0]);
                CircleRenderer(starSubSpheres_8_piSticks, starSubSphereColors, &viewportStuff[0]);
                CircleRenderer(starSubSpheres_9_saltBridges, starSubSphereColors, &viewportStuff[0]);

                // handle model rendering mode (single/allInCurPocket/allOnProtein)
                if (starSubSpheres_onSelect.Count() > 0) {
                    this->renVarsSMS.baseColor = {0.278f, 0.537f, 0.988f, 1.0};
                    this->renVarsSMS.selectedSphereID =
                        this->renVarsSMS.selectedSphereID == -1 ? 0 : this->renVarsSMS.selectedSphereID;
                    int sel = mouseRayIntersection.getFirstMouseRayIntersection(
                        this->starSubSpheres_onSelect, this->cam, this->mouseX, this->mouseY, this->isClicked);
                    if (sel != -1 && this->isClicked == 1) {
                        this->renVarsSMS.selectedSphereID = sel;
                        this->modelRenderModes.curMode = sel;
                        this->lmdcOut_newData = true;
                    }
                    if (sel != -1 && this->renVarsSMS.hoverSphereID != sel) {
                        this->renVarsSMS.hoverSphereID = sel;
                        this->renVarsSMS.hideOnHover = true;
                        this->renVarsSMS.hoverColor = this->renVarsSMS.baseColor;
                        if (this->renVarsSMS.selectedSphereID == sel) {
                            this->renVarsSMS.hoverColor = this->renVarsSMS.selectionColor;
                        }else{
                            this->renVarsSMS.hoverColor = this->renVarsSMS.baseColor;
                        }
                    } else if (this->renVarsSMS.hoverSphereID == sel) {
                        if (this->renVarsSMS.selectedSphereID == sel) {
                            this->renVarsSMS.hoverColor = this->renVarsSMS.selectionColor;
                        }
                    }
                     else {
                        this->renVarsSMS.hoverSphereID = -1;
                    }
                    // is rendered later because render-order plays a role (zoomed textured circles will have text in front)
                    //circleTextureRenderer( this->starSubSpheres_onSelect, this->tf.tf_functGroups, this->renVarsSMS, viewportStuff);
                }
            }

            ////////////////////////////
            ////////////////////////////
            ////////////////////////////
            ////////////////////////////
            ////////////////////////////
            ////////////////////////////
            ////////////////////////////

            glDisable(GL_LIGHTING);
            glEnable(GL_BLEND);
            glEnable(GL_POLYGON_SMOOTH);
            glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
            glHint(GL_POLYGON_SMOOTH, GL_NICEST);
            glColor4f(0.f, 0.f, 0.f, 1.0);
            glPolygonMode(GL_POLYGON, GL_LINE);
            // draw SVG ploygons (in this case trianlges)
            for (int i = 0; i < this->SVGPolygons.Count() / 3; i++) {
                glBegin(GL_POLYGON);
                // glColor4f(SVGLineColors[i].X(), SVGLineColors[i].Y(), SVGLineColors[i].Z(),
                // SVGLineColors[i].W());
                for (int j = 0; j < 3; j++) {
                    vec3 test = vec3(this->SVGPolygons[i * 3 + j].X(), this->SVGPolygons[i * 3 + j].Y(),
                        this->SVGPolygons[i * 3 + j].Z());
                    glVertex3f(this->SVGPolygons[i * 3 + j].X(), this->SVGPolygons[i * 3 + j].Y(),
                        this->SVGPolygons[i * 3 + j].Z());
                }
                glEnd();
            }

            // draw SVG lines
            for (int i = 0; i < SVGLineCoords.Count() / 4; i++) {

                glBegin(GL_TRIANGLE_STRIP);
                glColor4f(SVGLineColors[i].X(), SVGLineColors[i].Y(), SVGLineColors[i].Z(), SVGLineColors[i].W());
                for (int j = 0; j < 4; j++) {
                    glVertex3f(
                        SVGLineCoords[i * 4 + j].X(), SVGLineCoords[i * 4 + j].Y(), SVGLineCoords[i * 4 + j].Z());
                }
                glEnd();
            }
            glEnable(GL_LIGHTING);
            glDisable(GL_POLYGON_SMOOTH);
            glDisable(GL_BLEND);

            // draw SVG text
            for (int i = 0; i < this->SVGTextColors.Count(); i++) {

                auto col = this->SVGTextColors[i];

                float fontColor[4] = {col.GetX(), col.GetY(), col.GetZ(), col.GetW()};
                vec3 a =
                    vec3(this->SVGTextCoords[i].GetX(), this->SVGTextCoords[i].GetY(), this->SVGTextCoords[i].GetZ());

                font.SetRenderType(font.RENDERTYPE_FILL);
                this->font.DrawString(fontColor, a.x, a.y, a.z, this->SVGTextSizes[i] * this->SVGScales[i] * 0.022f,
                    false, this->SVGTextContent[i].PeekBuffer(), utility::SDFFont::Alignment::ALIGN_LEFT_BOTTOM);
            }

            ////////////////////////////

            // create star menu TEXT and COLOR
            for (int i = 0; i < this->starMenu.circleCnt; i++) {
                if (this->starMenu.curPageStartIDx + i < clusterSize) {

                    float textSize = this->starMenu.textSize * this->subCircleSizefactor * 1.43f * postsubTextScale;
                    for (int j = 0; j < this->subCirclePointers.size(); j++) {
                        auto subS = subCirclePointers[j];
                        if (this->font.IsInitialised()) {
                            if (this->renderDataControl.improveTextReadability) {
                                font.SetRenderType(font.RENDERTYPE_OUTLINE);
                                this->font.DrawString(&fontOutline_color.x, subS[i * 4 + 0], subS[i * 4 + 1],
                                    subS[i * 4 + 2], textSize, false,
                                    this->starMenu.curSubCircleContent[i][j].PeekBuffer(),
                                    utility::SDFFont::Alignment::ALIGN_CENTER_MIDDLE);
                            }

                            font.SetRenderType(font.RENDERTYPE_FILL);
                            this->font.DrawString(this->subCircleColorPointers[j], subS[i * 4 + 0], subS[i * 4 + 1],
                                subS[i * 4 + 2], textSize, false, this->starMenu.curSubCircleContent[i][j].PeekBuffer(),
                                utility::SDFFont::Alignment::ALIGN_CENTER_MIDDLE);
                        }
                    }
                }
            }
            vislib::StringA tmpStr;
            tmpStr.Format("×");
            font.SetRenderType(font.RENDERTYPE_FILL);
            this->font.DrawString(&fontOutline_color.x, starCircles[(indexOfLastStarCircle + 1) * 4 + 0],
                starCircles[(indexOfLastStarCircle + 1) * 4 + 1], starCircles[(indexOfLastStarCircle + 1) * 4 + 2],
                this->starMenu.textSize * postsubTextScale, false, tmpStr.PeekBuffer(),
                utility::SDFFont::Alignment::ALIGN_CENTER_MIDDLE);

            /********************************************
             ****** Render Functional Groups Lable ******
            /********************************************/
            if (this->fgsPerClusterFiltered[this->s.curClusterID].clusterCnt > 0) {
                float upScale = 1.25f;
                float textScale = this->starMenu.textSize * 0.15;
                int IDx = this->renVarsFGS.selectedSphereID;
                if (IDx != -1) {
                    // printf("IDx: %d FuncGroup: %d Name:%s\n", IDx, this->funcID_to_ClustCentroidID[IDx],
                    //     lmc->functionalGroupsIDsToWord[this->funcID_to_ClustCentroidID[IDx]].c_str());
                    vislib::StringA funcGroupString = "";
                    vislib::StringA openBracket = "(";
                    auto fgsMap = this->fgsPerClusterFiltered[this->s.curClusterID].fgsMaps_GblMDlIDs[IDx];


                    map<uint, std::vector<uint>>::iterator it;
                    for (it = fgsMap.begin(); it != fgsMap.end(); it++) {
                        int fgsTypeID = it->first;

                        int funcCentroidSize = it->second.size();
                        tmpStr.Format("%d", funcCentroidSize);
                        if (funcGroupString.Length() > 0) {

                            funcGroupString = funcGroupString + " \n " + openBracket + tmpStr + ") " +
                                              lmc->functionalGroupsIDsToWord[fgsTypeID].c_str();
                        } else {
                            funcGroupString =
                                openBracket + tmpStr + ") " + lmc->functionalGroupsIDsToWord[fgsTypeID].c_str();
                        }
                    }

                    if (this->renderDataControl.improveTextReadability) {
                        this->font.SetRenderType(font.RENDERTYPE_OUTLINE);
                        this->font.DrawString(&fontOutline_color.x,
                            this->fgsPerClusterFiltered[this->s.curClusterID].clusterCentroids[IDx * 4 + 0] +
                                UpVecStarMenu.x * upScale,
                            this->fgsPerClusterFiltered[this->s.curClusterID].clusterCentroids[IDx * 4 + 1] +
                                UpVecStarMenu.y * upScale,
                            this->fgsPerClusterFiltered[this->s.curClusterID].clusterCentroids[IDx * 4 + 2] +
                                UpVecStarMenu.z * upScale,
                            textScale, false, funcGroupString.PeekBuffer(),
                            utility::SDFFont::Alignment::ALIGN_CENTER_MIDDLE);
                    }

                    this->font.SetRenderType(font.RENDERTYPE_FILL);
                    this->font.DrawString(&fontFGSclustLables_color.x,
                        this->fgsPerClusterFiltered[this->s.curClusterID].clusterCentroids[IDx * 4 + 0] +
                            UpVecStarMenu.x * upScale,
                        this->fgsPerClusterFiltered[this->s.curClusterID].clusterCentroids[IDx * 4 + 1] +
                            UpVecStarMenu.y * upScale,
                        this->fgsPerClusterFiltered[this->s.curClusterID].clusterCentroids[IDx * 4 + 2] +
                            UpVecStarMenu.z * upScale,
                        textScale, false, funcGroupString.PeekBuffer(),
                        utility::SDFFont::Alignment::ALIGN_CENTER_MIDDLE);
                }
            }


            /**
             * End of Text Rendering
             */
        } // end if (this->selectedIDx >= 0)
          // render residue lables
    }     // end else (cluster sphere selected)


    if (starCircles.Count() > 1) {
    circleTextureRenderer(
        this->starSubSpheres_onSelect, this->tf.tf_functGroups, this->renVarsSMS, viewportStuff);
    }

    if (this->renderDataControl.residueLabels) {
        float correction = this->scaleWorld / 0.0300282277;
        for (int i = 0; i < this->residueLables.Count(); i++) {

            if (this->renderDataControl.improveTextReadability) {
                this->font.SetRenderType(font.RENDERTYPE_OUTLINE);
                this->font.DrawString(&fontOutline_color.x, this->residueLablesPositions[i].GetX(),
                    this->residueLablesPositions[i].GetY(), this->residueLablesPositions[i].GetZ(),
                    this->starMenu.textSize * 0.15 * correction, false, this->residueLables[i].PeekBuffer(),
                    utility::SDFFont::Alignment::ALIGN_CENTER_MIDDLE);
            }

            this->font.SetRenderType(font.RENDERTYPE_FILL);
            this->font.DrawString(&fontResidueLables_color.x, this->residueLablesPositions[i].GetX(),
                this->residueLablesPositions[i].GetY(), this->residueLablesPositions[i].GetZ(),
                this->starMenu.textSize * 0.15 * correction, false, this->residueLables[i].PeekBuffer(),
                utility::SDFFont::Alignment::ALIGN_CENTER_MIDDLE);
        }
    }

    glEnable(GL_DEPTH_TEST);


    // control set of protein residues for selected functional group centroid
    if (this->renVarsFGS.selectedSphereID != -1 &&
        this->renVarsFGS.old_selectedSphereID != this->renVarsFGS.selectedSphereID) {
        this->pmdcOut_newData = true;
        this->renVarsFGS.old_selectedSphereID = this->renVarsFGS.selectedSphereID;
    }

    this->isClicked = -1;
    this->s.old_HoveringIDx = this->s.hoveringIDx;
    this->s.old_SelectedIDx = this->s.selectedIDx;
    this->firstRenderCall = false;
    this->s.changedcurClusterID = 0;

    if (lmc->getCurrentGlobalModelID() != this->s.selectedGlobalMdlID && this->s.selectedGlobalMdlID != -1) {
        lmc->SetCurrentLigandAndModel_byGlobalMdlID(this->s.selectedGlobalMdlID);
        (*lmc)(LigandModelCall::CallForGetData);
        if (!this->webDataChanged) {
            lmc->getWebData()->setLigandID(lmc->getCurrentLigandID());
            lmc->getWebData()->setModelID(lmc->getCurrentModelID());
        }
    }

    // handle sending data to python server


    auto fgs = lmc->getClusteredFunctionalGroups();
    bool modelClustRenderer_busy = 0;
    bool socketCommManager_busy = 1;
    bool firstRenderCall_complete = 1;
    if ((fgs->getDataHash() >= 0 && fgs->getDataHash() != this->old_functionalGroup_dataHash || this->webDataChanged) &&
        !this->withoutSocket) {
        this->old_functionalGroup_dataHash = fgs->getDataHash();
        lmc->SetlmcStatus(modelClustRenderer_busy, socketCommManager_busy, firstRenderCall_complete);
        (*lmc)(LigandModelCall::CallForSetLMCStatus);
    }
    if ((lmc->getClusterAndCentroidData().DBSCANParams.paramsChanged == true) && !this->withoutSocket) {
        lmc->SetlmcStatus(modelClustRenderer_busy, socketCommManager_busy, firstRenderCall_complete);
        (*lmc)(LigandModelCall::CallForSetLMCStatus);
        lmc->ResetDBSCANChanged();
    }
    (*lmc)(LigandModelCall::CallForGetDataSilent);

    // default values
    lmc->SetlmcStatus(0, 0, 1);
    (*lmc)(LigandModelCall::CallForSetLMCStatus);

    this->fgcClusterParamsChanged = false;
    lmc->ResetDBSCANChanged();
    this->s.changedCurForceType = false;
    this->starMenu.changed_curSortingFroceType = false;

    // lmc->SetCurrentLigandAndModel_byGlobalMdlID(this->s.selectedGlobalMdlID);
    // int lig = lmc->getCurrentLigandID();
    // int model = lmc->getCurrentModelID();
    //printf("gblID: %i current lig: %i model: %i\n", this->s.selectedGlobalMdlID, lig, model);
    // this->webDataChanged = false;
    glPopMatrix();
    return true;
}

void ModelClusterRenderer::getLineCoords(vec3 startPoint, float length, float width, vec3 lineDirection,
    vec3 eyeDirection, vislib::Array<vislib::math::Vector<float, 3>>& lines, int offset = 0) {


    lineDirection = normalize(lineDirection);
    eyeDirection = normalize(eyeDirection);
    vec3 rotVec = glm::cross(lineDirection, eyeDirection);
    float angle = glm::acos(glm::dot(lineDirection, eyeDirection)) - (float)(90 * vislib::math::PI_DOUBLE / 180);
    lineDirection = glm::rotate(lineDirection, angle, rotVec);
    vec3 orthoDirection = glm::rotate(lineDirection, (float)(90 * vislib::math::PI_DOUBLE / 180), eyeDirection);
    orthoDirection = {orthoDirection.x * width, orthoDirection.y * width, orthoDirection.z * width};

    glm::vec3 starMenuInner = {startPoint.x + lineDirection.x * offset, startPoint.y + lineDirection.y * offset,
        startPoint.z + lineDirection.z * offset};
    glm::vec3 starMenuOuter = {startPoint.x + lineDirection.x * length, startPoint.y + lineDirection.y * length,
        startPoint.z + lineDirection.z * length};


    // [A] -> rectangle inner star
    lines.Append(
        {starMenuInner.x - orthoDirection.x, starMenuInner.y - orthoDirection.y, starMenuInner.z - orthoDirection.z});
    // [D] -> rectangle inner star
    lines.Append(
        {starMenuInner.x + orthoDirection.x, starMenuInner.y + orthoDirection.y, starMenuInner.z + orthoDirection.z});
    // [B] -> rectangle outer star
    lines.Append(
        {starMenuOuter.x - orthoDirection.x, starMenuOuter.y - orthoDirection.y, starMenuOuter.z - orthoDirection.z});
    // [C] -> rectangle outer star
    lines.Append(
        {starMenuOuter.x + orthoDirection.x, starMenuOuter.y + orthoDirection.y, starMenuOuter.z + orthoDirection.z});
}


bool ModelClusterRenderer::paramsStarMenuChanged(param::ParamSlot& param) {
    // set param as clean again
    param.ResetDirty();
    int starMenuTextSize = this->starMenuTextSizeParam.Param<param::IntParam>()->Value();
    int starMenuCnt = this->starMenuCntParam.Param<param::IntParam>()->Value();
    int starMenuSize = this->starMenuSizeParam.Param<param::IntParam>()->Value();
    // int starMenuInnerStart = this->starMenuInnerStartParam.Param<param::IntParam>()->Value();
    int subSphereSizefactor = this->subSphereSizefactorParam.Param<param::IntParam>()->Value();

    // MenuTextSize: handle out of bounds values
    if (starMenuTextSize < 0) {
        this->starMenu.textSize == 0;
        this->starMenuTextSizeParam.Param<param::IntParam>()->SetValue(0);
    } else {
        this->starMenu.textSize = (float)(starMenuTextSize) / 100;
    }

    if (subSphereSizefactor < 0) {
        this->subCircleSizefactor == 0;
        this->subSphereSizefactorParam.Param<param::IntParam>()->SetValue(0);
    } else {
        this->subCircleSizefactor = (float)(subSphereSizefactor) / 100;
        this->subSphereSizefactorParam.Param<param::IntParam>()->SetValue(subSphereSizefactor);
    }
    // starMenu.circleCnt: handle out of bounds values
    if (starMenuCnt < 1) {
        this->starMenu.circleCnt == 1;
        this->starMenuCntParam.Param<param::IntParam>()->SetValue(1);
    } else if (starMenuCnt > this->clusterData.DBSCANParams.minPts || starMenuCnt > this->starMenu.maxCircleCnt) {
        this->starMenu.circleCnt = starMenuCnt > this->starMenu.maxCircleCnt ? this->starMenu.maxCircleCnt
                                                                             : this->clusterData.DBSCANParams.minPts;
        this->starMenuCntParam.Param<param::IntParam>()->SetValue(this->starMenu.circleCnt);
    } else {
        this->starMenu.circleCnt = starMenuCnt;
    }
    // starMenuSize; handle out of bounds values
    if (starMenuSize < 0) {
        this->starMenu.size == 0;
        this->starMenuSizeParam.Param<param::IntParam>()->SetValue(0);
    } else {
        this->starMenu.size = starMenuSize;
    }

    // starMenuInnerStart; handle out of bounds and not valid values
    /*
    if (starMenuInnerStart < 0) {
        this->starMenuInnerStart == 0;
        this->starMenuInnerStartParam.Param<param::IntParam>()->SetValue(0);
    } else if (starMenuInnerStart > this->starMenuSize) {
        this->starMenuInnerStart = this->starMenuSize;
        this->starMenuInnerStartParam.Param<param::IntParam>()->SetValue(this->starMenuInnerStart);
    } else {
        this->starMenuInnerStart = starMenuInnerStart;
    }
    */
    return true;
};


float ModelClusterRenderer::getElementValue_svgLine(std::string line, std::string elementString) {

    elementString.append("=\"[0-9]*\.?[0-9]*\"");

    std::regex r1(elementString);
    std::regex r2("\"[0-9]+\.?[0-9]*\"");
    std::regex r3("[0-9]+\.?[0-9]*");

    std::smatch m1;
    std::smatch m2;
    std::smatch m3;
    std::regex_search(line, m1, r1);
    if (m1.size() <= 0) {
        printf("ERROR: XML-Parsing '%s' failed\n", elementString.c_str());
        return 0.0f;
    }


    std::string tmpString1 = m1.str();
    std::regex_search(tmpString1, m2, r2);


    std::string tmpString2 = m2.str();
    std::regex_search(tmpString2, m3, r3);

    return std::stof(m3.str());
}

vislib::math::Vector<int, 3> ModelClusterRenderer::getElementRGB_svgLine(std::string line, std::string elementString) {

    vislib::math::Vector<int, 3> RGBvector;
    elementString.append("=\"rgb{1}\\({1}[0-9]*\,[0-9]*\,[0-9]*\\){1}\"");

    std::regex r(elementString);
    std::smatch m1;
    std::smatch m2;
    std::regex_search(line, m1, r);
    if (m1.size() <= 0) {
        printf("ERROR: XML-Parsing '%s' failed\n", elementString.c_str());
        RGBvector.SetX(0);
        RGBvector.SetY(0);
        RGBvector.SetZ(0);
        return RGBvector;
    }

    std::regex r2("[0-9]*\,[0-9]*\,[0-9]*");
    std::string zu = m1.str();
    std::regex_search(zu, m2, r2);

    std::vector<std::string> tmpSplit;

    split(m2.str(), ",", tmpSplit);
    RGBvector.SetX(std::stoi(tmpSplit[0]));
    RGBvector.SetY(std::stoi(tmpSplit[1]));
    RGBvector.SetZ(std::stoi(tmpSplit[2]));

    return RGBvector;
}

void ModelClusterRenderer::getPoints_svgPolygon(
    std::string line, vislib::Array<vislib::math::Vector<float, 2>>& points) {
    std::string elementString;
    vislib::math::Vector<int, 3> RGBvector;
    elementString.append("(points){1}=\"([0-9]+(.[0-9]+)*\\s)+\"");

    std::regex r(elementString);
    std::smatch m1;
    std::smatch m2;
    std::regex_search(line, m1, r);
    if (m1.size() <= 0) {
        printf("ERROR: XML-Parsing '%s' failed\n", elementString.c_str());
        points.Append({0.f, 0.f});
        points.Append({0.f, 0.f});
        points.Append({0.f, 0.f});
    }
    std::regex r2("(([0-9]+(.[0-9]+)*)\\s*)+");
    std::string m1String = m1.str();
    std::regex_search(m1String, m2, r2);

    std::string tmpMatch1 = m2.str();
    std::vector<std::string> tmpSplit;
    split(tmpMatch1, " ", tmpSplit);
    int counter = 0;
    for (int i = 0; i < tmpSplit.size() - 1; i += 2) {
        points.Append({std::stof(tmpSplit[i + 0]), std::stof(tmpSplit[i + 1])});
        counter++;
    }
}

std::string ModelClusterRenderer::getElementText_svgLine(std::string line, std::string elementString) {

    if (elementString == "text") {


        std::regex r1("\>{1}[a-zA-Z0-9_+-.]*\<{1}");

        std::smatch m1;
        std::smatch m2;
        std::regex_search(line, m1, r1);
        if (m1.size() <= 0) {
            printf("ERROR: XML-Parsing '%s' failed\n", elementString.c_str());
            return "";
        }


        std::regex r2("[a-zA-Z0-9_+-.]+");
        std::string m1String = m1.str();
        std::regex_search(m1String, m2, r2);

        if (m2.str() == "." || m2.str() == "..") {
            return "";
        } else {
            return m2.str();
        }


    } else {

        elementString.append("=\"[a-z]*\"");

        std::regex r1(elementString);
        std::regex r2("\"[a-z]*\"");
        std::regex r3("[a-z]+");

        std::smatch m1;
        std::smatch m2;
        std::smatch m3;
        std::regex_search(line, m1, r1);
        if (m1.size() <= 0) {
            printf("ERROR: XML-Parsing '%s' failed\n", elementString.c_str());
            return "";
        }

        std::string m1String = m1.str();
        std::regex_search(m1String, m2, r2);

        std::string m2String = m2.str();
        std::regex_search(m2String, m3, r3);

        return m3.str();
    }
}

void ModelClusterRenderer::parseSVG(std::string filepath, class SVG1* svg1) {
    ifstream checkFile;
    checkFile.open(filepath.c_str());
    svg1->svgCentroid.SetX(0);
    svg1->svgCentroid.SetY(0);
    int centroidCnt = 0;

    if (checkFile.fail()) {
        printf("SVG file does not exist: %s", filepath.c_str());
    } else {
        int lCnt;
        int pCnt;
        int tCnt;
        std::string line;
        while (std::getline(checkFile, line)) {
            if (line._Starts_with("<line")) {
                // handle line
                Svg_line1 ll;
                svg1->svg_lineElements.push_back(ll);
                lCnt = svg1->lineCnt;
                svg1->svg_lineElements[lCnt].x1 = getElementValue_svgLine(line, "x1");
                svg1->svg_lineElements[lCnt].y1 = getElementValue_svgLine(line, "y1");
                svg1->svg_lineElements[lCnt].x2 = getElementValue_svgLine(line, "x2");
                svg1->svg_lineElements[lCnt].y2 = getElementValue_svgLine(line, "y2");
                svg1->svg_lineElements[lCnt].opacity = getElementValue_svgLine(line, "opacity");
                svg1->svg_lineElements[lCnt].strokeWidth = getElementValue_svgLine(line, "stroke-width");
                svg1->svg_lineElements[lCnt].color = getElementRGB_svgLine(line, "stroke");

                svg1->lineCnt++;

                svg1->svgCentroid.SetX(svg1->svgCentroid.GetX() + svg1->svg_lineElements[lCnt].x1);
                svg1->svgCentroid.SetX(svg1->svgCentroid.GetX() + svg1->svg_lineElements[lCnt].x2);
                svg1->svgMAX_xy.SetX(svg1->svgMAX_xy.GetX() < svg1->svg_lineElements[lCnt].x1
                                         ? svg1->svg_lineElements[lCnt].x1
                                         : svg1->svgMAX_xy.GetX());
                svg1->svgMAX_xy.SetX(svg1->svgMAX_xy.GetX() < svg1->svg_lineElements[lCnt].x2
                                         ? svg1->svg_lineElements[lCnt].x1
                                         : svg1->svgMAX_xy.GetX());
                svg1->svgMIN_xy.SetX(svg1->svgMIN_xy.GetX() > svg1->svg_lineElements[lCnt].x1
                                         ? svg1->svg_lineElements[lCnt].x1
                                         : svg1->svgMIN_xy.GetX());
                svg1->svgMIN_xy.SetX(svg1->svgMIN_xy.GetX() > svg1->svg_lineElements[lCnt].x2
                                         ? svg1->svg_lineElements[lCnt].x1
                                         : svg1->svgMIN_xy.GetX());
                centroidCnt++;
                svg1->svgCentroid.SetY(svg1->svgCentroid.GetY() + svg1->svg_lineElements[lCnt].y1);
                svg1->svgCentroid.SetY(svg1->svgCentroid.GetY() + svg1->svg_lineElements[lCnt].y2);
                svg1->svgMAX_xy.SetY(svg1->svgMAX_xy.GetY() < svg1->svg_lineElements[lCnt].y1
                                         ? svg1->svg_lineElements[lCnt].y1
                                         : svg1->svgMAX_xy.GetY());
                svg1->svgMAX_xy.SetY(svg1->svgMAX_xy.GetY() < svg1->svg_lineElements[lCnt].y2
                                         ? svg1->svg_lineElements[lCnt].y1
                                         : svg1->svgMAX_xy.GetY());
                svg1->svgMIN_xy.SetY(svg1->svgMIN_xy.GetY() > svg1->svg_lineElements[lCnt].y1
                                         ? svg1->svg_lineElements[lCnt].y1
                                         : svg1->svgMIN_xy.GetY());
                svg1->svgMIN_xy.SetY(svg1->svgMIN_xy.GetY() > svg1->svg_lineElements[lCnt].y2
                                         ? svg1->svg_lineElements[lCnt].y1
                                         : svg1->svgMIN_xy.GetY());

                centroidCnt++;
            }

            if (line._Starts_with("<polygon")) {
                // handle line
                Svg_ploygon1 pp;
                svg1->svg_polygons.push_back(pp);
                pCnt = svg1->polygonCnt;
                svg1->svg_polygons[pCnt].strokeWidth = getElementValue_svgLine(line, "stroke-width");
                svg1->svg_polygons[pCnt].color = getElementRGB_svgLine(line, "fill");
                getPoints_svgPolygon(line, svg1->svg_polygons[pCnt].points);

                // for SVG centroid calculation
                for (int i = 0; i < svg1->svg_polygons[pCnt].points.Count(); i++) {
                    auto points = svg1->svg_polygons[pCnt].points;
                    // svg->svgCentroid.SetX(svg->svgCentroid.GetX() + points[i].GetX());
                    // svg->svgCentroid.SetY(svg->svgCentroid.GetY() + points[i].GetY());
                    // centroidCnt++;
                    svg1->svgMAX_xy.SetX(
                        svg1->svgMAX_xy.GetX() < points[i].GetX() ? points[i].GetX() : svg1->svgMAX_xy.GetX());
                    svg1->svgMAX_xy.SetY(
                        svg1->svgMAX_xy.GetY() < points[i].GetY() ? points[i].GetY() : svg1->svgMAX_xy.GetY());
                    svg1->svgMIN_xy.SetX(
                        svg1->svgMIN_xy.GetX() > points[i].GetX() ? points[i].GetX() : svg1->svgMIN_xy.GetX());
                    svg1->svgMIN_xy.SetY(
                        svg1->svgMIN_xy.GetY() > points[i].GetY() ? points[i].GetY() : svg1->svgMIN_xy.GetY());
                }

                svg1->polygonCnt++;
            }


            if (line._Starts_with("<text")) {
                // handle text
                Svg_text1 tt;
                svg1->svg_textElements.push_back(tt);
                tCnt = svg1->textCnt;
                svg1->svg_textElements[tCnt].x = getElementValue_svgLine(line, "x");
                svg1->svg_textElements[tCnt].y = getElementValue_svgLine(line, "y");
                svg1->svg_textElements[tCnt].fontSize = getElementValue_svgLine(line, "font-size");
                svg1->svg_textElements[tCnt].fontWeight = getElementText_svgLine(line, "font-weight");
                svg1->svg_textElements[tCnt].strokeWidth = getElementValue_svgLine(line, "stroke-width");
                svg1->svg_textElements[tCnt].color = getElementRGB_svgLine(line, "fill");
                svg1->svg_textElements[tCnt].text = getElementText_svgLine(line, "text");

                // svg->svgCentroid.SetX(svg->svgCentroid.GetX() + svg->svg_textElements[tCnt].x);
                // svg->svgCentroid.SetY(svg->svgCentroid.GetY() + svg->svg_textElements[tCnt].y);
                // centroidCnt++;
                svg1->svgMAX_xy.SetX(svg1->svgMAX_xy.GetX() < svg1->svg_textElements[tCnt].x
                                         ? svg1->svg_textElements[tCnt].x
                                         : svg1->svgMAX_xy.GetX());
                svg1->svgMAX_xy.SetY(svg1->svgMAX_xy.GetY() < svg1->svg_textElements[tCnt].y
                                         ? svg1->svg_textElements[tCnt].y
                                         : svg1->svgMAX_xy.GetY());
                svg1->svgMIN_xy.SetX(svg1->svgMIN_xy.GetX() > svg1->svg_textElements[tCnt].x
                                         ? svg1->svg_textElements[tCnt].x
                                         : svg1->svgMIN_xy.GetX());
                svg1->svgMIN_xy.SetY(svg1->svgMIN_xy.GetY() > svg1->svg_textElements[tCnt].y
                                         ? svg1->svg_textElements[tCnt].y
                                         : svg1->svgMIN_xy.GetY());


                svg1->textCnt++;
            }
            if (line._Starts_with("</svg>")) {
                break;
            }
        }
        svg1->svgCentroid.SetX(svg1->svgCentroid.GetX() / centroidCnt);
        svg1->svgCentroid.SetY(svg1->svgCentroid.GetY() / centroidCnt);
        svg1->svgMiddle.SetX(svg1->svgMIN_xy.GetX() + (svg1->svgMAX_xy.GetX() - svg1->svgMIN_xy.GetX()) / 2);
        svg1->svgMiddle.SetY(svg1->svgMIN_xy.GetY() + (svg1->svgMAX_xy.GetY() - svg1->svgMIN_xy.GetY()) / 2);

        checkFile.close();
    }
};


void ModelClusterRenderer::MouseRayIntersection::release() {}

/**
 * returns indiex of first intersecting sphere
 *
 * @param spheres		array with spheres (local coordinates)
 * @param tmpCam		vislib::graphics::gl::CameraOpenGL&
 * @param rayOrigin		(world coordinates)
 * @param rayDirection	(world coordinates)
 * @param mouseAction	(0 or 1->action)
 */
int ModelClusterRenderer::MouseRayIntersection::getFirstMouseRayIntersection(vislib::Array<float>& spheres,
    vislib::graphics::gl::CameraOpenGL& tmpCam, float mouseX, float mouseY, int mouseAction) {
    this->firstRayIntersetion = -1;
    getMouseRayIntersections(spheres, tmpCam, mouseX, mouseY, mouseAction);
    return this->firstRayIntersetion;
};

/**
 * returns indices of intersecting spheres
 *
 * @param spheres		array with spheres (local coordinates)
 * @param tmpCam		vislib::graphics::gl::CameraOpenGL&
 * @param rayOrigin		(world coordinates)
 * @param rayDirection	(world coordinates)
 * @param mouseAction	(0 or 1->action)
 */
vislib::Array<int> ModelClusterRenderer::MouseRayIntersection::getMouseRayIntersections(vislib::Array<float>& spheres,
    vislib::graphics::gl::CameraOpenGL& tmpCam, float mouseX, float mouseY, int mouseAction) {

    this->mouseX = mouseX;
    this->mouseY = mouseY;
    this->spheres = spheres;


    /****** initiate first ray ******/
    /*
     * check for changes of projection matrix
     * info: using mux renderer has an effect on projection matrix after some frames
     */
    // float projectionMatrix[16];
    // glGetFloatv(GL_PROJECTION_MATRIX, projectionMatrix);
    // int changed_ProjectionMatrix;
    // roundf(projectionMatrix[0] * 1000) != roundf(this->transformMatrices.camera_to_canonical[0][0] * 1000)
    //    ? changed_ProjectionMatrix = 1
    //    : changed_ProjectionMatrix = 0;

    /****** check if update of transform matrices is necessary ******/
    /*
     * check if camera moved and if movement is finished
     */
    // int changedCamera;
    // this->cam.Parameters()->Position().GetX() != this->priviousCamPos.x ||
    //        this->cam.Parameters()->Position().GetY() != this->priviousCamPos.y ||
    //        this->cam.Parameters()->Position().GetZ() != this->priviousCamPos.z
    //    ? changedCamera = 1
    //    : changedCamera = 0;


    // if ((changedCamera && mouseAction) || changed_ProjectionMatrix) {
    this->cam = tmpCam;
    setTransformMatrices(tmpCam);
    setRayOrigin(tmpCam.Parameters()->Position());
    //}
    setRayDirection(mouseX, mouseY, getTransformMatrices(tmpCam));
    setRayIntersections(this->spheres, this->rayOrigin, this->rayDirection);

    return this->rayIntersections;
}


/**
 * sets indices of intersecting spheres to an array
 *
 * @param spheres		array with spheres (local coordinates)
 * @param rayOrigin		(world coordinates)
 * @param rayDirection	(world coordinates)
 */
void ModelClusterRenderer::MouseRayIntersection::setRayIntersections(
    vislib::Array<float> spheres, glm::vec3 rayOrigin, glm::vec3 rayDirection) {
    //  ********* intersection detection *********

    // https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection

    glm::vec3 camPos, L, center;
    glm::vec4 center_world, center_plus_r_world, camPos4;
    float radiusWorld, radius2, d2, tca, thc, t, t0, t1;

    vislib::Array<float> t_store;
    vislib::Array<int> t_storeID;

    camPos4 = glm::vec4(0, 0, 0, 1);
    camPos = glm::vec3(camPos4 * this->transformMatrices.camera_to_world);

    for (int i = 0; i < spheres.Count() / 4; i++) {
        t, t0, t1 = 0;
        center = glm::vec3(spheres[i * 4 + 0], spheres[i * 4 + 1], spheres[i * 4 + 2]);
        center_world = glm::vec4(center, 1);
        center_world = this->transformMatrices.local_to_world * center_world;

        glm::vec4 radiusLocal4 = glm::vec4(spheres[i * 4 + 3], 0, 0, 1);
        radiusWorld = (this->transformMatrices.local_to_world * radiusLocal4).x;

        L = vec3(center_world) - camPos;
        tca = dot(L, rayDirection);

        d2 = dot(L, L) - (tca * tca);
        radius2 = radiusWorld * radiusWorld;

        if ((d2 < radius2)) {

            thc = sqrt(radius2 - d2);
            t0 = tca - thc;
            t1 = tca + thc;

            if (t0 > t1) std::swap(t0, t1);

            if (t0 < 0) {
                t0 = t1; // if t0 is negative, let's use t1 instead
                // if (t0 < 0) return false; // both t0 and t1 are negative
            }
            t = t0;

            if (t != 0) {
                t_storeID.Append(i);
                t_store.Append(t);
            }
        }
    }

    float min_t = FLT_MAX;
    int min_tID = -1;
    // get minimum t and the stored ID for it to get the first intersection with the ray
    for (size_t i = 0; i < t_store.Count(); i++) {
        if (min_t > t_store[i]) {
            min_t = t_store[i];
            min_tID = t_storeID[i];
        }
    }
    this->firstRayIntersetion = min_tID;
    this->rayIntersections = t_storeID;
};

/**
 * sets the origin of the ray in world coordinates (e.g camera position)
 *
 * @param rayOrigin		(world coordinates)
 */
void ModelClusterRenderer::MouseRayIntersection::setRayOrigin(vislib::math::Vector<float, 3> rayOrigin) {
    this->rayOrigin = glm::vec3(rayOrigin.GetX(), rayOrigin.GetY(), rayOrigin.GetZ());
};


/**
 * sets the direction vector of the ray
 *
 * @param mouseX	(screen coordinates)
 * @param mouseY	(screen coordinates)
 * @param matrices	(object with transformation matrices for different 3D-spaces)
 */
void ModelClusterRenderer::MouseRayIntersection::setRayDirection(
    float mouseX, float mouseY, TransformMatrices matrices) {
    float vport[4];
    glGetFloatv(GL_VIEWPORT, vport);

    float normalised_x = 2 * this->mouseX / vport[2] - 1;
    float normalised_y = 1 - 2 * this->mouseY / vport[3];

    // https://community.khronos.org/t/mouse-picking-with-opengl/76055

    vec4 mouseCanonicalVec4;

    mouseCanonicalVec4[0] = normalised_x;
    mouseCanonicalVec4[1] = normalised_y;
    mouseCanonicalVec4[2] = -1;
    mouseCanonicalVec4[3] = 1;

    vec4 mouseCameraVec4 = mouseCanonicalVec4 * matrices.canonical_to_camera;
    mouseCameraVec4 = glm::vec4(mouseCameraVec4.x, mouseCameraVec4.y, mouseCanonicalVec4.z, 0.0f);
    vec4 rayWorld = mouseCameraVec4 * matrices.camera_to_world;
    vec3 rayDirection = glm::normalize(glm::vec3(rayWorld));
    this->rayDirection = rayDirection;
};


/**
 * sets all transformation matrices (from local-space until screen-space and the other way around)
 *
 * @param tmpCam (vislib::graphics::gl::CameraOpenGL& object)
 */
void ModelClusterRenderer::MouseRayIntersection::setTransformMatrices(vislib::graphics::gl::CameraOpenGL& tmpCam) {

    this->priviousCamPos = glm::vec3(this->cam.Parameters()->Position().GetX(),
        this->cam.Parameters()->Position().GetY(), this->cam.Parameters()->Position().GetZ());

    float32 local_to_camera1D[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, local_to_camera1D);

    float world_to_camera1D[16];
    float projectionMatrix1D[16];
    this->cam.ViewMatrix(&world_to_camera1D[0]);
    this->cam.ProjectionMatrix(&projectionMatrix1D[0]);

    mat4 projectionMatrix;

    for (size_t i = 0; i < 4; i++) {
        for (size_t j = 0; j < 4; j++) {
            this->transformMatrices.world_to_camera[j][i] = world_to_camera1D[j + i * 4];
            this->transformMatrices.camera_to_canonical[j][i] = projectionMatrix1D[j + i * 4];
            this->transformMatrices.local_to_camera[j][i] = local_to_camera1D[j + i * 4];
        }
    }

    this->transformMatrices.camera_to_local = glm::inverse(this->transformMatrices.local_to_camera);

    this->transformMatrices.canonical_to_camera = glm::inverse(this->transformMatrices.camera_to_canonical);
    this->transformMatrices.camera_to_world = glm::inverse(this->transformMatrices.world_to_camera);
    this->transformMatrices.local_to_world =
        this->transformMatrices.local_to_camera * this->transformMatrices.camera_to_world;
    this->transformMatrices.world_to_local = glm::inverse(this->transformMatrices.local_to_world);
}

/**
 * returns the transformation matrices object (for different 3D-sapces)
 *
 * @param tmpCam		vislib::graphics::gl::CameraOpenGL&
 */
ModelClusterRenderer::MouseRayIntersection::TransformMatrices
ModelClusterRenderer::MouseRayIntersection::getTransformMatrices(vislib::graphics::gl::CameraOpenGL& tmpCam) {
    setTransformMatrices(tmpCam);
    return this->transformMatrices;
}


void ModelClusterRenderer::planeRenderer(float* clipPoint, float* clipDat, float* clipCol) {

    this->simpleAmbientShader.Enable();
    glUniform4fv(this->simpleAmbientShader.ParameterLocation("baseColor"), 1, &clipCol[0]);
    glUniformMatrix4fvARB(
        this->simpleAmbientShader.ParameterLocation("modVIEW"), 1, false, this->modelMatrix->PeekComponents());
    glUniformMatrix4fvARB(
        this->simpleAmbientShader.ParameterLocation("proj"), 1, false, this->projMatrix->PeekComponents());
    float lightAmbient = 1.0f;
    glUniform1f(this->simpleAmbientShader.ParameterLocation("lightAmbient"), lightAmbient);

   
    vislib::Array<vislib::math::Vector<float, 3>> planeCoords;
    planeCoords.SetCount(4);
    float scale = 200.0f;
    glm::vec3 xPlaneDir;
    glm::vec3 yPlaneDir;
    float fact = 1.0000;
    glm::vec3 point(clipPoint[0] * fact, clipPoint[1] * fact, clipPoint[2] * fact);
    glm::vec3 normal(clipDat[0], clipDat[1], clipDat[2]);
    float w = clipDat[3];

    vec3 tmpP = normal;
    xPlaneDir = glm::normalize(glm::cross(vec3(1.0, 0.0, 0.0), normal));
    if (isnan(xPlaneDir.x)) {
        tmpP = glm::normalize(vec3(0.0, 1.0, 0.0));
        xPlaneDir = glm::normalize(glm::cross(tmpP, normal));
    }

    yPlaneDir = glm::normalize(glm::cross(xPlaneDir, normal));


    vec3 normAdd = glm::normalize(xPlaneDir + yPlaneDir);
    vec3 normSub = glm::normalize(xPlaneDir - yPlaneDir);

    glm::vec3 v2 = point - normAdd * scale;
    glm::vec3 v1 = point - normSub * scale;
    glm::vec3 v3 = point + normAdd * scale;
    glm::vec3 v4 = point + normSub * scale;


    glBegin(GL_TRIANGLE_STRIP);
    glColor4f(0.5, 0.5, 0.5, 1.0);

    glVertex3f(v1.x, v1.y, v1.z);
    glVertex3f(v2.x, v2.y, v2.z);
    glVertex3f(v3.x, v3.y, v3.z);
    glVertex3f(v4.x, v4.y, v4.z);

    glEnd();

     /*******************/
    this->simpleAmbientShader.Disable();
}


/**
 * renders the mesh data of the protein
 * (sorts vertices for transpararent rendering)
 *
 * @param tmc		CallTriMeshData
 * @param call		megamol::core::Call& (AbstractCall3D)
 */
bool ModelClusterRenderer::proteinMeshRenderer(core::view::AbstractCallRender3D* cr3d, CallTriMeshData* tmc,
    megamol::core::Call& call, core::view::CallGetTransferFunction* gtf) {

    if (tmc == NULL) return false;
    misc::VolumetricDataCall* vdc = this->volumetricData_callerSlot.CallAs<misc::VolumetricDataCall>();


    float opacity_noInteraction;
    if (this->colorControl.transparency) {
        opacity_noInteraction = 0.2f;
    } else {
        opacity_noInteraction = 1.0f;
    }

    vec4 interactionColor = {1.0, 0.5, 0.0, 0.81f};
    if (this->colorControl.opaquePocket || !this->colorControl.transparency) {
        interactionColor.w = 1.0f;
    }

    vec4 noInteractionColor = {49.0f / 255.0f, 104.0f / 255.0f, 145.0f / 255.0f, opacity_noInteraction};


    //***************************************************************


    const megamol::geocalls::CallTriMeshData::Mesh& obj = tmc->Objects()[0];

    // array allocations
    if (tmc->DataHash() != this->old_tmc_dataHash) {
        // TODO: at the moment MSMSloader dont detect change in PDB loader
        // so tmc->DataHash() is always the same!
        if (this->old_tmc_dataHash != -1) {
            printf("resize proteinMeshRenderer VBOs");
            cudaGraphicsUnregisterResource(this->gl_vertices_resource);
            cudaGraphicsUnregisterResource(this->gl_indices_resource);
            cudaGraphicsUnregisterResource(this->gl_sortedIndices_resource);
            glDeleteBuffers(1, &this->gl_meshVertices);
            glDeleteBuffers(1, &this->gl_meshIndices);
            glDeleteBuffers(1, &this->gl_meshNormals);
            // glDeleteBuffers(1, &this->gl_meshVertices_ColorSwitch);
            // glDeleteBuffers(1, &this->gl_meshVertices_ClusterAssignment);
            glDeleteBuffers(1, &this->gl_meshSortedIndices);
            cudaFree(this->d_TriangleToCamDistances);
            cudaFree(this->d_camPos);
        }


        uint sizeForAllVertices = obj.GetVertexCount() * 3 * sizeof(float);
        uint sizeForAllColors = obj.GetVertexCount() * 3 * sizeof(unsigned char);
        glGenBuffers(1, &this->gl_meshVertices);
        glBindBuffer(GL_ARRAY_BUFFER, this->gl_meshVertices);
        glBufferData(GL_ARRAY_BUFFER, sizeForAllVertices, obj.GetVertexPointerFloat(), GL_STATIC_DRAW);
        cudaGraphicsGLRegisterBuffer(
            &this->gl_vertices_resource, this->gl_meshVertices, cudaGraphicsMapFlagsWriteDiscard);

        glGenBuffers(1, &this->gl_meshNormals);
        glBindBuffer(GL_ARRAY_BUFFER, this->gl_meshNormals);
        glBufferData(GL_ARRAY_BUFFER, sizeForAllVertices, obj.GetNormalPointerFloat(), GL_STATIC_DRAW);

        glGenBuffers(1, &this->gl_meshColors);
        glBindBuffer(GL_ARRAY_BUFFER, this->gl_meshColors);
        glBufferData(GL_ARRAY_BUFFER, sizeForAllColors, obj.GetColourPointerByte(), GL_STATIC_DRAW);

        if (this->old_tmc_dataHash == -1) {
            glGenBuffers(1, &this->gl_meshVertices_ColorSwitch);
            glGenBuffers(1, &this->gl_meshVertices_ClusterAssignment);
        }

        glGenBuffers(1, &this->gl_meshIndices);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->gl_meshIndices);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, obj.GetTriCount() * 3 * sizeof(GLuint), obj.GetTriIndexPointerUInt32(),
            GL_STATIC_DRAW);
        cudaGraphicsGLRegisterBuffer(
            &this->gl_indices_resource, this->gl_meshIndices, cudaGraphicsMapFlagsWriteDiscard);
        // gl_meshSortedIndices: Data will be uploaded every draw call ()
        glGenBuffers(1, &this->gl_meshSortedIndices);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->gl_meshSortedIndices);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, obj.GetTriCount() * 3 * sizeof(GLuint), obj.GetTriIndexPointerUInt32(),
            GL_DYNAMIC_DRAW);
        cudaGraphicsGLRegisterBuffer(
            &this->gl_sortedIndices_resource, this->gl_meshSortedIndices, cudaGraphicsMapFlagsWriteDiscard);

        cudaMalloc((void**)&this->d_TriangleToCamDistances, obj.GetTriCount() * sizeof(float));
        cudaMalloc((void**)&this->d_camPos, 1 * sizeof(float3));

        this->old_tmc_dataHash = tmc->DataHash();
    }


    // interaction area detection and coloring
    if ((this->s.old_SelectedIDx != this->s.selectedIDx ||
            this->old_clusterDataHash != this->clusterData.clusterDataHash) &&
            this->clusterData.clusterDataHash > 0 ||
        this->s.changedCurForceType) {
        if (this->old_clusterDataHash != this->clusterData.clusterDataHash || this->s.changedCurForceType) {
            /********************************************
             ********** FIND INTERACTION-AREAS **********
             ********************************************/
            if (this->old_clusterDataHash != this->clusterData.clusterDataHash) {
                this->old_clusterDataHash = this->clusterData.clusterDataHash;
                getInteractionAreaTriangles();
                this->atomPostions.Resize(0);
            }


            glDeleteBuffers(1, &this->gl_meshVertices_ColorSwitch);
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, this->gl_meshVertices_ColorSwitch);
            glBufferStorage(GL_SHADER_STORAGE_BUFFER, obj.GetVertexCount() * 1 * sizeof(int),
                this->markedVerticesIndicesMap[s.curForceType_surfaceColoring].data(), GL_DYNAMIC_STORAGE_BIT);
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

            // get texture
            // (*gtf)(0) sets also range values to 0,1
            gtf->SetRange({0.f, (float)this->maxForceCntPerProtAtomMap[s.curForceType_surfaceColoring]});

            // setInteractionResiduesMDC(this->hBondAtomCnt);

            /********************************************/
        }
    }


    LigandModelCall* lmc = this->ligandModel_callerSlot.CallAs<LigandModelCall>();

    /**************************************
     ********** PREPARE AO VOLUME **********
     ***************************************/

    if (!(*vdc)(core::misc::VolumetricDataCall::IDX_GET_EXTENTS)) return false;
    if (!(*vdc)(core::misc::VolumetricDataCall::IDX_GET_METADATA)) return false;
    if (!(*vdc)(core::misc::VolumetricDataCall::IDX_GET_DATA)) return false;

    if (this->old_vdc_DataHash != vdc->DataHash()) {
        this->old_vdc_DataHash = vdc->DataHash();
        auto const rawData = vdc->GetData();
        this->volumeMetaData = vdc->GetMetadata();

        GLenum internal_format;
        GLenum format;
        GLenum type;
        format = GL_RED;
        type = GL_FLOAT;
        internal_format = GL_R32F;
        glowl::TextureLayout volume_layout(internal_format, this->volumeMetaData->Resolution[0],
            this->volumeMetaData->Resolution[1], this->volumeMetaData->Resolution[2], format, type, 1,
            {{GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER}, {GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER},
                {GL_TEXTURE_WRAP_R, GL_CLAMP_TO_BORDER}, {GL_TEXTURE_MIN_FILTER, GL_LINEAR},
                {GL_TEXTURE_MAG_FILTER, GL_LINEAR}},
            {});

        aoVolumeTexture = std::make_unique<glowl::Texture3D>("AO_volume_texture", volume_layout, rawData);
    }

    /*****************************
    ****** MESH RENDERING ********
    ******************************/

    /*
    float amb, diff, spec, expo;
    amb = 0.2f;
    diff = 0.8f;
    spec = 1.0f;
    expo = 2.0f;
    glm::vec4 ambient(amb, amb, amb, 1.f);
    glm::vec4 diffuse(diff, diff, diff, 1.f);
    glm::vec4 specular(spec, spec, spec, 1.f);
    */

    this->triangleMeshShader.Enable();

    gtf->BindConvenience(this->triangleMeshShader, GL_TEXTURE1, 1);

    glm::mat4 view;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            view[i][j] = this->modelMatrix->PeekComponents()[i * 4 + j];
        }
    }

    glm::mat4 proj;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            proj[i][j] = this->projMatrix->PeekComponents()[i * 4 + j];
        }
    }

    glm::mat4 MVinv = glm::inverse(view);

    glUniformMatrix4fvARB(
        this->triangleMeshShader.ParameterLocation("modelView"), 1, false, this->modelMatrix->PeekComponents());
    glUniformMatrix4fvARB(
        this->triangleMeshShader.ParameterLocation("projection"), 1, false, this->projMatrix->PeekComponents());
    /*
    GLfloat invMV[16];
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            invMV[i * 4 + j] = this->transformMatrices->camera_to_local[j][i];
            // invMV[i * 4 + j] = 1.0;
        }
    }
    */
    glUniformMatrix4fvARB(this->triangleMeshShader.ParameterLocation("MVinv"), 1, false, glm::value_ptr(MVinv));
    glUniform3fv(this->triangleMeshShader.ParameterLocation("ligthSource"), 1, value_ptr(this->camPos_local));


    ::glEnable(GL_DEPTH_TEST);
    ::glEnable(GL_LIGHT0);
    ::glEnable(GL_BLEND);
    ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    ::glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    //::glEnable(GL_COLOR_MATERIAL);
    ::glDisable(GL_TEXTURE_2D);
    // ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    //::glEnable(GL_NORMALIZE);


    GLint cfm;
    ::glGetIntegerv(GL_CULL_FACE_MODE, &cfm);
    GLint pm[2];
    ::glGetIntegerv(GL_POLYGON_MODE, pm);
    GLint twr;
    ::glGetIntegerv(GL_FRONT_FACE, &twr);
    ::glFrontFace(GL_CCW);

    int fpm, bpm;
    fpm = GL_FILL;
    bpm = GL_FILL;

    ::glPolygonMode(GL_FRONT, fpm);
    ::glPolygonMode(GL_BACK, bpm);
    ::glDisable(GL_CULL_FACE);


    cudaMemcpy(this->d_camPos, &this->camPos_local, 1 * sizeof(float3), cudaMemcpyHostToDevice);


    cudaGraphicsUnregisterResource(this->gl_sortedIndices_resource);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->gl_meshSortedIndices);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, obj.GetTriCount() * 3 * sizeof(GLuint), obj.GetTriIndexPointerUInt32(),
        GL_DYNAMIC_DRAW);
    cudaGraphicsGLRegisterBuffer(
        &this->gl_sortedIndices_resource, this->gl_meshSortedIndices, cudaGraphicsMapFlagsWriteDiscard);


    size_t resourceSize;
    float* vertices;
    uint* indices;
    uint3* sortedIndices;
    uint triaCnt = obj.GetTriCount();

    int threads = NUM_THREADS;
    int blocks;
    threads = std::min(threads, (int)triaCnt);
    blocks = (triaCnt % threads != 0) ? ((int)triaCnt / threads + 1) : ((int)triaCnt / threads);

    cudaGraphicsMapResources(1, &gl_vertices_resource, 0);
    cudaGraphicsResourceGetMappedPointer((void**)(&vertices), &resourceSize, gl_vertices_resource);

    cudaGraphicsMapResources(1, &gl_indices_resource, 0);
    cudaGraphicsResourceGetMappedPointer((void**)(&indices), &resourceSize, gl_indices_resource);

    cudaGraphicsMapResources(1, &gl_sortedIndices_resource, 0);
    cudaGraphicsResourceGetMappedPointer((void**)(&sortedIndices), &resourceSize, gl_sortedIndices_resource);


    /* Call CUDA kernels
     * f(): get distances to camera
     * f(): sort from back to front
     */
    getTriangleToCamDistancesCUDA(blocks, threads, triaCnt, vertices, indices, d_camPos, d_TriangleToCamDistances);
    cudaDeviceSynchronize();
    SortTrianglesDevice(triaCnt, (uint3*)sortedIndices, d_TriangleToCamDistances);

    // Unmap VBOs.
    cudaGraphicsUnmapResources(1, &this->gl_vertices_resource, 0);
    cudaGraphicsUnmapResources(1, &this->gl_indices_resource, 0);
    cudaGraphicsUnmapResources(1, &this->gl_sortedIndices_resource, 0);
    cudaDeviceSynchronize();

    // bind 3Dvolume as texture for AO
    glActiveTexture(GL_TEXTURE0);
    glEnable(GL_TEXTURE_3D);
    aoVolumeTexture->bindTexture();
    glUniform1i(this->triangleMeshShader.ParameterLocation("aoVol"), 0);

    int sx = this->volumeMetaData->Resolution[0];
    int sy = this->volumeMetaData->Resolution[1];
    int sz = this->volumeMetaData->Resolution[2];

    float minOSx = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().Left();
    float minOSy = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().Bottom();
    float minOSz = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().Back();
    float rangeOSx = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().Width();
    float rangeOSy = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().Height();
    float rangeOSz = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().Depth();

    rangeOSx /= (1.0f - 2.0f / static_cast<float>(sx));
    rangeOSy /= (1.0f - 2.0f / static_cast<float>(sy));
    rangeOSz /= (1.0f - 2.0f / static_cast<float>(sz));

    minOSx -= rangeOSx / static_cast<float>(sx);
    minOSy -= rangeOSy / static_cast<float>(sy);
    minOSz -= rangeOSz / static_cast<float>(sz);

    // get Params
    this->aoSampFact = this->aoSampFactParam.Param<param::FloatParam>()->Value();
    // this->aoSampFactParam.Param<param::FloatParam>()->SetValue(this->aoSampFact);

    this->aoSampDist = this->aoSampDistParam.Param<param::FloatParam>()->Value();
    // this->aoSampDistParam.Param<param::FloatParam>()->SetValue(this->aoSampDist);

    this->aoIntensityMin = this->aoIntensityMinParam.Param<param::FloatParam>()->Value();
    // this->aoIntensityMinParam.Param<param::FloatParam>()->SetValue(this->aoIntensityMin);

    this->aoIntensityMax = this->aoIntensityMaxParam.Param<param::FloatParam>()->Value();
    // this->aoIntensityMaxParam.Param<param::FloatParam>()->SetValue(this->aoIntensityMax);

    // upload uniforms for AO
    this->triangleMeshShader.SetParameter("aoIntensityMin", this->aoIntensityMin);
    this->triangleMeshShader.SetParameter("aoIntensityMax", this->aoIntensityMax);
    this->triangleMeshShader.SetParameter("aoSampFact", this->aoSampFact);
    this->triangleMeshShader.SetParameter("posOrigin", minOSx, minOSy, minOSz);
    this->triangleMeshShader.SetParameter("posExtents", rangeOSx, rangeOSy, rangeOSz);
    this->triangleMeshShader.SetParameter("aoSampDist",
        this->aoSampDist * (vdc->AccessBoundingBoxes().ObjectSpaceClipBox().Width() / sx),
        this->aoSampDist * (vdc->AccessBoundingBoxes().ObjectSpaceClipBox().Height() / sy),
        this->aoSampDist * (vdc->AccessBoundingBoxes().ObjectSpaceClipBox().Depth() / sz));

    // bind VBOs and render them
    // Vetrices
    glEnableClientState(GL_VERTEX_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, this->gl_meshVertices);
    glVertexPointer(3, GL_FLOAT, 0, (void*)(0));
    // Normals
    glEnableClientState(GL_NORMAL_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, this->gl_meshNormals);
    glNormalPointer(GL_FLOAT, 0, (void*)(0));
    // Colors
    glEnableClientState(GL_COLOR_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, this->gl_meshColors);
    glColorPointer(3, GL_UNSIGNED_BYTE, 0, (void*)(0));
    // get texture

    if (this->s.changedcurClusterID || !lmc->getLMCStatus().isFirstRenderCallComplete ||
        this->clusterData.DBSCANParams.paramsChanged) {
        glDeleteBuffers(1, &this->gl_meshVertices_ClusterAssignment);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, this->gl_meshVertices_ClusterAssignment);
        glBufferStorage(GL_SHADER_STORAGE_BUFFER, obj.GetVertexCount() * 1 * sizeof(int),
            clusterAreaVertices[this->s.curClusterID + 1].PeekElements(), GL_DYNAMIC_STORAGE_BIT);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    }


    glBindBufferBase(
        GL_SHADER_STORAGE_BUFFER, this->meshVertices_ColorSwitch_SSBObindingPoint, this->gl_meshVertices_ColorSwitch);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, this->meshVertices_ClusterAssignment_SSBObindingPoint,
        this->gl_meshVertices_ClusterAssignment);

    this->triangleMeshShader.SetParameter(
        "interactionColor", interactionColor.x, interactionColor.y, interactionColor.z, interactionColor.w);
    this->triangleMeshShader.SetParameter(
        "noInteractionColor", noInteractionColor.x, noInteractionColor.y, noInteractionColor.z, noInteractionColor.w);
    this->triangleMeshShader.SetParameter("curClusterID", this->s.curClusterID + 1);
    this->triangleMeshShader.SetParameter("interactionForce", this->colorControl.interactionForce);
    this->triangleMeshShader.SetParameter(
        "interactionForceAmountColoring", this->colorControl.interactionForceAmountColoring);
    this->triangleMeshShader.SetParameter("pocketArea", this->colorControl.pocketArea);
    this->triangleMeshShader.SetParameter("chemicalProperty", this->colorControl.chemicalProperty);

    glUniform4fv(this->triangleMeshShader.ParameterLocation("simpleColor"), 1,
        &this->forceColorsMap[s.curForceType_surfaceColoring].x);

    // clipping plane
    float clipPoint[4];
    float clipCol[4];
    float clipPlane[4];


    this->getClipData(&this->clipDat, clipCol, clipPoint);
    clipPlane[0] = clipDat.GetX();
    clipPlane[1] = clipDat.GetY();
    clipPlane[2] = clipDat.GetZ();
    clipPlane[3] = clipDat.GetW();


    glUniform4fv(this->triangleMeshShader.ParameterLocation("move"), 1, clipPoint);
    glUniform4fv(this->triangleMeshShader.ParameterLocation("clipDat"), 1, clipPlane);
    glUniform4fv(this->triangleMeshShader.ParameterLocation("clipCol"), 1, clipCol);

    // Indices
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->gl_meshSortedIndices);

    // glColor4f(1.0f, 0.5f, 0.0f, 0.5f);

    // Stencil Magic
    glEnable(GL_CULL_FACE);
    glEnable(GL_STENCIL_TEST);
    glClear(GL_STENCIL_BUFFER_BIT);
    glDisable(GL_DEPTH_TEST);
    glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

    // first pass
    glStencilFunc(GL_ALWAYS, 0, 0);
    glStencilOp(GL_KEEP, GL_KEEP, GL_INCR);
    glCullFace(GL_FRONT); // render back faces only
    glDrawElements(GL_TRIANGLES, triaCnt * 3, GL_UNSIGNED_INT, (void*)(0));

    // second pass: decrement stencil buffer value on front faces
    glStencilOp(GL_KEEP, GL_KEEP, GL_DECR);
    glCullFace(GL_BACK); // render front faces only
    glDrawElements(GL_TRIANGLES, triaCnt * 3, GL_UNSIGNED_INT, (void*)(0));


    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
    glEnable(GL_DEPTH_TEST);
    glStencilFunc(GL_NOTEQUAL, 0, ~0);
    this->triangleMeshShader.Disable();

    /**************************
     ***** PANE RENDERING *****
     **************************/
    planeRenderer(&clipPoint[0], &clipPlane[0], &clipCol[0]);
    /**************************************/


    this->triangleMeshShader.Enable();
    // render the plane here
    glDisable(GL_STENCIL_TEST);

    glDisable(GL_CULL_FACE);
    // glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glDrawElements(GL_TRIANGLES, triaCnt * 3, GL_UNSIGNED_INT, (void*)(0));


    ::glBindBuffer(GL_ARRAY_BUFFER, 0);
    ::glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    ::glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    ::glDisable(GL_TEXTURE0);
    ::glDisableClientState(GL_VERTEX_ARRAY);
    ::glDisableClientState(GL_NORMAL_ARRAY);
    ::glDisableClientState(GL_COLOR_ARRAY);
    ::glDisable(GL_BLEND);
    ::glDisable(GL_TEXTURE_3D);
    gtf->UnbindConvenience();


    this->triangleMeshShader.Disable();

    /*
    vislib::Array<float> oneSphere;
    oneSphere.SetCount(32);

    oneSphere[0] = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().GetLeftBottomFront().GetX();
    oneSphere[1] = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().GetLeftBottomFront().GetY();
    oneSphere[2] = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().GetLeftBottomFront().GetZ();
    oneSphere[3] = 1.0f;

    oneSphere[4] =  vdc->AccessBoundingBoxes().ObjectSpaceClipBox().GetLeftBottomBack().GetX();
    oneSphere[5] = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().GetLeftBottomBack().GetY();
    oneSphere[6] = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().GetLeftBottomBack().GetZ();
    oneSphere[7] = 1.0f;

    oneSphere[8] = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().GetLeftTopFront().GetX();
    oneSphere[9] = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().GetLeftTopFront().GetY();
    oneSphere[10] = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().GetLeftTopFront().GetZ();
    oneSphere[11] = 1.0f;

    oneSphere[12] = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().GetLeftTopBack().GetX();
    oneSphere[13] = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().GetLeftTopBack().GetY();
    oneSphere[14] = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().GetLeftTopBack().GetZ();
    oneSphere[15] = 1.0f;

    oneSphere[16] = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().GetRightBottomBack().GetX();
    oneSphere[17] = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().GetRightBottomBack().GetY();
    oneSphere[18] = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().GetRightBottomBack().GetZ();
    oneSphere[19] = 1.0f;

    oneSphere[20] = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().GetRightBottomFront().GetX();
    oneSphere[21] = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().GetRightBottomFront().GetY();
    oneSphere[22] = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().GetRightBottomFront().GetZ();
    oneSphere[23] = 1.0f;

    oneSphere[24] = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().GetRightTopBack().GetX();
    oneSphere[25] = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().GetRightTopBack().GetY();
    oneSphere[26] = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().GetRightTopBack().GetZ();
    oneSphere[27] = 1.0f;

    oneSphere[28] = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().GetRightTopFront().GetX();
    oneSphere[29] = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().GetRightTopFront().GetY();
    oneSphere[30] = vdc->AccessBoundingBoxes().ObjectSpaceClipBox().GetRightTopFront().GetZ();
    oneSphere[31] = 1.0f;
    printf("cam: %f %f %f\n", oneSphere[4], oneSphere[5], oneSphere[6]);


    float viewportStuff[4] = {this->cameraInfo->TileRect().Left(), this->cameraInfo->TileRect().Bottom(),
        this->cameraInfo->TileRect().Width(), this->cameraInfo->TileRect().Height()};
    if (viewportStuff[2] < 1.0f) viewportStuff[2] = 1.0f;
    if (viewportStuff[3] < 1.0f) viewportStuff[3] = 1.0f;
    viewportStuff[2] = 2.0f / viewportStuff[2];
    viewportStuff[3] = 2.0f / viewportStuff[3];
    this->sphereRenderer(oneSphere, gtf, this->renVarsCNS, viewportStuff);
    */

    return true;
}


bool ModelClusterRenderer::CircleRenderer(
    vislib::Array<float>& SphereCirles, vislib::Array<float>& colors, float* viewportStuff) {
    if (SphereCirles.Count() > 1) {

        const GLfloat* viewportStuff_ = viewportStuff;
        this->circleShader.Enable();
        // set shader variables
        glUniform4fv(
            this->circleShader.ParameterLocation("lightPos"), 1, this->cameraInfo->Position().PeekCoordinates());
        glUniform4fvARB(this->circleShader.ParameterLocation("viewAttr"), 1, viewportStuff_);
        glUniform3fvARB(this->circleShader.ParameterLocation("camIn"), 1, this->cameraInfo->Front().PeekComponents());
        glUniform3fvARB(
            this->circleShader.ParameterLocation("camRight"), 1, this->cameraInfo->Right().PeekComponents());
        glUniform3fvARB(this->circleShader.ParameterLocation("camUp"), 1, this->cameraInfo->Up().PeekComponents());

        glUniformMatrix4fvARB(
            this->circleShader.ParameterLocation("modelview"), 1, false, this->modelMatrix->PeekComponents());
        glUniformMatrix4fvARB(
            this->circleShader.ParameterLocation("proj"), 1, false, this->projMatrix->PeekComponents());

        GLint vertexPos = glGetAttribLocation(this->circleShader, "vertex");
        GLint vertexColor = glGetAttribLocationARB(this->circleShader, "color");

        // set vertex and color pointers and draw them
        glEnableVertexAttribArray(vertexPos);
        glEnableVertexAttribArrayARB(vertexColor);

        glVertexAttribPointerARB(vertexPos, 4, GL_FLOAT, GL_FALSE, 0, SphereCirles.PeekElements());
        glVertexAttribPointerARB(vertexColor, 3, GL_FLOAT, GL_FALSE, 0, colors.PeekElements());

        glDrawArrays(GL_POINTS, 0, (SphereCirles.Count() / 4));

        glDisableVertexAttribArray(vertexPos);
        glDisableVertexAttribArrayARB(vertexColor);

        this->circleShader.Disable();

        return true;
    } else {
        return false;
    }
}


namespace {
glm::vec3 GetNormal(glm::vec3 v1, glm::vec3 v2, glm::vec3 v3) {
    vec3 a = v1 - v2;
    vec3 b = v3 - v2;
    return normalize(cross(a, b));
}

void pushVertex(vislib::Array<float>* vertices, glm::vec4 v) {
    vertices->Append(v.x);
    vertices->Append(v.y);
    vertices->Append(v.z);
    vertices->Append(v.w);
};


void constructFace(vislib::Array<float>* vertices, vislib::Array<float>* normals, glm::vec3 v1, glm::vec3 v2,
    glm::vec3 v3, glm::vec3 v4, glm::vec3 v5) {

    vec3 tmpN = GetNormal(v1, v2, v5);


    pushVertex(vertices, vec4(v1, 1.0));
    pushVertex(vertices, vec4(v2, 1.0));
    pushVertex(vertices, vec4(v5, 1.0));

    pushVertex(vertices, vec4(v5, 1.0));
    pushVertex(vertices, vec4(v2, 1.0));
    pushVertex(vertices, vec4(v4, 1.0));

    pushVertex(vertices, vec4(v4, 1.0));
    pushVertex(vertices, vec4(v2, 1.0));
    pushVertex(vertices, vec4(v3, 1.0));

    for (int i = 0; i < 9; i++) {
        normals->Append(tmpN.x);
        normals->Append(tmpN.y);
        normals->Append(tmpN.z);
    }
}

glm::vec3 setAndScalePoint(glm::vec3 v1, glm::vec3 point, float scale) {
    vec3 tmpDir = normalize(v1 - point);
    return point + tmpDir * (scale);
}

void constructDodecahedron(vislib::Array<float>* vertices, vislib::Array<float>* normals, glm::vec4& inPoint) {
    vec4 inPos = inPoint;

    float a = 1.0;

    float phi = 0.5 * (a + sqrt(5));
    vec3 point = vec3(inPoint);

    // phi = sqrt(2*rm/a);
    // https://en.wikipedia.org/wiki/Regular_dodecahedron

    // orange vertices
    // bottom rectangle
    vec3 v16 = vec3(-a, -a, -a) + point; // 1. orange
    vec3 v1 = vec3(-a, a, -a) + point;   // 2. orange
    vec3 v3 = vec3(a, a, -a) + point;    // 3. orange
    vec3 v15 = vec3(a, -a, -a) + point;  // 4. orange
                                         // top rectangle
    vec3 v18 = vec3(-a, -a, a) + point;  // 5. orange
    vec3 v7 = vec3(-a, a, a) + point;    // 6. orange
    vec3 v9 = vec3(a, a, a) + point;     // 7. orange
    vec3 v13 = vec3(a, -a, a) + point;   // 8. orange

    // green vertices
    vec3 v20 = vec3(0, -phi, -a / phi) + point; // 1. green
    vec3 v19 = vec3(0, -phi, a / phi) + point;  // 2. green
    vec3 v2 = vec3(0, phi, -a / phi) + point;   // 3. green
    vec3 v8 = vec3(0, phi, a / phi) + point;    // 4. green

    // blue vertices
    vec3 v5 = vec3(-a / phi, 0, -phi) + point; // 1. blue
    vec3 v4 = vec3(a / phi, 0, -phi) + point;  // 2. blue
    vec3 v12 = vec3(a / phi, 0, phi) + point;  // 3. blue
    vec3 v11 = vec3(-a / phi, 0, phi) + point; // 4. blue

    // pink vertices
    vec3 v17 = vec3(-phi, -a / phi, 0) + point; // 1. pink
    vec3 v14 = vec3(phi, -a / phi, 0) + point;  // 2. pink
    vec3 v10 = vec3(phi, a / phi, 0) + point;   // 3. pink
    vec3 v6 = vec3(-phi, a / phi, 0) + point;   // 4. pink

    v1 = setAndScalePoint(v1, point, inPos.w);
    v2 = setAndScalePoint(v2, point, inPos.w);
    v3 = setAndScalePoint(v3, point, inPos.w);
    v4 = setAndScalePoint(v4, point, inPos.w);
    v5 = setAndScalePoint(v5, point, inPos.w);
    v6 = setAndScalePoint(v6, point, inPos.w);
    v7 = setAndScalePoint(v7, point, inPos.w);
    v8 = setAndScalePoint(v8, point, inPos.w);
    v9 = setAndScalePoint(v9, point, inPos.w);
    v10 = setAndScalePoint(v10, point, inPos.w);
    v11 = setAndScalePoint(v11, point, inPos.w);
    v12 = setAndScalePoint(v12, point, inPos.w);
    v13 = setAndScalePoint(v13, point, inPos.w);
    v14 = setAndScalePoint(v14, point, inPos.w);
    v15 = setAndScalePoint(v15, point, inPos.w);
    v16 = setAndScalePoint(v16, point, inPos.w);
    v17 = setAndScalePoint(v17, point, inPos.w);
    v18 = setAndScalePoint(v18, point, inPos.w);
    v19 = setAndScalePoint(v19, point, inPos.w);
    v20 = setAndScalePoint(v20, point, inPos.w);


    /// 12 pentagon patches : per triangle "CCW" -> count clock wise
    // 1 Patch 1,5,4,3,2
    constructFace(vertices, normals, v1, v5, v4, v3, v2);

    // 2 Patch 8,2,3,10,9
    constructFace(vertices, normals, v8, v2, v3, v10, v9);

    // 3 Patch 6,1,2,8,7
    constructFace(vertices, normals, v6, v1, v2, v8, v7);

    // 4 Patch 3,4,15,14,10
    constructFace(vertices, normals, v3, v4, v15, v14, v10);

    // 5 Patch 7,8,9,12,11
    constructFace(vertices, normals, v7, v8, v9, v12, v11);

    // 6 Patch 9,10,14,13,12
    constructFace(vertices, normals, v9, v10, v14, v13, v12);

    // 7 Patch 18,17,6,7,11
    constructFace(vertices, normals, v18, v17, v6, v7, v11);

    // 8 Patch 15,4,5,16,20
    constructFace(vertices, normals, v15, v4, v5, v16, v20);

    // 9 Patch 14,15,20,19,13
    constructFace(vertices, normals, v14, v15, v20, v19, v13);

    // 10 Patch 17,16,5,1,6
    constructFace(vertices, normals, v17, v16, v5, v1, v6);

    // 11 Patch 19,20,16,17,18
    constructFace(vertices, normals, v19, v20, v16, v17, v18);

    // 12 Patch 18,11,12,13,19
    constructFace(vertices, normals, v18, v11, v12, v13, v19);
}
} // namespace

bool ModelClusterRenderer::dodecahedronRenderer(vislib::Array<float>& spheres, vislib::Array<float>* vertices,
    vislib::Array<float>* normals, core::view::CallGetTransferFunction* gtf, SphereRendererVars& renVars,
    float* viewportStuff) {
    if (spheres.Count() > 1) {
        const GLfloat* viewportStuff_ = viewportStuff;


        if (!vertices->Count() > 0) {
            for (int i = 0; i < spheres.Count() / 4; i++) {
                vec4 tpmPoint = vec4(spheres[i * 4 + 0], spheres[i * 4 + 1], spheres[i * 4 + 2], spheres[i * 4 + 3]);
                constructDodecahedron(vertices, normals, tpmPoint);
            }
        }


        this->dodecahedron_Shader.Enable();

        // set shader variables
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, renVars.sphereValues_SSBObindingPoint, renVars.gl_sphereValuesSSBO);

        gtf->BindConvenience(this->dodecahedron_Shader, GL_TEXTURE2, 2);
        glUniform4fv(
            this->dodecahedron_Shader.ParameterLocation("lightPos"), 1, this->cameraInfo->Position().PeekCoordinates());
        glUniform4fvARB(this->dodecahedron_Shader.ParameterLocation("viewAttr"), 1, viewportStuff_);
        glUniform3fvARB(
            this->dodecahedron_Shader.ParameterLocation("camIn"), 1, this->cameraInfo->Front().PeekComponents());
        glUniform4fv(this->dodecahedron_Shader.ParameterLocation("hoverColor"), 1, renVars.hoverColor.PeekComponents());
        glUniform4fv(
            this->dodecahedron_Shader.ParameterLocation("selectionColor"), 1, renVars.selectionColor.PeekComponents());
        glUniform4fv(this->dodecahedron_Shader.ParameterLocation("baseColor"), 1, renVars.baseColor.PeekComponents());
        glUniform1i(this->dodecahedron_Shader.ParameterLocation("hoverSphereID"), renVars.hoverSphereID);
        glUniform1i(this->dodecahedron_Shader.ParameterLocation("selectedSphereID"), renVars.selectedSphereID);
        glUniform1i(this->dodecahedron_Shader.ParameterLocation("hideOnHover"), renVars.hideOnHover);
        glUniform1i(this->dodecahedron_Shader.ParameterLocation("hideOnSelection"), renVars.hideOnSelection);
        glUniform3fvARB(
            this->dodecahedron_Shader.ParameterLocation("camRight"), 1, this->cameraInfo->Right().PeekComponents());
        glUniform3fvARB(
            this->dodecahedron_Shader.ParameterLocation("camUp"), 1, this->cameraInfo->Up().PeekComponents());

        glUniformMatrix4fvARB(
            this->dodecahedron_Shader.ParameterLocation("modVIEW"), 1, false, this->modelMatrix->PeekComponents());
        glUniformMatrix4fvARB(
            this->dodecahedron_Shader.ParameterLocation("proj"), 1, false, this->projMatrix->PeekComponents());

        GLint vertexPos = glGetAttribLocation(this->dodecahedron_Shader, "vertex");
        GLint normalIn = glGetAttribLocation(this->dodecahedron_Shader, "normal");

        // set vertex and color pointers and draw them
        glEnableVertexAttribArray(vertexPos);
        glVertexAttribPointerARB(vertexPos, 4, GL_FLOAT, GL_FALSE, 0, vertices->PeekElements());

        glEnableVertexAttribArray(normalIn);
        glVertexAttribPointerARB(normalIn, 3, GL_FLOAT, GL_FALSE, 0, normals->PeekElements());

        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
        // glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDrawArrays(GL_TRIANGLES, 0, (vertices->Count() / 4));

        glDisableVertexAttribArray(vertexPos);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
        gtf->UnbindConvenience();

        this->dodecahedron_Shader.Disable();

        return true;
    } else {
        return false;
    }
}


bool ModelClusterRenderer::sphereRenderer(float* spheres, int sphereCnt, core::view::CallGetTransferFunction* gtf,
    SphereRendererVars& renVars, float* viewportStuff) {
    if (sphereCnt > 1) {
        const GLfloat* viewportStuff_ = viewportStuff;

        this->sphereShader.Enable();
        // set shader variables
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, renVars.sphereValues_SSBObindingPoint, renVars.gl_sphereValuesSSBO);

        gtf->BindConvenience(this->sphereShader, GL_TEXTURE2, 2);
        glUniform4fv(
            this->sphereShader.ParameterLocation("lightPos"), 1, this->cameraInfo->Position().PeekCoordinates());
        glUniform4fvARB(this->sphereShader.ParameterLocation("viewAttr"), 1, viewportStuff_);
        glUniform3fvARB(this->sphereShader.ParameterLocation("camIn"), 1, this->cameraInfo->Front().PeekComponents());
        glUniform4fv(this->sphereShader.ParameterLocation("hoverColor"), 1, renVars.hoverColor.PeekComponents());
        glUniform4fv(
            this->sphereShader.ParameterLocation("selectionColor"), 1, renVars.selectionColor.PeekComponents());
        glUniform4fv(this->sphereShader.ParameterLocation("baseColor"), 1, renVars.baseColor.PeekComponents());
        glUniform1i(this->sphereShader.ParameterLocation("hoverSphereID"), renVars.hoverSphereID);
        glUniform1i(this->sphereShader.ParameterLocation("selectedSphereID"), renVars.selectedSphereID);
        glUniform1i(this->sphereShader.ParameterLocation("hideOnHover"), renVars.hideOnHover);
        glUniform1i(this->sphereShader.ParameterLocation("hideOnSelection"), renVars.hideOnSelection);
        glUniform3fvARB(
            this->sphereShader.ParameterLocation("camRight"), 1, this->cameraInfo->Right().PeekComponents());
        glUniform3fvARB(this->sphereShader.ParameterLocation("camUp"), 1, this->cameraInfo->Up().PeekComponents());

        glUniformMatrix4fvARB(
            this->sphereShader.ParameterLocation("modelview"), 1, false, this->modelMatrix->PeekComponents());
        glUniformMatrix4fvARB(
            this->sphereShader.ParameterLocation("proj"), 1, false, this->projMatrix->PeekComponents());

        GLint vertexPos = glGetAttribLocation(this->sphereShader, "vertex");

        // set vertex and color pointers and draw them
        glEnableVertexAttribArray(vertexPos);

        glVertexAttribPointerARB(vertexPos, 4, GL_FLOAT, GL_FALSE, 0, spheres);

        glDrawArrays(GL_POINTS, 0, (sphereCnt));

        glDisableVertexAttribArray(vertexPos);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
        gtf->UnbindConvenience();
        this->sphereShader.Disable();

        return true;
    } else {
        return false;
    }
}

bool ModelClusterRenderer::circleTextureRenderer(vislib::Array<float>& spheres,
    core::view::CallGetTransferFunction* gtf, SphereRendererVars& renVars, float* viewportStuff) {
    if (spheres.Count() > 1) {
        const GLfloat* viewportStuff_ = viewportStuff;
        this->circleBillTextureShader.Enable();
        /////////////////////////////////////////////////////////////////

        // if not done load textures
        if (this->modelModeIcons.Count() <= 0) {
            loadTexture("pocket_single_gray.png");
            loadTexture("pocket_multiple_gray.png");
            loadTexture("pocket_multiple_plus_protein.png");
        }

        glEnable(GL_TEXTURE_2D);
        // max(i) should be 3
        for (int i = 0; i < this->modelModeIcons.Count(); i++) {
            glActiveTexture(GL_TEXTURE0 + i);
            this->modelModeIcons[i]->Bind();
        }

        // set texture binding points
        glUniform1i(this->circleBillTextureShader.ParameterLocation("icon0"), 0);
        glUniform1i(this->circleBillTextureShader.ParameterLocation("icon1"), 1);
        glUniform1i(this->circleBillTextureShader.ParameterLocation("icon2"), 2);

        /////////////////////////////////////////////////////////////////


        // set shader variables
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, renVars.sphereValues_SSBObindingPoint, renVars.gl_sphereValuesSSBO);
        gtf->BindConvenience(this->circleBillTextureShader, GL_TEXTURE3, 3);

        glUniform4fv(
            this->circleBillTextureShader.ParameterLocation("hoverColor"), 1, renVars.hoverColor.PeekComponents());
        glUniform4fv(this->circleBillTextureShader.ParameterLocation("selectionColor"), 1,
            renVars.selectionColor.PeekComponents());
        glUniform4fv(
            this->circleBillTextureShader.ParameterLocation("baseColor"), 1, renVars.baseColor.PeekComponents());
        glUniform1i(this->circleBillTextureShader.ParameterLocation("hoverSphereID"), renVars.hoverSphereID);
        glUniform1i(this->circleBillTextureShader.ParameterLocation("selectedSphereID"), renVars.selectedSphereID);
        glUniform1i(this->circleBillTextureShader.ParameterLocation("hideOnHover"), renVars.hideOnHover);
        glUniform1i(this->circleBillTextureShader.ParameterLocation("hideOnSelection"), renVars.hideOnSelection);

        glUniform4fv(this->circleBillTextureShader.ParameterLocation("camPosition"), 1, value_ptr(this->camPos_local));
        glUniform3fvARB(
            this->circleBillTextureShader.ParameterLocation("camIn"), 1, this->cameraInfo->Front().PeekComponents());
        glUniform3fvARB(
            this->circleBillTextureShader.ParameterLocation("camRight"), 1, this->cameraInfo->Right().PeekComponents());
        glUniform3fvARB(
            this->circleBillTextureShader.ParameterLocation("camUp"), 1, this->cameraInfo->Up().PeekComponents());
        glUniformMatrix4fvARB(this->circleBillTextureShader.ParameterLocation("modelview"), 1, false,
            this->modelMatrix->PeekComponents());
        glUniformMatrix4fvARB(
            this->circleBillTextureShader.ParameterLocation("proj"), 1, false, this->projMatrix->PeekComponents());

        GLint vertexPos = glGetAttribLocation(this->circleBillTextureShader, "vertex");

        // set vertex and color pointers and draw them
        glEnableVertexAttribArray(vertexPos);

        glVertexAttribPointerARB(vertexPos, 4, GL_FLOAT, GL_FALSE, 0, spheres.PeekElements());

        glDrawArrays(GL_POINTS, 0, (spheres.Count() / 4));

        glDisableVertexAttribArray(vertexPos);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
        gtf->UnbindConvenience();
        glDisable(GL_TEXTURE_2D);
        this->circleBillTextureShader.Disable();

        return true;
    } else {
        return false;
    }
}


bool ModelClusterRenderer::coneRenderer(vislib::Array<float>& positions, vislib::Array<float>& directions,
    vislib::Array<float>& coneLengths, vislib::Array<float>& coneRadii, glm::vec4 color,
    ModelClusterRenderer::SphereRendererVars& renVars, float* viewportStuff) {
    if (positions.Count() > 1) {
        const GLfloat* viewportStuff_ = viewportStuff;

        // Camera
        /*
        view::Camera_2 cam;
        call.GetCamera(cam);
        cam_type::snapshot_type snapshot;
        cam_type::matrix_type viewTemp, projTemp;
        cam.calc_matrices(snapshot, viewTemp, projTemp, thecam::snapshot_content::all);
        glm::vec4 cam_view = snapshot.view_vector;
        glm::vec4 cam_right = snapshot.right_vector;
        glm::vec4 cam_up = snapshot.up_vector;
        */

        // Matrices
        glm::mat4 view;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                view[i][j] = this->modelMatrix->PeekComponents()[i * 4 + j];
            }
        }

        glm::mat4 proj;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                proj[i][j] = this->projMatrix->PeekComponents()[i * 4 + j];
            }
        }

        glm::mat4 MVinv = glm::inverse(view);
        glm::mat4 MVtransp = glm::transpose(view);
        glm::mat4 MVP = proj * view;
        glm::mat4 MVPinv = glm::inverse(MVP);
        glm::mat4 MVPtransp = glm::transpose(MVP);

        float lengthFilter = 0.0f;

        glDisable(GL_BLEND);

        glEnable(GL_DEPTH_TEST);
        glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);


        this->coneShader.Enable();
        glUniformMatrix4fv(this->coneShader.ParameterLocation("MVinv"), 1, GL_FALSE, glm::value_ptr(MVinv));
        glUniformMatrix4fv(this->coneShader.ParameterLocation("MVtransp"), 1, GL_FALSE, glm::value_ptr(MVtransp));
        glUniformMatrix4fv(this->coneShader.ParameterLocation("MVP"), 1, GL_FALSE, glm::value_ptr(MVP));
        glUniformMatrix4fv(this->coneShader.ParameterLocation("MVPinv"), 1, GL_FALSE, glm::value_ptr(MVPinv));
        glUniformMatrix4fv(this->coneShader.ParameterLocation("MVPtransp"), 1, GL_FALSE, glm::value_ptr(MVPtransp));
        glUniform4fvARB(this->coneShader.ParameterLocation("viewAttr"), 1, viewportStuff_);
        glUniform4fv(this->coneShader.ParameterLocation("lightDir"), 1, value_ptr(this->camPos_local));
        this->coneShader.SetParameter("lengthFilter", lengthFilter);

        // glUniform1i(this->coneShader.ParameterLocation("drawCylinder"), renVars.hideOnHover);
        // glUniform4fv(this->coneShader.ParameterLocation("clipDat"), 1, clipDat);
        // glUniform3fv(this->coneShader.ParameterLocation("clipCol"), 1, clipCol);


        unsigned int cial = glGetAttribLocationARB(this->coneShader, "colIdx");
        unsigned int tpal = glGetAttribLocationARB(this->coneShader, "dir");
        unsigned int coneLength = glGetAttribLocationARB(this->coneShader, "coneLength");
        unsigned int coneRadius = glGetAttribLocationARB(this->coneShader, "coneRadius");
        bool useFlags = false;
        glUniform1ui(this->coneShader.ParameterLocation("flagsAvailable"), useFlags ? 1 : 0);


        float minC = 0.0f, maxC = 0.0f;
        unsigned int colTabSize = 0;

        // colour
        glColor3f(color[0], color[1], color[2]);


        /*
        glEnable(GL_TEXTURE_1D);
        view::CallGetTransferFunction* cgtf = this->getTFSlot.CallAs<view::CallGetTransferFunction>();
        if ((cgtf != nullptr) && ((*cgtf)())) {
            glBindTexture(GL_TEXTURE_1D, cgtf->OpenGLTexture());
            colTabSize = cgtf->TextureSize();
        } else {
            glBindTexture(GL_TEXTURE_1D, this->greyTF);
            colTabSize = 2;
        }
        */


        // radius and position
        glEnableClientState(GL_VERTEX_ARRAY);
        glUniform4f(this->coneShader.ParameterLocation("inConsts1"), -1.0f, minC, maxC, float(colTabSize));
        glVertexPointer(3, GL_FLOAT, 0, positions.PeekElements());

        // direction
        glEnableVertexAttribArrayARB(tpal);
        glVertexAttribPointerARB(tpal, 3, GL_FLOAT, GL_FALSE, 0, directions.PeekElements());

        // coneLength
        glEnableVertexAttribArrayARB(coneLength);
        glVertexAttribPointerARB(coneLength, 1, GL_FLOAT, GL_FALSE, 0, coneLengths.PeekElements());

        // conrRadius
        glEnableVertexAttribArrayARB(coneRadius);
        glVertexAttribPointerARB(coneRadius, 1, GL_FLOAT, GL_FALSE, 0, coneRadii.PeekElements());


        // draw call
        glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(positions.Count() / 3));


        // reset/unbind things
        glColorPointer(3, GL_FLOAT, 0, nullptr);
        glVertexPointer(5, GL_FLOAT, 0, nullptr);
        glDisableClientState(GL_COLOR_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);

        glVertexAttribPointerARB(tpal, 4, GL_FLOAT, GL_FALSE, 0, nullptr);
        glDisableVertexAttribArrayARB(tpal);
        glDisable(GL_TEXTURE_1D);
    }


    this->coneShader.Disable();


    glDisable(GL_DEPTH_TEST);
    glDisable(GL_VERTEX_PROGRAM_POINT_SIZE);

    return true;
}

bool ModelClusterRenderer::cylinderRenderer(float* positions, float* directions, float* coneLengths, int cylinderCnt,
    ModelClusterRenderer::SphereRendererVars& renVars, float* viewportStuff) {
    if (cylinderCnt > 0) {
        const GLfloat* viewportStuff_ = viewportStuff;

        // Matrices
        glm::mat4 view;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                view[i][j] = this->modelMatrix->PeekComponents()[i * 4 + j];
            }
        }

        glm::mat4 proj;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                proj[i][j] = this->projMatrix->PeekComponents()[i * 4 + j];
            }
        }

        glm::mat4 MVinv = glm::inverse(view);
        glm::mat4 MVtransp = glm::transpose(view);
        glm::mat4 MVP = proj * view;
        glm::mat4 MVPinv = glm::inverse(MVP);
        glm::mat4 MVPtransp = glm::transpose(MVP);

        float lengthFilter = 0.0f;

        glDisable(GL_BLEND);

        glEnable(GL_DEPTH_TEST);
        glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);


        this->cylinderShader.Enable();
        glUniformMatrix4fv(this->cylinderShader.ParameterLocation("MVinv"), 1, GL_FALSE, glm::value_ptr(MVinv));
        glUniformMatrix4fv(this->cylinderShader.ParameterLocation("MVtransp"), 1, GL_FALSE, glm::value_ptr(MVtransp));
        glUniformMatrix4fv(this->cylinderShader.ParameterLocation("MVP"), 1, GL_FALSE, glm::value_ptr(MVP));
        glUniformMatrix4fv(this->cylinderShader.ParameterLocation("MVPinv"), 1, GL_FALSE, glm::value_ptr(MVPinv));
        glUniformMatrix4fv(this->cylinderShader.ParameterLocation("MVPtransp"), 1, GL_FALSE, glm::value_ptr(MVPtransp));
        glUniform4fvARB(this->cylinderShader.ParameterLocation("viewAttr"), 1, viewportStuff_);
        glUniform4fv(
            this->cylinderShader.ParameterLocation("lightDir"), 1, this->cameraInfo->EyeDirection().PeekComponents());
        this->cylinderShader.SetParameter("lengthFilter", lengthFilter);

        unsigned int cial = glGetAttribLocationARB(this->cylinderShader, "colIdx");
        unsigned int tpal = glGetAttribLocationARB(this->cylinderShader, "dir");
        unsigned int coneLength = glGetAttribLocationARB(this->cylinderShader, "coneLength");
        unsigned int coneRadius = glGetAttribLocationARB(this->cylinderShader, "coneRadius");
        bool useFlags = false;
        glUniform1ui(this->cylinderShader.ParameterLocation("flagsAvailable"), useFlags ? 1 : 0);


        float minC = 0.0f, maxC = 0.0f;
        unsigned int colTabSize = 0;

        // colour
        glColor3f(renVars.baseColor.GetX(), renVars.baseColor.GetY(), renVars.baseColor.GetZ());


        // radius and position
        glEnableClientState(GL_VERTEX_ARRAY);
        glUniform4f(this->cylinderShader.ParameterLocation("inConsts1"), -1.0f, minC, maxC, float(colTabSize));
        glVertexPointer(4, GL_FLOAT, 0, positions);
        //(x,y,z,radius,x,y,z,radius,x,y,z,radius)


        // direction
        glEnableVertexAttribArrayARB(tpal);
        glVertexAttribPointerARB(tpal, 3, GL_FLOAT, GL_FALSE, 0, directions);

        // coneLength
        glEnableVertexAttribArrayARB(coneLength);
        glVertexAttribPointerARB(coneLength, 1, GL_FLOAT, GL_FALSE, 0, coneLengths);

        // conrRadius
        // glEnableVertexAttribArrayARB(coneRadius);
        // glVertexAttribPointerARB(coneRadius, 1, GL_FLOAT, GL_FALSE, 16, positions);


        // draw call
        glDrawArrays(GL_POINTS, 0, cylinderCnt);


        // reset/unbind things
        glColorPointer(3, GL_FLOAT, 0, nullptr);
        glVertexPointer(5, GL_FLOAT, 0, nullptr);
        glDisableClientState(GL_COLOR_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);

        glVertexAttribPointerARB(tpal, 4, GL_FLOAT, GL_FALSE, 0, nullptr);
        glDisableVertexAttribArrayARB(tpal);
        glDisable(GL_TEXTURE_1D);
    }


    this->cylinderShader.Disable();


    glDisable(GL_DEPTH_TEST);
    glDisable(GL_VERTEX_PROGRAM_POINT_SIZE);

    return true;
}

bool ModelClusterRenderer::barChartRenderer(std::vector<std::map<string, float>>& clusterForceAreas,
    std::vector<std::map<string, std::map<uint, uint>>>& clusterForceResidueCnts, std::vector<float>& clusterAreas,
    megamol::core::Call* ligandModelCall, int ID, float liftSize) {

    float sizeFactor = 0.85f;
    float liftLength = liftSize;

    vislib::Array<vislib::math::Vector<float, 3>> bars;
    std::vector<float*> barsColors;
    LigandModelCall* lmc = dynamic_cast<LigandModelCall*>(ligandModelCall);
    if (lmc == NULL) return false;

    // top data
    if (clusterForceAreas.size() == 0) return false;
    if (clusterAreas.size() == 0) return false;

    // bottom data
    if (clusterForceResidueCnts.size() == 0) return false;

    // get the camera Up vector
    auto c = this->cameraInfo;
    glm::vec3 rightVec = glm::vec3(c->Right().GetX(), c->Right().GetY(), c->Right().GetZ());
    glm::vec3 eyeDirection = glm::vec3(c->EyeDirection().GetX(), c->EyeDirection().GetY(), c->EyeDirection().GetZ());
    glm::vec3 CameraUp = vec3(c->Up().GetX(), c->Up().GetY(), c->Up().GetZ());
    glm::vec3 upVec = glm::cross(rightVec, eyeDirection);
    upVec = normalize(upVec);

    float barWidth = 0.65f * sizeFactor;
    float axisHight = 0.15f * sizeFactor;
    glm::vec4 axisColor = glm::vec4(0.f, 0.f, 0.f, 1.0f);
    glm::vec4 connectorColor = glm::vec4(0.302f, 0.686f, 0.29f, 1.0f);

    // get the bar-chart center point
    // for (int i = 0; i < lmc->getClusterAndCentroidData().clusterCentroids->Count() / 4; i++) {
    int barCnt = clusterForceResidueCnts[ID].size();
    float spacer = barWidth * 2.0 + (0.2f * sizeFactor);
    auto a = *lmc->getClusterAndCentroidData().clusterCentroids;
    glm::vec3 origPos = glm::vec3(a[ID * 4 + 0], a[ID * 4 + 1], a[ID * 4 + 2]);

    barsColors.push_back(&connectorColor.x);
    getLineCoords(origPos, liftLength - 6.5f, axisHight, upVec, eyeDirection, bars, 0);

    glm::vec3 centerPoint = glm::vec3(a[ID * 4 + 0], a[ID * 4 + 1], a[ID * 4 + 2]) + upVec * liftLength;

    int cnter = 0;
    if (barCnt % 2 == 0) {
        centerPoint = centerPoint - rightVec * (spacer * (float)std::ceil(barCnt / 2) - barWidth / 2.f);
    } else {
        centerPoint = centerPoint - rightVec * (spacer * (float)std::ceil(barCnt / 2));
    }

    // create x-axis
    float axeisLength = spacer * (barCnt);
    glm::vec3 axisStartPoint = (centerPoint - rightVec * barWidth);
    getLineCoords(axisStartPoint, axeisLength, axisHight, rightVec, eyeDirection, bars, 0);
    barsColors.push_back(&axisColor.x);

    //	int ID = drawOrder[i].GetX();


    glm::vec3 pointA;
    // get max top
    float maxTop = 0;
    float minTop = 0;
    std::map<string, std::map<uint, uint>>::iterator itForce = clusterForceResidueCnts[ID].begin();
    while (itForce != clusterForceResidueCnts[ID].end()) {

        maxTop = itForce->second.size() > maxTop ? itForce->second.size() : maxTop;
        minTop = itForce->second.size() < minTop ? itForce->second.size() : minTop;
        itForce++;
    }

    float rangeTop = (maxTop - minTop);
    float scale = 5.f;
    if (rangeTop > 0) {
        scale = 5.f / rangeTop;
    } else {
        minTop = 0;
    }


    // create top bars
    glm::vec3 centerPointTop = centerPoint + upVec * axisHight * 1.05f;
    itForce = clusterForceResidueCnts[ID].begin();
    while (itForce != clusterForceResidueCnts[ID].end()) {
        pointA = glm::vec3(centerPointTop + rightVec * spacer * (float)cnter);
        float length = ((float)itForce->second.size() - minTop) * scale * sizeFactor;

        barsColors.push_back(&this->forceColorsMap[itForce->first][0]);
        getLineCoords(pointA, length, barWidth, upVec, eyeDirection, bars, 0);

        itForce++;
        cnter++;
    }

    // get max bot
    float maxBot = 0;
    float minBot = 0;
    std::map<string, float>::iterator itForceAreaS = clusterForceAreas[ID].begin();
    while (itForceAreaS != clusterForceAreas[ID].end()) {

        maxBot = itForceAreaS->second > maxBot ? itForceAreaS->second : maxBot;
        minBot = itForceAreaS->second < minBot ? itForceAreaS->second : minBot;
        itForceAreaS++;
    }

    // create bottom bars
    cnter = 0;
    std::map<string, float>::iterator itForceArea = clusterForceAreas[ID].begin();
    float rangeBot = (maxBot - minBot);
    if (rangeBot > 0) {
        scale = 5.f / rangeBot;
    } else {
        minBot = 0;
    }

    glm::vec3 centerPointBot = centerPoint - upVec * axisHight * 1.05f;
    while (itForceArea != clusterForceAreas[ID].end()) {
        pointA = glm::vec3(centerPointBot + rightVec * spacer * (float)cnter);
        float length = (itForceArea->second - minBot) * scale * sizeFactor;

        barsColors.push_back(&this->forceColorsMap[itForceArea->first][0]);
        getLineCoords(pointA, length, barWidth, -upVec, eyeDirection, bars, 0);

        itForceArea++;
        cnter++;
    }
    //}


    // draw bar-chart
    glDisable(GL_LIGHTING);
    glEnable(GL_POLYGON_SMOOTH);
    glDisable(GL_DEPTH_TEST);
    for (int i = 0; i < bars.Count() / 4; i++) {
        glDisable(GL_DEPTH_TEST);

        // first element (i=0) is the leading line (this will be drawn as real 3D-Object)
        if (i == 0) {
            glEnable(GL_DEPTH_TEST);
        }
        glBegin(GL_TRIANGLE_STRIP);
        glColor4f(barsColors[i][0], barsColors[i][1], barsColors[i][2], barsColors[i][3]);
        for (int j = 0; j < 4; j++) {
            glVertex3f(bars[i * 4 + j].X(), bars[i * 4 + j].Y(), bars[i * 4 + j].Z());
        }
        glEnd();
    }
    glEnable(GL_LIGHTING);
    glDisable(GL_POLYGON_SMOOTH);
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);

    return true;
}


/**
 * determines which triangles of the protein mesh a part of interaction areas with models/ligands
 * (uses nanoflann for different neighbour searches)
 *
 * @param ??		CallTriMeshData
 * @param call		??
 */
bool ModelClusterRenderer::getInteractionAreaTriangles() {
    auto start_c = std::chrono::high_resolution_clock::now();
    printf("Calculating protein mesh colors...\n");
    if (this->clusterAreaVertices.Count() > 0) {
        for (int i = 0; i < this->clusterAreaVertices.Count(); i++) {
            this->clusterAreaVertices[i].Resize(0);
        }
        this->clusterAreaVertices.Resize(0);
    }
    this->clusterAreas.clear();
    this->clusterAreaAtomsToClusterIDs.clear();

    MolecularDataCall* pmdc = this->proteinMolecularData_callerSlot.CallAs<MolecularDataCall>();
    LigandModelCall* lmc = this->ligandModel_callerSlot.CallAs<LigandModelCall>();
    MolecularDataCall* mdc = this->ligandMolecularData_callerSlot.CallAs<MolecularDataCall>();
    CallTriMeshData* tmc = this->triMeshData_callerSlot.CallAs<CallTriMeshData>();
    const megamol::geocalls::CallTriMeshData::Mesh& obj = tmc->Objects()[0];

    if (pmdc == NULL) {
        return 0;
    }
    if (!(*pmdc)(MolecularDataCall::CallForGetData)) {
        return 0;
    }

    // handle change of protein data
    if (old_pmdc_DataHash != pmdc->DataHash()) {
        this->atomToVerticesList.Resize(0);
        this->atomToVerticesList.SetCount(pmdc->AtomCount() + 1);

        this->atomToTriangleList.Resize(0);
        this->atomToTriangleList.SetCount(pmdc->AtomCount() + 1);

        this->atomToArea.Resize(0);
        this->atomToArea.SetCount(pmdc->AtomCount() + 1);

        auto vertAtomIDs = obj.GetVertexAttribPointerUInt32(0);
        auto triangleVertIDs = obj.GetTriIndexPointerUInt32();
        auto vertPos = obj.GetVertexPointerFloat();
        for (int i = 0; i < obj.GetVertexCount(); i++) {
            this->atomToVerticesList[vertAtomIDs[i]].AssertCapacity(15);
            this->atomToVerticesList[vertAtomIDs[i]].Append(i);
        }

        for (int i = 0; i < obj.GetTriCount(); i++) {
            auto atomID0 = vertAtomIDs[triangleVertIDs[i * 3 + 0]];
            auto atomID1 = vertAtomIDs[triangleVertIDs[i * 3 + 1]];
            auto atomID2 = vertAtomIDs[triangleVertIDs[i * 3 + 2]];
            // this->atomToTriangleList[vertAtomIDs[triangleVertIDs[i*3+0]]].AssertCapacity(15);
            if (this->atomToTriangleList[atomID0].Find(i) == NULL) this->atomToTriangleList[atomID0].Append(i);
            if (this->atomToTriangleList[atomID1].Find(i) == NULL) this->atomToTriangleList[atomID1].Append(i);
            if (this->atomToTriangleList[atomID2].Find(i) == NULL) this->atomToTriangleList[atomID2].Append(i);
        }

        old_pmdc_DataHash = pmdc->DataHash();

        // get area for each atom
        for (int i = 0; i < this->atomToTriangleList.Count(); i++) {
            float area = 0;
            int ID = i;
            for (int j = 0; j < this->atomToTriangleList[ID].Count(); j++) {
                int vert_1 = triangleVertIDs[this->atomToTriangleList[ID][j] * 3 + 0];
                int vert_2 = triangleVertIDs[this->atomToTriangleList[ID][j] * 3 + 1];
                int vert_3 = triangleVertIDs[this->atomToTriangleList[ID][j] * 3 + 2];
                glm::vec3 a = vec3(vertPos[vert_1 * 3 + 0], vertPos[vert_1 * 3 + 1], vertPos[vert_1 * 3 + 2]);
                glm::vec3 b = vec3(vertPos[vert_2 * 3 + 0], vertPos[vert_2 * 3 + 1], vertPos[vert_2 * 3 + 2]);
                glm::vec3 c = vec3(vertPos[vert_3 * 3 + 0], vertPos[vert_3 * 3 + 1], vertPos[vert_3 * 3 + 2]);
                glm::vec3 ab = b - a;
                glm::vec3 ac = c - a;
                area += glm::length(glm::cross(ab, ac)) * 0.5;
            }
            atomToArea[i] = area;
        }
    }

    this->vertexCount = obj.GetVertexCount();


    /***************************
     *** 1. NEIGHBOUR SEARCH ***
     ***************************/
    FFRNNS rns;
    std::vector<std::vector<uint>> tmpNeighbours;
    // collect 1. data
    visFloat4 pointCloud;
    pointCloud.SetCount(pmdc->AtomCount());

#pragma omp parallel for
    for (int i = 0; i < pmdc->AtomCount(); i++) {
        const megamol::protein_calls::MolecularDataCall::AtomType atom0 = pmdc->AtomTypes()[pmdc->AtomTypeIndices()[i]];
        pointCloud[i].SetX(pmdc->AtomPositions()[i * 3 + 0]);
        pointCloud[i].SetY(pmdc->AtomPositions()[i * 3 + 1]);
        pointCloud[i].SetZ(pmdc->AtomPositions()[i * 3 + 2]);
        pointCloud[i].SetW(atom0.Radius());
    }


    // create 1. search structure
    auto start = std::chrono::high_resolution_clock::now();


    vislib::Array<float> tmpClusterCentroid;
    tmpClusterCentroid.SetCount(4);
    float3 tmpClusterCentroid2;

    auto clusterCentroids = this->clusterData.clusterCentroids->PeekElements();
    auto clusterCentroidsCuboids = this->clusterData.clusterCentroidsCuboids->PeekElements();
    this->clusterAreaVertices.SetCount((this->clusterData.clusterCentroids->Count() / 4) + 1);
    this->clusterAreas.resize(this->clusterAreaVertices.Count());
    // allocate clusterAreaVertices[0] containing all protein atoms involved in at least on of all clsuters

    for (int i = 0; i <= this->clusterData.clusterCentroids->Count() / 4; i++) {
        this->clusterAreaVertices[i].SetCount(this->vertexCount);
        this->clusterAreas[i] = 0;
        for (int z = 0; z < this->vertexCount; z++) {
            this->clusterAreaVertices[i][z] = 0;
        }
    }


    if (this->proteinPocketAtoms.size() > 0) {
        for (int i = 0; i < this->proteinPocketAtoms.size(); i++) {
            this->proteinPocketAtoms[i].resize(0);
        }
        this->proteinPocketAtoms.resize(0);
    }

    // initialize
    this->clusterAreaAtomsToClusterIDs.resize(pmdc->AtomCount(), 0);


    if (true) {
        int a = lmc->getCurrentGlobalModelID();
        // slow exact pocket surface determination
        visFloat4 clusterAtomData;
        this->proteinPocketAtoms.resize(0);
        this->proteinPocketAtoms.resize(this->clusterData.clusterSizes->Count() + 1);

        auto clusterIndexSum = this->clusterData.clusterIndexSum->PeekElements();
        for (int i = 0; i < this->clusterData.clusterCentroids->Count() / 4; i++) {
            clusterAtomData.Resize(0);
            int clusterIDx_start = i > 0 ? clusterIndexSum[i] : clusterIndexSum[0];
            int clusterIDx_end = clusterIndexSum[i + 1];

            for (int j = clusterIDx_start; j < clusterIDx_end; j++) {
                lmc->SetCurrentLigandAndModel_byGlobalMdlID(
                    this->clusterData.assign_modelID_sphereID->PeekElements()[j].GetY());
                (*lmc)(LigandModelCall::CallForGetData);
                (*mdc)(MolecularDataCall::CallForGetData);

                for (int z = 0; z < mdc->AtomCount(); z++) {
                    const megamol::protein_calls::MolecularDataCall::AtomType atom0 =
                        mdc->AtomTypes()[mdc->AtomTypeIndices()[z]];
                    vislib::math::Vector<float, 4> tmp;
                    tmp.SetX(mdc->AtomPositions()[z * 3 + 0]);
                    tmp.SetY(mdc->AtomPositions()[z * 3 + 1]);
                    tmp.SetZ(mdc->AtomPositions()[z * 3 + 2]);
                    tmp.SetW(atom0.Radius());

                    clusterAtomData.Append(tmp);
                }
            }
            rns.setData(clusterAtomData, this->searchRadius_PocketPatch);
            tmpNeighbours = rns.getNeigbours(pointCloud);

            if (tmpNeighbours.size() > 0) {
                for (int j = 0; j < pmdc->AtomCount(); j++) {
                    if (tmpNeighbours[j].size() > 0) {
                        int neighbourAtom = j;


                        this->clusterAreaAtomsToClusterIDs[neighbourAtom] = i + 1;
                        this->clusterAreas[i + 1] += this->atomToArea[neighbourAtom];
                        this->clusterAreas[0] += this->atomToArea[neighbourAtom];
                        this->proteinPocketAtoms[i].push_back(j);
                        for (int p = 0; p < this->atomToVerticesList[neighbourAtom].Count(); p++) {
                            clusterAreaVertices[i + 1][this->atomToVerticesList[neighbourAtom][p]] = 1;

                            // index = 0 contains all atoms of all clusters
                            clusterAreaVertices[0][this->atomToVerticesList[neighbourAtom][p]] = 1;
                        }
                    }
                }
            }
        }
        // printf("totalArea: %f\n", this->clusterAreas[0]);

    } else {
        // fast approximated pocket surface determination
        for (int i = 0; i < this->clusterData.clusterCentroids->Count() / 4; i++) {

            tmpClusterCentroid[0] = clusterCentroids[i * 4 + 0];
            tmpClusterCentroid[1] = clusterCentroids[i * 4 + 1];
            tmpClusterCentroid[2] = clusterCentroids[i * 4 + 2];
            tmpClusterCentroid[3] = 0.0f;

            tmpClusterCentroid2 = {0.0f, 0.0f, 0.0f};
            float3 tmpPoint = {clusterCentroidsCuboids[i].Width(), clusterCentroidsCuboids[i].Height(),
                clusterCentroidsCuboids[i].Depth()};
            float searchDistance = distanceTwoVectors(tmpPoint, tmpClusterCentroid2) / 6.0f; // +5.0f;
            rns.setData(pointCloud, searchDistance);
            tmpNeighbours = rns.getNeigbours(tmpClusterCentroid);
#pragma omp parallel for
            for (int j = 0; j < tmpNeighbours.size(); j++) {
                for (int k = 0; k < tmpNeighbours[j].size(); k++) {
                    int neighbourAtom = tmpNeighbours[j][k];
                    for (int p = 0; p < this->atomToVerticesList[neighbourAtom].Count(); p++) {
                        clusterAreaVertices[i + 1][this->atomToVerticesList[neighbourAtom][p]] = 1;
                        // index = 0 contains all atoms of all clusters
                        clusterAreaVertices[0][this->atomToVerticesList[neighbourAtom][p]] = 1;
                    }
                }
            }
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time 1. neighbour search - getInteractionAreaTriangles: " << elapsed.count() << " s\n";

    // mark the vertices of the protein surface regrding the forces
    markForceVertices_getBarchatDat(this->atomToVerticesList, this->vertexCount, pmdc->AtomCount());


    auto finish_c = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_c = finish_c - start_c;
    std::cout << "Elapsed time 1. complete getInteractionAreaTriangles: " << elapsed_c.count() << " s\n";

    return 1;
}


bool ModelClusterRenderer::markForceVertices_getBarchatDat(
    vislib::Array<vislib::Array<uint>>& atomToVerticesList, int vertexCount, int protAtomCnt) {

    LigandModelCall* lmc = this->ligandModel_callerSlot.CallAs<LigandModelCall>();
    if (lmc == NULL) {
        return false;
    }

    std::vector<uint*> pointerStore;
    pointerStore.push_back(
        this->forcesLineStorage.getCurForceCntPerProtAtomArray(interForceTypes.saltBridges, protAtomCnt + 1, lmc));
    pointerStore.push_back(
        this->forcesLineStorage.getCurForceCntPerProtAtomArray(interForceTypes.halogenBonds, protAtomCnt + 1, lmc));
    pointerStore.push_back(this->forcesLineStorage.getCurForceCntPerProtAtomArray(
        interForceTypes.hydrophobicInteractions, protAtomCnt + 1, lmc));
    pointerStore.push_back(
        this->forcesLineStorage.getCurForceCntPerProtAtomArray(interForceTypes.metalComplexes, protAtomCnt + 1, lmc));
    pointerStore.push_back(this->forcesLineStorage.getCurForceCntPerProtAtomArray(
        interForceTypes.piCationInteractions, protAtomCnt + 1, lmc));
    pointerStore.push_back(
        this->forcesLineStorage.getCurForceCntPerProtAtomArray(interForceTypes.piStacks, protAtomCnt + 1, lmc));
    pointerStore.push_back(
        this->forcesLineStorage.getCurForceCntPerProtAtomArray(interForceTypes.HBonds, protAtomCnt + 1, lmc));


    MolecularDataCall* pmdc = this->proteinMolecularData_callerSlot.CallAs<MolecularDataCall>();
    if (pmdc == NULL) {
        return 0;
    }

    // initialize tmpPerClustCounter
    std::vector<std::map<string, std::vector<uint>>> tmpPerClustCounter;

    this->clusterForceAreas.clear();
    this->clusterForceResidueCnts.clear();
    this->clusterForceAreas.resize(this->clusterData.clusterCentroids->Count() / 4);
    this->clusterForceResidueCnts.resize(this->clusterData.clusterCentroids->Count() / 4);

    tmpPerClustCounter.resize(this->clusterData.clusterCentroids->Count() / 4);
    for (int i = 0; i < tmpPerClustCounter.size(); i++) {
        for (int j = 0; j < this->interForceTypes.forceTypeVec.size(); j++) {
            tmpPerClustCounter[i][this->interForceTypes.forceTypeVec[j]].resize(pmdc->AtomCount() + 1, 0);
        }
    }


    std::string curFType = "";
    // mark protein acceptors involved in H-bonds
    for (int z = 0; z < pointerStore.size(); z++) {
        auto marked = &this->markedVerticesIndicesMap[this->interForceTypes.saltBridges];
        auto maxMarked = &this->maxForceCntPerProtAtomMap[this->interForceTypes.saltBridges];
        curFType = this->interForceTypes.saltBridges;
        if (z == 1) {
            marked = &this->markedVerticesIndicesMap[this->interForceTypes.halogenBonds];
            maxMarked = &this->maxForceCntPerProtAtomMap[this->interForceTypes.halogenBonds];
            curFType = this->interForceTypes.halogenBonds;
        } else if (z == 2) {
            marked = &this->markedVerticesIndicesMap[this->interForceTypes.hydrophobicInteractions];
            maxMarked = &this->maxForceCntPerProtAtomMap[this->interForceTypes.hydrophobicInteractions];
            curFType = this->interForceTypes.hydrophobicInteractions;
        } else if (z == 3) {
            marked = &this->markedVerticesIndicesMap[this->interForceTypes.metalComplexes];
            maxMarked = &this->maxForceCntPerProtAtomMap[this->interForceTypes.metalComplexes];
            curFType = this->interForceTypes.metalComplexes;
        } else if (z == 4) {
            marked = &this->markedVerticesIndicesMap[this->interForceTypes.piCationInteractions];
            maxMarked = &this->maxForceCntPerProtAtomMap[this->interForceTypes.piCationInteractions];
            curFType = this->interForceTypes.piCationInteractions;
        } else if (z == 5) {
            marked = &this->markedVerticesIndicesMap[this->interForceTypes.piStacks];
            maxMarked = &this->maxForceCntPerProtAtomMap[this->interForceTypes.piStacks];
            curFType = this->interForceTypes.piStacks;
        } else if (z == 6) {
            marked = &this->markedVerticesIndicesMap[this->interForceTypes.HBonds];
            maxMarked = &this->maxForceCntPerProtAtomMap[this->interForceTypes.HBonds];
            curFType = this->interForceTypes.HBonds;
        }

        marked[0].resize(vertexCount);
#pragma omp parallel for
        for (int i = 0; i < vertexCount; i++) {
            marked[0][i] = 0;
        }
        maxMarked[0] = 0;
        if (pointerStore[z] == NULL) {
            continue;
        }

        // init
        for (int d = 0; d < lmc->getClusterAndCentroidData().clusterCentroids->Count() / 4; d++) {
            this->clusterForceResidueCnts[d][curFType];
            this->clusterForceAreas[d][curFType] = 0;
        }

        //#pragma omp parallel for
        for (int i = 0; i < protAtomCnt; i++) {
            // only mark vertices if they are part of the determined cluster areas;
            if ((pointerStore[z][i + 1] != 0) && this->clusterAreaAtomsToClusterIDs[i] > 0) {

                // get the surface areas for the forces
                int clusterID = this->clusterAreaAtomsToClusterIDs[i] - 1;
                float area = this->atomToArea[i];
                clusterForceAreas[clusterID][curFType] += area;

                if (tmpPerClustCounter[clusterID][curFType][i] == 0) {
                    tmpPerClustCounter[clusterID][curFType][i] = 1;

                    float area = this->atomToArea[i];
                    clusterForceAreas[clusterID][curFType] += area;
                }

                // get the residues for the forces
                int residueID = pmdc->AtomResidueIndices()[i];
                this->clusterForceResidueCnts[clusterID][curFType][residueID] = 1;

                for (int j = 0; j < atomToVerticesList[i].Count(); j++) {
                    marked[0][atomToVerticesList[i][j]] = pointerStore[z][i + 1];
                }
                maxMarked[0] = maxMarked[0] < pointerStore[z][i + 1] ? pointerStore[z][i + 1] : maxMarked[0];
            }
        }
    }

    // this->max_hBondAtomCnt = lmc->getInteractionForces().H_bonds.max_hBondAtomCnt_protein;
    printf("marked force involved vertices");


    return true;
}


/**********************************************************************
********* PMDC MolecularDataCall for Protein Residues Model(s) ********
***********************************************************************/
bool ModelClusterRenderer::getExtentPMDC(core::Call& call) {
    MolecularDataCall* protMdcOut = dynamic_cast<MolecularDataCall*>(&call);
    if (protMdcOut == NULL) return false;

    MolecularDataCall* protMdcIn = this->proteinMolecularData_callerSlot.CallAs<MolecularDataCall>();
    if (protMdcIn == NULL) return false;
    (*protMdcIn)(MolecularDataCall::CallForGetExtent);

    protMdcOut->AccessBoundingBoxes().Clear();
    protMdcOut->AccessBoundingBoxes().SetObjectSpaceBBox(protMdcIn->AccessBoundingBoxes().ObjectSpaceBBox());
    protMdcOut->AccessBoundingBoxes().SetObjectSpaceClipBox(protMdcIn->AccessBoundingBoxes().ObjectSpaceClipBox());

    protMdcOut->SetFrameCount(1);
    if (this->pmdcOut_newData) {
        int oldHash = protMdcOut->DataHash();
        protMdcOut->SetDataHash(protMdcOut->DataHash() + 1);
    }

    return true;
}

// call copy the ingoing MolecularDataCall
bool ModelClusterRenderer::getDataPMDC(core::Call& call) {

    if (this->pmdcOut_newData && this->renderDataControl.interactionResidues) {
        if (this->atomPostions.Count() > 0) {
            this->atomPostions.Clear();
            this->atomTypeIdx.Clear();
            this->atomCharges.Clear();
            this->connections.Clear();
            this->atomRadii.Clear();
            setInteractionResiduesPMDC(call);
            this->pmdcOut_newData = false;
            return true;
        } else {
            setInteractionResiduesPMDC(call);
            return true;
        }
    } else if (this->renderDataControl.interactionResidues) {
        return true;
        // set MDC to empty (case renderDataControl.interactionResidues == false &&  renderDataControl.residues
        // == false)
    } else if (!this->renderDataControl.residues) {
        setInteractionResiduesPMDC(call);

        /***************************************************
         * if whole protein --> just pass through MDC data *
         ***************************************************/
    } else if (this->pmdcOut_newData) {


        using megamol::core::utility::log::Log;

        MolecularDataCall* protMdcOut = dynamic_cast<MolecularDataCall*>(&call);
        if (protMdcOut == NULL) return false;

        MolecularDataCall* protMdcIn = this->proteinMolecularData_callerSlot.CallAs<MolecularDataCall>();
        if (protMdcIn == NULL) return false;
        (*protMdcIn)(MolecularDataCall::CallForGetData);

        protMdcOut->SetDataHash(protMdcIn->DataHash());
        protMdcOut->SetFrameCount(1);

        protMdcOut->SetAtoms(protMdcIn->AtomCount(), protMdcIn->AtomTypeCount(), protMdcIn->AtomTypeIndices(),
            protMdcIn->AtomPositions(), protMdcIn->AtomTypes(), 0, 0, protMdcIn->AtomCharges(), 0);

        protMdcOut->SetBFactorRange(0.0f, 0.0f);
        protMdcOut->SetChargeRange(protMdcIn->MinimumCharge(), protMdcIn->MinimumCharge());
        protMdcOut->SetOccupancyRange(0.0f, 0.0f);
        protMdcOut->SetFormerAtomIndices(0);

        protMdcOut->SetConnections(protMdcIn->ConnectionCount(), protMdcIn->Connection());
        protMdcOut->SetResidues(protMdcIn->ResidueCount(), (const MolecularDataCall::Residue**)protMdcIn->Residues());
        protMdcOut->SetSolventResidueIndices(0, 0);
        protMdcOut->SetResidueTypeNames(
            protMdcIn->ResidueTypeNameCount(), (vislib::StringA*)protMdcIn->ResidueTypeNames());
        protMdcOut->SetMolecules(protMdcIn->MoleculeCount(), (MolecularDataCall::Molecule*)protMdcIn->Molecules());
        protMdcOut->SetChains(protMdcIn->ChainCount(), (MolecularDataCall::Chain*)protMdcIn->Chains());

        // Set the filter array for the molecular data call
        protMdcOut->SetFilter(0);
        this->pmdcOut_newData = false;
    }
    return true;
}

bool ModelClusterRenderer::setInteractionResiduesPMDC(core::Call& call) {

    int residueCnt = 0;

    residueLables.Resize(0);
    residueLablesPositions.Resize(0);

    MolecularDataCall* protMdcIn = this->proteinMolecularData_callerSlot.CallAs<MolecularDataCall>();
    if (protMdcIn == NULL) return false;
    (*protMdcIn)(MolecularDataCall::CallForGetData);

    MolecularDataCall* protMdcOut = dynamic_cast<MolecularDataCall*>(&call);
    if (protMdcOut == NULL) return false;

    LigandModelCall* lmc = this->ligandModel_callerSlot.CallAs<LigandModelCall>();
    MolecularDataCall* mdc = this->ligandMolecularData_callerSlot.CallAs<MolecularDataCall>();
    if (mdc == NULL) {
        return false;
    }
    if (lmc == NULL) {
        return false;
    }
    (*mdc)(MolecularDataCall::CallForGetData);
    (*lmc)(MolecularDataCall::CallForGetExtent);
    (*lmc)(MolecularDataCall::CallForGetData);
    // protMdcOut->SetDataHash(protMdcIn->DataHash());
    // protMdcOut->SetFrameCount(1);

    vislib::Array<int> involvedResiduesIndices;
    involvedResiduesIndices.SetCount(protMdcIn->ResidueCount());

    // initialize with zero
    for (int i = 0; i < involvedResiduesIndices.Count(); i++) {
        involvedResiduesIndices[i] = 0;
    }

    // sets MDC to empty if renderDataControl.residues == false
    if (this->renderDataControl.residues) {
        // show residues of selected cluster
        if (this->s.curClusterID > -1) {
            if (this->renVarsFGS.selectedSphereID != -1 || false) {
                getInteactionResidues_FuncGroups(this->fgsPerClusterFiltered[this->s.curClusterID],
                    involvedResiduesIndices, protMdcIn, this->renVarsFGS.selectedSphereID);
            } else {
                if (this->proteinPocketAtoms.size() == 0) {
                    return false;
                }
                for (int i = 0; i < (int)(this->proteinPocketAtoms[this->s.curClusterID].size() - 1); i++) {

                    int protID = this->proteinPocketAtoms[this->s.curClusterID][i];
                    int resIDx = protMdcIn->AtomResidueIndices()[protID];
                    involvedResiduesIndices[resIDx]++;
                }
            }
            // show residues of all clusters
        } else {
            // check if clusterData is already present (DataHash != 0)
            if (lmc->getClusterAndCentroidData().clusterDataHash != 0) {
                if (lmc->getClusterAndCentroidData().clusterSizes->Count() + 1 == this->proteinPocketAtoms.size()) {
                    for (int i = 0; i < lmc->getClusterAndCentroidData().clusterCentroids->Count() / 4; i++) {
                        for (int j = 0; j < this->proteinPocketAtoms[i].size(); j++) {
                            int protID = this->proteinPocketAtoms[i][j];
                            int resIDx = protMdcIn->AtomResidueIndices()[protID];
                            involvedResiduesIndices[resIDx]++;
                        }
                    }
                }
            }
        }
    }

    int old_resIDx = -1;
    std::map<int, vislib::math::Cuboid<float>> residueAtomsCuboid;
    std::map<int, int> residueAtomsCnt;
    int cuboidCnt = 0;
    for (int i = 1; i < protMdcIn->AtomCount(); i++) {
        int resIDx = protMdcIn->AtomResidueIndices()[i];
        if (involvedResiduesIndices[resIDx] > 0) {
            int test = protMdcIn->Residues()[resIDx]->Identifier();
            MolecularDataCall::AminoAcid* amino;

            // discard atoms which are part of an anminoacid
            if (test == MolecularDataCall::Residue::AMINOACID) {
                amino = (MolecularDataCall::AminoAcid*)(protMdcIn->Residues()[resIDx]);
            } else {
                continue;
            }


            // discard the backbone
            if (!this->renderDataControl.residuesBackbone) {
                if (amino->CAlphaIndex() == i || amino->CCarbIndex() == i || amino->NIndex() == i ||
                    amino->OIndex() == i) {
                    continue;
                }
            }

            // discard hydrogens
            if (protMdcIn->AtomTypes()[protMdcIn->AtomTypeIndices()[i]].Name().StartsWith("H") ||
                protMdcIn->AtomTypes()[protMdcIn->AtomTypeIndices()[i]].Name().StartsWith("1") ||
                protMdcIn->AtomTypes()[protMdcIn->AtomTypeIndices()[i]].Name().StartsWith("2") ||
                protMdcIn->AtomTypes()[protMdcIn->AtomTypeIndices()[i]].Name().StartsWith("3")) {
                continue;
            }


            if (residueAtomsCuboid[resIDx].IsEmpty()) {
                residueAtomsCuboid[resIDx].Set(protMdcIn->AtomPositions()[i * 3 + 0],
                    protMdcIn->AtomPositions()[i * 3 + 1], protMdcIn->AtomPositions()[i * 3 + 2],
                    protMdcIn->AtomPositions()[i * 3 + 0], protMdcIn->AtomPositions()[i * 3 + 1],
                    protMdcIn->AtomPositions()[i * 3 + 2] + 0.01f);
            } else {
                residueAtomsCuboid[resIDx].GrowToPoint(protMdcIn->AtomPositions()[i * 3 + 0],
                    protMdcIn->AtomPositions()[i * 3 + 1], protMdcIn->AtomPositions()[i * 3 + 2]);
            }


            this->atomPostions.Append(protMdcIn->AtomPositions()[i * 3 + 0]);
            this->atomPostions.Append(protMdcIn->AtomPositions()[i * 3 + 1]);
            this->atomPostions.Append(protMdcIn->AtomPositions()[i * 3 + 2]);
            this->atomTypeIdx.Append(protMdcIn->AtomTypeIndices()[i]);
            this->atomCharges.Append(protMdcIn->AtomCharges()[i]);
            this->atomRadii.Append(protMdcIn->AtomTypes()[protMdcIn->AtomTypeIndices()[i]].Radius());

            if (old_resIDx != resIDx) {
                residueCnt++;
                old_resIDx = resIDx;
                std::string chainName;
                vislib::StringA type, resIDx_str, lable;

                int firstAtomID = protMdcIn->Residues()[resIDx]->FirstAtomIndex();
                int indetifier = protMdcIn->Residues()[resIDx]->Type();
                int moleculeID = protMdcIn->Residues()[resIDx]->MoleculeIndex();

                chainName = protMdcIn->Chains()[protMdcIn->Molecules()[moleculeID].ChainIndex()].Name();
                type = protMdcIn->ResidueTypeNames()[indetifier];
                resIDx_str.Format("%d", resIDx + 1);
                lable = type + ":" + resIDx_str + ":" + chainName.c_str();
                this->residueLables.Append(lable);
            }
        }
    }

    // set the text label locations
    map<int, vislib::math::Cuboid<float>>::iterator it;
    for (it = residueAtomsCuboid.begin(); it != residueAtomsCuboid.end(); it++) {
        auto cPoint = it->second.CalcCenter();
        this->residueLablesPositions.Append({cPoint.GetX(), cPoint.GetY(), cPoint.GetZ()});
    }

    // TODO: fix one out of range message when switching back from all residues to only involved residues!
    vislib::math::Vector<float, 3> atomPos0, atomPos1;
    for (int i = 0; i < this->atomTypeIdx.Count(); i++) {

        // get first atom positions
        atomPos0.Set(this->atomPostions[3 * i + 0], this->atomPostions[3 * i + 1], this->atomPostions[3 * i + 2]);
        for (int j = i + 1; j < this->atomTypeIdx.Count(); j++) {
            // get second atom positions
            atomPos1.Set(this->atomPostions[3 * j + 0], this->atomPostions[3 * j + 1], this->atomPostions[3 * j + 2]);

            // check distance
            if ((atomPos0 - atomPos1).Length() < 0.58f * (this->atomRadii[i] + this->atomRadii[j])) {
                // add connection
                this->connections.Append(i);
                this->connections.Append(j);
            }
        }
    }

    // protMdcOut->SetDataHash(protMdcIn->DataHash() + 1);

    protMdcOut->SetAtoms(this->atomPostions.Count() / 3, protMdcIn->AtomTypeCount(), this->atomTypeIdx.PeekElements(),
        this->atomPostions.PeekElements(), protMdcIn->AtomTypes(), 0, 0, this->atomCharges.PeekElements(), 0);

    protMdcOut->SetConnections(this->connections.Count() / 2, this->connections.PeekElements());


    return true;
}


/************************************************************
********* LMDC MolecaulrDataCall for Ligand Model(s) ********
*************************************************************/

bool ModelClusterRenderer::getExtentLMDC(core::Call& call) {

    if (this->firstRenderCall) return false;

    MolecularDataCall* ligandMdcOut = dynamic_cast<MolecularDataCall*>(&call);
    if (ligandMdcOut == NULL) return false;

    MolecularDataCall* ligMdcIn = this->ligandMolecularData_callerSlot.CallAs<MolecularDataCall>();
    if (ligMdcIn == NULL) return false;
    (*ligMdcIn)(MolecularDataCall::CallForGetExtent);

    ligandMdcOut->AccessBoundingBoxes().Clear();
    ligandMdcOut->AccessBoundingBoxes().SetObjectSpaceBBox(ligMdcIn->AccessBoundingBoxes().ObjectSpaceBBox());
    ligandMdcOut->AccessBoundingBoxes().SetObjectSpaceClipBox(ligMdcIn->AccessBoundingBoxes().ObjectSpaceClipBox());

    ligandMdcOut->SetFrameCount(1);
    if (this->lmdcOut_newData) {
        int oldHash = ligandMdcOut->DataHash();
        ligandMdcOut->SetDataHash(ligandMdcOut->DataHash() + 1);
    }

    return true;
}

// call copy the ingoing MolecularDataCall

bool ModelClusterRenderer::getDataLMDC(core::Call& call) {

    if (this->firstRenderCall) return false;

    LigandModelCall* lmc = this->ligandModel_callerSlot.CallAs<LigandModelCall>();
    if (lmc == NULL) return false;
    (*lmc)(MolecularDataCall::CallForGetExtent);
    (*lmc)(MolecularDataCall::CallForGetData);
    if (lmc->getClusterAndCentroidData().DBSCANParams.paramsChanged) return false;


    if (this->lmdcOut_newData && this->modelRenderModes.curMode != this->modelRenderModes.single) {
        if (this->lmdc_positions.Count() > 0) {
            this->lmdc_positions.Clear();
            this->lmdc_charges.Clear();
            this->lmdc_connections.Clear();
            this->lmdc_atomTypeIdx.Clear();
            this->lmdc_chains.Clear();
            this->lmdc_residues.Clear();
            this->lmdc_moleclues.Clear();
            this->lmdc_residuesStore.Clear();
            this->lmdc_residuesIndices.Clear();
            setMultipleLigandModelsLMDC(call);
            this->lmdcOut_newData = false;
            return true;
        } else {
            setMultipleLigandModelsLMDC(call);
            this->lmdcOut_newData = false;
            return true;
        }
    } else if (this->modelRenderModes.curMode != this->modelRenderModes.single) {
        return true;

        /*****************************************************
         * if one lingad pose --> just pass through MDC data *
         *****************************************************/
    } else if (this->lmdcOut_newData && this->modelRenderModes.curMode == this->modelRenderModes.single) {
        this->lmdc_positions.Clear();
        using megamol::core::utility::log::Log;

        MolecularDataCall* ligandMdcOut = dynamic_cast<MolecularDataCall*>(&call);
        if (ligandMdcOut == NULL) return false;

        MolecularDataCall* ligMdcIn = this->ligandMolecularData_callerSlot.CallAs<MolecularDataCall>();
        if (ligMdcIn == NULL) return false;
        (*ligMdcIn)(MolecularDataCall::CallForGetData);

        ligandMdcOut->SetDataHash(ligMdcIn->DataHash());
        ligandMdcOut->SetFrameCount(1);

        auto a = ligMdcIn->AtomCount();
        auto b = ligMdcIn->AtomTypeCount();
        auto c = ligMdcIn->AtomTypeIndices();
        auto d = ligMdcIn->AtomPositions();
        auto e = ligMdcIn->AtomTypes();
        auto f = ligMdcIn->AtomCharges();

        ligandMdcOut->SetAtoms(ligMdcIn->AtomCount(), ligMdcIn->AtomTypeCount(), ligMdcIn->AtomTypeIndices(),
            ligMdcIn->AtomPositions(), ligMdcIn->AtomTypes(), 0, 0, ligMdcIn->AtomCharges(), 0);

        ligandMdcOut->SetBFactorRange(0.0f, 0.0f);
        ligandMdcOut->SetChargeRange(ligMdcIn->MinimumCharge(), ligMdcIn->MinimumCharge());
        ligandMdcOut->SetOccupancyRange(0.0f, 0.0f);
        ligandMdcOut->SetFormerAtomIndices(0);

        ligandMdcOut->SetConnections(ligMdcIn->ConnectionCount(), ligMdcIn->Connection());
        ligandMdcOut->SetResidues(ligMdcIn->ResidueCount(), (const MolecularDataCall::Residue**)ligMdcIn->Residues());
        ligandMdcOut->SetSolventResidueIndices(0, 0);
        ligandMdcOut->SetResidueTypeNames(
            ligMdcIn->ResidueTypeNameCount(), (vislib::StringA*)ligMdcIn->ResidueTypeNames());
        ligandMdcOut->SetMolecules(ligMdcIn->MoleculeCount(), (MolecularDataCall::Molecule*)ligMdcIn->Molecules());
        ligandMdcOut->SetChains(ligMdcIn->ChainCount(), (MolecularDataCall::Chain*)ligMdcIn->Chains());

        // Set the filter array for the molecular data call
        ligandMdcOut->SetFilter(0);
    }
    return true;
}


bool ModelClusterRenderer::loadTexture(std::string filename) {
    using megamol::core::utility::log::Log;

    static sg::graphics::PngBitmapCodec pbc;
    pbc.Image() = &this->img;
    ::glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    void* buf = NULL;
    SIZE_T size = 0;
    vislib::StringA name = filename.c_str();
    if ((size = megamol::core::utility::ResourceWrapper::LoadResource(
             this->GetCoreInstance()->Configuration(), name, &buf)) > 0) {
        if (pbc.Load(buf, size)) {
            this->img.Convert(vislib::graphics::BitmapImage::TemplateByteRGBA);
            this->modelModeIcons.Add(vislib::SmartPtr<vislib::graphics::gl::OpenGLTexture2D>());
            this->modelModeIcons.Last() = new vislib::graphics::gl::OpenGLTexture2D();
            if (this->modelModeIcons.Last()->Create(this->img.Width(), this->img.Height(), false,
                    this->img.PeekDataAs<BYTE>(), GL_RGBA) != GL_NO_ERROR) {
                Log::DefaultLog.WriteError("could not load %s texture.", filename.c_str());
                ARY_SAFE_DELETE(buf);
                return false;
            }
            this->modelModeIcons.Last()->Bind();
            this->modelModeIcons.Last()->SetFilter(GL_LINEAR, GL_LINEAR);
            glGenerateMipmap(GL_TEXTURE_2D);
            glBindTexture(GL_TEXTURE_2D, 0);
            ARY_SAFE_DELETE(buf);
            return true;
        } else {
            Log::DefaultLog.WriteError("could not read %s texture.", filename.c_str());
        }
    } else {
        Log::DefaultLog.WriteError("could not find %s texture.", filename.c_str());
    }
    return false;
}


bool ModelClusterRenderer::setMultipleLigandModelsLMDC(core::Call& call) {

    MolecularDataCall* ligandMdcOut = dynamic_cast<MolecularDataCall*>(&call);
    if (ligandMdcOut == NULL) return false;

    MolecularDataCall* ligMdcIn = this->ligandMolecularData_callerSlot.CallAs<MolecularDataCall>();
    if (ligMdcIn == NULL) return false;
    (*ligMdcIn)(MolecularDataCall::CallForGetExtent);
    (*ligMdcIn)(MolecularDataCall::CallForGetData);

    LigandModelCall* lmc = this->ligandModel_callerSlot.CallAs<LigandModelCall>();
    if (lmc == NULL) return false;
    (*lmc)(MolecularDataCall::CallForGetExtent);
    (*lmc)(MolecularDataCall::CallForGetData);


    int curLigID = lmc->getCurrentLigandID();
    vislib::Array<uint> models;
    models.Resize(0);
    if (this->modelRenderModes.curMode == this->modelRenderModes.allOnProtein) {
        for (int i = 0; i < *lmc->getModelCount(); i++) {
            models.Append(i);
        }
    } else if (this->modelRenderModes.curMode == this->modelRenderModes.allInCurPocket) {
        auto clustDat = lmc->getClusterAndCentroidData();
        // get start and end index for iterating over ligands of the cluster (clusterID)
        int clusterIDx_start = s.curClusterID > 0 ? clustDat.clusterIndexSum->PeekElements()[s.curClusterID]
                                                  : clustDat.clusterIndexSum->PeekElements()[0];
        int clusterIDx_end = clustDat.clusterIndexSum->PeekElements()[s.curClusterID + 1];

        // collect func groups of the cluster/pocket
        // iterating over ligands
        for (int j = clusterIDx_start; j < clusterIDx_end; j++) {

            // set current ligand
            lmc->SetCurrentLigandAndModel_byGlobalMdlID(clustDat.assign_modelID_sphereID->PeekElements()[j].GetY());
            // get data of current ligand
            (*lmc)(LigandModelCall::CallForGetData);
            if (lmc->getCurrentLigandID() == curLigID) {
                models.Append(lmc->getCurrentModelID());
            }
        }
    }
    uint mdlCnt = models.Count();


    int lmdc_atomCnt = ligMdcIn->AtomCount() * mdlCnt;
    int lmdc_connectCnt = ligMdcIn->ConnectionCount() * mdlCnt;

    this->lmdc_positions.Resize(lmdc_atomCnt * 3);
    this->lmdc_positions.SetCount(lmdc_atomCnt * 3);

    this->lmdc_charges.Resize(lmdc_atomCnt);
    this->lmdc_charges.SetCount(lmdc_atomCnt);

    this->lmdc_atomTypeIdx.Resize(lmdc_atomCnt);
    this->lmdc_atomTypeIdx.SetCount(lmdc_atomCnt);

    this->lmdc_residuesIndices.Resize(lmdc_atomCnt);
    this->lmdc_residuesIndices.SetCount(lmdc_atomCnt);

    this->lmdc_connections.Resize(lmdc_connectCnt * 2);
    this->lmdc_connections.SetCount(lmdc_connectCnt * 2);


    // SET THE DATA
    for (int i = 0; i < models.Count(); i++) {

        // call the ligand data
        lmc->SetCurrentLigandAndModel(curLigID, models[i]);
        (*lmc)(MolecularDataCall::CallForGetData);
        (*ligMdcIn)(MolecularDataCall::CallForGetData);
        ligandMdcOut->SetDataHash(ligMdcIn->DataHash());
        ligandMdcOut->SetFrameCount(1);

        // set per atom data
#pragma omp parallel for
        for (int j = 0; j < ligMdcIn->AtomCount(); j++) {
            int idx = i * ligMdcIn->AtomCount() + j;
            this->lmdc_positions[idx * 3 + 0] = ligMdcIn->AtomPositions()[j * 3 + 0];
            this->lmdc_positions[idx * 3 + 1] = ligMdcIn->AtomPositions()[j * 3 + 1];
            this->lmdc_positions[idx * 3 + 2] = ligMdcIn->AtomPositions()[j * 3 + 2];
            this->lmdc_charges[idx] = ligMdcIn->AtomCharges()[j];
            this->lmdc_atomTypeIdx[idx] = ligMdcIn->AtomTypeIndices()[j];
            this->lmdc_residuesIndices[idx] = i;
        }

        // set connections between atoms (x->y)
#pragma omp parallel for
        for (int j = 0; j < ligMdcIn->ConnectionCount(); j++) {
            int idx = (i * ligMdcIn->ConnectionCount() + j) * 2;
            int add = ligMdcIn->AtomCount() * i;
            this->lmdc_connections[idx + 0] = ligMdcIn->Connection()[j * 2 + 0] + add;
            this->lmdc_connections[idx + 1] = ligMdcIn->Connection()[j * 2 + 1] + add;
        }

        // set atom->residue information
        MolecularDataCall::Residue residue = *ligMdcIn->Residues()[0];
        residue.SetPosition(ligMdcIn->AtomCount() * i, ligMdcIn->AtomCount() * (i + 1));
        residue.SetMoleculeIndex(i);
        this->lmdc_residuesStore.Append(residue);
        this->lmdc_residues.Append(&this->lmdc_residuesStore.Last());

        // set residue->molecule information
        MolecularDataCall::Molecule mol = ligMdcIn->Molecules()[0];
        mol.SetPosition(i, 1);
        this->lmdc_moleclues.Append(mol);

        // set molecule->chain information
        MolecularDataCall::Chain chain = ligMdcIn->Chains()[0];
        chain.SetPosition(i, 1);
        this->lmdc_chains.Append(chain);
    }


    ligandMdcOut->SetAtoms(this->lmdc_positions.Count() / 3, ligMdcIn->AtomTypeCount(),
        this->lmdc_atomTypeIdx.PeekElements(), this->lmdc_positions.PeekElements(), ligMdcIn->AtomTypes(),
        this->lmdc_residuesIndices.PeekElements(), 0, this->lmdc_charges.PeekElements(), 0);
    ligandMdcOut->SetBFactorRange(0.0f, 0.0f);
    ligandMdcOut->SetChargeRange(ligMdcIn->MinimumCharge(), ligMdcIn->MinimumCharge());
    ligandMdcOut->SetOccupancyRange(0.0f, 0.0f);
    ligandMdcOut->SetFormerAtomIndices(0);

    ligandMdcOut->SetResidues(
        this->lmdc_residues.Count(), (const MolecularDataCall::Residue**)this->lmdc_residues.PeekElements());
    ligandMdcOut->SetSolventResidueIndices(0, 0);
    ligandMdcOut->SetResidueTypeNames(ligMdcIn->ResidueTypeNameCount(), (vislib::StringA*)ligMdcIn->ResidueTypeNames());
    ligandMdcOut->SetMolecules(mdlCnt, this->lmdc_moleclues.PeekElements());
    ligandMdcOut->SetChains(mdlCnt, this->lmdc_chains.PeekElements());

    // Set the filter array for the molecular data call
    ligandMdcOut->SetFilter(0);
    ligandMdcOut->SetConnections(
        static_cast<unsigned int>(this->lmdc_connections.Count() / 2), this->lmdc_connections.PeekElements());


    return true;
}


bool ModelClusterRenderer::getInteactionResidues_FuncGroups(LigandModelCall::FGS_structN& functGH_struct,
    vislib::Array<int>& involvedResiduesIndices, protein_calls::MolecularDataCall* protMdcIn, int funcGrpID) {
    visFloat4 protAtoms;
    protAtoms.SetCount(protMdcIn->AtomCount() * 4);

    vislib::Array<float> centroid;
    centroid.SetCount(1 * 4);

    // check if func group data is present
    if (functGH_struct.clusterCnt > 0) {
        centroid[0] = functGH_struct.clusterCentroids[funcGrpID * 4 + 0];
        centroid[1] = functGH_struct.clusterCentroids[funcGrpID * 4 + 1];
        centroid[2] = functGH_struct.clusterCentroids[funcGrpID * 4 + 2];
        centroid[3] = 0.0f;
    }
    for (int i = 0; i < protMdcIn->AtomCount(); i++) {
        int typeIDx = protMdcIn->AtomTypeIndices()[i];
        auto atomType = protMdcIn->AtomTypes()[typeIDx];
        protAtoms[i] = {protMdcIn->AtomPositions()[3 * i + 0], protMdcIn->AtomPositions()[3 * i + 1],
            protMdcIn->AtomPositions()[3 * i + 2], atomType.Radius()};
    }

    std::vector<std::vector<uint>> neighbourProtAtoms;
    FFRNNS fns;
    fns.setData(protAtoms, this->searchRadius_InteractionResidues);
    neighbourProtAtoms = fns.getNeigbours(centroid);

    for (int i = 0; i < neighbourProtAtoms.size(); i++) {
        for (int j = 0; j < neighbourProtAtoms[i].size(); j++) {
            int resIDx = protMdcIn->AtomResidueIndices()[neighbourProtAtoms[i][j]];
            involvedResiduesIndices[resIDx]++;
        }
    }

    return true;
};


/*
 * ModelClusterRenderer::getClipData
 */
void ModelClusterRenderer::getClipData(vislib::math::Vector<double, 4>* clipDat, float* clipCol, float* clipPoint) {
    view::CallClipPlane* ccp = this->getClipPlane_callerSlot.CallAs<view::CallClipPlane>();
    if ((ccp != NULL) && (*ccp)()) {
        clipDat->SetX(ccp->GetPlane().Normal().X());
        clipDat->SetY(ccp->GetPlane().Normal().Y());
        clipDat->SetZ(ccp->GetPlane().Normal().Z());
        vislib::math::Vector<float, 3> grr(ccp->GetPlane().Point().PeekCoordinates());
        clipDat->SetW(grr.Dot(ccp->GetPlane().Normal()));
        clipCol[0] = static_cast<float>(ccp->GetColour()[0]) / 255.0f;
        clipCol[1] = static_cast<float>(ccp->GetColour()[1]) / 255.0f;
        clipCol[2] = static_cast<float>(ccp->GetColour()[2]) / 255.0f;
        clipCol[3] = static_cast<float>(ccp->GetColour()[3]) / 255.0f;
        clipPoint[0] = grr.GetX();
        clipPoint[1] = grr.GetY();
        clipPoint[2] = grr.GetZ();
        //  clipPoint[3] = d.Dot(ccp->GetPlane().Normal());

    } else {

        clipDat->SetX(0.0f);
        clipDat->SetY(0.0f);
        clipDat->SetZ(0.0f);
        clipDat->SetW(0.0f);
        clipCol[0] = clipCol[1] = clipCol[2] = 0.75f;
        clipCol[3] = 1.0f;
    }
}


bool ModelClusterRenderer::prepareHBond_Cones() {
    LigandModelCall* lmc = this->ligandModel_callerSlot.CallAs<LigandModelCall>();
    if (lmc == NULL) {
        return false;
    }
    class protH_storage {
    public:
        glm::vec3 fromPos = glm::vec3(0.0, 0.0, 0.0);
        vislib::Array<glm::vec3> toPositions;
        float length = 0.0;
        uint posCnt = 0;
        glm::vec3 direction = glm::vec3(0.0, 0.0, 0.0);
    };


    if (lmc->getClusterAndCentroidData().DBSCANParams.paramsChanged ||
        this->s.selectedGlobalMdlID != this->s.old_selectedGlobalMdlID) {
        this->s.old_selectedGlobalMdlID = this->s.selectedGlobalMdlID;

        this->protHBonds_HDirection.Resize(0);
        this->protHBonds_HPos.Resize(0);
        this->protHBonds_coneLengths.Resize(0);
        this->protHBonds_coneRadii.Resize(0);

        this->modelHBonds_HDirection.Resize(0);
        this->modelHBonds_HPos.Resize(0);
        this->modelHBonds_coneLengths.Resize(0);
        this->modelHBonds_coneRadii.Resize(0);

        std::vector<uint> mldStor;
        std::map<uint, protH_storage> protMap;
        std::map<uint, protH_storage> mdlMap;

        int start = 0;
        int end = 0;

        // TODO: this is to expensive when deselcting a cluster
        // add one array which stores the data for a single model (which content is changed often)
        // add another array which permanently hold all the data for all models
        // => to get rid of interating over all models again and again!
        if (this->s.selectedGlobalMdlID > -1) {
            start = this->s.selectedGlobalMdlID;
            end = start + 1;
        } else {
            start = 0;
            end = lmc->getTotalModelCount();
        }
        for (int x = start; x < end; x++) {
            lmc->SetCurrentLigandAndModel_byGlobalMdlID(x);
            (*lmc)(LigandModelCall::CallForGetData);
            if (lmc->getCurrentBindingenergy() > lmc->getClusterAndCentroidData().DBSCANParams.minBindEnergy) {
                continue;
            }
            auto hbondsDON = lmc->getInteractionForces().hProtAcceptors;
            auto hbondsACC = lmc->getInteractionForces().hProtDonors;
            if (this->s.selectedGlobalMdlID == x || this->s.selectedGlobalMdlID == -1) {

                /************************************
                 ****** ligand model acceptor *******
                 ************************************/
                for (int i = 0; i < hbondsACC->cnt; i++) {

                    // model
                    glm::vec3 mdlAtom = hbondsACC->getLigPos(i);
                    // protein
                    glm::vec3 protAtomH = hbondsACC->getProtPos(i);
                    glm::vec3 direction = glm::normalize(mdlAtom - protAtomH);
                    float length = glm::length(mdlAtom - protAtomH);

                    int protID_Donor = hbondsACC->getProtID(i);
                    // H-Bond direction
                    // printf("pos%i-%i: %f %f %f--%f %f %f\n", mdlID, protID_H, mdlAtom.x, mdlAtom.y,
                    // mdlAtom.z,
                    //    protAtomH.x, protAtomH.y, protAtomH.z);
                    if (protMap.find(protID_Donor) != protMap.end()) {
                        // add to existing object

                        auto p = protMap[protID_Donor];
                        p.direction += direction;
                        p.toPositions.Append(mdlAtom);
                        p.posCnt++;
                        p.length += length;
                        protMap[protID_Donor] = p;
                    } else {
                        // create object
                        protH_storage t;
                        protMap[protID_Donor] = t;

                        auto p = protMap[protID_Donor];
                        p.direction += direction;
                        p.toPositions.Append(mdlAtom);
                        p.fromPos = protAtomH;
                        p.length += length;
                        p.posCnt = 1;
                        protMap[protID_Donor] = p;
                    }
                }


                /************************************
                 ******** ligand model donor ********
                 ************************************/
                bool donorToAcceptor = false;
                for (int i = 0; i < hbondsDON->cnt; i++) {


                    glm::vec3 mdlAtomH = hbondsDON->getLigPos(i);
                    glm::vec3 protAtomAccept = hbondsDON->getProtPos(i);
                    int mdlID_Hdonor = hbondsDON->getLigID(i);
                    int protID_Acceptor = hbondsDON->getProtID(i);
                    // printf("%d prot:%d\n", i, protID_Acceptor);
                    if (this->s.selectedGlobalMdlID == x || this->s.selectedGlobalMdlID == -1) {
                        // model
                        float length = glm::length(protAtomAccept - mdlAtomH);
                        if (donorToAcceptor) {
                            int mdlID_donorReal = mdlID_Hdonor + x;
                            glm::vec3 direction = glm::normalize(protAtomAccept - mdlAtomH);
                            if (mdlMap.find(mdlID_donorReal) != mdlMap.end()) {
                                // add to existing object

                                auto p = mdlMap[mdlID_donorReal];
                                p.direction += direction;
                                p.toPositions.Append(protAtomAccept);
                                p.posCnt++;
                                p.length += length;
                                mdlMap[mdlID_donorReal] = p;
                            } else {
                                // create object
                                protH_storage t;
                                mdlMap[mdlID_donorReal] = t;

                                auto p = mdlMap[mdlID_donorReal];
                                p.direction += direction;
                                p.toPositions.Append(protAtomAccept);
                                p.fromPos = mdlAtomH;
                                p.length += length;
                                p.posCnt = 1;
                                mdlMap[mdlID_donorReal] = p;
                            }
                        } else {
                            int mdlID_donorReal = protID_Acceptor;
                            glm::vec3 direction = glm::normalize(mdlAtomH - protAtomAccept);
                            if (mdlMap.find(mdlID_donorReal) != mdlMap.end()) {
                                // add to existing object

                                auto p = mdlMap[mdlID_donorReal];
                                p.direction += direction;
                                p.toPositions.Append(mdlAtomH);
                                p.posCnt++;
                                p.length += length;
                                mdlMap[mdlID_donorReal] = p;
                            } else {
                                // create object
                                protH_storage t;
                                mdlMap[mdlID_donorReal] = t;

                                auto p = mdlMap[mdlID_donorReal];
                                p.direction += direction;
                                p.toPositions.Append(mdlAtomH);
                                p.fromPos = protAtomAccept;
                                p.length += length;
                                p.posCnt = 1;
                                mdlMap[mdlID_donorReal] = p;
                            }
                        }
                    }
                }
            }
        }

        float minConeRad = 0.1f;
        float maxConeRad = 2.0f;
        /************************************
         ****** ligand model acceptor *******
         ************************************/
        // printf("ligand model acceptor:\n");
        float radius = 0;
        int HBondCnt = 0;
        int itCnt = 0;
        map<uint, protH_storage>::iterator it;
        for (it = protMap.begin(); it != protMap.end(); it++) {
            radius = 0;
            HBondCnt = 0;

            // set H-Bond start pos
            // determine average H-Bond direction
            // determine average H-Bond length
            protHBonds_HDirection.Append(it->second.direction.x / it->second.posCnt);
            protHBonds_HDirection.Append(it->second.direction.y / it->second.posCnt);
            protHBonds_HDirection.Append(it->second.direction.z / it->second.posCnt);
            protHBonds_HPos.Append(it->second.fromPos.x);
            protHBonds_HPos.Append(it->second.fromPos.y);
            protHBonds_HPos.Append(it->second.fromPos.z);
            protHBonds_coneLengths.Append(it->second.length / it->second.posCnt);

            // determine average radius of H-Bond-Cone
            // iterate over all H-bonds
            glm::vec3 startPos(
                protHBonds_HPos[3 * itCnt + 0], protHBonds_HPos[3 * itCnt + 1], protHBonds_HPos[3 * itCnt + 2]);
            glm::vec3 dir(protHBonds_HDirection[3 * itCnt + 0], protHBonds_HDirection[3 * itCnt + 1],
                protHBonds_HDirection[3 * itCnt + 2]);
            float avgLength = this->protHBonds_coneLengths[itCnt];
            dir = glm::normalize(dir);
            for (int j = 0; j < it->second.toPositions.Count(); j++) {
                glm::vec3 toPosDir = glm::normalize(it->second.toPositions[j] - startPos);
                float dot = glm::dot(toPosDir, dir);
                float radians = glm::tan(glm::acos(dot));
                radius += glm::max(radians * avgLength, minConeRad);
                HBondCnt = it->second.posCnt;
            }
            float tmpRad = radius / (float)HBondCnt;
            if (tmpRad > maxConeRad) {
                tmpRad = maxConeRad;
            }
            protHBonds_coneRadii.Append(tmpRad);
            //  printf("Radius: %f \t Hcnt:%i rad:%f length:%f\n", protHBonds_coneRadii[itCnt], HBondCnt,
            //  radius,
            //      avgLength);
            itCnt++;
        }
        // calculate radius with pyhtagoras => tan(toPosDir,dir) * ankathete(avgLength);


        /************************************
         ******** ligand model donor ********
         ************************************/
        // printf("ligand model donor:\n");
        itCnt = 0;
        for (it = mdlMap.begin(); it != mdlMap.end(); it++) {
            radius = 0;
            HBondCnt = 0;

            // set H-Bond start pos
            // determine average H-Bond direction
            // determine average H-Bond length
            modelHBonds_HDirection.Append(it->second.direction.x / it->second.posCnt);
            modelHBonds_HDirection.Append(it->second.direction.y / it->second.posCnt);
            modelHBonds_HDirection.Append(it->second.direction.z / it->second.posCnt);
            modelHBonds_HPos.Append(it->second.fromPos.x);
            modelHBonds_HPos.Append(it->second.fromPos.y);
            modelHBonds_HPos.Append(it->second.fromPos.z);
            modelHBonds_coneLengths.Append(it->second.length / it->second.posCnt);

            // determine average radius of H-Bond-Cone
            // iterate over all H-bonds
            glm::vec3 startPos(
                modelHBonds_HPos[3 * itCnt + 0], modelHBonds_HPos[3 * itCnt + 1], modelHBonds_HPos[3 * itCnt + 2]);
            glm::vec3 dir(modelHBonds_HDirection[3 * itCnt + 0], modelHBonds_HDirection[3 * itCnt + 1],
                modelHBonds_HDirection[3 * itCnt + 2]);
            float avgLength = this->modelHBonds_coneLengths[itCnt];
            dir = glm::normalize(dir);
            for (int j = 0; j < it->second.toPositions.Count(); j++) {
                glm::vec3 toPosDir = glm::normalize(it->second.toPositions[j] - startPos);
                float dot = glm::dot(toPosDir, dir);
                float radians = glm::tan(glm::acos(dot));
                radius += glm::max(radians * avgLength, minConeRad);
                HBondCnt = it->second.posCnt;
            }
            float tmpRad = radius / (float)HBondCnt;
            if (tmpRad > maxConeRad) {
                tmpRad = maxConeRad;
            }
            modelHBonds_coneRadii.Append(tmpRad);
            //  printf("Radius: %f \t Hcnt:%i rad:%f length:%f\n", modelHBonds_coneRadii[itCnt], HBondCnt,
            //  radius,
            //       avgLength);
            itCnt++;
        }
        // calculate radius with pyhtagoras => tan(toPosDir,dir) * ankathete(avgLength);
    }
}


int ModelClusterRenderer::ForcesLineStorage::getPosCnt(std::string forceName, megamol::core::Call* ligandModelCall) {
    // handle inputs and set class states
    LigandModelCall* lmc = dynamic_cast<LigandModelCall*>(ligandModelCall);
    if (lmc->getWebData()->getClusterID() == -1) {
        this->singleModel = false;
        checkDataExists(forceName, lmc);
        return this->posArraysMap[forceName].size() / 4;
    }
    this->singleModel = true;
    // get the data
    checkDataExists(forceName, lmc);
    return this->posArraysMap_singleModel[forceName].size() / 4;
};


float* ModelClusterRenderer::ForcesLineStorage::getCurPosArray(
    std::string forceName, megamol::core::Call* ligandModelCall) {
    // handle inputs and set class states
    LigandModelCall* lmc = dynamic_cast<LigandModelCall*>(ligandModelCall);
    if (lmc->getWebData()->getClusterID() == -1) {
        this->singleModel = false;
        checkDataExists(forceName, lmc);
        return this->posArraysMap[forceName].data();
    }
    this->singleModel = true;
    // get the data
    checkDataExists(forceName, lmc);
    return this->posArraysMap_singleModel[forceName].data();
};

float* ModelClusterRenderer::ForcesLineStorage::getCurStartPosArray(
    std::string forceName, megamol::core::Call* ligandModelCall) {
    // handle inputs and set class states
    LigandModelCall* lmc = dynamic_cast<LigandModelCall*>(ligandModelCall);
    if (lmc->getWebData()->getClusterID() == -1) {
        this->singleModel = false;
        checkDataExists(forceName, lmc);
        return this->startPositionsMap[forceName].data();
    }
    this->singleModel = true;
    // get the data
    checkDataExists(forceName, lmc);
    return this->startPositionsMap_singleModel[forceName].data();
};


float* ModelClusterRenderer::ForcesLineStorage::getCurDirArray(
    std::string forceName, megamol::core::Call* ligandModelCall) {
    // handle inputs and set class states
    LigandModelCall* lmc = dynamic_cast<LigandModelCall*>(ligandModelCall);
    if (lmc->getWebData()->getClusterID() == -1) {
        this->singleModel = false;
        checkDataExists(forceName, lmc);
        return this->dircetionsMap[forceName].data();
    }
    this->singleModel = true;
    // get the data
    checkDataExists(forceName, lmc);
    return this->dircetionsMap_singleModel[forceName].data();
};


float* ModelClusterRenderer::ForcesLineStorage::getCurLengthArray(
    std::string forceName, megamol::core::Call* ligandModelCall) {
    // handle inputs and set class states
    LigandModelCall* lmc = dynamic_cast<LigandModelCall*>(ligandModelCall);
    if (lmc->getWebData()->getClusterID() == -1) {
        this->singleModel = false;
        checkDataExists(forceName, lmc);
        return this->lengthsMap[forceName].data();
    }
    this->singleModel = true;
    // get the data
    checkDataExists(forceName, lmc);
    return this->lengthsMap_singleModel[forceName].data();
};


uint* ModelClusterRenderer::ForcesLineStorage::getCurForceCntPerProtAtomArray(
    std::string forceName, int protAtomCnt, megamol::core::Call* ligandModelCall) {
    // handle inputs and set class states
    LigandModelCall* lmc = dynamic_cast<LigandModelCall*>(ligandModelCall);
    checkDataExists(forceName, lmc, protAtomCnt);
    return this->forceCntPerProtAtomMap[forceName].data();
};


bool ModelClusterRenderer::ForcesLineStorage::checkDataExists(
    std::string forceName, megamol::core::Call* call, int protAtomCnt) {
    LigandModelCall* lmc = dynamic_cast<LigandModelCall*>(call);
    if (lmc == NULL) return false;

    auto posMap = &this->posArraysMap;
    auto startPosMap = &this->startPositionsMap;
    auto dirMap = &this->dircetionsMap;
    auto lengthMap = &this->lengthsMap;


    if (this->singleModel) {
        if (this->old_glblMdlID != lmc->getCurrentGlobalModelID()) {
            this->old_glblMdlID = lmc->getCurrentGlobalModelID();
            printf("gblMdlID: %d\n", lmc->getCurrentGlobalModelID());
            auto ligForce = lmc->getInteractionForces();
            this->posArraysMap_singleModel.clear();
            this->startPositionsMap_singleModel.clear();
            this->dircetionsMap_singleModel.clear();
            this->lengthsMap_singleModel.clear();
        }
        posMap = &this->posArraysMap_singleModel;
        startPosMap = &this->startPositionsMap_singleModel;
        dirMap = &this->dircetionsMap_singleModel;
        lengthMap = &this->lengthsMap_singleModel;
    }


    if (posMap->find(forceName) == posMap->end()) {
        // froce key not found
        int rangeStart, rangeEnd;
        auto curFPostions = &posMap[0][forceName];
        auto curStartPos = &startPosMap[0][forceName];
        auto curDir = &dirMap[0][forceName];
        auto curLenght = &lengthMap[0][forceName];


        if (lmc->getWebData()->getClusterID() != -1) {
            rangeStart = lmc->getCurrentGlobalModelID();
            rangeEnd = rangeStart + 1;
        } else {
            rangeStart = 0;
            rangeEnd = lmc->getTotalModelCount();
        }
        for (int i = rangeStart; i < rangeEnd; i++) {
            lmc->SetCurrentLigandAndModel_byGlobalMdlID(i);
            (*lmc)(LigandModelCall::CallForGetData);
            (*lmc)(LigandModelCall::CallForGetClusterAndCentroidData);
            if (lmc->getCurrentBindingenergy() <= lmc->getClusterAndCentroidData().DBSCANParams.minBindEnergy &&
                lmc->getClusterAndCentroidData().assignedClusters[0][lmc->getCurrentGlobalModelID()] != -1) {
                auto fcc = lmc->getInteractionForces();

                glm::vec3 a, b;
                // saltBridges
                if (forceName == this->interForceTypes.saltBridges) {
                    auto fc = fcc.saltBridges;
                    for (int j = 0; j < fc->cnt; j++) {
                        a = fc->getProtPos(j);
                        b = fc->getLigPos(j);
                        curFPostions->push_back(a.x);
                        curFPostions->push_back(a.y);
                        curFPostions->push_back(a.z);
                        curFPostions->push_back(this->radius);
                        curFPostions->push_back(b.x);
                        curFPostions->push_back(b.y);
                        curFPostions->push_back(b.z);
                        curFPostions->push_back(this->radius);

                        glm::vec3 dir = b - a;
                        curStartPos->push_back(a.x);
                        curStartPos->push_back(a.y);
                        curStartPos->push_back(a.z);
                        curStartPos->push_back(this->radius);
                        float length = glm::length(dir);
                        dir = glm::normalize(dir);
                        curDir->push_back(dir.x);
                        curDir->push_back(dir.y);
                        curDir->push_back(dir.z);
                        curLenght->push_back(length);
                    }
                    // halogenBonds
                } else if (forceName == this->interForceTypes.halogenBonds) {
                    auto fc = fcc.halogenBonds;
                    for (int j = 0; j < fc->cnt; j++) {
                        a = fc->getProtPos(j);
                        b = fc->getLigPos(j);
                        curFPostions->push_back(a.x);
                        curFPostions->push_back(a.y);
                        curFPostions->push_back(a.z);
                        curFPostions->push_back(this->radius);
                        curFPostions->push_back(b.x);
                        curFPostions->push_back(b.y);
                        curFPostions->push_back(b.z);
                        curFPostions->push_back(this->radius);

                        glm::vec3 dir = b - a;
                        curStartPos->push_back(a.x);
                        curStartPos->push_back(a.y);
                        curStartPos->push_back(a.z);
                        curStartPos->push_back(this->radius);
                        float length = glm::length(dir);
                        dir = glm::normalize(dir);
                        curDir->push_back(dir.x);
                        curDir->push_back(dir.y);
                        curDir->push_back(dir.z);
                        curLenght->push_back(length);
                    }
                    // hydrophobicInteractions
                } else if (forceName == this->interForceTypes.hydrophobicInteractions) {
                    auto fc = fcc.hydrophobicInteractions;
                    for (int j = 0; j < fc->cnt; j++) {
                        a = fc->getProtPos(j);
                        b = fc->getLigPos(j);
                        curFPostions->push_back(a.x);
                        curFPostions->push_back(a.y);
                        curFPostions->push_back(a.z);
                        curFPostions->push_back(this->radius);
                        curFPostions->push_back(b.x);
                        curFPostions->push_back(b.y);
                        curFPostions->push_back(b.z);
                        curFPostions->push_back(this->radius);

                        glm::vec3 dir = b - a;
                        curStartPos->push_back(a.x);
                        curStartPos->push_back(a.y);
                        curStartPos->push_back(a.z);
                        curStartPos->push_back(this->radius);
                        float length = glm::length(dir);
                        dir = glm::normalize(dir);
                        curDir->push_back(dir.x);
                        curDir->push_back(dir.y);
                        curDir->push_back(dir.z);
                        curLenght->push_back(length);
                    }
                    // metalComplexes
                } else if (forceName == this->interForceTypes.metalComplexes) {
                    auto fc = fcc.metalComplexes;
                    for (int j = 0; j < fc->cnt; j++) {
                        a = fc->getMetalPos(j);
                        b = fc->getTargetPos(j);
                        curFPostions->push_back(a.x);
                        curFPostions->push_back(a.y);
                        curFPostions->push_back(a.z);
                        curFPostions->push_back(this->radius);
                        curFPostions->push_back(b.x);
                        curFPostions->push_back(b.y);
                        curFPostions->push_back(b.z);
                        curFPostions->push_back(this->radius);

                        glm::vec3 dir = b - a;
                        curStartPos->push_back(a.x);
                        curStartPos->push_back(a.y);
                        curStartPos->push_back(a.z);
                        curStartPos->push_back(this->radius);
                        float length = glm::length(dir);
                        dir = glm::normalize(dir);
                        curDir->push_back(dir.x);
                        curDir->push_back(dir.y);
                        curDir->push_back(dir.z);
                        curLenght->push_back(length);
                    }
                    // piCationInteractions
                } else if (forceName == this->interForceTypes.piCationInteractions) {
                    auto fc = fcc.piCationInteractions;
                    for (int j = 0; j < fc->cnt; j++) {
                        a = fc->getProtPos(j);
                        b = fc->getLigPos(j);
                        curFPostions->push_back(a.x);
                        curFPostions->push_back(a.y);
                        curFPostions->push_back(a.z);
                        curFPostions->push_back(this->radius);
                        curFPostions->push_back(b.x);
                        curFPostions->push_back(b.y);
                        curFPostions->push_back(b.z);
                        curFPostions->push_back(this->radius);

                        glm::vec3 dir = b - a;
                        curStartPos->push_back(a.x);
                        curStartPos->push_back(a.y);
                        curStartPos->push_back(a.z);
                        curStartPos->push_back(this->radius);
                        float length = glm::length(dir);
                        dir = glm::normalize(dir);
                        curDir->push_back(dir.x);
                        curDir->push_back(dir.y);
                        curDir->push_back(dir.z);
                        curLenght->push_back(length);
                    }
                    // piStacks
                } else if (forceName == this->interForceTypes.piStacks) {
                    auto fc = fcc.piStacks;
                    for (int j = 0; j < fc->cnt; j++) {
                        a = fc->getProtPos(j);
                        b = fc->getLigPos(j);
                        curFPostions->push_back(a.x);
                        curFPostions->push_back(a.y);
                        curFPostions->push_back(a.z);
                        curFPostions->push_back(this->radius);
                        curFPostions->push_back(b.x);
                        curFPostions->push_back(b.y);
                        curFPostions->push_back(b.z);
                        curFPostions->push_back(this->radius);

                        glm::vec3 dir = b - a;
                        curStartPos->push_back(a.x);
                        curStartPos->push_back(a.y);
                        curStartPos->push_back(a.z);
                        curStartPos->push_back(this->radius);
                        float length = glm::length(dir);
                        dir = glm::normalize(dir);
                        curDir->push_back(dir.x);
                        curDir->push_back(dir.y);
                        curDir->push_back(dir.z);
                        curLenght->push_back(length);
                    }
                    // Hbonds_protAcceptor
                } else if (forceName == this->interForceTypes.Hbonds_protAcceptor) {

                    auto fc = fcc.hProtAcceptors;

                    for (int j = 0; j < fc->cnt; j++) {
                        a = fc->getProtPos(j);
                        b = fc->getLigPos(j);
                        curFPostions->push_back(a.x);
                        curFPostions->push_back(a.y);
                        curFPostions->push_back(a.z);
                        curFPostions->push_back(this->radius);
                        curFPostions->push_back(b.x);
                        curFPostions->push_back(b.y);
                        curFPostions->push_back(b.z);
                        curFPostions->push_back(this->radius);


                        glm::vec3 dir = b - a;
                        curStartPos->push_back(a.x);
                        curStartPos->push_back(a.y);
                        curStartPos->push_back(a.z);
                        curStartPos->push_back(this->radius);
                        float length = glm::length(dir);
                        dir = glm::normalize(dir);
                        curDir->push_back(dir.x);
                        curDir->push_back(dir.y);
                        curDir->push_back(dir.z);
                        curLenght->push_back(length);
                    }
                    // Hbonds_protDonor
                } else if (forceName == this->interForceTypes.Hbonds_protDonor) {

                    auto fc = fcc.hProtDonors;

                    for (int j = 0; j < fc->cnt; j++) {
                        a = fc->getProtPos(j);
                        b = fc->getLigPos(j);
                        curFPostions->push_back(a.x);
                        curFPostions->push_back(a.y);
                        curFPostions->push_back(a.z);
                        curFPostions->push_back(this->radius);
                        curFPostions->push_back(b.x);
                        curFPostions->push_back(b.y);
                        curFPostions->push_back(b.z);
                        curFPostions->push_back(this->radius);


                        glm::vec3 dir = b - a;
                        curStartPos->push_back(a.x);
                        curStartPos->push_back(a.y);
                        curStartPos->push_back(a.z);
                        curStartPos->push_back(this->radius);
                        float length = glm::length(dir);
                        dir = glm::normalize(dir);
                        curDir->push_back(dir.x);
                        curDir->push_back(dir.y);
                        curDir->push_back(dir.z);
                        curLenght->push_back(length);
                    }
                } else if (forceName == this->interForceTypes.HBonds) {
                    // do nothing
                } else {
                    throw printf("\nUnkown force name!: %s\n", forceName.c_str());
                }
            }
        }
    }
    auto protAtomCntMap = &this->forceCntPerProtAtomMap;
    if (protAtomCntMap->find(forceName) == protAtomCntMap->end() && protAtomCnt > 1) {
        // froce key not found

        auto curProtAtomCntMap = &protAtomCntMap[0][forceName];
        int rangeStart = 0;
        int rangeEnd = lmc->getTotalModelCount();

        for (int i = rangeStart; i < rangeEnd; i++) {
            lmc->SetCurrentLigandAndModel_byGlobalMdlID(i);
            (*lmc)(LigandModelCall::CallForGetData);
            (*lmc)(LigandModelCall::CallForGetClusterAndCentroidData);
            const auto& clusterData = lmc->getClusterAndCentroidData();
            if (lmc->getCurrentBindingenergy() < clusterData.DBSCANParams.minBindEnergy 
				&& clusterData.assignedClusters[0][lmc->getCurrentGlobalModelID()] != -1) 
			{
                auto fcc = lmc->getInteractionForces();

                glm::vec3 a, b;
                if (forceName == "saltBridges") {
                    auto fc = fcc.saltBridges;
                    curProtAtomCntMap[0].resize(protAtomCnt, 0);
                    for (int j = 0; j < fc->cnt; j++) {
                        for (int g = 0; g < fc->getProtIDsCnt(j); g++) {
                            curProtAtomCntMap[0][(uint)fc->getProtIDs(j)[g]]++;
                        }
                    }
                } else if (forceName == "halogenBonds") {
                    auto fc = fcc.halogenBonds;
                    curProtAtomCntMap[0].resize(protAtomCnt, 0);
                    for (int j = 0; j < fc->cnt; j++) {
                        curProtAtomCntMap[0][(uint)fc->getProtID(j)]++;
                    }
                } else if (forceName == "hydrophobicInteractions") {
                    auto fc = fcc.hydrophobicInteractions;
                    curProtAtomCntMap[0].resize(protAtomCnt, 0);
                    for (int j = 0; j < fc->cnt; j++) {
                        curProtAtomCntMap[0][(uint)fc->getProtID(j)]++;
                    }
                } else if (forceName == "metalComplexes") {
                    auto fc = fcc.metalComplexes;
                    curProtAtomCntMap[0].resize(protAtomCnt, 0);
                    for (int j = 0; j < fc->cnt; j++) {
                        curProtAtomCntMap[0][(uint)fc->getMetalID(j)]++;
                    }
                } else if (forceName == "piCationInteractions") {
                    auto fc = fcc.piCationInteractions;
                    curProtAtomCntMap[0].resize(protAtomCnt, 0);
                    for (int j = 0; j < fc->cnt; j++) {
                        for (int g = 0; g < fc->getProtIDsCnt(j); g++) {
                            curProtAtomCntMap[0][(uint)fc->getProtIDs(j)[g]]++;
                        }
                    }
                } else if (forceName == "piStacks") {
                    auto fc = fcc.piStacks;
                    curProtAtomCntMap[0].resize(protAtomCnt, 0);
                    for (int j = 0; j < fc->cnt; j++) {
                        for (int g = 0; g < fc->getProtIDsCnt(j); g++) {
                            curProtAtomCntMap[0][(uint)fc->getProtIDs(j)[g]]++;
                        }
                    }
                } else if (forceName == this->interForceTypes.HBonds) {
                    auto fcA = fcc.hProtAcceptors;
                    auto fcD = fcc.hProtDonors;
                    curProtAtomCntMap[0].resize(protAtomCnt, 0);
                    for (int j = 0; j < fcA->cnt; j++) {
                        curProtAtomCntMap[0][(uint)fcA->getProtID(j)]++;
                    }
                    for (int j = 0; j < fcD->cnt; j++) {
                        curProtAtomCntMap[0][(uint)fcD->getProtID(j)]++;
                    }

                } else {
                    throw printf("\nUnkown force name!: %s\n", forceName.c_str());
                }
            }
        }
    }
    return true;
};

bool ModelClusterRenderer::updateBoolParams(
    core::param::ParamSlot* param, std::string paramName, LigandModelCall::WebData* webData, int* renderBool) {
    auto webGUI = (*webData->GUI);
    bool changedParam = false;
    if (param->IsDirty() || this->webDataChanged) {
        auto p = param->Param<param::BoolParam>();
        if (webGUI[paramName] != (float)p->Value() && this->webDataChanged) {
            p->SetValue((bool)webGUI[paramName]);
            param->ForceSetDirty();
        } else if (param->IsDirty()) {
            *renderBool = p->Value();
            GUI[paramName] = (float)p->Value();
            changedParam = true;
            param->ResetDirty();
        }
    }
    return changedParam;
}


bool ModelClusterRenderer::updateIntParams(
    core::param::ParamSlot* param, std::string paramName, LigandModelCall::WebData* webData, int* renderInt) {
    auto webGUI = (*webData->GUI);
    bool changedParam = false;
    if (param->IsDirty() || this->webDataChanged) {
        auto p = param->Param<param::IntParam>();
        if (webGUI[paramName] != (float)p->Value() && this->webDataChanged) {
            p->SetValue((int)webGUI[paramName]);
            param->ForceSetDirty();
        } else if (param->IsDirty()) {
            *renderInt = p->Value();
            GUI[paramName] = (float)p->Value();
            changedParam = true;
            param->ResetDirty();
        }
    }
    return changedParam;
}

bool ModelClusterRenderer::updateEnumParams(
    core::param::ParamSlot* param, std::string paramName, LigandModelCall::WebData* webData, int* renderInt) {
    auto webGUI = (*webData->GUI);
    bool changedParam = false;
    if (param->IsDirty() || this->webDataChanged) {
        auto p = param->Param<param::EnumParam>();
        if (webGUI[paramName] != (float)p->Value() && this->webDataChanged) {
            p->SetValue((int)webGUI[paramName]);
            param->ForceSetDirty();
        } else if (param->IsDirty()) {
            *renderInt = p->Value();
            GUI[paramName] = (float)p->Value();
            changedParam = true;
            param->ResetDirty();
        }
    }
    return changedParam;
}

///////////////////////////////////////////


/**
 * sort function for a int vec3 by Y()
 */
int sortVecf3ByY(const vislib::math::Vector<float, 3U>& lhs, const vislib::math::Vector<float, 3U>& rhs) {
    if (lhs.GetY() < rhs.GetY()) {
        return -1;
    } else if (lhs.GetY() > rhs.GetY()) {
        return 1;
    } else {
        return 0;
    }
}

/**
 * sort function for a int vec3 by Y     if Z is 0 for lhs&rhs       or if Z > 0 for lhs&rhs
 */
int sortVecf3ByYifZisNoneZero(const vislib::math::Vector<float, 3U>& lhs, const vislib::math::Vector<float, 3U>& rhs) {
    if (lhs.GetZ() > 0 && rhs.GetZ() > 0) {
        return sortVecf3ByY(lhs, rhs);
    } else if (lhs.GetZ() == 0 && rhs.GetZ() > 0) {
        return 1;
    } else if (lhs.GetZ() > 0 && rhs.GetZ() == 0) {
        return -1;
    } else {
        return sortVecf3ByY(lhs, rhs);
    }
}
