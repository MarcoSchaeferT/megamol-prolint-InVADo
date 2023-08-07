/*
 * ModelClusterRenderer.h
 *
 * Copyright (C) 2019 by University of Tuebingen (BDVA).
 * All rights reserved.
 */

#ifndef PROLINT_PLUGIN_MODELCLUSTERRENDERER_H_INCLUDED
#define PROLINT_PLUGIN_MODELCLUSTERRENDERER_H_INCLUDED
typedef unsigned int uint;

#include "LigandModelCall.h"
#include "geometry_calls/CallTriMeshData.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"
#include "mmcore/misc/VolumetricDataCall.h"
#include "mmcore/param//EnumParam.h"
#include "mmcore/param/BoolParam.h"
#include "mmcore/param/FloatParam.h"
#include "mmcore/param/IntParam.h"
#include "mmcore/param/ParamSlot.h"
#include "mmcore/param/Vector3fParam.h"
#include "mmcore/utility/log/Log.h"
#include "mmcore/view/AbstractCallRender3D.h"
#include "mmcore/view/Renderer3DModule.h"
#include "protein_calls/MolecularDataCall.h"
#include "vislib/Array.h"
#include "vislib/graphics/Cursor2D.h"
#include "vislib/graphics/gl/GLSLGeometryShader.h"
#include "vislib/graphics/gl/GLSLShader.h"
#include "vislib/graphics/gl/ShaderSource.h"
#include "vislib/math/Vector.h"

#include "glowl/Texture.hpp"
#include "glowl/Texture2D.hpp"
#include "glowl/Texture3D.hpp"

#include <GL/glu.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include "ModelClusterRenderer.cuh"
#include "MouseInputCall.h"
#include "mmcore/utility/SDFFont.h"
#include "mmcore/view/AbstractView.h"
#include "mmcore/view/CallClipPlane.h"
#include "mmcore/view/CallGetTransferFunction.h"
#include "mmcore/view/MouseFlags.h"

#include <filesystem>
#include <map>
#include "functionalGroupHandler.h"
#include "png.h"


namespace megamol {
namespace prolint {
typedef vislib::Array<vislib::math::Vector<float, 4>> visFloat4;

/************************************
 **** ModelClusterRenderer CLASS ****
 ************************************/

class ModelClusterRenderer : public megamol::core::view::Renderer3DModule {


public:
    /** Ctor. */
    ModelClusterRenderer();
    /** Dtor. */
    virtual ~ModelClusterRenderer();

    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName(void) { return "ModelClusterRenderer"; }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description(void) {
        return "Selects docking data, calculates centroids, cluster these, and renders the center of the clsuters as "
               "spheres.";
    }

    /**
     * Answers whether this module is available on the current system.
     *
     * @return 'true' if the module is available, 'false' otherwise.
     */
    static bool IsAvailable(void) { return true; }

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
     * Callback for mouse events (move, press, and release)
     *
     * @param[in] x The x coordinate of the mouse in screen space
     * @param[in] y The y coordinate of the mouse in screen space
     * @param[in] flags The mouse flags
     * @return 'true' on success
     */
    virtual bool MouseEvent(int x, int y, core::view::MouseFlags flags);

    /**
     * Callback function for StarMenu params
     */
    bool paramsStarMenuChanged(megamol::core::param::ParamSlot& param);

    /**
     * calculates coordinates for lines (2 trianlges-> 4 points)
     *
     * @param startPoint	defines the point where the line starts
     * @param length		float factor that will be apllied on a normalized line direction vector to grow line
     * @param width			float factor that will be apllied on a orthogonal (to line direction) vector // + and - from
     * the start- and end-point
     * @param lineDirection	the direction of the line
     * @param eyeDirection	the look at direction of the camera
     * @param lines			a buffer for storing coordinates of the lines (each 4 points) // will be expanded
     * dynamically
     */
    void getLineCoords(glm::vec3 startPoint, float length, float width, glm::vec3 lineDirection, glm::vec3 eyeDirection,
        vislib::Array<vislib::math::Vector<float, 3>>& lines, int offset);


    /*
     * variables
     */
    float scaleWorld;
    bool verbose = false;

    // a buffer for storing coordinates of the lines
    vislib::Array<vislib::math::Vector<float, 3>> lines;
    vislib::Array<vislib::math::Vector<float, 3>> SVGLineCoords;
    vislib::Array<vislib::math::Vector<float, 4>> SVGLineColors;
    vislib::Array<vislib::math::Vector<float, 3>> SVGTextCoords;
    vislib::Array<vislib::math::Vector<float, 4>> SVGTextColors;
    vislib::Array<vislib::StringA> SVGTextContent;
    vislib::Array<int> SVGTextSizes;
    vislib::Array<float> SVGScales;
    vislib::Array<vislib::math::Vector<float, 3>> SVGPolygons;

    /**************************************
     **** MOUSE RAY INTERSECTION CLASS ****
     **************************************/ /* ModuleClusterRenderer-->here */

    class MouseRayIntersection {
    public:
        /** Ctor. */
        MouseRayIntersection() { firstRayIntersetion = -2; };
        /** Dtor. */
        virtual ~MouseRayIntersection() { release(); }
        /**
         * Answer the name of this module.
         *
         * @return The name of this module.
         */
        static const char* ClassName(void) { return "MouseRayIntersection"; }

        /**
         * Answer a human readable description of this module.
         *
         * @return A human readable description of this module.
         */
        static const char* Description(void) {
            return "calculates mouse ray and check for intersections for a sphere list";
        }


        /*
         * public functions
         */


    protected:
        /**
         * Implementation of 'MouseRayIntersection'.
         */
        void release(void);


        //
        /** variables for MouseRayIntersection
         *
         *
         * mouseX, mouseY:		mouse position (screen coordinates)
         * firstRayIntersetion:	indiex of first intersecting sphere
         * rayIntersections:	indices of intersecting spheres
         * rayOrigin:			origin of the ray (world coordinates)
         * rayDirection:		direction of the ray (world coordinates)
         * priviousCamPos:		privious postion of camera (world coordinates)
         * cam:					cam object (used to get the view matrix [world_to_camera])
         **/

        float mouseX, mouseY;
        int firstRayIntersetion;
        vislib::Array<int> rayIntersections;
        vislib::Array<float> spheres;
        glm::vec3 rayOrigin;
        glm::vec3 rayDirection;
        glm::vec3 priviousCamPos;
        vislib::graphics::gl::CameraOpenGL cam;

        /**********************************
         **** TRANSFROM MATRICES CLASS ****
         **********************************/ /* ModuleClusterRenderer-->MouseRayIntersection-->here */
    public:
        class TransformMatrices {
        public:
            /** Ctor. */
            TransformMatrices(){};
            /** Dtor. */
            virtual ~TransformMatrices(){};

            /*
             * Pscreen = Mviewport * (Mprojection * (Mview * (Mmodel * Pobject)))
             *
             * GL_MODELVIEW_MATRIX	local -> camera
             * model matrix	 =	local -> world
             * camera matrix =	world -> camera
             *
             *
             * GL_PROJECTION_MATRIX
             * projection matrix	=	camera -> canonical view
             *
             * MEGAMOL viewport
             * viewport matrix = canonical view -> screen
             */

            //***** FORWARD *****
            glm::mat4 local_to_camera;
            glm::mat4 local_to_world;
            glm::mat4 world_to_camera;
            glm::mat4 camera_to_canonical;
            // viewport canonical_to_screen

            // ***** BACKWARD *****
            // viewport screen_to_canoniocal
            glm::mat4 canonical_to_camera;
            glm::mat4 camera_to_local;
            glm::mat4 camera_to_world;
            glm::mat4 world_to_local;
        };

    private:
        TransformMatrices transformMatrices;

        /**
         * sets all transformation matrices (from local-space until screen-space and the othe way around)
         *
         * @param tmpCam	(vislib::graphics::gl::CameraOpenGL& object)
         */
        void setTransformMatrices(vislib::graphics::gl::CameraOpenGL& tmpCam);

        /**
         * sets the origin of the ray in world coordinates (e.g camera position)
         *
         * @param rayOrigin		(world coordinates)
         */
        void setRayOrigin(vislib::math::Vector<float, 3> rayOrigin);

        /**
         * sets the direction vector of the ray (world coordinates)
         *
         * @param mouseX	(screen coordinates)
         * @param mouseY	(screen coordinates)
         * @param matrices	(object with transformation matrices for different 3D-spaces)
         */
        void setRayDirection(float mouseX, float mouseY, TransformMatrices matrices);

        /**
         * sets indices of intersecting spheres to an array
         *
         * @param spheres		array with spheres (local coordinates)
         * @param rayOrigin		(world coordinates)
         * @param rayDirection	(world coordinates)
         */
        void setRayIntersections(vislib::Array<float> spheres, glm::vec3 rayOrigin, glm::vec3 rayDirection);

    public:
        /**
         * returns the transformation matrices object (for different 3D-sapces)
         *
         * @param tmpCam		vislib::graphics::gl::CameraOpenGL&
         */
        TransformMatrices getTransformMatrices(vislib::graphics::gl::CameraOpenGL& tmpCam);

        /**
         * returns indices of intersecting spheres
         *
         * @param spheres		array with spheres (local coordinates)
         * @param tmpCam		vislib::graphics::gl::CameraOpenGL&
         * @param rayOrigin		(world coordinates)
         * @param rayDirection	(world coordinates)
         * @param mouseAction	(0 or 1->action)
         */
        vislib::Array<int> getMouseRayIntersections(vislib::Array<float>& spheres,
            vislib::graphics::gl::CameraOpenGL& tmpCam, float mouseX, float mouseY, int mouseAction);

        /**
         * returns indiex of first intersecting sphere
         *
         * @param spheres		array with spheres (local coordinates)
         * @param tmpCam		vislib::graphics::gl::CameraOpenGL&
         * @param rayOrigin		(world coordinates)
         * @param rayDirection	(world coordinates)
         * @param mouseAction	(0 or 1->action)
         */
        int getFirstMouseRayIntersection(vislib::Array<float>& spheres, vislib::graphics::gl::CameraOpenGL& tmpCam,
            float mouseX, float mouseY, int mouseAction);


    }; // namespace prolint
    MouseRayIntersection mouseRayIntersection;
    MouseRayIntersection::TransformMatrices transformMatrices;


private:
    /**********************************************************************
     * 'render'-functions
     **********************************************************************/

    /**
     * The get extents callback. The module should set the members of
     * 'call' to tell the caller the extents of its data (bounding boxes
     * and times).
     *
     * @param call The calling call.
     *
     * @return The return value of the function.
     */
    virtual bool GetExtents(megamol::core::Call& call);


    /**
     * The Open GL Render callback.
     *
     * @param call The calling call.
     * @return The return value of the function.
     */
    virtual bool Render(megamol::core::Call& call);


    /*******************************************
    ********* COMMON RENDERER VARIABLES ********
    ********************************************/
    GLfloat modelMatrix_column[16];
    vislib::math::Matrix<GLfloat, 4, vislib::math::COLUMN_MAJOR>* modelMatrix;
    vislib::math::Matrix<GLfloat, 4, vislib::math::COLUMN_MAJOR>* projMatrix;
    GLfloat projMatrix_column[16];
    glm::vec3 camPos_local;

    /*******************************
    ********* MESH RENDERER ********
    ********************************/
    // binding points
    const GLuint meshVertices_ColorSwitch_SSBObindingPoint = 0;
    const GLuint meshVertices_ClusterAssignment_SSBObindingPoint = 1;

    // variables
    GLuint gl_meshVertices = 0;
    GLuint gl_meshIndices = 0;
    GLuint gl_meshSortedIndices = 0;
    cudaGraphicsResource* gl_vertices_resource;
    cudaGraphicsResource* gl_indices_resource;
    cudaGraphicsResource* gl_sortedIndices_resource;

    GLuint gl_meshNormals = 0;
    GLuint gl_meshColors = 0;
    GLuint gl_meshVertices_ColorSwitch = 0;
    GLuint gl_meshVertices_ClusterAssignment = 0;


    size_t old_tmc_dataHash;
    float* d_TriangleToCamDistances;
    float3* d_camPos;

    // Ambient Occlusion
    std::unique_ptr<glowl::Texture3D> aoVolumeTexture;
    const core::misc::VolumetricDataCall::Metadata* volumeMetaData;
    float aoSampFact;
    float aoSampDist;
    float aoIntensityMin;
    float aoIntensityMax;
    int old_vdc_DataHash;
    float searchRadius_PocketPatch;


    class InterForceTypes {
    public:
        InterForceTypes() {
            this->hydrophobicInteractions = "hydrophobicInteractions";
            this->saltBridges = "saltBridges";
            this->piStacks = "piStacks";
            this->piCationInteractions = "piCationInteractions";
            this->halogenBonds = "halogenBonds";
            this->metalComplexes = "metalComplexes";
            this->Hbonds_protAcceptor = "Hbonds_protAcceptor";
            this->Hbonds_protDonor = "Hbonds_protDonor";
            this->HBonds = "HBonds";
            this->bestScores = "bestScores";
        };

        std::string hydrophobicInteractions;
        std::string saltBridges;
        std::string piStacks;
        std::string piCationInteractions;
        std::string halogenBonds;
        std::string metalComplexes;
        std::string Hbonds_protAcceptor;
        std::string Hbonds_protDonor;
        std::string HBonds;
        std::string bestScores;
        std::vector<string> forceTypeVec;
    };
    InterForceTypes interForceTypes;

    // forces line pos storage
    // @info: all forces without hBonbds
    class ForcesLineStorage {
    public:
        ForcesLineStorage() {
            this->curForce = "";
            this->singleModel = false;
            this->old_glblMdlID = -1;
            this->radius = 0.05;
        };

        // variables
        bool singleModel;
        int old_glblMdlID;
        std::string curForce;
        float radius;

        InterForceTypes interForceTypes;

        std::map<std::string, std::vector<uint>> forceCntPerProtAtomMap;

        std::map<std::string, std::vector<float>> posArraysMap;
        std::map<std::string, std::vector<float>> startPositionsMap;
        std::map<std::string, std::vector<float>> dircetionsMap;
        std::map<std::string, std::vector<float>> lengthsMap;

        std::map<std::string, std::vector<float>> startPositionsMap_singleModel;
        std::map<std::string, std::vector<float>> dircetionsMap_singleModel;
        std::map<std::string, std::vector<float>> lengthsMap_singleModel;
        std::map<std::string, std::vector<float>> posArraysMap_singleModel;

        // functions
        // access the data
        int getPosCnt(std::string forceName, megamol::core::Call* lmc);
        int getForceCnt(std::string forceName, megamol::core::Call* lmc) { return getPosCnt(forceName, lmc) / 2; };
        float* getCurPosArray(std::string forceName, megamol::core::Call* lmc);
        float* getCurStartPosArray(std::string forceName, megamol::core::Call* lmc);
        float* getCurDirArray(std::string forceName, megamol::core::Call* lmc);
        float* getCurLengthArray(std::string forceName, megamol::core::Call* lmc);

        uint* getCurForceCntPerProtAtomArray(std::string forceName, int protAtomCnt, megamol::core::Call* lmc);

        // get the last added force name
        std::string getCurForce() { return this->curForce; }

        // empty the map
        bool resetAll() {
            posArraysMap.clear();
            startPositionsMap.clear();
            dircetionsMap.clear();
            lengthsMap.clear();
            forceCntPerProtAtomMap.clear();
            return true;
        }

    private:
        bool checkDataExists(std::string forceName, megamol::core::Call* lmc, int protAtomCnt = 1);
    };
    ForcesLineStorage forcesLineStorage;


    class SphereRendererVars {
    public:
        SphereRendererVars() {
            this->hoverSphereID = -1;
            this->hideOnHover = true;
            this->selectedSphereID = -1;
            this->old_selectedSphereID = -1;
            this->hideOnSelection = false;
            this->hoverColor = {0.99f, 0, 0, 1};
            this->selectionColor = {1, 0.6f, 0, 1};
            this->baseColor = {0, 1, 1, 1};
        };
        ~SphereRendererVars(){};

        int hoverSphereID;
        bool hideOnHover;
        int selectedSphereID;
        int old_selectedSphereID;
        bool hideOnSelection;
        vislib::math::Vector<float, 4> hoverColor;
        vislib::math::Vector<float, 4> selectionColor;
        vislib::math::Vector<float, 4> baseColor;

        const GLuint sphereValues_SSBObindingPoint = 1;
        GLuint gl_sphereValuesSSBO = 0;
    };

    // CNS = cluster noise spheres
    SphereRendererVars renVarsCNS;
    // FGS = functional group spheres
    SphereRendererVars renVarsFGS;
    // ACAS = all cluster atoms
    SphereRendererVars renVarsACAS;
    // SMS = star menu circles
    SphereRendererVars renVarsSMS;
    // IFC = interaction force cylinders
    SphereRendererVars renVarsIFC;


    /******************************************
    ******** FUNCTINAL GROUP CLUSTERS *********
    *******************************************/
    // param variables for clustering
    int minPTS_fgc;
    float searchRadius_fgc;
    bool fgcClusterParamsChanged = true;

    LigandModelCall::ClusteredFunctionalGroups* clusteredFunctionalGroups;

    // new
    megamol::prolint::FGS_Handler fgsHandler;
    // complete fgsHierarchy
    std::vector<LigandModelCall::FGS_structN> fgsPerClusterComplete;
    // filterd regarding its fgsHierarchy (only roots of fgs)
    std::vector<LigandModelCall::FGS_structN> fgsPerClusterFiltered;

    // Render Functions
    bool proteinMeshRenderer(core::view::AbstractCallRender3D* cr3d, megamol::geocalls::CallTriMeshData* tmc,
        megamol::core::Call& call, core::view::CallGetTransferFunction* gtf);

    bool CircleRenderer(vislib::Array<float>& SphereCirles, vislib::Array<float>& colors, float* viewportStuff);

    bool sphereRenderer(float* spheres, int sphereCnt, core::view::CallGetTransferFunction* gtf,
        ModelClusterRenderer::SphereRendererVars& renVars, float* viewportStuff);

    bool circleTextureRenderer(vislib::Array<float>& spheres, core::view::CallGetTransferFunction* gtf,
        ModelClusterRenderer::SphereRendererVars& renVars, float* viewportStuff);

    bool coneRenderer(vislib::Array<float>& positions, vislib::Array<float>& directions,
        vislib::Array<float>& coneLengths, vislib::Array<float>& coneRadii, glm::vec4 color,
        ModelClusterRenderer::SphereRendererVars& renVars, float* viewportStuff);

    bool cylinderRenderer(float* positions, float* directions, float* coneLengths, int cylinderCnt,
        ModelClusterRenderer::SphereRendererVars& renVars, float* viewportStuff);

    bool dodecahedronRenderer(vislib::Array<float>& spheres, vislib::Array<float>* vertices,
        vislib::Array<float>* normals, core::view::CallGetTransferFunction* gtf, SphereRendererVars& renVars,
        float* viewportStuff);

    bool barChartRenderer(std::vector<std::map<string, float>>& clusterForceAreas,
        std::vector<std::map<string, std::map<uint, uint>>>& clusterForceResidueCnts, std::vector<float>& clusterAreas,
        megamol::core::Call* lmc, int ID, float liftSize);


    vislib::Array<float> clusterDodeca_verts;
    vislib::Array<float> clusterDodeca_normals;

    vislib::Array<float> funcGroupDodeca_verts;
    vislib::Array<float> funcGroupDodeca_normals;


    /********************************************
    ********* INTERACTION AREA TRIANGLES ********
    *********************************************/
    /** getInteractionAreaTriangles()
     *
     * interactionAreaAtoms:		stores all atomIDs which belongs to found pockets
     * interactionAreaAtomsClusterAssignment:	stores the assignment of which atomID bleogns to which cluster
     * atomToVerticesList:			2D array for assignment of atomIDs to all verticeIDs of an atom
     * (atomToVerticesList[atomID][i] = vertexID)
     * markedVerticesIndices:		stores all vertices which belogns to found clusters
     * interactionAreaVertices:		stores all vertices from the found clusters where a ligandAtom to protein
     *surface contact was calculated
     * old_MeshSelectedIDx:			check for change of selected cluster (performance: control recalculation)
     * old_checkForDataChange:		-||-
     * old_pmdcDataHash:			-||-
     */


    // functions
    bool getInteractionAreaTriangles();

    bool markForceVertices_getBarchatDat(
        vislib::Array<vislib::Array<uint>>& atomToVerticesList, int vertexCount, int protAtomCnt);

    // variables
    vislib::Array<vislib::Array<uint>> atomToVerticesList;
    vislib::Array<vislib::Array<uint>> atomToTriangleList;
    std::map<std::string, std::vector<uint>> markedVerticesIndicesMap;
    std::map<std::string, uint> maxForceCntPerProtAtomMap;
    vislib::Array<vislib::Array<uint>> clusterAreaVertices;

    // assigns the surface area to the atoms (which they contibute to the whole protein-SES)
    vislib::Array<float> atomToArea;
    // assigns the clusterID to the atoms
    std::vector<uint> clusterAreaAtomsToClusterIDs;
    // stores the surface-area / pocket-area for the clsuters
    std::vector<float> clusterAreas;

    // force specific variables/maps
   
	// stores of each pocket/cluster: the area of a certain interaction force which take part in one specific interaction force (only once)
    std::vector<std::map<string, float>> clusterForceAreas;
    
	// stores of each pocket/cluster: the count (size of the map) of individual residues which take part in one specific interaction force (only once)
    std::vector<std::map<string, std::map<uint, uint>>> clusterForceResidueCnts;
    std::map<std::string, glm::vec4> forceColorsMap;

    std::map<std::string, float> GUI;

    int vertexCount;

    vislib::Array<float> atomSphere;

    bool webDataChanged = false;
    int old_getWebData_dataHash = -1;
    int old_clusterDataHash;
    int old_pmdc_DataHash;
    int old_functionalGroup_dataHash;


    // h-bond line and cone rendering
    bool prepareHBond_Cones();

    vislib::Array<float> protHBonds_HDirection;
    vislib::Array<float> protHBonds_HPos;
    vislib::Array<float> protHBonds_coneLengths;
    vislib::Array<float> protHBonds_coneRadii;

    vislib::Array<float> modelHBonds_HDirection;
    vislib::Array<float> modelHBonds_HPos;
    vislib::Array<float> modelHBonds_coneLengths;
    vislib::Array<float> modelHBonds_coneRadii;


    /***********************************************************************
     ********* PMDC MolecularDataCall for Protein Residues Model(s) ********
     ***********************************************************************/
    // functions
    bool getDataPMDC(core::Call& call);
    bool getExtentPMDC(core::Call& call);

    bool setInteractionResiduesPMDC(core::Call& call);

    bool getInteactionResidues_FuncGroups(LigandModelCall::FGS_structN& functGH_struct,
        vislib::Array<int>& involvedResiduesIndices, protein_calls::MolecularDataCall* protMdcIn, int funcGrpID);

    // variables

    // protein residues and labels
    vislib::Array<float> atomPostions;
    vislib::Array<uint> atomTypeIdx;
    vislib::Array<float> atomCharges;
    vislib::Array<unsigned int> connections;
    vislib::Array<float> atomRadii;

    vislib::Array<vislib::StringA> residueLables;
    vislib::Array<vislib::math::Vector<float, 3>> residueLablesPositions;
    std::vector<std::vector<int>> proteinPocketAtoms;

    // controls if new data has to be set for MolecularDataCall -> out
    bool pmdcOut_newData = true;
    float searchRadius_InteractionResidues;


    /*************************************************************
     ********* LMDC MolecualrDataCall for Ligand Model(s) ********
     *************************************************************/
    // functions
    bool getDataLMDC(core::Call& call);
    bool getExtentLMDC(core::Call& call);
    bool setMultipleLigandModelsLMDC(core::Call& call);
    bool loadTexture(std::string filename);

    vislib::Array<float> lmdc_charges;
    vislib::Array<uint> lmdc_connections;
    vislib::Array<float> lmdc_positions;
    vislib::Array<uint> lmdc_atomTypeIdx;
    vislib::Array<protein_calls::MolecularDataCall::Chain> lmdc_chains;
    vislib::Array<megamol::protein_calls::MolecularDataCall::Residue> lmdc_residuesStore;
    vislib::Array<megamol::protein_calls::MolecularDataCall::Residue*> lmdc_residues;
    vislib::Array<int> lmdc_residuesIndices;
    vislib::Array<protein_calls::MolecularDataCall::Molecule> lmdc_moleclues;
    vislib::Array<vislib::SmartPtr<vislib::graphics::gl::OpenGLTexture2D>> modelModeIcons;
    vislib::graphics::BitmapImage img;
    GLuint bindID;

    bool lmdcOut_newData = true;

    class ModelRenderModes {
    public:
        ModelRenderModes() {
            this->curMode = 0;
            this->single = 0;
            this->allInCurPocket = 1;
            this->allOnProtein = 2;
        }
        uint curMode;

        uint single;
        uint allInCurPocket;
        uint allOnProtein;
    };
    ModelRenderModes modelRenderModes;


    // variables

    /****************************
    ********* CLIP PLANE ********
    *****************************/

    // variables
    vislib::math::Vector<double, 4> clipDat;


    // functions
    void getClipData(vislib::math::Vector<double, 4>* clipDat, float* clipCol, float* clipPoint);

    void planeRenderer(float* clipPoint, float* clipDat, float* clipCol);


    // control protein mesh coloring
    class ColoringControl {
    public:
        int interactionForce;
        int interactionForceAmountColoring;
        int pocketArea;
        int opaquePocket;
        int transparency;
        int chemicalProperty;
        int surfaceColoring_dropdown;
    };
    ColoringControl colorControl;


    // control render data
    class RenderDataControl {
    public:
        int clusterNoise;
        int clusterSpheres;
        int interactionResidues;
        int residues;
        int residuesBackbone;
        int residueLabels;
        int clusterAtoms;
        int funcGroupCentroids;
        int improveTextReadability;
        int barchartGylphe;
        class InteractionForces {
        public:
            int hydrophobicInteractions;
            int saltBridges;
            int piStacks;
            int piCationInteractions;
            int halogenBonds;
            int metalComplexes;
            int hydrogenBonds;
            int hydrogenCones;
        };
        InteractionForces intForces;
    };
    RenderDataControl renderDataControl;


    // control hoverings and selections
    /** infos:
     *
     * selectedIDx: IDx of clicked/selected sphere (combined_noiseModelCentroids_clusterCentroids[IDx])
     * hoveringIDx:	IDx of current hovered sphere (combined_noiseModelCentroids_clusterCentroids[IDx])
     */
    class SelectionControl {
    public:
        int selectedIDx;
        int old_SelectedIDx;

        int hoveringIDx;
        int old_HoveringIDx;

        int selectedGlobalMdlID;
        int old_selectedGlobalMdlID;

        int curClusterID;
        int old_curClusterID;
        bool changedcurClusterID;

        int selectedFgsClusterID;

        std::string curForceType_surfaceColoring;

        bool changedCurForceType;
    };
    SelectionControl s;
    vislib::Array<float> fgsClustSizes;

    // transferFunctions
    class TransferFunctions {
    public:
        core::view::CallGetTransferFunction* tf_interactionForce;
        core::view::CallGetTransferFunction* tf_clusterSpheres;
        core::view::CallGetTransferFunction* tf_functGroups;
    };
    TransferFunctions tf;


    // text colors
    glm::vec4 fontOutline_color;
    glm::vec4 fontStarMenu_color;
    glm::vec4 fontResidueLables_color;
    glm::vec4 fontFGSclustLables_color;


    /***********************************/

    /** camera information */
    vislib::SmartPtr<vislib::graphics::CameraParameters> cameraInfo;

    /** callee slots */
    //** MDC protein residues Callee Slot */
    core::CalleeSlot dataOutSlot_MDC_protein;

    //** MDC ligand poses Callee Slot */
    core::CalleeSlot dataOutSlot_MDC_ligand;

    //** Mouse Callee Slot */
    core::CalleeSlot dataSlot_MousePos;

    /** caller slots */
    megamol::core::CallerSlot ligandMolecularData_callerSlot;
    megamol::core::CallerSlot ligandModel_callerSlot;
    megamol::core::CallerSlot triMeshData_callerSlot;
    megamol::core::CallerSlot proteinMolecularData_callerSlot;
    megamol::core::CallerSlot volumetricData_callerSlot;
    megamol::core::CallerSlot transferFunction_interactionForce_callerSlot;
    megamol::core::CallerSlot transferFunction_clusterSpheres_callerSlot;
    megamol::core::CallerSlot transferFunction_functGroups_callerSlot;

    /** The call for clipping plane */
    megamol::core::CallerSlot getClipPlane_callerSlot;


    /** params slots **/
    // control render data
    core::param::ParamSlot clusterNoiseParam;
    core::param::ParamSlot clusterSpheresParam;
    core::param::ParamSlot hydrogenBondsParam;
    core::param::ParamSlot hydrogenConesParam;
    core::param::ParamSlot interactionResiduesParam;
    core::param::ParamSlot residuesParam;
    core::param::ParamSlot residuesBackboneParam;
    core::param::ParamSlot residuesLabelsParam;
    core::param::ParamSlot clusterAtomsParam;
    core::param::ParamSlot funcGroupCentroidsParam;
    core::param::ParamSlot improveTextReadabilityParam;
    core::param::ParamSlot barchartGylpheParam;

    // control render data => plipForces
    core::param::ParamSlot halogenBondsParam;
    core::param::ParamSlot hydrophobicInteractionsParam;
    core::param::ParamSlot metalComplexesParam;
    core::param::ParamSlot piCationInteractionsParam;
    core::param::ParamSlot piStacksParam;
    core::param::ParamSlot saltBridgesParam;

    // StarMenu
    // core::param::ParamSlot starMenuInnerStartParam;
    core::param::ParamSlot starMenuTextSizeParam;
    core::param::ParamSlot starMenuCntParam;
    core::param::ParamSlot starMenuSizeParam;
    core::param::ParamSlot subSphereSizefactorParam;
    core::param::ParamSlot starMenuSorting_dropdownParam;

    // control protein mesh coloring
    core::param::ParamSlot pocketAreaParam;
    core::param::ParamSlot opaquePocketParam;
    core::param::ParamSlot transparencyParam;
    core::param::ParamSlot chemicalPropertyParam;
    // control protein mesh coloring => inter force types
    // sc = surface coloring
    core::param::ParamSlot interactionForceParam;
    core::param::ParamSlot interactionForceAmountColoringParam;
    core::param::ParamSlot surfaceColoring_dropdownParam;


    // triangleMesh Rendering AO
    core::param::ParamSlot aoSampFactParam;
    core::param::ParamSlot aoSampDistParam;
    core::param::ParamSlot aoIntensityMinParam;
    core::param::ParamSlot aoIntensityMaxParam;

    // control functional group clustering
    core::param::ParamSlot minPTS_fgcParam;
    core::param::ParamSlot searchRadius_fgcParam;

    // extern param slots (from other modules)
    // MSMS mesh loader :: Module
    vislib::StringA paramName_proteinColoring0;
    vislib::StringA paramName_proteinColoring1;
    core::param::ParamSlot* proteinColoringMode0_Param;
    core::param::ParamSlot* proteinColoringMode1_Param;
    // MultiPDBQTLoader :: Module
    vislib::StringA paramName_modelClustering_eps;
    vislib::StringA paramName_modelClustering_minBindEnergy;
    vislib::StringA paramName_modelClustering_minPts;
    core::param::ParamSlot* modelClustering_eps_Param;
    core::param::ParamSlot* modelClustering_minBindEnergy_Param;
    core::param::ParamSlot* modelClustering_minPts_Param;
    // ClipPlane :: Module
    vislib::StringA paramName_clipPlaneEnable;
    vislib::StringA paramName_clipPlaneNormal;
    vislib::StringA paramName_clipPlaneDist;
    core::param::ParamSlot* clipPlaneEnable_Param;
    core::param::ParamSlot* clipPlaneNormal_Param;
    core::param::ParamSlot* clipPlaneDist_Param;
    std::map<string, glm::vec3> planeOrientations;
    int old_clipPlaneOrientation;
    float old_clipPlaneDist;
    float old_clipPlane_distVal;
    // protein cartoon renderer (Mux4Renderer3D :: Module)
    vislib::StringA paramName_cartoonRnederer;
    core::param::ParamSlot* cartoonRenderer_Param;

   

    // start without server
    core::param::ParamSlot withoutSocketParam;
    bool withoutSocket;


    std::map<string, uint> subCircleContentIDxes;

    // StarMenu
    /**variables */
    class StarMenu {
    public:
        StarMenu(){ 
			this->changed_curSortingFroceType = false;
		};

        // the number of current main circles
        int circleCnt;
        // the number of current main circles
        int curDrawCnt;
        // the old number of current main circles
        int old_circleCnt;
        // the size of the star menu on screen
        int size;
        // max number of current main circles which is allowed
        int maxCircleCnt;
        // the main circleID which is currently hovered
        int hoverIDx;
        // the main circleID which was clicked/selected
        int selectedIDx;
        // the old main cirlcleID which was clicked/selected before
        int old_selectedIDx;
        // the text size of the star menu
        float textSize;
        // the font color of the star menu
        glm::vec4 fontColor;
        // starMenuSorting_dropdown
        std::string curSortingFroceType;
        // changed cur force Type
        bool changed_curSortingFroceType;
        // int innerStart;

        std::vector<std::vector<vislib::StringA>> curSubCircleContent;

        // contains the current IDx of the cluster where the star menu gets the data for its first value from
        int curPageStartIDx;
        int curPage;
        int old_curPage;
    };
    StarMenu starMenu;
    megamol::core::utility::SDFFont font;


    /**********************************************************************
     * 'other'-functions
     **********************************************************************/
    bool getMousePos(core::Call& call);

    bool updateBoolParams(
        core::param::ParamSlot* param, std::string paramName, LigandModelCall::WebData* webData, int* renderBoo);
    bool updateIntParams(
        core::param::ParamSlot* param, std::string paramName, LigandModelCall::WebData* webData, int* renderInt);
    bool updateEnumParams(
        core::param::ParamSlot* param, std::string paramName, LigandModelCall::WebData* webData, int* renderInt);

    vislib::Array<vislib::math::Vector<float, 3>> oneClusterData_sortedWithEnergy;
    vislib::Array<float> starCircles;
    std::vector<float*> subCirclePointers;
    std::vector<float*> subCircleColorPointers;
    // the order numnber
    vislib::Array<float> starSubSpheres_1_orderNumber;
    // the score
    vislib::Array<float> starSubSpheres_2_dockScore;
    // subCircles H-bonds
    vislib::Array<float> starSubSpheres_3_Hbond;
    // subCircles halogen bond
    vislib::Array<float> starSubSpheres_4_halogenBond;
    // subCircles hydrophobicInters
    vislib::Array<float> starSubSpheres_5_hydrophobicInters;
    // subCircles for metalComplexes
    vislib::Array<float> starSubSpheres_6_metalComplexes;
    // subCircles for piCationInters
    vislib::Array<float> starSubSpheres_7_piCationInters;
    // subCircles for piSticks
    vislib::Array<float> starSubSpheres_8_piSticks;
    // subCircles for saltBridges
    vislib::Array<float> starSubSpheres_9_saltBridges;
    // the lignad model control (ball-stick rendering)
    vislib::Array<float> starSubSpheres_onSelect;

    int subCircleCnt;
    float subCircleSizefactor;
    vislib::Array<float> arrowSpheres;
    vislib::Array<float> textSphere;


    int curSelectedStarMenuModelIDx;
    int arrowUpDown;
    int indexOfLastStarCircle;
    vislib::Array<float> allClusterAtoms;

    /** Variables */

    /** Cluster Data
     * noiseModelCentroidsIDXs: stores indices of models which not belongs to any cluster
     * combined_noiseModelCentroids_clusterCentroids:	 final render data (name = description)
     */

    LigandModelCall::ClusterAndCentroidData clusterData;
    vislib::Array<float> tf_clusterSizes;
    vislib::Array<int> noiseModelCentroidsIDXs;
    vislib::Array<float> combined_noiseModelCentroids_clusterCentroids;


    /** filled with MouseInputCall
     *
     * mouseX, mouseY:	mouse position (screen coordinates)
     * cam:				vislib::graphics::gl::CameraOpenGL object
     * mouseFlag:		(e.g button x,y,z is pressed)
     * isClicked:		when the mousbutton releases it is set to 1
     */
    float mouseX, mouseY;
    vislib::graphics::gl::CameraOpenGL cam;
    int mouseFlag;
    int isClicked;


    vislib::Array<float> sphereColorsScales;
    GLuint gl_sphereColorsScales = 0;
    GLuint gl_funczGroupColorsScales = 0;

    float noiseModel_CentroidRadius;
    float cluster_CentroidRadius;
    float higlight_noiseModel_CentroidRadius;
    float higlight_cluster_CentroidRadius;

    /** shader for triangle meshes */
    vislib::graphics::gl::GLSLShader triangleMeshShader;

    /** shader for the circles (sphereShader without light - raycasting view) */
    vislib::graphics::gl::GLSLGeometryShader circleShader;

    /** shader for the spheres (raycasting view) */
    vislib::graphics::gl::GLSLGeometryShader sphereShader;

    /** shader for the textured cricles (raycasting view) */
    vislib::graphics::gl::GLSLGeometryShader circleBillTextureShader;

    /** shader for the dodecahedron */
    vislib::graphics::gl::GLSLGeometryShader dodecahedron_Shader;

    /** shader for the simpleAmbient */
    vislib::graphics::gl::GLSLGeometryShader simpleAmbientShader;

    /** shader for the cones */
    vislib::graphics::gl::GLSLGeometryShader coneShader;

    /** shader for the cones */
    vislib::graphics::gl::GLSLGeometryShader cylinderShader;

    /*********************************************
     ************** SVG XML PARSING **************
     *********************************************/

    //********* SVG DATA STURCTURE/CLASS ***********//
    class Svg_line1 {
    public:
        float x1;
        float y1;
        float x2;
        float y2;
        float opacity;
        vislib::math::Vector<int, 3> color;
        float strokeWidth;
    };

    class Svg_ploygon1 {
    public:
        vislib::Array<vislib::math::Vector<float, 2>> points;
        float strokeWidth;
        vislib::math::Vector<int, 3> color;
    };

    class Svg_text1 {
    public:
        float x;
        float y;
        vislib::math::Vector<float, 3> color;
        float strokeWidth;
        std::string fontWeight;
        float fontSize;
        std::string text;
    };

    class SVG1 {
    public:
        vislib::math::Vector<float, 2> svgCentroid;
        vislib::math::Vector<float, 2> svgMAX_xy;
        vislib::math::Vector<float, 2> svgMIN_xy;
        vislib::math::Vector<float, 2> svgMiddle;
        std::vector<Svg_line1> svg_lineElements;
        std::vector<Svg_text1> svg_textElements;
        std::vector<Svg_ploygon1> svg_polygons;
        int lineCnt = 0;
        int polygonCnt = 0;
        int textCnt = 0;
    };
    std::map<int, int> int2String;
    std::vector<SVG1> svg1;
    int SVGcnt;

    bool firstRenderCall;
    std::vector<std::string> SVGPaths;

    //********* SVG XML PARSING FUNCTIONS ***********//

    float getElementValue_svgLine(std::string line, std::string elementString);

    vislib::math::Vector<int, 3> getElementRGB_svgLine(std::string line, std::string elementString);

    std::string getElementText_svgLine(std::string line, std::string elementString);

    void getPoints_svgPolygon(std::string line, vislib::Array<vislib::math::Vector<float, 2>>& points);

    void parseSVG(std::string filepath, class SVG1* svg1);
};


} // namespace prolint
} /* end namespace megamol */

int sortVec2ByX(const vislib::math::Vector<int, 2U>& lhs, const vislib::math::Vector<int, 2U>& rhs);
int sortVecf3ByY(const vislib::math::Vector<float, 3U>& lhs, const vislib::math::Vector<float, 3U>& rhs);
int sortVecf3ByYifZisNoneZero(const vislib::math::Vector<float, 3U>& lhs, const vislib::math::Vector<float, 3U>& rhs);

#endif // PROLINT_PLUGIN_ModelClusterRenderer_H_INCLUDED
