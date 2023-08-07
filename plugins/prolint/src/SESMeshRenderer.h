/*
 * SESMeshRenderer.h
 *
 * Copyright (C) 2019 by Universitaet Tuebingen (BDVA).
 * All rights reserved.
 */

#ifndef MMPRLINTPLUGIN_SESMESHREN_H_INCLUDED
#define MMPRLINTPLUGIN_SESMESHREN_H_INCLUDED
#if (defined(_MSC_VER) && (_MSC_VER > 1000))
#    pragma once
#endif /* (defined(_MSC_VER) && (_MSC_VER > 1000)) */

#include "mmcore/CallerSlot.h"
#include "mmcore/param/ParamSlot.h"
#include "mmcore/view/CallRender3D_2.h"
#include "mmcore/view/Renderer3DModule_2.h"
#include "protein/Color.h"
#include "protein_calls/BindingSiteCall.h"
#include "protein_calls/MolecularDataCall.h"

#include "vislib/graphics/gl/GLSLGeometryShader.h"
#include "vislib/graphics/gl/GLSLShader.h"
#include "vislib/graphics/gl/GLSLTesselationShader.h"

#include <math.h>
#include "Icosphere.h"
#include "SESTrilaterationCUDA.cuh"
#include "glm/gtc/type_ptr.hpp"


namespace megamol {
namespace prolint {


/*
 * Test Renderer class
 */

class SESMeshRenderer : public megamol::core::view::Renderer3DModule_2 {
public:
    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName(void) { return "SESMeshRenderer"; }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description(void) { return "Offers SES renderings for molecular data."; }

    /**
     * Answers whether this module is available on the current system.
     *
     * @return 'true' if the module is available, 'false' otherwise.
     */
    static bool IsAvailable(void) { return true; }

    /** Ctor. */
    SESMeshRenderer(void);

    /** Dtor. */
    virtual ~SESMeshRenderer(void);

protected:
    /**
     * Implementation of 'Create'.
     *
     * @return 'true' on success, 'false' otherwise.
     */
    virtual bool create(void);

    /**
     * Implementation of 'release'.
     */
    virtual void release(void);

private:
    /**********************************************************************
     * 'render'-functions
     **********************************************************************/

    /**
     * The get capabilities callback. The module should set the members
     * of 'call' to tell the caller its capabilities.
     *
     * @param call The calling call.
     *
     * @return The return value of the function.
     */
    // virtual bool GetCapabilities( megamol::core::Call& call);

    /**
     * The get extents callback. The module should set the members of
     * 'call' to tell the caller the extents of its data (bounding boxes
     * and times).
     *
     * @param call The calling call.
     *
     * @return The return value of the function.
     */
    virtual bool GetExtents(core::view::CallRender3D_2& call);


    /**
     * The Open GL Render callback.
     *
     * @param call The calling call.
     * @return The return value of the function.
     */
    virtual bool Render(core::view::CallRender3D_2& call);


    /**
     * Render the molecular data in SES Trilateration mode.
     *
     * @param mol        Pointer to the data call.
     * @param atomPos    Pointer to the interpolated atom positions.
     */
    void RenderSES(megamol::protein_calls::MolecularDataCall* mol, const float* atomPos, core::view::Camera_2& cam);

    /**
     * Update all parameter slots.
     *
     * @param mol   Pointer to the data call.
     */
    void UpdateParameters(
        const megamol::protein_calls::MolecularDataCall* mol, const protein_calls::BindingSiteCall* bs = 0);


    /**********************************************************************
     * variables
     **********************************************************************/

    /** MolecularDataCall caller slot */
    megamol::core::CallerSlot molDataCallerSlot;
    /** BindingSiteCall caller slot */
    megamol::core::CallerSlot bsDataCallerSlot;

    /** parameter slot for SAS probe radius */
    megamol::core::param::ParamSlot probeRadiusParam;
    /** parameter slot for color table filename */
    megamol::core::param::ParamSlot colorTableFileParam;
    /** parameter slot for coloring mode */
    megamol::core::param::ParamSlot coloringModeParam0;
    /** parameter slot for coloring mode */
    megamol::core::param::ParamSlot coloringModeParam1;
    /** parameter slot for coloring mode weighting*/
    megamol::core::param::ParamSlot cmWeightParam;
    /** parameter slot for min color of gradient color mode */
    megamol::core::param::ParamSlot minGradColorParam;
    /** parameter slot for mid color of gradient color mode */
    megamol::core::param::ParamSlot midGradColorParam;
    /** parameter slot for max color of gradient color mode */
    megamol::core::param::ParamSlot maxGradColorParam;
    /** list of molecule indices */
    megamol::core::param::ParamSlot molIdxListParam;
    /** parameter slot for special color */
    megamol::core::param::ParamSlot specialColorParam;

    /** shader for the spheres (raycasting view) */
    vislib::graphics::gl::GLSLShader sphereShader;
    /** shader for the spherical triangles (raycasting view) */
    vislib::graphics::gl::GLSLShader sphericalTriangleShader;
    /** shader for torus (raycasting view) */
    vislib::graphics::gl::GLSLShader torusShader;
    /** shader for triangles */
    vislib::graphics::gl::GLSLShader triangleShader;

    /** The atom color table for rendering */
    vislib::Array<float> atomColorTable;
    vislib::Array<float> atomColorTableRGB;
    bool uploadColors;

    /** vertex array for spheres */
    vislib::Array<float> vertSpheres;

    // the list of molecular indices
    vislib::Array<vislib::StringA> molIdxList;

    /** The hash of the lastly rendered molecular data call*/
    float lastTimestamp;

    /** defintions for neighbour-search and 3 spheres intersections*/
    unsigned int gridsize;  // number of cells that was used for last memory allocation
    unsigned int atomCnt;   // number of atoms that was used for last memory allocation
    float3 move_origin_pdb; // x,y,z for moving positions out of negative values
    float4* d_atomPos;      // array of atom positions x,y,z * atom_count (device)
    int2* d_insert_grid;    // index = atomID PDB; x=gridcell, y=intern index of atoms in this cell
    int* d_grid;            // index = gridcell; values = #atoms in grid cell
    int2* d_sorted; // x=atomID PDB; y=gridcell => sorted  d_grid (index = atomID PDB; value=gridcell) now sorted after
                    // gridcell little > big
    int2* d_cell_start_end; // index+1=gridcell; x= start index (from d_sorted) of gridcell; y=end index (from d_sorted)
                            // of gridcell => range of a gridcell cell in d_sorted (50=>7;50=>12;50=>9) [50 x=0 y=2]
    unsigned int* d_neighbours; // 1D array of neighbours; size 100 * atom_count; all neighbours of atomID=1 in 0-99;
                                // atomID=2 100-199...
    unsigned int* d_neighbourCounts; // index = atomID; value=number of neighbours; size=atom_count (device)
    float probeRad;
    float grid_step_size;

    float3* d_sphereVertices;
    size_t sphereVerticesCnt;
    uint3* d_sphereTriangles;
    size_t sphereTriangleCnt;
    uint3* d_outerSphereTriangles;
    uint3* d_cutSphereTriangles;
    size_t totalSphereTriangleCnt;

    // OpenGL buffers (will be filled by CUDA)

    // Basic SSBOs for atomPos and color (for spheres renderer)
    GLuint gl_atomPos = 0;
    cudaGraphicsResource* gl_atomPosResource;
    size_t gl_atomPosSize;

    GLuint gl_atomPosIdx = 0;
    cudaGraphicsResource* gl_atomPosIdxResource;
    size_t gl_atomPosIdxSize;
    int* d_atomPosIdx;
    unsigned int size_AtomPosIdx;

    GLuint gl_atomColorTable = 0;

    /** The current coloring mode */
    protein::Color::ColoringMode currentColoringMode0;
    protein::Color::ColoringMode currentColoringMode1;

    /** The color lookup table (for chains, amino acids,...) */
    vislib::Array<vislib::math::Vector<float, 3>> colorLookupTable;
    /** The color lookup table which stores the rainbow colors */
    vislib::Array<vislib::math::Vector<float, 3>> rainbowColors;

    // isosphere/geosphere
    Icosphere* sphere;

    float* pos0;
    float* pos1;
    size_t posArraySize;
    vislib::Array<float> posInter;
};


} /* end namespace prolint */
} /* end namespace megamol */

#endif // MMPRLINTPLUGIN_SESMESHREN_H_INCLUDED
