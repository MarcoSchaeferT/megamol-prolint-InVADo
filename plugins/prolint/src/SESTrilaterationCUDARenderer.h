/*
 * SESTrilaterationCUDARenderer.h
 *
 * Copyright (C) 2019 by Universitaet Tuebingen (BDVA).
 * All rights reserved.
 */

#ifndef MMPRLINTPLUGIN_SESTRILATCUDAREN_H_INCLUDED
#define MMPRLINTPLUGIN_SESTRILATCUDAREN_H_INCLUDED
#if (defined(_MSC_VER) && (_MSC_VER > 1000))
#    pragma once
#endif /* (defined(_MSC_VER) && (_MSC_VER > 1000)) */

#include "protein_calls/BindingSiteCall.h"
#include "protein_calls/MolecularDataCall.h"
#include "protein/Color.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/param/ParamSlot.h"
#include "mmcore/view/CallRender3D_2.h"
#include "mmcore/view/Renderer3DModule_2.h"


#include "vislib/graphics/gl/GLSLTesselationShader.h"
#include "vislib/graphics/gl/GLSLGeometryShader.h"
#include "vislib/graphics/gl/GLSLShader.h"

#include <math.h>
#include "glm/gtc/type_ptr.hpp"
#include "SESTrilaterationCUDA.cuh"


namespace megamol {
namespace prolint {


/*
 * Test Renderer class
 */

class SESTrilaterationCUDARenderer : public megamol::core::view::Renderer3DModule_2 {
public:
    /** The names of the render modes */
    enum RenderMode {
        LINES = 0,
        STICK = 1,
        BALL_AND_STICK = 2,
        SPACEFILLING = 3,
        SAS = 4,
        LINES_FILTER = 5,
        STICK_FILTER = 6,
        SPACEFILL_FILTER = 7,
		SES_TRILATERATION = 8
    };

    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName(void) { return "SESTrilaterationCUDARenderer"; }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description(void) { return "Offers molecule renderings."; }

    /**
     * Answers whether this module is available on the current system.
     *
     * @return 'true' if the module is available, 'false' otherwise.
     */
    static bool IsAvailable(void) { return true; }

    /** Ctor. */
    SESTrilaterationCUDARenderer(void);

    /** Dtor. */
    virtual ~SESTrilaterationCUDARenderer(void);

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
	void RenderSES_Trilateration(megamol::protein_calls::MolecularDataCall* mol, const float* atomPos, core::view::Camera_2 &cam);

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

    /** parameter slot for rendering mode */
    megamol::core::param::ParamSlot renderModeParam;
    /** parameter slot for stick radius */
    megamol::core::param::ParamSlot stickRadiusParam;
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
    /** DEBUG */
    //megamol::core::param::ParamSlot debugNeighbors;

    /** shader for the spheres (raycasting view) */
    vislib::graphics::gl::GLSLShader sphereShader;
    /** shader for the spherical triangles (raycasting view) */
    vislib::graphics::gl::GLSLShader sphericalTriangleShader;
    /** shader for torus (raycasting view) */
    vislib::graphics::gl::GLSLShader torusShader;
    /** shader for the cylinders (raycasting view) */
    vislib::graphics::gl::GLSLShader cylinderShader;
	
    /** The atom color table for rendering */
    vislib::Array<float> atomColorTable;
    vislib::Array<float> atomColorTableRGB;
    bool uploadColors;

    /** The current rendering mode */
    RenderMode currentRenderMode;

    /** vertex array for spheres */
    vislib::Array<float> vertSpheres;

    // the list of molecular indices
    vislib::Array<vislib::StringA> molIdxList;

    /** The hash of the lastly rendered molecular data call*/
    float lastTimestamp;

    /** defintions for neighbour-search and 3 spheres intersections*/
    unsigned int gridsize;						// number of cells that was used for last memory allocation
    unsigned int atomCnt;						// number of atoms that was used for last memory allocation
    float3 move_origin_pdb;						// x,y,z for moving positions out of negative values
    float4* d_atomPos;							// array of atom positions x,y,z * atom_count (device)
    int2* d_insert_grid;						// index = atomID PDB; x=gridcell, y=intern index of atoms in this cell
    int* d_grid;								// index = gridcell; values = #atoms in grid cell
    int2* d_sorted;								// x=atomID PDB; y=gridcell => sorted  d_grid (index = atomID PDB; value=gridcell) now sorted after gridcell little > big
	int2* d_cell_start_end;						// index+1=gridcell; x= start index (from d_sorted) of gridcell; y=end index (from d_sorted) of gridcell => range of a gridcell cell in d_sorted (50=>7;50=>12;50=>9) [50 x=0 y=2]
    unsigned int* d_neighbours;					// 1D array of neighbours; size 100 * atom_count; all neighbours of atomID=1 in 0-99; atomID=2 100-199...
    unsigned int* d_reduced_neighbours;			// like d_neighbour but with only neighbours greater than the actual atom-index => used for calc of all unique 3er Combos
    unsigned int* d_neighbourCounts;			// index = atomID; value=number of neighbours; size=atom_count (device)
    unsigned int* d_reduced_neighbourCounts;	// like d_neighbourCounts but with only neighboursCounts of neighbours greater than the actual atom-index => used for calc of all unique 3er Combos
	int* d_neighbour_combination_numbers;		// stores the number of max 3er combos per atom
	unsigned int d_3atomComboSize;				// stortes the toatal number of all 3er combos
    unsigned int d_intersectionsSize;           // size that should theoretically be used for memory allocation
    unsigned int d_intersectionMemSize;         // actual size that was used for last memory allocation
    float4* d_intersections;					// stores all Intersections; index = S1 // index+1 = S2;  size=sum of all overlaps * 2
    int3* d_SpheresComb3;						// stores all 3Spheres combinations (atom1 atom2 atom3) where Intersections exists; size=sum of all combinations with overlaps => used for filtering
    int4* d_filtered_spheresComb3;				// stores all the 3Spheres combinations for which an intersection exists and for this intersection not intersects with other atoms
    int* d_filtered_dim;						// number of intersections which not lie in any other sphere (device)		#
    int h_filtered_dim;							// number of intersections which not lie in any other sphere (host) => used for pushing to sphere renderer
	int* d_filtered_prepos;						// index=one of all interections (s1), index+1=one of all interections (s2), value=is value of an incremented counter of intersections which not lie in a sphere (otherwise=0), size = sum of all overlaps * 2		#
    unsigned int d_filtered_IntersectionsSize;  // size that was used for last memory allocation
    int3* d_filtered_3spheresComb;				// stores the corresponding 3spheresComb to the intersections in d_filtered_Intersections  (device)		#
    int2* d_torusAxes;							// stores the torus axes (int2 data but int4 for getting SSBO correct into shader)
    unsigned int numTorusAxes;
    unsigned int run_number;
    int* d_atomicIdx;
    float probeRad; 
	float grid_step_size;

    float3* d_intersecting_intersections;		// stores all intersections which intersecting with each other
    int2 *d_intersection_insert_grid;			// index = intersection => x=gridcell y=intern index of the intersection in this gridcell
    int2* d_intersecting_sorted;				// x=intersection; y=gridcell => sorted  d_grid after gridcell little > big
    int2* d_intersecting_cell_start_end;		// index=gridcell; x= start index (from d_intersecting_sorted) of gridcell; y=end index (from d_intersecting_sorted) of gridcell 
    float4* d_intersecting_neighbours;	// 1D array of neighbours; size 50 * d_filtered_IntersectionsSize (stores probe positions)
    int* d_intersecting_neighbourCounts;		// index = intersectionID; value=number of neighbours;  (device)


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

	// SSBOs for spherical triangles (renderer) 
	GLuint gl_intersecting_neighbours = 1;
	cudaGraphicsResource *gl_intersecting_neighbours_resource;
	size_t gl_intersecting_neighbours_size;

	GLuint gl_filtered_spheresComb3 = 0;
	cudaGraphicsResource* gl_filtered_spheresComb3_resource;
	size_t gl_intersectionsSize;

	GLuint gl_intersections = 0;
	cudaGraphicsResource* gl_intersectionsResource;

	//SSBOs for torus parameter (renderer)
    GLuint gl_torusAxes = 0;
    cudaGraphicsResource* gl_torusAxesResource;
    size_t gl_torusAxesSize;

    /** The current coloring mode */
    protein::Color::ColoringMode currentColoringMode0;
    protein::Color::ColoringMode currentColoringMode1;

    /** The color lookup table (for chains, amino acids,...) */
    vislib::Array<vislib::math::Vector<float, 3>> colorLookupTable;
    /** The color lookup table which stores the rainbow colors */
    vislib::Array<vislib::math::Vector<float, 3>> rainbowColors;

    float* pos0;
    float* pos1;
    size_t posArraySize;
    vislib::Array<float> posInter;
};


} /* end namespace prolint */
} /* end namespace megamol */

#endif // MMPRLINTPLUGIN_SESTRILATCUDAREN_H_INCLUDED
