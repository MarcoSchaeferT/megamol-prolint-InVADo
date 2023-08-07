/*
 * SESMeshRenderer.cpp
 *
 * Copyright (C) 2019 by Universitaet Tuebingen (BDVA).
 * All rights reserved.
 */

#include "stdafx.h"

#define _USE_MATH_DEFINES 1
#define mem_manage 1
#define NUM_THREADS 256

#include "vislib/graphics/gl/IncludeAllGL.h"

#include <GL/glu.h>
#include <omp.h>
#include "SESMeshRenderer.h"
#include "mmcore/CoreInstance.h"
#include "mmcore/param/BoolParam.h"
#include "mmcore/param/EnumParam.h"
#include "mmcore/param/FloatParam.h"
#include "mmcore/param/IntParam.h"
#include "mmcore/param/StringParam.h"
#include "mmcore/utility/ColourParser.h"
#include "mmcore/utility/ShaderSourceFactory.h"
#include "vislib/OutOfRangeException.h"
#include "vislib/String.h"
#include "vislib/StringConverter.h"
#include "vislib/Trace.h"
#include "vislib/assert.h"
#include "vislib/graphics/gl/AbstractOpenGLShader.h"
#include "vislib/graphics/gl/ShaderSource.h"
#include "vislib/math/Matrix.h"
#include "vislib/math/Quaternion.h"
#include "mmcore/utility/sys/ASCIIFileBuffer.h"

#include <Windows.h>

//#define TEST
#ifdef TEST
#    include <fstream>
#    include <iostream>
#    include <time.h>
clock_t tStart = clock();
#endif

#include "mmcore/CoreInstance.h"
#include "vislib/graphics/gl/ShaderSource.h"

#include <cuda_gl_interop.h>


#define CHECK_FOR_OGL_ERROR()                                                                                          \
    do {                                                                                                               \
        GLenum err;                                                                                                    \
        err = glGetError();                                                                                            \
        if (err != GL_NO_ERROR) {                                                                                      \
            fprintf(stderr, "%s(%d) glError: %s\n", __FILE__, __LINE__, gluErrorString(err));                          \
        }                                                                                                              \
    } while (0)

using namespace megamol;
using namespace megamol::core;
using namespace megamol::prolint;
using namespace megamol::protein_calls;


const GLuint probePosSSBObindingPoint = 2;
const GLuint atomPosSSBObindingPoint = 3;
const GLuint atomColorTableSSBObindingPoint = 4;
const GLuint intersectionsSSBObindingPoint = 5;
const GLuint filtered_spheresComb3SSBObindingPoint = 6;
const GLuint torusAxesSSBObindingPoint = 7;
const GLuint atomPosIdxSSBObindingPoint = 8;


/*
 * prolint::SESMeshRenderer::SESMeshRenderer (CTOR)
 */
SESMeshRenderer::SESMeshRenderer(void)
    : Renderer3DModule_2()
    , molDataCallerSlot("getData", "Connects the molecule rendering with molecule data storage")
    , bsDataCallerSlot("getBindingSites", "Connects the molecule rendering with binding site data storage")
    , probeRadiusParam("probeRadius", "The probe radius for SAS rendering")
    , colorTableFileParam("color::colorTableFilename", "The filename of the color table.")
    , coloringModeParam0("color::coloringMode0", "The first coloring mode.")
    , coloringModeParam1("color::coloringMode1", "The second coloring mode.")
    , cmWeightParam("color::colorWeighting", "The weighting of the two coloring modes.")
    , minGradColorParam("color::minGradColor", "The color for the minimum value for gradient coloring")
    , midGradColorParam("color::midGradColor", "The color for the middle value for gradient coloring")
    , maxGradColorParam("color::maxGradColor", "The color for the maximum value for gradient coloring")
    , molIdxListParam("molIdxList", "The list of molecule indices for RS computation:")
    , specialColorParam("color::specialColor", "The color for the specified molecules")
    , gridsize(0)
    , atomCnt(0)
    , d_sphereVertices(0)
    , sphereVerticesCnt(0)
	, d_sphereTriangles(0)
	, sphereTriangleCnt(0)
    , totalSphereTriangleCnt(0)
    , probeRad(0)
    , lastTimestamp(-1.0f)
    , posArraySize(0)
    , uploadColors(true)
	, sphere(0) {
    this->molDataCallerSlot.SetCompatibleCall<MolecularDataCallDescription>();
    this->MakeSlotAvailable(&this->molDataCallerSlot);
    this->bsDataCallerSlot.SetCompatibleCall<BindingSiteCallDescription>();
    this->MakeSlotAvailable(&this->bsDataCallerSlot);

    // molecular indices list param
    this->molIdxList.Clear();
    this->molIdxListParam.SetParameter(new param::StringParam(""));
    this->MakeSlotAvailable(&this->molIdxListParam);

    // ----- color setup -----

    // fill color table with default values and set the filename param
    this->probeRadiusParam.SetParameter(new param::FloatParam(1.4f));
    this->MakeSlotAvailable(&this->probeRadiusParam);

    // fill color table with default values and set the filename param
    vislib::StringA filename("colors.txt");
    protein::Color::ReadColorTableFromFile(filename, this->colorLookupTable);
    this->colorTableFileParam.SetParameter(new param::StringParam(A2T(filename)));
    this->MakeSlotAvailable(&this->colorTableFileParam);

    // coloring modes
    this->currentColoringMode0 = protein::Color::CHAIN;
    this->currentColoringMode1 = protein::Color::ELEMENT;
    param::EnumParam* cm0 = new param::EnumParam(int(this->currentColoringMode0));
    param::EnumParam* cm1 = new param::EnumParam(int(this->currentColoringMode1));
    MolecularDataCall* mol = new MolecularDataCall();
    BindingSiteCall* bs = new BindingSiteCall();
    unsigned int cCnt;
    protein::Color::ColoringMode cMode;
    for (cCnt = 0; cCnt < protein::Color::GetNumOfColoringModes(mol, bs); ++cCnt) {
        cMode = protein::Color::GetModeByIndex(mol, bs, cCnt);
        cm0->SetTypePair(cMode, protein::Color::GetName(cMode).c_str());
        cm1->SetTypePair(cMode, protein::Color::GetName(cMode).c_str());
    }
    delete mol;
    delete bs;
    this->coloringModeParam0 << cm0;
    this->coloringModeParam1 << cm1;
    this->MakeSlotAvailable(&this->coloringModeParam0);
    this->MakeSlotAvailable(&this->coloringModeParam1);

    // Color weighting parameter
    this->cmWeightParam.SetParameter(new param::FloatParam(0.5f, 0.0f, 1.0f));
    this->MakeSlotAvailable(&this->cmWeightParam);

    // the color for the minimum value (gradient coloring
    this->minGradColorParam.SetParameter(new param::StringParam("#146496"));
    this->MakeSlotAvailable(&this->minGradColorParam);

    // the color for the middle value (gradient coloring
    this->midGradColorParam.SetParameter(new param::StringParam("#f0f0f0"));
    this->MakeSlotAvailable(&this->midGradColorParam);

    // the color for the maximum value (gradient coloring
    this->maxGradColorParam.SetParameter(new param::StringParam("#ae3b32"));
    this->MakeSlotAvailable(&this->maxGradColorParam);

    // make the rainbow color table
    protein::Color::MakeRainbowColorTable(100, this->rainbowColors);

    // the color for the maximum value (gradient coloring
    this->specialColorParam.SetParameter(new param::StringParam("#228B22"));
    this->MakeSlotAvailable(&this->specialColorParam);

}

/*
 * prolint::SESMeshRenderer::~SESMeshRenderer (DTOR)
 */
SESMeshRenderer::~SESMeshRenderer(void) { this->Release(); }

/*
 * prolint::SESMeshRenderer::released_intersecting_neighbours
 */
void SESMeshRenderer::release(void) {
    if (this->posArraySize > 0) {
        delete[] this->pos0;
        delete[] this->pos1;
    }
}

/*
 * prolint::SESMeshRenderer::create
 */
bool SESMeshRenderer::create(void) {
    if (!isExtAvailable("GL_ARB_vertex_program") || !ogl_IsVersionGEQ(2, 0)) return false;

    if (!vislib::graphics::gl::GLSLShader::InitialiseExtensions()) return false;
    if (!vislib::graphics::gl::GLSLTesselationShader::InitialiseExtensions()) return false;

    // make default sphere
    if (!this->sphere) {
        this->sphere = new Icosphere(2.0f, 1, true);
    }
    using megamol::core::utility::log::Log;
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glEnable(GL_VERTEX_PROGRAM_TWO_SIDE);
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_ARB);

    using namespace vislib::sys;
    using namespace vislib::graphics::gl;

    ShaderSource vertSrc;
    ShaderSource fragSrc;
    ShaderSource geomSrc;

    // Load sphere shader
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource("prolint-ses::ses::sphereVertex", vertSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to load vertex shader source for sphere shader");
        return false;
    }
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource("prolint-ses::ses::sphereFragment", fragSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to load vertex shader source for sphere shader");
        return false;
    }
    try {
        if (!this->sphereShader.Create(vertSrc.Code(), vertSrc.Count(), fragSrc.Code(), fragSrc.Count())) {
            throw vislib::Exception("Generic creation failure", __FILE__, __LINE__);
        }
    } catch (vislib::Exception e) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to create sphere shader: %s\n", e.GetMsgA());
        return false;
    }

    // Load spherical triangle shader
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource(
            "prolint-ses::ses::sphericaltriangleVertex", vertSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to load vertex shader source for spherical triangle shader");
        return false;
    }
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource(
            "prolint-ses::ses::sphericaltriangleFragment", fragSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to load vertex shader source for spherial triangle shader");
        return false;
    }
    try {
        if (!this->sphericalTriangleShader.Create(vertSrc.Code(), vertSrc.Count(), fragSrc.Code(), fragSrc.Count())) {
            throw vislib::Exception("Generic creation failure", __FILE__, __LINE__);
        }
    } catch (vislib::Exception e) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to create spherical triangle shader: %s\n", e.GetMsgA());
        return false;
    }

    // Load torus shader
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource("prolint-ses::ses::torusVertex", vertSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to load vertex shader source for torus shader");
        return false;
    }
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource("prolint-ses::ses::torusFragment", fragSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to load vertex shader source for torus shader");
        return false;
    }
    try {
        if (!this->torusShader.Create(vertSrc.Code(), vertSrc.Count(), fragSrc.Code(), fragSrc.Count())) {
            throw vislib::Exception("Generic creation failure", __FILE__, __LINE__);
        }
    } catch (vislib::Exception e) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to create torus shader: %s\n", e.GetMsgA());
        return false;
    }

    // Load triangle shader
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource("prolint-ses::triangle::vertex", vertSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to load vertex shader source for triangle rendering");
        return false;
    }
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource("prolint-ses::triangle::fragment", fragSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to load vertex shader source for triangle rendering");
        return false;
    }
    try {
        if (!this->triangleShader.Create(vertSrc.Code(), vertSrc.Count(), fragSrc.Code(), fragSrc.Count())) {
            throw vislib::Exception("Generic creation failure", __FILE__, __LINE__);
        }
    } catch (vislib::Exception e) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to create triangle rendering shader: %s\n", e.GetMsgA());
        return false;
    }

    return true;
}


/*
 * prolint::SESMeshRenderer::GetExtents
 */
bool SESMeshRenderer::GetExtents(core::view::CallRender3D_2& call) {
    MolecularDataCall* mol = this->molDataCallerSlot.CallAs<MolecularDataCall>();
    if (mol == NULL) return false;
    if (!(*mol)(MolecularDataCall::CallForGetExtent)) return false;

    call.AccessBoundingBoxes().SetBoundingBox(mol->AccessBoundingBoxes().ObjectSpaceBBox());
    call.AccessBoundingBoxes().SetClipBox(mol->AccessBoundingBoxes().ObjectSpaceClipBox());
    call.SetTimeFramesCount(mol->FrameCount());

    return true;
}

/**********************************************************************
 * 'render'-functions
 **********************************************************************/


/*
 * prolint::SESMeshRenderer::Render
 */
bool SESMeshRenderer::Render(core::view::CallRender3D_2& call) {

    view::Camera_2 cam;
    call.GetCamera(cam);

    float callTime = call.Time();

    // get pointer to MolecularDataCall
    MolecularDataCall* mol = this->molDataCallerSlot.CallAs<MolecularDataCall>();
    if (mol == NULL) return false;

    // get pointer to BindingSiteCall
    BindingSiteCall* bs = this->bsDataCallerSlot.CallAs<BindingSiteCall>();
    if (bs) {
        (*bs)(BindingSiteCall::CallForGetData);
    }

    // set call time
    mol->SetCalltime(callTime);
    // set frame ID and call data
    mol->SetFrameID(static_cast<int>(callTime));

    if (!(*mol)(MolecularDataCall::CallForGetData)) return false;
    // check if atom count is zero
    if (mol->AtomCount() == 0) return true;

    // get positions of the first frame
    if (this->posArraySize < mol->AtomCount()) {
        if (this->posArraySize > 0) {
            delete[] this->pos0;
            delete[] this->pos1;
        }
        pos0 = new float[mol->AtomCount() * 3];
        pos1 = new float[mol->AtomCount() * 3];
        this->posArraySize = mol->AtomCount();
    }
    memcpy(pos0, mol->AtomPositions(), mol->AtomCount() * 3 * sizeof(float));

    // set next frame ID and get positions of the second frame
    if ((static_cast<int>(callTime) + 1) < int(mol->FrameCount()))
        mol->SetFrameID(static_cast<int>(callTime) + 1);
    else
        mol->SetFrameID(static_cast<int>(callTime));
    // return false if calling for data fails
    if (!(*mol)(MolecularDataCall::CallForGetData)) {
        if (this->posArraySize > 0) {
            delete[] this->pos0;
            delete[] this->pos1;
        }
        this->posArraySize = 0;
        return false;
    }
    // time steps must have the same number of atoms
    if (this->posArraySize != mol->AtomCount()) {
        if (this->posArraySize > 0) {
            delete[] this->pos0;
            delete[] this->pos1;
        }
        this->posArraySize = 0;
        return false;
    }
    memcpy(pos1, mol->AtomPositions(), mol->AtomCount() * 3 * sizeof(float));


    // ---------- update parameters ----------
    this->UpdateParameters(mol, bs);

    // recompute color table, if necessary
    if (this->atomColorTableRGB.Count() / 3 < mol->AtomCount()) {
        // Mix two coloring modes
        protein::Color::MakeColorTable(mol, this->currentColoringMode0, this->currentColoringMode1,
            cmWeightParam.Param<param::FloatParam>()->Value(),        // weight for the first cm
            1.0f - cmWeightParam.Param<param::FloatParam>()->Value(), // weight for the second cm
            this->atomColorTableRGB, this->colorLookupTable, this->rainbowColors,
            this->minGradColorParam.Param<param::StringParam>()->Value(),
            this->midGradColorParam.Param<param::StringParam>()->Value(),
            this->maxGradColorParam.Param<param::StringParam>()->Value(), true);

        // copy to RGBA array
        this->atomColorTable.SetCount(4 * mol->AtomCount());
#pragma omp parallel for
        for (int i = 0; i < mol->AtomCount(); ++i) {
            this->atomColorTable[4 * i + 0] = this->atomColorTableRGB[3 * i + 0];
            // static_cast<float>(mol->AtomTypes()[mol->AtomTypeIndices()[i]].Colour()[0]) / 255.0f;
            this->atomColorTable[4 * i + 1] = this->atomColorTableRGB[3 * i + 1];
            // static_cast<float>(mol->AtomTypes()[mol->AtomTypeIndices()[i]].Colour()[1]) / 255.0f;
            this->atomColorTable[4 * i + 2] = this->atomColorTableRGB[3 * i + 2];
            // static_cast<float>(mol->AtomTypes()[mol->AtomTypeIndices()[i]].Colour()[2]) / 255.0f;
            this->atomColorTable[4 * i + 3] = 1.0f;
            this->uploadColors = true;
        }
    }

    // interpolate atom positions between frames
    this->posInter.SetCount(mol->AtomCount() * 4);
    float inter = callTime - static_cast<float>(static_cast<int>(callTime));
    float threshold = vislib::math::Min(mol->AccessBoundingBoxes().ObjectSpaceBBox().Width(),
                          vislib::math::Min(mol->AccessBoundingBoxes().ObjectSpaceBBox().Height(),
                              mol->AccessBoundingBoxes().ObjectSpaceBBox().Depth())) *
                      0.75f;
    int cnt;
#pragma omp parallel for
    for (cnt = 0; cnt < int(mol->AtomCount()); ++cnt) {
        if (std::sqrt(std::pow(pos0[3 * cnt + 0] - pos1[3 * cnt + 0], 2) +
                      std::pow(pos0[3 * cnt + 1] - pos1[3 * cnt + 1], 2) +
                      std::pow(pos0[3 * cnt + 2] - pos1[3 * cnt + 2], 2)) < threshold) {
            this->posInter[4 * cnt + 0] = (1.0f - inter) * pos0[3 * cnt + 0] + inter * pos1[3 * cnt + 0];
            this->posInter[4 * cnt + 1] = (1.0f - inter) * pos0[3 * cnt + 1] + inter * pos1[3 * cnt + 1];
            this->posInter[4 * cnt + 2] = (1.0f - inter) * pos0[3 * cnt + 2] + inter * pos1[3 * cnt + 2];
        } else if (inter < 0.5f) {
            this->posInter[4 * cnt + 0] = pos0[3 * cnt + 0];
            this->posInter[4 * cnt + 1] = pos0[3 * cnt + 1];
            this->posInter[4 * cnt + 2] = pos0[3 * cnt + 2];
        } else {
            this->posInter[4 * cnt + 0] = pos1[3 * cnt + 0];
            this->posInter[4 * cnt + 1] = pos1[3 * cnt + 1];
            this->posInter[4 * cnt + 2] = pos1[3 * cnt + 2];
        }
        this->posInter[4 * cnt + 3] = mol->AtomTypes()[mol->AtomTypeIndices()[cnt]].Radius() * 0.2f;
    }

    // ---------- render ----------
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glEnable(GL_VERTEX_PROGRAM_TWO_SIDE);
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_ARB);

    this->RenderSES(mol, this->posInter.PeekElements(), cam);

    // unlock the current frame
    mol->Unlock();

    return true;
}


/*
 * Render the molecular data in SES Trilateration mode.
 */
void SESMeshRenderer::RenderSES(MolecularDataCall* mol, const float* atomPos, view::Camera_2& cam) {

    cudaError_t err;

    /*Marco...--------------------------------------------------------------------------------------------*/
#ifdef TEST
    float time0, time1, time2, time3, time4, time5, time6, time7, time8;
    cudaEvent_t start1, stop1;
    cudaEventCreate(&start1);
    cudaEventCreate(&stop1);
    cudaEventRecord(start1, 0);
    cudaDeviceSynchronize();
#endif

    // get minimum of every x,y,z coordinate for moving the origin of the grid to positive coordinates
    auto bbox = mol->AccessBoundingBoxes().ObjectSpaceBBox();
    bbox.EnforcePositiveSize();
    float grid_x_max = bbox.GetRightTopBack().GetX();
    float grid_x_min = bbox.GetLeftBottomFront().GetX();

    float grid_y_max = bbox.GetRightTopBack().GetY();
    float grid_y_min = bbox.GetLeftBottomFront().GetY();

    float grid_z_min = bbox.GetRightTopBack().GetZ();
    float grid_z_max = bbox.GetLeftBottomFront().GetZ();

    float grid_x_range = grid_x_max - grid_x_min;
    float grid_y_range = grid_y_max - grid_y_min;
    float grid_z_range = grid_z_max - grid_z_min;


    //***************************************************************************************************//
    //***************************** GRID INSERT - ATOMS INTO BINS AND COUNT *****************************//
    //***************************************************************************************************//
    int atom_count = mol->AtomCount();

    int threads = NUM_THREADS;
    int blocks;
    threads = std::min(threads, atom_count);
    blocks = (atom_count % threads != 0) ? (atom_count / threads + 1) : (atom_count / threads);

    // get gird-demension - considering grid_step_size

    // move_origin_pdb: x,y,z for moving positions out of negative values (host)
    move_origin_pdb.x = grid_x_min;
    move_origin_pdb.y = grid_y_min;
    move_origin_pdb.z = grid_z_min;

    // TODO remove "|| true" -- this is just for testing the performance!
    if (mol->Calltime() != this->lastTimestamp || mem_manage == 0 ||
        this->probeRad != this->probeRadiusParam.Param<megamol::core::param::FloatParam>()->Value() || true) {

        this->probeRad = this->probeRadiusParam.Param<megamol::core::param::FloatParam>()->Value();
        this->grid_step_size = (this->probeRad * 2.0f + 4.0f);

        int3 griddims = {(int)ceilf(grid_x_range / grid_step_size), (int)ceilf(grid_y_range / grid_step_size),
            (int)ceilf(grid_z_range / grid_step_size)};

        this->lastTimestamp = mol->Calltime();

        int tmpGridsize = griddims.x * griddims.y * griddims.z;

        // (re-)allocate all arrays that depend in the number of grid cells
        if (this->gridsize < tmpGridsize || mem_manage == 0) {
            // d_grid: index+1 = atomID PDB; value=gridcell
            if (this->gridsize > 0) {
                cudaFree(d_grid);
                cudaFree(this->d_cell_start_end);
            }

            this->gridsize = static_cast<int>(static_cast<float>(tmpGridsize) * mem_reserve);

            // update gridsize: number of cells
            cudaMalloc((void**)&d_grid, this->gridsize * sizeof(int));
            cudaMalloc((void**)&d_cell_start_end, this->gridsize * sizeof(int2));
        }

        // (re-)allocate all arrays that depend in the number of atoms
        if (this->atomCnt < atom_count || mem_manage == 0) {
            if (this->atomCnt > 0) {
                // cudaFree(d_atomPos);
                cudaGraphicsUnregisterResource(this->gl_atomPosResource);
                glDeleteBuffers(1, &this->gl_atomPos);
                cudaGraphicsUnregisterResource(this->gl_atomPosIdxResource);
                glDeleteBuffers(1, &this->gl_atomPosIdx);
                glDeleteBuffers(1, &this->gl_atomColorTable);
                cudaFree(d_insert_grid);
                cudaFree(this->d_sorted);
                cudaFree(this->d_neighbours);
                cudaFree(this->d_neighbourCounts);
            }

            // d_insert_grid: index+1 = atomID PDB; x=gridcell, y= intern index of atoms in this cell
            cudaMalloc((void**)&d_insert_grid, (atom_count) * sizeof(int2));

            cudaMalloc((void**)&this->d_sorted, (atom_count) * sizeof(int2));
            // d_neighbours: 1D array of neighbours
            cudaMalloc((void**)&this->d_neighbours, (atom_count * max_num_neigbors) * sizeof(unsigned int));
            // d_neighbourCounts: index+1 = atomID; value=number of neighbours; size=atom_count (device)
            cudaMalloc((void**)&this->d_neighbourCounts, (atom_count) * sizeof(unsigned int));

            this->atomCnt = atom_count;

            // set SSBO for atomPos
            if (!glIsBuffer(this->gl_atomPos)) {
                glGenBuffers(1, &this->gl_atomPos);
            }
            CHECK_FOR_OGL_ERROR();
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, this->gl_atomPos);
            CHECK_FOR_OGL_ERROR();
            auto memsize = (this->atomCnt) * 4 * sizeof(float);
            glBufferStorage(GL_SHADER_STORAGE_BUFFER, memsize, 0, GL_DYNAMIC_STORAGE_BIT);
            CHECK_FOR_OGL_ERROR();
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
            CHECK_FOR_OGL_ERROR();

            err = cudaGraphicsGLRegisterBuffer(
                &this->gl_atomPosResource, this->gl_atomPos, cudaGraphicsRegisterFlagsWriteDiscard);

            // set SSBO for atomColorTable
            if (!glIsBuffer(this->gl_atomColorTable)) {
                glGenBuffers(1, &this->gl_atomColorTable);
            }
            CHECK_FOR_OGL_ERROR();
        }

        err = cudaGraphicsMapResources(1, &this->gl_atomPosResource, 0);
        err = cudaGraphicsResourceGetMappedPointer(
            (void**)(&this->d_atomPos), &this->gl_atomPosSize, this->gl_atomPosResource);

#ifdef reduceAtoms
        err = cudaGraphicsMapResources(1, &this->gl_atomPosIdxResource, 0);
        err = cudaGraphicsResourceGetMappedPointer(
            (void**)(&this->d_atomPosIdx), &this->gl_atomPosIdxSize, this->gl_atomPosIdxResource);
#endif
        err = cudaMemcpy(d_atomPos, atomPos, (atom_count) * sizeof(float4), cudaMemcpyHostToDevice);

        err = cudaMemset(d_insert_grid, 0, (atom_count) * sizeof(int2)); // must be set to 0 because of using atomicAdd

        // CUDA: insert atoms into bins and count atoms per bin (atomicAdd)
        err = cudaMemset(d_grid, 0, tmpGridsize * sizeof(int));
        gridInsertCuda(
            blocks, threads, d_grid, d_insert_grid, d_atomPos, move_origin_pdb, griddims, atom_count, grid_step_size);

#ifdef TEST
        FILE* myfile;
        int write = 0;
        // if ((double)(clock() - tStart) / CLOCKS_PER_SEC >= 20) write = 1;

        if (write) {

            myfile = fopen("prolint_time_measure.csv", "a");
        }
        cudaDeviceSynchronize();
        cudaEventRecord(stop1, 0);
        cudaEventSynchronize(stop1);
        cudaEventElapsedTime(&time1, start1, stop1);
        if (write) {
            fprintf(myfile, "GRID INSERT; %f \n", time1);
        } else
            printf("Time for GRID INSERT: \t\t %f ms\n", time1);

        cudaEvent_t start2, stop2;
        cudaEventCreate(&start2);
        cudaEventCreate(&stop2);
        cudaEventRecord(start2, 0);
        cudaDeviceSynchronize();
#endif

        //*************************************************************************//
        //***************************** COUNTING SORT *****************************//
        //*************************************************************************//

        // calculate Prefix Sum
        ScanCuda(d_grid, d_grid + tmpGridsize, d_grid); // in-place scan

        // TODO this seems to be unnecessary
        cudaMemset(this->d_sorted, 0, (atom_count) * sizeof(int2));
        // CUDA: counting sort on GPU
        countingSortCuda(blocks, threads, d_insert_grid, this->d_sorted, atom_count, d_grid);

#ifdef TEST
        cudaDeviceSynchronize();
        cudaEventRecord(stop2, 0);
        cudaEventSynchronize(stop2);
        cudaEventElapsedTime(&time2, start2, stop2);
        if (write) {
            fprintf(myfile, "COUNTING SORT; %f \n", time2);
        } else
            printf("Time for COUNTING SORT: \t %f ms\n", time2);

        cudaEvent_t start3, stop3;
        cudaEventCreate(&start3);
        cudaEventCreate(&stop3);
        cudaEventRecord(start3, 0);
        cudaDeviceSynchronize();

#endif
        //**********************************************************************//
        //***************************** REINDEXING *****************************//
        //**********************************************************************//

        // CUDA reindexing on GPU (cell_start, cell_end)
        cudaMemset(d_cell_start_end, 0, tmpGridsize * sizeof(int2));
        reindexCuda(blocks, threads, d_sorted, d_cell_start_end, atom_count);

#ifdef TEST
        cudaDeviceSynchronize();
        cudaEventRecord(stop3, 0);
        cudaEventSynchronize(stop3);
        cudaEventElapsedTime(&time3, start3, stop3);
        if (write) {
            fprintf(myfile, "REINDEXING; %f \n", time3);
        } else
            printf("Time for REINDEXING: \t\t %f ms\n", time3);

        cudaEvent_t start4, stop4;
        cudaEventCreate(&start4);
        cudaEventCreate(&stop4);
        cudaEventRecord(start4, 0);
        cudaDeviceSynchronize();

#endif

        //*****************************************************************************//
        //***************************** NEIGHBOUR SEARCH ******************************//
        //*****************************************************************************//

        // CUDA: searching neighbours for each atom on GPU (cell_start, cell_end)
        NeighbourSearchCuda(blocks, threads, d_atomPos, atom_count, move_origin_pdb, griddims, d_cell_start_end,
            d_sorted, d_neighbours, 0, d_neighbourCounts, 0, this->probeRad,
            this->grid_step_size);

#ifdef TEST
        cudaDeviceSynchronize();
        cudaEventRecord(stop4, 0);
        cudaEventSynchronize(stop4);
        cudaEventElapsedTime(&time4, start4, stop4);
        if (write) {
            fprintf(myfile, "NEIGBHOUR SERACH; %f \n", time4);
        } else
            printf("Time for NEIGBHOUR SERACH: \t %f ms\n", time4);

        cudaEvent_t start5, stop5;
        cudaEventCreate(&start5);
        cudaEventCreate(&stop5);
        cudaEventRecord(start5, 0);
        cudaDeviceSynchronize();

#endif

        //*****************************************************************************//
        //***************************** SPHERES ******************************//
        //*****************************************************************************//

		if (this->sphereVerticesCnt != this->sphere->getVertexCount()) {
            if (this->sphereTriangleCnt > 0) {
                err = cudaFree(this->d_sphereVertices);
                err = cudaFree(this->d_sphereTriangles);
            }
            this->sphereVerticesCnt = this->sphere->getVertexCount();
            this->sphereTriangleCnt = this->sphere->getTriangleCount();
            err = cudaMalloc((void**)&d_sphereVertices, this->sphereVerticesCnt * sizeof(float3));
            err = cudaMalloc((void**)&d_sphereTriangles, this->sphereTriangleCnt * sizeof(uint3));
        }
		// copy sphere data (for one default sphere) to GPU
        err = cudaMemcpy(this->d_sphereVertices, this->sphere->getVertices(), this->sphereVerticesCnt * sizeof(float3),
            cudaMemcpyHostToDevice);
        err = cudaMemcpy(this->d_sphereTriangles, this->sphere->getIndices(), this->sphereTriangleCnt * sizeof(uint3),
            cudaMemcpyHostToDevice);

		// allocate memory for sphere vertices for all atoms
        if (this->totalSphereTriangleCnt != this->sphere->getVertexCount() * this->atomCnt ) {
            if (this->totalSphereTriangleCnt > 0) {
                err = cudaFree(this->d_outerSphereTriangles);
                err = cudaFree(this->d_cutSphereTriangles);
            }
            this->totalSphereTriangleCnt = this->sphere->getVertexCount() * this->atomCnt;
            err = cudaMalloc((void**)&d_outerSphereTriangles, this->totalSphereTriangleCnt * sizeof(uint3));
            err = cudaMalloc((void**)&d_cutSphereTriangles, this->totalSphereTriangleCnt * sizeof(uint3));
        }


        // BASIC ARRAYS -> atomPos and color
        err = cudaGraphicsUnmapResources(1, &gl_atomPosResource, 0);
    }

    //######################################################################################################################

    //************************************************************************************//
    //************************************ RENDERING *************************************//
    //************************************************************************************//

    cam_type::snapshot_type snapshot;
    cam_type::matrix_type viewTemp, projTemp;

    // Generate complete snapshot and calculate matrices
    cam.calc_matrices(snapshot, viewTemp, projTemp, thecam::snapshot_content::all);

    glm::mat4 view = viewTemp;
    glm::mat4 proj = projTemp;
    glm::mat4 MVinv = glm::inverse(view);
    glm::mat4 MVP = proj * view;
    glm::mat4 MVPinv = glm::inverse(MVP);
    glm::mat4 MVPtransp = glm::transpose(MVP);

    if (this->uploadColors) {
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, this->gl_atomColorTable);
        CHECK_FOR_OGL_ERROR();
        glBufferData(GL_SHADER_STORAGE_BUFFER, this->atomColorTable.Count() * sizeof(float),
            this->atomColorTable.PeekElements(), GL_DYNAMIC_COPY);
        CHECK_FOR_OGL_ERROR();
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
        CHECK_FOR_OGL_ERROR();
        this->uploadColors = false;
    }

    using namespace vislib::sys;
    float fogStart = 0.5f;
    float* clearColor = new float[4];
    vislib::math::Vector<float, 3> fogCol(clearColor[0], clearColor[1], clearColor[2]);

    // get viewpoint parameters for raycasting
    float viewportStuff[4] = {0.0f, 0.0f, cam.resolution_gate().width(), cam.resolution_gate().height()};
    if (viewportStuff[2] < 1.0f) viewportStuff[2] = 1.0f;
    if (viewportStuff[3] < 1.0f) viewportStuff[3] = 1.0f;
    viewportStuff[2] = 2.0f / viewportStuff[2];
    viewportStuff[3] = 2.0f / viewportStuff[3];

    //************************************************************************************//
    //********************************* SPHERE RENDERER **********************************//
    //************************************************************************************//

    this->sphereShader.Enable();

    // set shader variables
    glUniform4fvARB(this->sphereShader.ParameterLocation("viewAttr"), 1, viewportStuff);
    glUniform3fvARB(
        this->sphereShader.ParameterLocation("camIn"), 1, glm::value_ptr(static_cast<glm::vec4>(snapshot.view_vector)));
    glUniform3fvARB(this->sphereShader.ParameterLocation("camRight"), 1,
        glm::value_ptr(static_cast<glm::vec4>(snapshot.right_vector)));
    glUniform3fvARB(
        this->sphereShader.ParameterLocation("camUp"), 1, glm::value_ptr(static_cast<glm::vec4>(snapshot.up_vector)));

    glUniform1fARB(this->sphereShader.ParameterLocation("alpha"), 1.0f);
    glUniformMatrix4fv(this->sphereShader.ParameterLocation("MVP"), 1, false, glm::value_ptr(MVP));
    glUniformMatrix4fv(this->sphereShader.ParameterLocation("MVinv"), 1, false, glm::value_ptr(MVinv));
    glUniformMatrix4fv(this->sphereShader.ParameterLocation("MVPinv"), 1, false, glm::value_ptr(MVPinv));
    glUniformMatrix4fv(this->sphereShader.ParameterLocation("MVPtransp"), 1, false, glm::value_ptr(MVPtransp));

    // bind vertex and color SSBOs and draw them
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, atomColorTableSSBObindingPoint, this->gl_atomColorTable);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, atomPosSSBObindingPoint, this->gl_atomPos);

    glDrawArrays(GL_POINTS, 0, atomCnt);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, atomPosSSBObindingPoint, 0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, atomColorTableSSBObindingPoint, 0);
	
    this->sphereShader.Disable();

	// ---------- sphere mesh rendering ----------

	glColor3f(1.0f, 0.5f, 0.0f);
    //unsigned int atomId = 30;
    this->triangleShader.Enable();
    glUniformMatrix4fv(this->triangleShader.ParameterLocation("MVP"), 1, false, glm::value_ptr(MVP));
    glUniformMatrix4fv(this->triangleShader.ParameterLocation("MVinv"), 1, false, glm::value_ptr(MVinv));
    for (unsigned int atomId = 0; atomId < mol->AtomCount(); atomId++) {
		glUniform3fv(this->triangleShader.ParameterLocation("origin"), 1, &this->posInter.PeekElements()[4*atomId]);
		this->sphere->draw();
    }
    this->triangleShader.Disable();
}


/*
 * update parameters
 */
void SESMeshRenderer::UpdateParameters(const MolecularDataCall* mol, const protein_calls::BindingSiteCall* bs) {

    // get molecule lust
    if (this->molIdxListParam.IsDirty()) {
        vislib::StringA tmpStr(this->molIdxListParam.Param<param::StringParam>()->Value());
        this->molIdxList = vislib::StringTokeniser<vislib::CharTraitsA>::Split(tmpStr, ';', true);
        this->molIdxListParam.ResetDirty();
    }

    bool colorTableChanged = false;
    // color table param
    if (this->colorTableFileParam.IsDirty()) {
        protein::Color::ReadColorTableFromFile(
            this->colorTableFileParam.Param<param::StringParam>()->Value(), this->colorLookupTable);
        this->colorTableFileParam.ResetDirty();
        colorTableChanged = true;
    }
    // Recompute color table
    if ((this->coloringModeParam0.IsDirty()) || (this->coloringModeParam1.IsDirty()) ||
        (this->cmWeightParam.IsDirty()) || colorTableChanged) {

        this->currentColoringMode0 =
            static_cast<protein::Color::ColoringMode>(int(this->coloringModeParam0.Param<param::EnumParam>()->Value()));

        this->currentColoringMode1 =
            static_cast<protein::Color::ColoringMode>(int(this->coloringModeParam1.Param<param::EnumParam>()->Value()));

        // Mix two coloring modes
        protein::Color::MakeColorTable(mol, this->currentColoringMode0, this->currentColoringMode1,
            cmWeightParam.Param<param::FloatParam>()->Value(),        // weight for the first cm
            1.0f - cmWeightParam.Param<param::FloatParam>()->Value(), // weight for the second cm
            this->atomColorTableRGB, this->colorLookupTable, this->rainbowColors,
            this->minGradColorParam.Param<param::StringParam>()->Value(),
            this->midGradColorParam.Param<param::StringParam>()->Value(),
            this->maxGradColorParam.Param<param::StringParam>()->Value(), true);

        // copy to RGBA array
        this->atomColorTable.SetCount(4 * mol->AtomCount());
#pragma omp parallel for
        for (int i = 0; i < mol->AtomCount(); ++i) {
            this->atomColorTable[4 * i + 0] = this->atomColorTableRGB[3 * i + 0];
            // static_cast<float>(mol->AtomTypes()[mol->AtomTypeIndices()[i]].Colour()[0]) / 255.0f;
            this->atomColorTable[4 * i + 1] = this->atomColorTableRGB[3 * i + 1];
            // static_cast<float>(mol->AtomTypes()[mol->AtomTypeIndices()[i]].Colour()[1]) / 255.0f;
            this->atomColorTable[4 * i + 2] = this->atomColorTableRGB[3 * i + 2];
            // static_cast<float>(mol->AtomTypes()[mol->AtomTypeIndices()[i]].Colour()[2]) / 255.0f;
            this->atomColorTable[4 * i + 3] = 1.0f;
        }

        this->coloringModeParam0.ResetDirty();
        this->coloringModeParam1.ResetDirty();
        this->cmWeightParam.ResetDirty();
        this->uploadColors = true;
    }
}
