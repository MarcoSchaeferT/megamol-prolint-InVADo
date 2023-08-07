/*
 * LigandModelSelector.cpp
 *
 * Copyright (C) 2019 by University of Tuebingen (BDVA).
 * All rights reserved.
 */


#include "stdafx.h"
#include "LigandModelSelector.h"
#include "mmcore/param/IntParam.h"


using namespace megamol;
using namespace megamol::core;
using namespace megamol::prolint;

/**
 * Constructor
 */
LigandModelSelector::LigandModelSelector(void)
    : Module()
	, ligandModelCall("getData", "Call for a specific model of a specific ligand")
    , ligandIdxParam("LigandIdx", "Index of a ligand")
    , modelIdxParam("ModelIdx", "Specifies which model of the current ligand will be choosen") {
    // set up caller slot
    this->ligandModelCall.SetCompatibleCall<LigandModelCallDescription>();
    this->MakeSlotAvailable(&this->ligandModelCall);

	// set up parameters
    this->ligandIdxParam << new param::IntParam(0);
    this->ligandIdxParam.SetUpdateCallback(&LigandModelSelector::currentLigandOrModelChanged);
    this->MakeSlotAvailable(&this->ligandIdxParam);
    this->modelIdxParam << new param::IntParam(0);
    this->modelIdxParam.SetUpdateCallback(&LigandModelSelector::currentLigandOrModelChanged);
    this->MakeSlotAvailable(&this->modelIdxParam);

}

/**
 * Destructor
 */
LigandModelSelector::~LigandModelSelector(void) { this->Release(); }

bool LigandModelSelector::create(void) {
    // TODO allocate variables etc.
    return true;
}

void LigandModelSelector::release(void) {

}


bool LigandModelSelector::currentLigandOrModelChanged(param::ParamSlot& param) {
    // set param as clean again
    param.ResetDirty();
	// get pointer to LigandModelCall
    LigandModelCall* lmc = this->ligandModelCall.CallAs<LigandModelCall>();
    if (lmc == NULL) return false;
	// call for extent (number of ligands and models per ligand)
    (*lmc)(LigandModelCall::CallForGetExtent);

	// check if the values were out of bounds 
    int ligandIdx = this->ligandIdxParam.Param<param::IntParam>()->Value();
    int modelIdx = this->modelIdxParam.Param<param::IntParam>()->Value();
    if (ligandIdx >= (int)(lmc->getLigandCount())) {
        ligandIdx = lmc->getLigandCount()-1;
        this->ligandIdxParam.Param<param::IntParam>()->SetValue(ligandIdx);
    } else if (ligandIdx < 0) {
        ligandIdx = 0;
        this->ligandIdxParam.Param<param::IntParam>()->SetValue(ligandIdx);
	}

	// set ligand and CALL for data to be able to check the number of models for this ligand
    lmc->SetCurrentLigandAndModel(ligandIdx, 0);
    (*lmc)(LigandModelCall::CallForGetData);

    if (modelIdx > (int)(lmc->getModelCount()[ligandIdx])) {
        modelIdx = lmc->getModelCount()[ligandIdx];
        this->modelIdxParam.Param<param::IntParam>()->SetValue(modelIdx);
    } else if (modelIdx < 0) {
        modelIdx = 0;
        this->modelIdxParam.Param<param::IntParam>()->SetValue(modelIdx);
    }

	// set current ligand and model idx to the call
    lmc->SetCurrentLigandAndModel(ligandIdx, modelIdx);
	// call for get data (i.e., send requested ligand and model idx to called object)
    (*lmc)(LigandModelCall::CallForGetData);
    return true;
}
