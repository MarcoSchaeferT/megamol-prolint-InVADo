/*
 * LigandModelSelector.h
 *
 * Copyright (C) 2019 by University of Tuebingen (BDVA).
 * All rights reserved.
 */

#ifndef PROLINT_PLUGIN_LIGANDMODELSELECTOR_H_INCLUDED
#define PROLINT_PLUGIN_LIGANDMODELSELECTOR_H_INCLUDED
typedef unsigned int uint;

#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"
#include "mmcore/param/ParamSlot.h"
#include "LigandModelCall.h"
#include "mmcore/utility/log/Log.h"

using namespace std;

namespace megamol {
namespace prolint {


/****************************
 **** PDBQT LOADER CLASS ****
 ****************************/
class LigandModelSelector : public megamol::core::Module {


public:
    /** Ctor. */
    LigandModelSelector();
    /** Dtor. */
    virtual ~LigandModelSelector();

    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName(void) { return "LigandModelSelector"; }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description(void) { return "Selects docking data."; }

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
     * Callback function for current ligand parameter slot
     */
    bool currentLigandOrModelChanged(megamol::core::param::ParamSlot& param);

private:
    /** LigandModelCall caller slot */
    megamol::core::CallerSlot ligandModelCall;

    core::param::ParamSlot ligandIdxParam;
    core::param::ParamSlot modelIdxParam;

};


} // namespace prolint
} /* end namespace megamol */

#endif // PROLINT_PLUGIN_LIGANDMODELSELECTOR_H_INCLUDED
