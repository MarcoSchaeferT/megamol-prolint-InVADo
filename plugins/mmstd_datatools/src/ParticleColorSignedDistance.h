/*
 * ParticleColorSignedDistance.h
 *
 * Copyright (C) 2015 by S. Grottel
 * Alle Rechte vorbehalten.
 */

#ifndef MEGAMOLCORE_PARTICLECOLORSIGNEDDISTANCE_H_INCLUDED
#define MEGAMOLCORE_PARTICLECOLORSIGNEDDISTANCE_H_INCLUDED
#pragma once

#include "mmstd_datatools/AbstractParticleManipulator.h"
#include "mmcore/param/ParamSlot.h"
#include <vector>


namespace megamol {
namespace stdplugin {
namespace datatools {

    /**
     * Module overriding global attributes of particles
     */
    class ParticleColorSignedDistance : public AbstractParticleManipulator {
    public:

        /** Return module class name */
        static const char *ClassName(void) {
            return "ParticleColorSignedDistance";
        }

        /** Return module class description */
        static const char *Description(void) {
            return "Computes signed distances of particles to the closest particle with a color value of zero. The sign is preserved from the original color.";
        }

        /** Module is always available */
        static bool IsAvailable(void) {
            return true;
        }

        /** Ctor */
        ParticleColorSignedDistance(void);

        /** Dtor */
        virtual ~ParticleColorSignedDistance(void);

    protected:

        /**
         * Manipulates the particle data
         *
         * @remarks the default implementation does not changed the data
         *
         * @param outData The call receiving the manipulated data
         * @param inData The call holding the original data
         *
         * @return True on success
         */
        virtual bool manipulateData(
            megamol::core::moldyn::MultiParticleDataCall& outData,
            megamol::core::moldyn::MultiParticleDataCall& inData);

    private:

        void compute_colors(megamol::core::moldyn::MultiParticleDataCall& dat);
        void set_colors(megamol::core::moldyn::MultiParticleDataCall& dat);

        core::param::ParamSlot enableSlot;
        core::param::ParamSlot cyclXSlot;
        core::param::ParamSlot cyclYSlot;
        core::param::ParamSlot cyclZSlot;
        size_t datahash;
        unsigned int time;
        std::vector<float> newColors;
        float minCol, maxCol;

    };

} /* end namespace datatools */
} /* end namespace stdplugin */
} /* end namespace megamol */

#endif /* MEGAMOLCORE_PARTICLECOLORSIGNEDDISTANCE_H_INCLUDED */
