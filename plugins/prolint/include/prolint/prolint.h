/*
 * prolint.h
 * Copyright (C) 2009-2015 by MegaMol Team
 * Alle Rechte vorbehalten.
 */
#ifndef PROLINT_H_INCLUDED
#define PROLINT_H_INCLUDED
#if (defined(_MSC_VER) && (_MSC_VER > 1000))
#pragma once
#endif /* (defined(_MSC_VER) && (_MSC_VER > 1000)) */


#ifdef _WIN32
// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the TRISOUPPLUGIN_EXPORTS
// symbol defined on the command line. this symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// TRISOUPPLUGIN_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
#ifdef PROLINT_EXPORTS
#define PROLINT_API __declspec(dllexport)
#else
#define PROLINT_API __declspec(dllimport)
#endif
#else /* _WIN32 */
#define PROLINT_API
#endif /* _WIN32 */

#include "mmcore/utility/plugins/Plugin200Instance.h"

#ifdef PROLINT_EXPORTS
#ifdef __cplusplus
extern "C" {
#endif

/**
 * Returns the version of the MegaMol� plugin api used by this plugin.
 *
 * @return 200 -- (ver.: 2.00)
 */
PROLINT_API int mmplgPluginAPIVersion(void);

/**
 * Provides compatibility information
 *
 * @param onError Callback function pointer used when an error occures
 *
 * @return The compatibility information struct, or nullptr in case of an
 *         error.
 *
 * @remarks Always use 'mmplgReleasePluginCompatibilityInfo' to release the
 *          memory of the returned struct.
 */
PROLINT_API
::megamol::core::utility::plugins::PluginCompatibilityInfo *
mmplgGetPluginCompatibilityInfo(
    ::megamol::core::utility::plugins::ErrorCallback onError);

/**
 * Releases the memory of a compatibility information struct previously
 * returned by 'mmplgGetPluginCompatibilityInfo'
 *
 * @param ci The compatibility information struct to be released
 */
PROLINT_API void mmplgReleasePluginCompatibilityInfo(
    ::megamol::core::utility::plugins::PluginCompatibilityInfo* ci);

/**
 * Creates a new instance of this plugin
 *
 * @param onError Callback function pointer used when an error occures
 *
 * @return A new instance of this plugin, or nullptr in case of an error
 *
 * @remarks Always use 'mmplgReleasePluginInstance' to release the memory of
 *          the returned object.
 */
PROLINT_API
::megamol::core::utility::plugins::AbstractPluginInstance*
mmplgGetPluginInstance
    (::megamol::core::utility::plugins::ErrorCallback onError);

/**
 * Releases the memory of the plugin instance previously returned by
 * 'mmplgGetPluginInstance'
 *
 * @param pi The plugin instance to be released
 */
PROLINT_API void mmplgReleasePluginInstance(
    ::megamol::core::utility::plugins::AbstractPluginInstance* pi);

#ifdef __cplusplus
} /* extern "C" */
#endif
#endif

#endif /* PROLINT_H_INCLUDED */
