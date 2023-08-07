/*
 * MegaMolPlugin.cpp
 *
 * Copyright (C) 2010 by VISUS (Universitaet Stuttgart)
 * Alle Rechte vorbehalten.
 */
#include "stdafx.h"

/*
 * This file is used under Windows only!
 */
#ifdef _WIN32


BOOL APIENTRY DllMain(HMODULE hModule, DWORD  ul_reason_for_call,
        LPVOID lpReserved) {
    switch (ul_reason_for_call) {
    case DLL_PROCESS_ATTACH:
    case DLL_THREAD_ATTACH:
    case DLL_THREAD_DETACH:
    case DLL_PROCESS_DETACH:
        break;
    }
    return TRUE;
}

#endif /* _WIN32 */
