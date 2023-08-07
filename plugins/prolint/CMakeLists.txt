#
# MegaMol™ prolint Plugin
# Copyright 2017, by MegaMol Team
# Alle Rechte vorbehalten. All rights reserved.
#
include(CheckLanguage)

check_language(CUDA)

if(MSVC AND CMAKE_CUDA_COMPILER)
  option(BUILD_PROLINT_PLUGIN "Option to build protein_cuda" OFF)
else()
  set(BUILD_PROLINT_PLUGIN OFF)
endif()

set(Python3_FIND_ABI "ANY" "ANY" "ANY")
find_package(Python3 COMPONENTS Interpreter Development)

if(!Python3_FOUND)
  message(WARNING "unable to find Python 3...")
endif()

execute_process(COMMAND pip install -r ${CMAKE_CURRENT_SOURCE_DIR}\\requirements.txt)

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}\\InVADo_setup)
else()
	execute_process(COMMAND CMD /c cd ${CMAKE_CURRENT_SOURCE_DIR} && git clone https://github.com/MarcoSchaeferT/InVADo_setup.git)
	
	execute_process(COMMAND CMD /c curl -L https://github.com/Orbmu2k/nvidiaProfileInspector/releases/download/2.4.0.4/nvidiaProfileInspector.zip  --output ${CMAKE_CURRENT_SOURCE_DIR}\\InVADo_setup\\graphicSettings\\Nprofile.zip)
	execute_process(COMMAND CMD /c cd ${CMAKE_CURRENT_SOURCE_DIR}\\InVADo_setup\\graphicSettings && tar -xf Nprofile.zip)
	execute_process(COMMAND CMD /c cd ${CMAKE_CURRENT_SOURCE_DIR}\\InVADo_setup\\graphicSettings && del Nprofile.zip)
	
	execute_process(COMMAND CMD /c curl -L https://homepage.univie.ac.at/norbert.haider/download/chemistry/checkmol/bin/checkmol-latest-win32.zip  --output ${CMAKE_CURRENT_SOURCE_DIR}\\InVADo_setup\\checkmol.zip)
	execute_process(COMMAND CMD /c cd ${CMAKE_CURRENT_SOURCE_DIR}\\InVADo_setup && tar -xf checkmol.zip)
	execute_process(COMMAND CMD /c cd ${CMAKE_CURRENT_SOURCE_DIR}\\InVADo_setup && del checkmol.zip)
	
	execute_process(COMMAND CMD /c curl -L https://github.com/mozilla/geckodriver/releases/download/v0.33.0/geckodriver-v0.33.0-win64.zip --output ${CMAKE_CURRENT_SOURCE_DIR}\\InVADo_setup\\server\\geckodriver.zip)
	execute_process(COMMAND CMD /c cd ${CMAKE_CURRENT_SOURCE_DIR}\\InVADo_setup\\server && tar -xf geckodriver.zip)
	execute_process(COMMAND CMD /c cd ${CMAKE_CURRENT_SOURCE_DIR}\\InVADo_setup\\server && del geckodriver.zip)
	
	execute_process(COMMAND CMD /c cd ${CMAKE_CURRENT_SOURCE_DIR} && xcopy InVADo_setup\\checkmol.exe  /H /C /I /Y)
	execute_process(COMMAND CMD /c cd ${CMAKE_CURRENT_SOURCE_DIR} && xcopy InVADo_setup\\server server /E /H /C /I /Y)
endif()


cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# install npm packages for vue (used for website of the docking project)
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}\\server\\vue_app\\node_modules)
else()
	message("COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR}\\server\\vue_app && npm install && npm run build")
	execute_process(COMMAND CMD /c cd ${CMAKE_CURRENT_SOURCE_DIR}\\server\\vue_app && npm install && npm run build && echo OUTPUT_VARIABLE FOO)
endif()

if(BUILD_PROLINT_PLUGIN)
  project(prolint LANGUAGES CXX CUDA)

  string(REPLACE "." ";" CUDA_VERSION_LIST ${CMAKE_CUDA_COMPILER_VERSION})
  list(GET CUDA_VERSION_LIST 0 CUDA_VERSION_MAJOR)

  if(NOT CUDA_VERSION_MAJOR GREATER 7)
    message(FATAL_ERROR "PROLINT currently only supports CUDA Version 8 or above. The version found was CUDA ${CUDA_VERSION_MAJOR} (${CMAKE_CUDA_COMPILER_VERSION}). Please turn off the Protein_CUDA plugin or install a correct version of CUDA.")
  endif()

  string(TOUPPER prolint EXPORT_NAME)
  option(BUILD_${EXPORT_NAME}_PLUGIN "Option to build prolint" ON)

  if(BUILD_${EXPORT_NAME}_PLUGIN)
    project(prolint)

    # Check option dependencies: First list item is always the current build option,
    # afterwards list the dependencies.
    set(DEP_LIST "${DEP_LIST};BUILD_${EXPORT_NAME}_PLUGIN BUILD_PROTEIN_CALLS_PLUGIN BUILD_PROTEIN_PLUGIN BUILD_GEOMETRY_CALLS_PLUGIN BUILD_CORE" CACHE INTERNAL "")

    # NVCC options
    set(CMAKE_CUDA_STANDARD 11)
    set(CMAKE_CUDA_FLAGS_RELEASE "${CMAKE_CUDA_FLAGS} -O3")
    set(CMAKE_CUDA_ARCHITECTURES "50")

    get_filename_component(CUDA_COMPILER_DIRECTORY "${CMAKE_CUDA_COMPILER}" DIRECTORY)

    set(LIBS vislib core protein_calls geometry_calls protein glad)

    if(MSVC)
      include_directories("include" "src")
    else()
      get_filename_component(CUDA_BIN_DIR ${CMAKE_CUDA_COMPILER} DIRECTORY)
      include_directories("include" "src" "src/helper_includes" ${CUDA_BIN_DIR}/../include)
      add_definitions(-D__CUDACC__)
    endif()

    # Build MegaMol Plugin
    # public header files for installation
    file(GLOB_RECURSE public_header_files RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "include/*.h")

    # source files for compilation
    file(GLOB_RECURSE source_files RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "src/*.cpp")
    file(GLOB_RECURSE header_files RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "src/*.h")
    file(GLOB_RECURSE shader_files RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "Shaders/*")
    file(GLOB_RECURSE resource_files RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "Resources/*")

    # CUDA related files
    file(GLOB_RECURSE cuda_source_files RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "src/*.cu")
    file(GLOB_RECURSE cuda_header_files RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "src/*.cuh")

    file(GLOB_RECURSE cuda_std_header_files RELATIVE "${CUDA_COMPILER_DIRECTORY}" "../include/*.h")

    set_source_files_properties(cuda_source_files PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(cuda_header_files PROPERTIES LANGUAGE CUDA)

    # remove special files from the build
    list(REMOVE_ITEM cuda_source_files "src/cuenergy.cu")
    list(REMOVE_ITEM cuda_source_files "src/gpu_poisson_solver.cu")
    list(REMOVE_ITEM cuda_source_files "src/particles_kernel.cu")

    if(UNIX)
      list(REMOVE_ITEM source_files
        "src/dllmain.cpp"
      )
    endif()

    # Target definition
    add_library(${PROJECT_NAME} SHARED ${public_header_files} ${header_files} ${shader_files} ${source_files} ${cuda_header_files} ${cuda_source_files})

    # Set target naming conventions for different build types
    set_target_properties(${PROJECT_NAME} PROPERTIES SUFFIX ".mmplg")

    # Dependent projects only need to link against the core.so itself.
    target_link_libraries(${PROJECT_NAME} PRIVATE core protein_calls geometry_calls protein ${Python3_LIBRARIES})
    target_include_directories(${PROJECT_NAME} PUBLIC ${Python3_INCLUDE_DIRS} $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> "include" "src")

    string(TOUPPER ${PROJECT_NAME} EXPORT_NAME)
    target_compile_definitions(${PROJECT_NAME} PRIVATE ${EXPORT_NAME}_EXPORTS)

    # target_include_directories(${PROJECT_NAME} PUBLIC
    # $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    # "include" "src")
    # target_link_libraries(${PROJECT_NAME} PRIVATE core)

    # Installation rules for generated files
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION "include")
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Shaders/ DESTINATION "share/shaders")
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Resources/ DESTINATION "share/resources")

    if(WIN32)
      install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION "bin")
      install(TARGETS ${PROJECT_NAME} ARCHIVE DESTINATION "lib")
    else()
      install(TARGETS ${PROJECT_NAME} DESTINATION "lib" EXPORT ${PROJECT_NAME}-target)
    endif()

    # Grouping in Visual Studio
    set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER plugins)
    source_group("Public Header Files" FILES ${public_header_files})
    source_group("Header Files" FILES ${header_files})
    source_group("CUDA Header Files" FILES ${cuda_header_files})
    source_group("CUDA Source Files" FILES ${cuda_source_files})
    source_group("Source Files" FILES ${source_files})
    source_group("Shaders" FILES ${shader_files})
  endif()
endif(BUILD_PROLINT_PLUGIN)
