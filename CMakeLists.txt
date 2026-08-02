cmake_minimum_required(VERSION 3.20)

project(nmopt LANGUAGES CXX)

add_library(nmopt_contract INTERFACE)
target_include_directories(nmopt_contract INTERFACE include)
target_compile_features(nmopt_contract INTERFACE cxx_std_17)

include(CTest)

if(BUILD_TESTING)
  add_executable(nmopt_contract_tests tests/reduced_dto_contract.cc)
  target_link_libraries(nmopt_contract_tests PRIVATE nmopt_contract)
  add_test(NAME nmopt_contract_tests COMMAND nmopt_contract_tests)

  add_executable(nmopt_semantic_v1_contract_test tests/semantic_v1_contract.cc)
  target_link_libraries(nmopt_semantic_v1_contract_test PRIVATE nmopt_contract)
  add_test(NAME nmopt_semantic_v1_contract_test
    COMMAND nmopt_semantic_v1_contract_test)
endif()

option(NMOPT_ENABLE_DEAL_II "Build the serial deal.II reference lowerer" ON)

if(NMOPT_ENABLE_DEAL_II)
  find_package(deal.II CONFIG QUIET)

  if(deal.II_FOUND)
    deal_ii_initialize_cached_variables()
    add_library(nmopt_dealii_contract INTERFACE)
    target_link_libraries(nmopt_dealii_contract INTERFACE nmopt_contract)
    set(NMOPT_USING_DEAL_II_CMAKE_PACKAGE TRUE)
  else()
    find_path(NMOPT_DEAL_II_INCLUDE_DIR deal.II/base/config.h)
    find_library(NMOPT_DEAL_II_LIBRARY NAMES deal.ii)

    if(NMOPT_DEAL_II_INCLUDE_DIR AND NMOPT_DEAL_II_LIBRARY)
      add_library(nmopt_dealii_contract INTERFACE)
      target_link_libraries(nmopt_dealii_contract
        INTERFACE
          nmopt_contract
          "${NMOPT_DEAL_II_LIBRARY}")
      target_include_directories(nmopt_dealii_contract
        INTERFACE
          "${NMOPT_DEAL_II_INCLUDE_DIR}")
    endif()
  endif()

  if(TARGET nmopt_dealii_contract AND BUILD_TESTING)
    add_executable(nmopt_dealii_diffusion_contract_test
      tests/dealii_diffusion_contract.cc)
    target_link_libraries(nmopt_dealii_diffusion_contract_test
      nmopt_dealii_contract)
    if(NMOPT_USING_DEAL_II_CMAKE_PACKAGE)
      deal_ii_setup_target(nmopt_dealii_diffusion_contract_test)
    endif()
    add_test(NAME nmopt_dealii_diffusion_contract_test
      COMMAND nmopt_dealii_diffusion_contract_test)
  endif()
endif()
