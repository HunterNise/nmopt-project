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
endif()
