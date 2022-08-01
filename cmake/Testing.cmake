# enable testing with CTest (but prevent CTest from adding additional CDash targets)
set_property(GLOBAL PROPERTY CTEST_TARGETS_ADDED 1)
include(CTest)

# Catch2 v3
find_package(Catch2 3 CONFIG REQUIRED)
include(Catch)

# FakeIt
find_path(FAKEIT_INCLUDE_DIRS "boost/fakeit.hpp")
