TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    array_alloc.cpp \
    unit_tests.cpp \
    poisson.cpp \
    rossby_solver_1d.cpp \
    periodic_solver_1d.cpp \
    basin_solver_1d.cpp

HEADERS += \
    array_alloc.hpp \
    unit_tests.hpp \
    poisson.hpp \
    rossby_solver_1d.hpp \
    periodic_solver_1d.hpp \
    basin_solver_1d.hpp
