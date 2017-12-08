TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    array_alloc.cpp \
    unit_tests.cpp \
    poisson.cpp \
    periodic_solver_1d.cpp \
    basin_solver_1d.cpp \
    rossby_solver.cpp \
    basin_solver_2d.cpp \
    periodic_solver_2d.cpp

HEADERS += \
    array_alloc.hpp \
    unit_tests.hpp \
    poisson.hpp \
    periodic_solver_1d.hpp \
    basin_solver_1d.hpp \
    rossby_solver.hpp \
    basin_solver_2d.hpp \
    periodic_solver_2d.hpp

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib
LIBS += -larmadillo -llapack -lblas
