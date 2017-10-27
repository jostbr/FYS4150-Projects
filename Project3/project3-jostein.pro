TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    earth_sun.cpp \
    planet.cpp \
    nbody_solver.cpp

HEADERS += \
    planet.h \
    earth_sun.h \
    constants.h \
    nbody_solver.h

