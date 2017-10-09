TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    earth_sun.cpp \
    planet.cpp \
    ode_solver.cpp

HEADERS += \
    earth_sun.h \
    planet.h \
    ode_solver.h
