TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    ising.cpp \
    unit_tests.cpp

HEADERS += \
    ising.h \
    unit_tests.h
