TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    array_alloc.cpp \
    unit_tests.cpp \
    poisson.cpp

HEADERS += \
    array_alloc.hpp \
    unit_tests.hpp \
    poisson.hpp
