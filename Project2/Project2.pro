TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    jacobirota.cpp

INCLUDEPATH = \home\trude\Downloads

LIBS += -llapack -lblas -larmadillo

