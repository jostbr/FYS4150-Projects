TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    jacobi.cpp \
    unit_tests.cpp \
    initialize.cpp

HEADERS += \
    jacobi.h \
    unit_tests.h \
    initialize.h

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib
LIBS += -larmadillo -llapack -lblas
