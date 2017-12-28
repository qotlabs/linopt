TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

INCLUDEPATH += /usr/include/eigen3
INCLUDEPATH += ./lib

LIBS += -L../linopt/lib -llinopt

#QMAKE_CXXFLAGS_RELEASE += -g -pg
#QMAKE_LFLAGS_RELEASE += -g -pg
QMAKE_LFLAGS += -fopenmp
QMAKE_CXXFLAGS += -fopenmp
QMAKE_CXXFLAGS_RELEASE += -march=native -ffast-math -O3
