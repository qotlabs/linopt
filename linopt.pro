TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

#QMAKE_CXXFLAGS_RELEASE += -g -pg
#QMAKE_LFLAGS_RELEASE += -g -pg
QMAKE_CXXFLAGS_RELEASE += -march=native -ffast-math -O3
