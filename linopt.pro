TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    states.cpp \
    matrix.cpp \
    chip.cpp \
    cost_functor.cpp

HEADERS += \
    states.h \
    linopt.h \
    matrix.h \
    chip.h \
    cost_functor.h \
    bfgs.h \
    hco.h

#QMAKE_CXXFLAGS_RELEASE += -g -pg
#QMAKE_LFLAGS_RELEASE += -g -pg
QMAKE_CXXFLAGS_RELEASE += -march=native -ffast-math -O3
