TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    states.cpp \
    matrix.cpp

HEADERS += \
    states.h \
    linopt.h \
    matrix.h

QMAKE_CXXFLAGS_RELEASE += -march=native -ffast-math -O3
