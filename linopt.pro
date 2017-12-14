TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    states.cpp \
    matrix.cpp \
    chip.cpp

HEADERS += \
    states.h \
    linopt.h \
    matrix.h \
    chip.h

QMAKE_CXXFLAGS_RELEASE += -march=native -ffast-math -O3
