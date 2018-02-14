TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    variationalmontecarlo.cpp \
    hamiltonian.cpp \
    sampling.cpp \
    wavefunction.cpp

HEADERS += \
    variationalmontecarlo.h \
    hamiltonian.h \
    sampling.h \
    wavefunction.h
