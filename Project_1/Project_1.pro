TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    src/hamiltonian.cpp \
    src/sampling.cpp \
    src/variationalmontecarlo.cpp \
    src/wavefunction.cpp \

HEADERS += inc/hamiltonian.h \
    inc/sampling.h \
    inc/variationalmontecarlo.h \
    inc/wavefunction.h \
