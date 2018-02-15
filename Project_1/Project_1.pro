TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    src/hamiltonian.cpp \
    src/sampling.cpp \
    src/variationalmontecarlo.cpp \
    src/wavefunction.cpp \
    src/tests.cpp \
    src/matrix.cpp

HEADERS += inc/catch.hpp \
    inc/hamiltonian.h \
    inc/sampling.h \
    inc/variationalmontecarlo.h \
    inc/wavefunction.h \
    inc/matrix.h
