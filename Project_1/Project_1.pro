TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    variationalmontecarlo.cpp \
    hamiltonian.cpp \
    sampling.cpp \
    wavefunction.cpp \
    src/hamiltonian.cpp \
    src/sampling.cpp \
    src/variationalmontecarlo.cpp \
    src/wavefunction.cpp \
    main.cpp

HEADERS += \
    variationalmontecarlo.h \
    hamiltonian.h \
    sampling.h \
    wavefunction.h \
    inc/hamiltonian.h \
    inc/sampling.h \
    inc/variationalmontecarlo.h \
    inc/wavefunction.h
