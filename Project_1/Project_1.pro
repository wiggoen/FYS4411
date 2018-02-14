TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
#DEFINES += CATCH_CONFIG_MAIN # uncomment this when using unit_test.cpp. Remember to rebuild when uncommenting and commenting!

SOURCES += main.cpp \
    src/hamiltonian.cpp \
    src/sampling.cpp \
    src/unit_test.cpp \
    src/variationalmontecarlo.cpp \
    src/wavefunction.cpp

HEADERS += inc/catch.hpp \
    inc/hamiltonian.h \
    inc/sampling.h \
    inc/variationalmontecarlo.h \
    inc/wavefunction.h
