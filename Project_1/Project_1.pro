TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    src/hamiltonian.cpp \
    src/variationalmontecarlo.cpp \
    src/wavefunction.cpp \
    src/tests.cpp

HEADERS += inc/catch.hpp \
    inc/hamiltonian.h \
    inc/variationalmontecarlo.h \
    inc/wavefunction.h

INCLUDEPATH += /usr/local/Cellar/armadillo/8.400.0/include/

LIBS += -L/usr/local/Cellar/armadillo/8.400.0/lib/ -larmadillo
