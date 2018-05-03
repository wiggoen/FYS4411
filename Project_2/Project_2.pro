TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt


# Remove default optimization flag used by QT
QMAKE_CXXFLAGS_RELEASE -= -O2

# Add optimization flags
#QMAKE_CXXFLAGS_RELEASE *= -O0  # Run this when developing
QMAKE_CXXFLAGS_RELEASE *= -O3  # Run this on release


SOURCES += \
    main.cpp \
    src/hamiltonian.cpp \
    src/tests.cpp \
    src/variationalmontecarlo.cpp \
    src/wavefunction.cpp

HEADERS += \
    inc/json.hpp \
    inc/catch.hpp \
    inc/hamiltonian.h \
    inc/variationalmontecarlo.h \
    inc/wavefunction.h
