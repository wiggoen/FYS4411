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


    # MPI-flag for compiling
    #QMAKE_CXXFLAGS += -DMPI_ON         # NOT IMPLEMENTED YET


    # MPI Settings
MPI_ON {
    QMAKE_CXX = mpicxx
    QMAKE_CXX_RELEASE = $$QMAKE_CXX
    QMAKE_CXX_DEBUG = $$QMAKE_CXX
    QMAKE_LINK = $$QMAKE_CXX
    QMAKE_CC = mpicc
    QMAKE_CFLAGS += -O3 -std=c++11 $$system(mpicc --showme:compile)
    QMAKE_LFLAGS += $$system(mpicxx --showme:link)
    QMAKE_CXXFLAGS += -O3 -std=c++11 $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
    QMAKE_CXXFLAGS_RELEASE += -O3 -std=c++11 $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
}

# Linux specific
unix:!macx {
    LIBS += -llapack -lblas -larmadillo
}

# Mac specific
macx: {
    INCLUDEPATH += /usr/local/Cellar/armadillo/8.500.1/include/

    LIBS += -L/usr/local/Cellar/armadillo/8.500.1/lib/ -larmadillo
}
