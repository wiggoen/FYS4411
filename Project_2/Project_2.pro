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



# Linux specific
unix:!macx {

}

# Mac specific
macx: {
    # MPI header
    #HEADERS += /usr/local/Cellar/open-mpi/3.1.0/include/mpi.h

    # MPI Settings
    #QMAKE_CC = $$system(/usr/local/bin/mpicc)    # C compiler
    #QMAKE_CXX = $$system(/usr/local/bin/mpicxx)  # C++ compiler
    #QMAKE_CXX_RELEASE = $$QMAKE_CXX
    #QMAKE_CXX_DEBUG = $$QMAKE_CXX
    #QMAKE_LINK = $$QMAKE_CXX

    #QMAKE_CFLAGS += $$system(/usr/local/bin/mpicc --showme:compile) -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
    #QMAKE_LFLAGS += $$system(/usr/local/bin/mpicxx --showme:link)
    ##QMAKE_CXXFLAGS += $$system(/usr/local/bin/mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
    ##QMAKE_CXXFLAGS_RELEASE += $$system(/usr/local/bin/mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX

    INCLUDEPATH += /usr/local/Cellar/armadillo/8.500.1/include/ #\
                   #/usr/local/Cellar/open-mpi/3.1.0/include

    LIBS += -L/usr/local/Cellar/armadillo/8.500.1/lib/ -larmadillo #\
            #-L/usr/local/opt/libevent/lib -L/usr/local/Cellar/open-mpi/3.1.0/lib -lmpi
}
