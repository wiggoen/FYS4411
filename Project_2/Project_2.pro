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
        main.cpp

HEADERS += \
    inc/json.hpp
