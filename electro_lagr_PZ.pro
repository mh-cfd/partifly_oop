TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    Main.cpp \
    tools.cpp \
    particlessolver.cpp \
    fieldsloader.cpp

HEADERS += \
    tools.h \
    particlessolver.h \
    fieldsloader.h
#QMAKE_CXXFLAGS += -O2
QMAKE_CXXFLAGS_RELEASE += -O3 -ffast-math  -msse -std=c++11
#LIBS += -lopenGL32 -lGLU32 -lm
#LIBS += -L$$PWD/my_lib -lglut32

QMAKE_LFLAGS += -O3 -ffast-math  -msse -std=c++11

INCLUDEPATH += /usr/include/paraview
#INCLUDEPATH += /usr/lib/paraview

LIBS+=  -lGL -lGLU -lglut -lm -lpthread
LIBS += -L/usr/lib/paraview -lvtkCommonCore -lvtkFiltersCore -lvtkIOImage -lvtkIOCore -lvtkIOXMLParser -lvtkCommonDataModel -lvtkIOXML -lvtkIOCore -lvtkIOExport -lvtkIOLegacy

