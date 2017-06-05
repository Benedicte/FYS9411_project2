TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    basis.cpp \
    hartree_fock.cpp \
    coulomb_functions.cpp \
    ccd.cpp \
    test.cpp \
    examples.cpp \
    coulomb_function2.cpp \
    ccd_copy.cpp \
    ccd_matrix.cpp

HEADERS += \
    basis.h \
    hartree_fock.h \
    coulomb_functions.h \
    ccd.h \
    test.h \
    examples.h \
    coulomb_function2.h \
    abstract_coulomb.h \
    ccd_copy.h \
    ccd_matrix.h

LIBS += -L/usr/local/Cellar/armadillo/7.500.0/lib -larmadillo

INCLUDEPATH += /usr/local/Cellar/armadillo/7.500.0/include
