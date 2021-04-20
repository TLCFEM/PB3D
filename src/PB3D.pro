QT       += core gui opengl svg

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

greaterThan(QT_MAJOR_VERSION, 5): QT += openglwidgets svgwidgets

CONFIG += c++17

RC_ICONS = ../res/UC.ico

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

DEFINES += NDEBUG

SOURCES += \
    Database.cpp \
    Knock.cpp \
    ModelRenderer.cpp \
    ModelBuilder.cpp \
    PlotSetting.cpp

HEADERS += \
    Database.h \
    ModelBuilder.h \
    ModelRenderer.h \
    PlotSetting.h

FORMS += \
    ModelBuilder.ui

win32{
LIBS += -lopengl32
}

unix{
LIBS += -lGL -lglut
}

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

RESOURCES += \
    PB3D.qrc
