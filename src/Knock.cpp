////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2021 Theodore Chang, Minghao Li
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
////////////////////////////////////////////////////////////////////////////////

#include <QApplication>
#include <QScreen>
#include <QSettings>
#include <QStyleFactory>
#include <QSurfaceFormat>
#include "ModelBuilder.h"

int main(int argc, char* argv[]) {
    QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);

    QApplication app(argc, argv);
    QApplication::setApplicationName("PB3D");
    QApplication::setApplicationDisplayName("PB3D");
    QApplication::setOrganizationName("University of Canterbury");
    QApplication::setWindowIcon(QIcon(":/../res/UC.ico"));

    //    auto font = QApplication::font();
    //    const auto rec = QGuiApplication::primaryScreen()->availableGeometry();
    //    if(std::max(rec.height(), rec.width()) > 2000)
    //        font.setPointSize(14);
    //    QApplication::setFont(font);

    qApp->setStyle(QStyleFactory::create("fusion"));

    ModelBuilder win;
    win.setWindowTitle("PB3D");

    QSurfaceFormat format;
    format.setDepthBufferSize(24);
    format.setStencilBufferSize(8);
    format.setVersion(2, 0);
    format.setProfile(QSurfaceFormat::CompatibilityProfile);
    QSurfaceFormat::setDefaultFormat(format);

    win.show();

    return QApplication::exec();
}
