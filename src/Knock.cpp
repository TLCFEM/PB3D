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

#include "ModelBuilder.h"

#include <QApplication>
#include <QScreen>
#include <QSettings>
#include <QStyleFactory>
#include <QSurfaceFormat>

bool PB3D_DARK = false;

int main(int argc, char* argv[]) {
	QApplication app(argc, argv);
    QApplication::setApplicationName("PB3D");
    QApplication::setApplicationDisplayName("PB3D");
	QApplication::setOrganizationName("University of Canterbury");
	QApplication::setWindowIcon(QIcon(":/../res/UC.ico"));

	auto font = QApplication::font();
	const auto rec = QGuiApplication::primaryScreen()->availableGeometry();
	if(std::max(rec.height(), rec.width()) > 2000) font.setPointSize(12);
	else font.setPointSize(9);
	QApplication::setFont(font);

	qApp->setStyle(QStyleFactory::create("fusion"));

#ifdef Q_OS_WIN
	const QSettings settings(R"(HKEY_CURRENT_USER\Software\Microsoft\Windows\CurrentVersion\Themes\Personalize)", QSettings::NativeFormat);
	if(settings.value("AppsUseLightTheme") == 0) {
		QPalette darkPalette;
		const auto darkColor = QColor(45, 45, 45);
		const auto disabledColor = QColor(209, 17, 65);
		darkPalette.setColor(QPalette::Window, darkColor);
		darkPalette.setColor(QPalette::WindowText, Qt::white);
		darkPalette.setColor(QPalette::Base, QColor(18, 18, 18));
		darkPalette.setColor(QPalette::AlternateBase, darkColor);
		darkPalette.setColor(QPalette::ToolTipBase, Qt::white);
		darkPalette.setColor(QPalette::ToolTipText, Qt::white);
		darkPalette.setColor(QPalette::Text, Qt::white);
		darkPalette.setColor(QPalette::Disabled, QPalette::Text, disabledColor);
		darkPalette.setColor(QPalette::Button, darkColor);
		darkPalette.setColor(QPalette::ButtonText, Qt::white);
		darkPalette.setColor(QPalette::Disabled, QPalette::ButtonText, disabledColor);
		darkPalette.setColor(QPalette::BrightText, Qt::red);
		darkPalette.setColor(QPalette::Link, QColor(42, 130, 218));

		darkPalette.setColor(QPalette::Highlight, QColor(42, 130, 218));
		darkPalette.setColor(QPalette::HighlightedText, Qt::black);
		darkPalette.setColor(QPalette::Disabled, QPalette::HighlightedText, disabledColor);

		qApp->setPalette(darkPalette);

		qApp->setStyleSheet("QToolTip { color: #ffffff; background-color: #2a82da; border: 1px solid white; }");

        PB3D_DARK = true;
	}
#endif

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
