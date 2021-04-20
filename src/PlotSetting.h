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

#ifndef PLOTSETTING_H
#define PLOTSETTING_H

#include <QOpenGLWidget>
#include <QtWidgets>
#include <cmath>

#ifdef __unix
#define powf pow
#endif

extern bool PB3D_DARK;

class PlotSetting : public QOpenGLWidget, protected QOpenGLFunctions {
Q_OBJECT
public:
	using QOpenGLWidget::QOpenGLWidget;

	struct PlotColor {
		QColor BG = PB3D_DARK ? QColor(40, 40, 40) : QColor(250, 250, 250);
		QColor NODE = QColor(166, 206, 227);
		QColor MASS = QColor(31, 120, 180);
		QColor FRAME = QColor(51, 160, 44);
		QColor HL = QColor(227, 26, 28);
		QColor GRID = QColor(30, 30, 30);
		QColor BC = QColor(106, 61, 154);
		QColor LOAD = QColor(255, 127, 0);
		QColor BRACE = QColor(177, 89, 40);
		QColor WALL = QColor(255, 255, 153);
	};

	struct PlotSize {
		int GRID_NUM = 10;
		float AXIS = 100;
		float GRID = 1;
		float PT = std::powf(10.f, 1.2f);
		float LINE_WIDTH = 1;
		float BC = std::powf(.5f, 1.2f);
		float LOAD = std::powf(.5f, 1.2f);
		float XSHIFT = 0;
		float YSHIFT = 0;
		float ZSHIFT = 0;
		float MASS = std::powf(15.f, 1.2f);
	};

	struct PlotSwitch {
		bool AXIS = true;
		bool GRID = true;
		bool NODE_LABEL = true;
		bool ELEMENT_LABEL = true;
		bool BC = true;
		bool LOAD = true;
		bool FRAME = true;
		bool BRACE = true;
		bool WALL = true;
	};

	struct PlotView {
		float XR = 0;
		float YR = 0;
		float ZR = 0;
		float XT = 0;
		float YT = 0;
		float ZT = -20;

		float FOV = 15;
		float NEAR_PLANE = 0;
		float FAR_PLANE = 0;
	};

	PlotSwitch Switch;
	PlotColor Color;
	PlotSize Size;
	PlotView View;
public slots:
	void setColorBG();
	void setColorNode();
	void setColorElement();
	void setColorTruss();
	void setColorHighlight();
	void setColorGrid();
	void setColorBC();
	void setColorLoad();
	void setColorMass();
	void setColorWall();

	void setSizeGridNumber(int);
	void setSizeAxis(int);
	void setSizeGrid(int);
	void setSizePoint(int);
	void setSizeLineWidth(int);
	void setSizeBC(int);
	void setSizeLoad(int);
	void setSizeXShift(int);
	void setSizeYShift(int);
	void setSizeZShift(int);
	void setSizeMass(int);

	void setSwitchAxis(bool);
	void setSwitchGrid(bool);
	void setSwitchNodeLabel(bool);
	void setSwitchElementLabel(bool);
	void setSwitchBC(bool);
	void setSwitchLoad(bool);
	void setSwitchFrame(bool);
	void setSwitchBrace(bool);
	void setSwitchWall(bool);

	void setViewXR(float);
	void setViewYR(float);
	void setViewZR(float);
	void setViewXT(float);
	void setViewYT(float);
	void setViewZT(float);

	void setViewFOV(int);
	void setViewNearPlane(float);
	void setViewFarPlane(float);

private:
	QVector<float> axis_ref = QVector<float>{0.f, .1f, .2f, .3f, .5f, 1.f, 2.f, 3.f, 5.f, 10.f, 20.f, 30.f, 50.f, 100.f, 200.f, 300.f, 500.f, 1000.f, 2000.f, 3000.f, 5000.f, 10000.f};

	QColor getColor();

	static float normaliseAngle(float);
	static float scaleSize(int);

protected:
	QMatrix4x4 getTransformation() const;

	QMatrix4x4 current_trans;
};

#endif // PLOTSETTING_H
