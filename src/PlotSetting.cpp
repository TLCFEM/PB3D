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

#include "PlotSetting.h"

void PlotSetting::setColorBG() {
    Color.BG = getColor();
    update();
}

void PlotSetting::setColorNode() {
    Color.NODE = getColor();
    update();
}

void PlotSetting::setColorElement() {
    Color.FRAME = getColor();
    update();
}

void PlotSetting::setColorTruss() {
    Color.BRACE = getColor();
    update();
}

void PlotSetting::setColorHighlight() {
    Color.HL = getColor();
    update();
}

void PlotSetting::setColorGrid() {
    Color.GRID = getColor();
    update();
}

void PlotSetting::setColorBC() {
    Color.BC = getColor();
    update();
}

void PlotSetting::setColorLoad() {
    Color.LOAD = getColor();
    update();
}

void PlotSetting::setColorMass() {
    Color.MASS = getColor();
    update();
}

void PlotSetting::setColorWall() {
    Color.WALL = getColor();
    update();
}

void PlotSetting::setSizeGridNumber(const int F) {
    Size.GRID_NUM = std::max(0, F);
    update();
}

void PlotSetting::setSizeAxis(const int F) {
    Size.AXIS = scaleSize(F);
    update();
}

void PlotSetting::setSizeGrid(const int F) {
    if(F >= axis_ref.size())
        Size.GRID = axis_ref.back();
    else
        Size.GRID = axis_ref.at(F);
    update();
}

void PlotSetting::setSizePoint(const int F) {
    Size.PT = scaleSize(F);
    update();
}

void PlotSetting::setSizeLineWidth(const int F) {
    Size.LINE_WIDTH = std::min(10.f, std::max(static_cast<float>(F), 0.f));
    update();
}

void PlotSetting::setSizeBC(const int F) {
    Size.BC = scaleSize(F);
    update();
}

void PlotSetting::setSizeLoad(const int F) {
    Size.LOAD = scaleSize(F);
    update();
}

void PlotSetting::setSizeXShift(const int F) {
    Size.XSHIFT = scaleSize(F);
    update();
}

void PlotSetting::setSizeYShift(const int F) {
    Size.YSHIFT = scaleSize(F);
    update();
}

void PlotSetting::setSizeZShift(const int F) {
    Size.ZSHIFT = scaleSize(F);
    update();
}

void PlotSetting::setSizeMass(const int F) {
    Size.MASS = scaleSize(F);
    update();
}

void PlotSetting::setSwitchAxis(const bool F) {
    Switch.AXIS = F;
    update();
}

void PlotSetting::setSwitchGrid(const bool F) {
    Switch.GRID = F;
    update();
}

void PlotSetting::setSwitchNodeLabel(const bool F) {
    Switch.NODE_LABEL = F;
    update();
}

void PlotSetting::setSwitchElementLabel(const bool F) {
    Switch.ELEMENT_LABEL = F;
    update();
}

void PlotSetting::setSwitchBC(const bool F) {
    Switch.BC = F;
    update();
}

void PlotSetting::setSwitchLoad(const bool F) {
    Switch.LOAD = F;
    update();
}

void PlotSetting::setSwitchFrame(const bool F) {
    Switch.FRAME = F;
    update();
}

void PlotSetting::setSwitchBrace(const bool F) {
    Switch.BRACE = F;
    update();
}

void PlotSetting::setSwitchWall(const bool F) {
    Switch.WALL = F;
    update();
}

void PlotSetting::setViewXR(const float F) {
    View.XR = normaliseAngle(F);
    update();
}

void PlotSetting::setViewYR(const float F) {
    View.YR = normaliseAngle(F);
    update();
}

void PlotSetting::setViewZR(const float F) {
    View.ZR = normaliseAngle(F);
    update();
}

void PlotSetting::setViewXT(const float F) {
    View.XT = F;
    update();
}

void PlotSetting::setViewYT(const float F) {
    View.YT = F;
    update();
}

void PlotSetting::setViewZT(const float F) {
    View.ZT = F;
    update();
}

void PlotSetting::setViewFOV(const int F) {
    View.FOV = std::min(std::max(0.f, static_cast<float>(F)), 150.f);
    update();
}

void PlotSetting::setViewNearPlane(const float F) {
    View.NEAR_PLANE = std::max(.01f, F);
    update();
}

void PlotSetting::setViewFarPlane(const float F) {
    View.FAR_PLANE = std::max(.01f, F);
    update();
}

QColor PlotSetting::getColor() {
    auto dialog = QColorDialog(this);
    dialog.adjustSize();
    dialog.exec();
    return dialog.selectedColor();
}

float PlotSetting::normaliseAngle(float angle) {
    while(true) {
        if(angle < 0.f)
            angle += 360.f;
        else if(angle > 360.f)
            angle -= 360.f;
        else
            return angle;
    }
}

float PlotSetting::scaleSize(const int F) { return (F >= 0 ? 1.f : -1.f) * std::powf(.01f * static_cast<float>(F >= 0 ? F : -F), 1.2f); }

QMatrix4x4 PlotSetting::getTransformation() const {
    QMatrix4x4 trans_mat;

    trans_mat.setToIdentity();
    trans_mat.perspective(View.FOV, static_cast<GLfloat>(width()) / static_cast<GLfloat>(height()), View.NEAR_PLANE, View.FAR_PLANE);
    trans_mat.translate(View.XT, View.YT, View.ZT);
    trans_mat.rotate(View.XR, 1, 0, 0);
    trans_mat.rotate(View.YR, 0, 1, 0);
    trans_mat.rotate(View.ZR, 0, 0, 1);

    return trans_mat;
}
