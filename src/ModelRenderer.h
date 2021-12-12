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

#ifndef MODELRENDERER_H
#define MODELRENDERER_H

#include <PlotSetting.h>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>

class Database;

class ModelRenderer final : public PlotSetting {
    Q_OBJECT
public:
    using PlotSetting::PlotSetting;

    void setModel(Database*);

public slots:
    void resetView();

protected:
    void initializeGL() override;
    void paintGL() override;
    void mousePressEvent(QMouseEvent*) override;
    void mouseMoveEvent(QMouseEvent*) override;
    void wheelEvent(QWheelEvent*) override;

private:
    static const char* vertexSource;
    static const char* fragmentSource;

    Database* model_ptr = nullptr;

    QOpenGLBuffer m_buffer = QOpenGLBuffer(QOpenGLBuffer::Type::VertexBuffer);
    std::unique_ptr<QOpenGLShaderProgram> m_program = nullptr;

    QPoint m_last_pos;
    int m_trans_mat = 0;

    void setPlane();

    void paintAxis();
    void paintNode();
    void paintNodeLabel();
    void paintElementLabel();
    void paintElement();
    void paintBC();
    void paintLoad();
    void paintMass();

    void appendFixX(std::vector<GLfloat>&, const QVector3D&) const;
    void appendFixY(std::vector<GLfloat>&, const QVector3D&) const;
    void appendFixZ(std::vector<GLfloat>&, const QVector3D&) const;
    void appendFixRX(std::vector<GLfloat>&, const QVector3D&) const;
    void appendFixRY(std::vector<GLfloat>&, const QVector3D&) const;
    void appendFixRZ(std::vector<GLfloat>&, const QVector3D&) const;
    void appendFixColor(std::vector<GLfloat>&) const;

    void appendLoadX(std::vector<GLfloat>&, const QVector3D&, float) const;
    void appendLoadY(std::vector<GLfloat>&, const QVector3D&, float) const;
    void appendLoadZ(std::vector<GLfloat>&, const QVector3D&, float) const;
    void appendLoadRX(std::vector<GLfloat>&, const QVector3D&, float);
    void appendLoadRY(std::vector<GLfloat>&, const QVector3D&, float);
    void appendLoadRZ(std::vector<GLfloat>&, const QVector3D&, float);
    void appendLoadColor(std::vector<GLfloat>&) const;

    void renderLabel(QPainter&, const QVector3D&, const QString&) const;
};

#endif // MODELRENDERER_H
