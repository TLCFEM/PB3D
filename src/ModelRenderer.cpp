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

// ReSharper disable CppClangTidyBugproneNarrowingConversions
#include "ModelRenderer.h"
#include <Database.h>
#include <QMouseEvent>
#include <cmath>

void ModelRenderer::renderLabel(QPainter& painter, const QVector3D& position, const QString& string) const {
    const auto canvas_pos = current_trans * QVector4D(position.x(), position.y(), position.z(), 1.);

    const auto norm_z = canvas_pos.z() / canvas_pos.w();

    if(norm_z > 1.f || norm_z < -1.f) return;

    painter.drawText((.5 + .5 * canvas_pos.x() / canvas_pos.w()) * width(), (.5 - .5 * canvas_pos.y() / canvas_pos.w()) * height(), string);
}

const char* ModelRenderer::vertexSource =
    "#version 130\n"
    "in vec3 position;"
    "in vec3 i_color;"
    "out vec3 m_color;"
    "uniform mat4 trans_mat;"
    "void main(){"
    "gl_Position=trans_mat*vec4(position,1.0);"
    "m_color=i_color;"
    "}";

const char* ModelRenderer::fragmentSource =
    "#version 130\n"
    "in vec3 m_color;"
    "out vec4 o_color;"
    "void main(){"
    "o_color=vec4(m_color,1.);"
    "}";

void ModelRenderer::setModel(Database* ptr) { model_ptr = ptr; }

void ModelRenderer::resetView() {
    View = PlotView();

    repaint();
}

void ModelRenderer::initializeGL() {
    initializeOpenGLFunctions();
    glEnable(GL_DEPTH_TEST | GL_CULL_FACE | GL_LINE_SMOOTH);

    m_program = std::make_unique<QOpenGLShaderProgram>();
    m_program->addShaderFromSourceCode(QOpenGLShader::Vertex, vertexSource);
    m_program->addShaderFromSourceCode(QOpenGLShader::Fragment, fragmentSource);
    m_program->bindAttributeLocation("position", 0);
    m_program->bindAttributeLocation("i_color", 1);
    m_program->link();
    m_program->bind();

    m_trans_mat = m_program->uniformLocation("trans_mat");

    m_program->release();

    m_buffer.create();
}

void ModelRenderer::paintGL() {
    glClearColor(Color.BG.redF(), Color.BG.greenF(), Color.BG.blueF(), 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    setPlane();

    m_program->bind();
    m_program->setUniformValue(m_trans_mat, current_trans = getTransformation());

    glPointSize(1);
    glLineWidth(1);

    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);

    if(Switch.AXIS) paintAxis();

    glPointSize(Size.PT);
    glLineWidth(Size.LINE_WIDTH);

    paintNode();
    paintElement();
    paintBC();
    paintLoad();
    paintMass();

    if(Switch.NODE_LABEL) paintNodeLabel();
    if(Switch.ELEMENT_LABEL) paintElementLabel();

    m_program->release();
}

void ModelRenderer::paintAxis() {
    auto idx = 36;

    const auto limit = Size.GRID_NUM / 2;

    std::vector<GLfloat> verts(idx + (2llu * limit + 1llu) * 24, 0.);

    verts[0] = verts[13] = verts[26] = -0.f * (verts[6] = verts[19] = verts[32] = Size.AXIS);
    verts[3] = verts[9] = verts[16] = verts[22] = verts[29] = verts[35] = 1;

    const auto bound = static_cast<float>(limit) * Size.GRID;

    for(auto I = -limit; I <= limit; ++I) {
        const auto edge = static_cast<float>(I) * Size.GRID;
        verts[idx] = edge;
        verts[idx + 1ll] = -bound;
        verts[idx + 3ll] = Color.GRID.redF();
        verts[idx + 4ll] = Color.GRID.greenF();
        verts[idx + 5ll] = Color.GRID.blueF();
        idx += 6;
        verts[idx] = edge;
        verts[idx + 1ll] = 0 == I ? 0 : bound;
        verts[idx + 3ll] = Color.GRID.redF();
        verts[idx + 4ll] = Color.GRID.greenF();
        verts[idx + 5ll] = Color.GRID.blueF();
        idx += 6;
        verts[idx] = -bound;
        verts[idx + 1ll] = edge;
        verts[idx + 3ll] = Color.GRID.redF();
        verts[idx + 4ll] = Color.GRID.greenF();
        verts[idx + 5ll] = Color.GRID.blueF();
        idx += 6;
        verts[idx] = 0 == I ? 0 : bound;
        verts[idx + 1ll] = edge;
        verts[idx + 3ll] = Color.GRID.redF();
        verts[idx + 4ll] = Color.GRID.greenF();
        verts[idx + 5ll] = Color.GRID.blueF();
        idx += 6;
    }

    m_buffer.bind();

    m_buffer.allocate(verts.data(), sizeof(GLfloat) * static_cast<int>(verts.size()));

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), nullptr);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), reinterpret_cast<void*>(3 * sizeof(GLfloat)));

    m_buffer.release();

    glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(verts.size() / 6));
}

void ModelRenderer::paintNode() {
    std::vector<GLfloat> node_data;

    const auto& node_pool = model_ptr->getNodePool();

    node_data.reserve(6 * node_pool.size());

    for(auto& [fst, snd] : node_pool) {
        node_data.emplace_back(snd.position.x());
        node_data.emplace_back(snd.position.y());
        node_data.emplace_back(snd.position.z());
        if(snd.highlighted) {
            node_data.emplace_back(Color.HL.redF());
            node_data.emplace_back(Color.HL.greenF());
            node_data.emplace_back(Color.HL.blueF());
        }
        else {
            node_data.emplace_back(Color.NODE.redF());
            node_data.emplace_back(Color.NODE.greenF());
            node_data.emplace_back(Color.NODE.blueF());
        }
    }

    m_buffer.bind();

    m_buffer.allocate(node_data.data(), sizeof(GLfloat) * static_cast<int>(node_data.size()));

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), nullptr);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), reinterpret_cast<void*>(3 * sizeof(GLfloat)));

    glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(node_pool.size()));

    m_buffer.release();
}

void ModelRenderer::paintNodeLabel() {
    const auto& node_pool = model_ptr->getNodePool();

    QPainter painter(this);
    painter.setPen(Color.HL);
    QFont font;
    font.setPointSize(14);
    painter.setFont(font);

    for(auto& [fst, snd] : node_pool) {
        const auto text = QString::number(fst);
        renderLabel(painter, snd.position + QVector3D{Size.XSHIFT, Size.YSHIFT, Size.ZSHIFT}, text);
    }

    painter.end();
}

void ModelRenderer::paintElement() {
    std::vector<GLfloat> element_data;

    auto ele_color = [&](const Database::Element& snd) {
        if(snd.highlighted) {
            element_data.emplace_back(Color.HL.redF());
            element_data.emplace_back(Color.HL.greenF());
            element_data.emplace_back(Color.HL.blueF());
        }
        else if(snd.type == Database::Element::Type::Frame) {
            element_data.emplace_back(Color.FRAME.redF());
            element_data.emplace_back(Color.FRAME.greenF());
            element_data.emplace_back(Color.FRAME.blueF());
        }
        else if(snd.type == Database::Element::Type::Wall) {
            element_data.emplace_back(Color.WALL.redF());
            element_data.emplace_back(Color.WALL.greenF());
            element_data.emplace_back(Color.WALL.blueF());
        }
        else {
            element_data.emplace_back(Color.BRACE.redF());
            element_data.emplace_back(Color.BRACE.greenF());
            element_data.emplace_back(Color.BRACE.blueF());
        }
    };

    const auto& node_pool = model_ptr->getNodePool();
    const auto& element_pool = model_ptr->getElementPool();

    element_data.reserve(12 * element_pool.size());

    for(const auto& [fst, snd] : element_pool) {
        if(snd.type == Database::Element::Type::Wall && !Switch.WALL) continue;
        if(snd.type == Database::Element::Type::Frame && !Switch.FRAME) continue;
        if(snd.type == Database::Element::Type::Brace && !Switch.BRACE) continue;

        element_data.emplace_back(node_pool.at(snd.encoding.at(0)).position.x());
        element_data.emplace_back(node_pool.at(snd.encoding.at(0)).position.y());
        element_data.emplace_back(node_pool.at(snd.encoding.at(0)).position.z());
        ele_color(snd);
        element_data.emplace_back(node_pool.at(snd.encoding.at(1)).position.x());
        element_data.emplace_back(node_pool.at(snd.encoding.at(1)).position.y());
        element_data.emplace_back(node_pool.at(snd.encoding.at(1)).position.z());
        ele_color(snd);
    }

    m_buffer.bind();

    m_buffer.allocate(element_data.data(), sizeof(GLfloat) * static_cast<int>(element_data.size()));

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), nullptr);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), reinterpret_cast<void*>(3 * sizeof(GLfloat)));

    glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(element_data.size() / 6));

    m_buffer.release();
}

void ModelRenderer::paintElementLabel() {
    const auto& node_pool = model_ptr->getNodePool();
    const auto& element_pool = model_ptr->getElementPool();

    QPainter painter(this);
    painter.setPen(Color.HL);
    QFont font;
    font.setPointSize(14);
    painter.setFont(font);

    for(auto& [fst, snd] : element_pool) {
        if(snd.type == Database::Element::Type::Wall && !Switch.WALL) continue;
        if(snd.type == Database::Element::Type::Frame && !Switch.FRAME) continue;
        if(snd.type == Database::Element::Type::Brace && !Switch.BRACE) continue;

        QVector3D position(0, 0, 0);
        auto counter = 0.f;
        for(auto& I : snd.encoding) {
            position += node_pool.at(I).position;
            counter += 1.f;
        }
        const auto text = QString::number(fst);
        renderLabel(painter, position / counter + QVector3D{Size.XSHIFT, Size.YSHIFT, Size.ZSHIFT}, text);
    }

    painter.end();
}

void ModelRenderer::paintBC() {
    std::vector<GLfloat> data;
    std::vector<GLenum> type;

    const auto& node_pool = model_ptr->getNodePool();

    auto bc_num = 0;

    for(const auto& [fst, snd] : node_pool) {
        if(snd.fixity.at(0) || snd.fixity.at(3)) ++bc_num;
        if(snd.fixity.at(1) || snd.fixity.at(4)) ++bc_num;
        if(snd.fixity.at(2) || snd.fixity.at(5)) ++bc_num;
    }

    data.reserve(24llu * bc_num);
    type.reserve(bc_num);

    for(const auto& [fst, snd] : node_pool) {
        if(snd.fixity.at(0) && snd.fixity.at(3)) {
            appendFixX(data, snd.position);
            type.emplace_back(GL_QUADS);
        }
        else if(snd.fixity.at(0)) {
            appendFixX(data, snd.position);
            type.emplace_back(GL_LINE_LOOP);
        }
        else if(snd.fixity.at(3)) {
            appendFixRX(data, snd.position);
            type.emplace_back(GL_LINE_LOOP);
        }

        if(snd.fixity.at(1) && snd.fixity.at(4)) {
            appendFixY(data, snd.position);
            type.emplace_back(GL_QUADS);
        }
        else if(snd.fixity.at(1)) {
            appendFixY(data, snd.position);
            type.emplace_back(GL_LINE_LOOP);
        }
        else if(snd.fixity.at(4)) {
            appendFixRY(data, snd.position);
            type.emplace_back(GL_LINE_LOOP);
        }

        if(snd.fixity.at(2) && snd.fixity.at(5)) {
            appendFixZ(data, snd.position);
            type.emplace_back(GL_QUADS);
        }
        else if(snd.fixity.at(2)) {
            appendFixZ(data, snd.position);
            type.emplace_back(GL_LINE_LOOP);
        }
        else if(snd.fixity.at(5)) {
            appendFixRZ(data, snd.position);
            type.emplace_back(GL_LINE_LOOP);
        }
    }

    m_buffer.bind();

    m_buffer.allocate(data.data(), sizeof(GLfloat) * static_cast<int>(data.size()));

    for(auto I = 0llu, J = 0llu; I < data.size(); I += 24, ++J) {
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), reinterpret_cast<void*>(I * sizeof(GLfloat)));
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), reinterpret_cast<void*>((3 + I) * sizeof(GLfloat)));

        glDrawArrays(type[J], 0, static_cast<GLsizei>(4));
    }

    m_buffer.release();
}

void ModelRenderer::paintLoad() {
    std::vector<GLfloat> data;

    auto load_num = 0;

    const auto& node_pool = model_ptr->getNodePool();

    for(const auto& [fst, snd] : node_pool) {
        if(snd.load.at(0) != 0.f) ++load_num;
        if(snd.load.at(1) != 0.f) ++load_num;
        if(snd.load.at(2) != 0.f) ++load_num;
    }

    data.reserve(48llu * load_num);

    for(const auto& [fst, snd] : node_pool) {
        if(snd.load.at(0) != 0.f) appendLoadX(data, snd.position, snd.load.at(0));
        if(snd.load.at(1) != 0.f) appendLoadY(data, snd.position, snd.load.at(1));
        if(snd.load.at(2) != 0.f) appendLoadZ(data, snd.position, snd.load.at(2));
    }

    m_buffer.bind();

    m_buffer.allocate(data.data(), sizeof(GLfloat) * static_cast<int>(data.size()));

    for(auto I = 0llu, J = 0llu; I < data.size(); I += 48, ++J) {
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), reinterpret_cast<void*>(I * sizeof(GLfloat)));
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), reinterpret_cast<void*>((3 + I) * sizeof(GLfloat)));

        glDrawArrays(GL_LINE_STRIP, 0, static_cast<GLsizei>(8));
    }

    m_buffer.release();
}

void ModelRenderer::paintMass() {
    std::vector<GLfloat> node_data;

    const auto& node_pool = model_ptr->getNodePool();

    node_data.reserve(6 * node_pool.size());

    for(auto& [fst, snd] : node_pool)
        if(snd.mass > 0.) {
            node_data.emplace_back(snd.position.x());
            node_data.emplace_back(snd.position.y());
            node_data.emplace_back(snd.position.z());
            node_data.emplace_back(Color.MASS.redF());
            node_data.emplace_back(Color.MASS.greenF());
            node_data.emplace_back(Color.MASS.blueF());
        }

    if(node_data.size() == 0) return;

    glPointSize(Size.MASS);

    m_buffer.bind();

    m_buffer.allocate(node_data.data(), sizeof(GLfloat) * static_cast<int>(node_data.size()));

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), nullptr);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), reinterpret_cast<void*>(3 * sizeof(GLfloat)));

    glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(node_data.size() / 6));

    m_buffer.release();

    glPointSize(Size.PT);
}

void ModelRenderer::appendFixX(std::vector<GLfloat>& data, const QVector3D& position) const {
    const auto& x = position.x();
    const auto& y = position.y();
    const auto& z = position.z();

    data.emplace_back(x);
    data.emplace_back(y - Size.BC);
    data.emplace_back(z - Size.BC);
    appendFixColor(data);

    data.emplace_back(x);
    data.emplace_back(y + Size.BC);
    data.emplace_back(z - Size.BC);
    appendFixColor(data);

    data.emplace_back(x);
    data.emplace_back(y + Size.BC);
    data.emplace_back(z + Size.BC);
    appendFixColor(data);

    data.emplace_back(x);
    data.emplace_back(y - Size.BC);
    data.emplace_back(z + Size.BC);
    appendFixColor(data);
}

void ModelRenderer::appendFixY(std::vector<GLfloat>& data, const QVector3D& position) const {
    const auto& x = position.x();
    const auto& y = position.y();
    const auto& z = position.z();

    data.emplace_back(x - Size.BC);
    data.emplace_back(y);
    data.emplace_back(z - Size.BC);
    appendFixColor(data);

    data.emplace_back(x + Size.BC);
    data.emplace_back(y);
    data.emplace_back(z - Size.BC);
    appendFixColor(data);

    data.emplace_back(x + Size.BC);
    data.emplace_back(y);
    data.emplace_back(z + Size.BC);
    appendFixColor(data);

    data.emplace_back(x - Size.BC);
    data.emplace_back(y);
    data.emplace_back(z + Size.BC);
    appendFixColor(data);
}

void ModelRenderer::appendFixZ(std::vector<GLfloat>& data, const QVector3D& position) const {
    const auto& x = position.x();
    const auto& y = position.y();
    const auto& z = position.z();

    data.emplace_back(x - Size.BC);
    data.emplace_back(y - Size.BC);
    data.emplace_back(z);
    appendFixColor(data);

    data.emplace_back(x - Size.BC);
    data.emplace_back(y + Size.BC);
    data.emplace_back(z);
    appendFixColor(data);

    data.emplace_back(x + Size.BC);
    data.emplace_back(y + Size.BC);
    data.emplace_back(z);
    appendFixColor(data);

    data.emplace_back(x + Size.BC);
    data.emplace_back(y - Size.BC);
    data.emplace_back(z);
    appendFixColor(data);
}

void ModelRenderer::appendFixRX(std::vector<GLfloat>& data, const QVector3D& position) const {
    const auto& x = position.x();
    const auto& y = position.y();
    const auto& z = position.z();

    data.emplace_back(x);
    data.emplace_back(y - Size.BC);
    data.emplace_back(z);
    appendFixColor(data);

    data.emplace_back(x);
    data.emplace_back(y);
    data.emplace_back(z + Size.BC);
    appendFixColor(data);

    data.emplace_back(x);
    data.emplace_back(y + Size.BC);
    data.emplace_back(z);
    appendFixColor(data);

    data.emplace_back(x);
    data.emplace_back(y);
    data.emplace_back(z - Size.BC);
    appendFixColor(data);
}

void ModelRenderer::appendFixRY(std::vector<GLfloat>& data, const QVector3D& position) const {
    const auto& x = position.x();
    const auto& y = position.y();
    const auto& z = position.z();

    data.emplace_back(x - Size.BC);
    data.emplace_back(y);
    data.emplace_back(z);
    appendFixColor(data);

    data.emplace_back(x);
    data.emplace_back(y);
    data.emplace_back(z + Size.BC);
    appendFixColor(data);

    data.emplace_back(x + Size.BC);
    data.emplace_back(y);
    data.emplace_back(z);
    appendFixColor(data);

    data.emplace_back(x);
    data.emplace_back(y);
    data.emplace_back(z - Size.BC);
    appendFixColor(data);
}

void ModelRenderer::appendFixRZ(std::vector<GLfloat>& data, const QVector3D& position) const {
    const auto& x = position.x();
    const auto& y = position.y();
    const auto& z = position.z();

    data.emplace_back(x - Size.BC);
    data.emplace_back(y);
    data.emplace_back(z);
    appendFixColor(data);

    data.emplace_back(x);
    data.emplace_back(y + Size.BC);
    data.emplace_back(z);
    appendFixColor(data);

    data.emplace_back(x + Size.BC);
    data.emplace_back(y);
    data.emplace_back(z);
    appendFixColor(data);

    data.emplace_back(x);
    data.emplace_back(y - Size.BC);
    data.emplace_back(z);
    appendFixColor(data);
}

void ModelRenderer::appendFixColor(std::vector<GLfloat>& data) const {
    data.emplace_back(Color.BC.redF());
    data.emplace_back(Color.BC.greenF());
    data.emplace_back(Color.BC.blueF());
}

void ModelRenderer::appendLoadX(std::vector<GLfloat>& data, const QVector3D& position, const float load) const {
    const auto& x = position.x();
    const auto& y = position.y();
    const auto& z = position.z();

    const auto size = load > 0.f ? -Size.LOAD : Size.LOAD;
    const auto size_2 = 2 * size;
    const auto size_6 = 6 * size;

    data.emplace_back(x + size_6);
    data.emplace_back(y);
    data.emplace_back(z);
    appendLoadColor(data);
    data.emplace_back(x);
    data.emplace_back(y);
    data.emplace_back(z);
    appendLoadColor(data);
    data.emplace_back(x + size_2);
    data.emplace_back(y + size);
    data.emplace_back(z);
    appendLoadColor(data);
    data.emplace_back(x + size_2);
    data.emplace_back(y - size);
    data.emplace_back(z);
    appendLoadColor(data);
    data.emplace_back(x);
    data.emplace_back(y);
    data.emplace_back(z);
    appendLoadColor(data);
    data.emplace_back(x + size_2);
    data.emplace_back(y);
    data.emplace_back(z + size);
    appendLoadColor(data);
    data.emplace_back(x + size_2);
    data.emplace_back(y);
    data.emplace_back(z - size);
    appendLoadColor(data);
    data.emplace_back(x);
    data.emplace_back(y);
    data.emplace_back(z);
    appendLoadColor(data);
}

void ModelRenderer::appendLoadY(std::vector<GLfloat>& data, const QVector3D& position, const float load) const {
    const auto& x = position.x();
    const auto& y = position.y();
    const auto& z = position.z();

    const auto size = load > 0.f ? -Size.LOAD : Size.LOAD;
    const auto size_2 = 2 * size;
    const auto size_6 = 6 * size;

    data.emplace_back(x);
    data.emplace_back(y + size_6);
    data.emplace_back(z);
    appendLoadColor(data);
    data.emplace_back(x);
    data.emplace_back(y);
    data.emplace_back(z);
    appendLoadColor(data);
    data.emplace_back(x + size);
    data.emplace_back(y + size_2);
    data.emplace_back(z);
    appendLoadColor(data);
    data.emplace_back(x - size);
    data.emplace_back(y + size_2);
    data.emplace_back(z);
    appendLoadColor(data);
    data.emplace_back(x);
    data.emplace_back(y);
    data.emplace_back(z);
    appendLoadColor(data);
    data.emplace_back(x);
    data.emplace_back(y + size_2);
    data.emplace_back(z + size);
    appendLoadColor(data);
    data.emplace_back(x);
    data.emplace_back(y + size_2);
    data.emplace_back(z - size);
    appendLoadColor(data);
    data.emplace_back(x);
    data.emplace_back(y);
    data.emplace_back(z);
    appendLoadColor(data);
}

void ModelRenderer::appendLoadZ(std::vector<GLfloat>& data, const QVector3D& position, const float load) const {
    const auto& x = position.x();
    const auto& y = position.y();
    const auto& z = position.z();

    const auto size = load > 0.f ? -Size.LOAD : Size.LOAD;
    const auto size_2 = 2 * size;
    const auto size_6 = 6 * size;

    data.emplace_back(x);
    data.emplace_back(y);
    data.emplace_back(z + size_6);
    appendLoadColor(data);
    data.emplace_back(x);
    data.emplace_back(y);
    data.emplace_back(z);
    appendLoadColor(data);
    data.emplace_back(x);
    data.emplace_back(y + size);
    data.emplace_back(z + size_2);
    appendLoadColor(data);
    data.emplace_back(x);
    data.emplace_back(y - size);
    data.emplace_back(z + size_2);
    appendLoadColor(data);
    data.emplace_back(x);
    data.emplace_back(y);
    data.emplace_back(z);
    appendLoadColor(data);
    data.emplace_back(x + size);
    data.emplace_back(y);
    data.emplace_back(z + size_2);
    appendLoadColor(data);
    data.emplace_back(x - size);
    data.emplace_back(y);
    data.emplace_back(z + size_2);
    appendLoadColor(data);
    data.emplace_back(x);
    data.emplace_back(y);
    data.emplace_back(z);
    appendLoadColor(data);
}

void ModelRenderer::appendLoadRX(std::vector<GLfloat>&, const QVector3D&, float) {}

void ModelRenderer::appendLoadRY(std::vector<GLfloat>&, const QVector3D&, float) {}

void ModelRenderer::appendLoadRZ(std::vector<GLfloat>&, const QVector3D&, float) {}

void ModelRenderer::appendLoadColor(std::vector<GLfloat>& data) const {
    data.emplace_back(Color.LOAD.redF());
    data.emplace_back(Color.LOAD.greenF());
    data.emplace_back(Color.LOAD.blueF());
}

void ModelRenderer::mousePressEvent(QMouseEvent* event) { m_last_pos = event->pos(); }

void ModelRenderer::mouseMoveEvent(QMouseEvent* event) {
    const auto dx = static_cast<float>(event->pos().x() - m_last_pos.x());
    const auto dy = static_cast<float>(event->pos().y() - m_last_pos.y());

    if(event->buttons() & Qt::LeftButton) {
        setViewXR(View.XR + .5f * dy);
        setViewYR(View.YR + .5f * dx);
    }
    else if(event->buttons() & Qt::RightButton) {
        setViewXR(View.XR + .5f * dy);
        setViewZR(View.ZR + .5f * dx);
    }
    else if(event->buttons() & Qt::MiddleButton) {
        setViewXT(View.XT + .002f * dx * View.ZT);
        setViewYT(View.YT + .002f * dy * View.ZT);
    }

    m_last_pos = event->pos();

    repaint();
}

void ModelRenderer::wheelEvent(QWheelEvent* event) {
    setViewZT(View.ZT + View.ZT * static_cast<float>(event->angleDelta().y()) * .002f);

    repaint();
}

void ModelRenderer::setPlane() {
    const auto t_vec = QVector3D{View.XT, View.YT, View.ZT};

    View.NEAR_PLANE = View.FAR_PLANE = t_vec.lengthSquared();

    for(auto& [fst, snd] : model_ptr->getNodePool()) {
        const auto distance = (t_vec - snd.position).lengthSquared();
        if(distance < View.NEAR_PLANE)
            View.NEAR_PLANE = distance;
        else if(distance > View.FAR_PLANE)
            View.FAR_PLANE = distance;
    }

    View.NEAR_PLANE = .5f * std::sqrt(View.NEAR_PLANE);
    View.FAR_PLANE = 1.5f * std::sqrt(View.FAR_PLANE);
}
