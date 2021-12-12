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
#include <QSvgRenderer>
#include <QSvgWidget>
#include "ui_ModelBuilder.h"

ModelBuilder::ModelBuilder(QWidget* parent)
    : QMainWindow(parent)
    , ui(new Ui::ModelBuilder) {
    ui->setupUi(this);

    ui->canvas->setModel(&model);
}

ModelBuilder::~ModelBuilder() { delete ui; }

void ModelBuilder::highlightNode(const QString& text, const int index) {
    const auto tag = text.toInt();

    model.highlight<Database::Node>(highlighted_node.at(index), false);
    highlighted_node[index] = 0;

    if(highlighted_node.count(tag) == 0 || tag == 0) {
        model.highlight<Database::Node>(tag, true);
        highlighted_node[index] = tag;
    }

    ui->canvas->repaint();
}

void ModelBuilder::highlightNodeA(QString text) {
    highlightNode(text, 0);

    if(ui->box_modify_type->currentIndex() == 3 && text.toInt() != 0) {
        const auto coor = model.get<Database::Node>(text.toInt()).position;
        ui->input_modify_node_a->setText(QString::number(coor.x()));
        ui->input_modify_node_b->setText(QString::number(coor.y()));
        ui->input_modify_node_c->setText(QString::number(coor.z()));
    }
}

void ModelBuilder::highlightNodeB(QString text) { highlightNode(text, 1); }

void ModelBuilder::highlightNodeC(QString text) { highlightNode(text, 2); }

void ModelBuilder::highlightNodeD(QString text) { highlightNode(text, 3); }

void ModelBuilder::highlightNodeE() {
    const auto index = ui->box_modify_type->currentIndex();

    const auto text_a = ui->input_modify_node_a->text();
    const auto text_b = ui->input_modify_node_b->text();
    const auto text_c = ui->input_modify_node_c->text();

    for(const auto& I : highlighted_group) model.highlight<Database::Node>(I, false);

    highlighted_group.clear();

    if(1 == index) {
        auto a = text_a.toInt();
        auto b = text_b.toInt();
        auto c = text_c.toInt();
        if(b < 0) b = -b;
        if(c == 0) c = a;
        if(a == 0) a = c;
        if(c <= a) std::swap(a, c);
        if(b == 0) b = std::max(c - a, 1);
        while(a <= c) {
            model.highlight<Database::Node>(a, true);
            highlighted_group.append(a);
            a += b;
        }
    }
    else if(2 == index) {
        if(text_a.size() + text_b.size() + text_c.size() == 0) return;

        const auto& node_pool = model.getNodePool();

        auto I = node_pool.begin();
        while(I != node_pool.end()) {
            const auto flag_x = text_a.size() == 0 || std::fabs(text_a.toFloat() - I->second.position.x()) < 1E-8;
            const auto flag_y = text_b.size() == 0 || std::fabs(text_b.toFloat() - I->second.position.y()) < 1E-8;
            const auto flag_z = text_c.size() == 0 || std::fabs(text_c.toFloat() - I->second.position.z()) < 1E-8;

            if(flag_x && flag_y && flag_z) {
                model.highlight<Database::Node>(I->first, true);
                highlighted_group.append(I->first);
            }

            ++I;
        }
    }

    ui->canvas->repaint();
}

void ModelBuilder::highlightElement(const QString& text, const int index) {
    const auto tag = text.toInt();

    model.highlight<Database::Element>(highlighted_element.at(index), false);
    highlighted_element[index] = 0;

    if(highlighted_element.count(tag) == 0 || tag == 0) {
        model.highlight<Database::Element>(tag, true);
        highlighted_element[index] = tag;
    }

    ui->canvas->repaint();
}

void ModelBuilder::highlightElementA(QString text) { highlightElement(text, 0); }

void ModelBuilder::writeOutput() {
    QFileDialog dialog(this);
    dialog.setFileMode(QFileDialog::AnyFile);
    dialog.setAcceptMode(QFileDialog::AcceptSave);
    if(dialog.exec()) {
        const auto filename = dialog.selectedFiles();
        if(1 == filename.size() && model.saveModel(filename.at(0))) {
            saved_file_name = filename.at(0);

            saved = true;
        }
    }

    ui->canvas->repaint();
}

void ModelBuilder::saveScreenshot() {
    QFileDialog dialog(this);
    dialog.setFileMode(QFileDialog::AnyFile);
    dialog.setAcceptMode(QFileDialog::AcceptSave);
    if(dialog.exec()) {
        const auto filename = dialog.selectedFiles();
        if(1 == filename.size()) {
            auto path = filename.at(0);
            if(!path.endsWith(".png") && !path.endsWith(".PNG")) path.append(".png");
            if(!ui->canvas->grabFramebuffer().save(path, "PNG")) {
                QMessageBox msg(QMessageBox::Critical, tr("Error"), tr("Fail to save file."), QMessageBox::Ok, this);
                msg.exec();
            }
        }
    }
}

void ModelBuilder::showAbout() {
    QDialog about(this);

    about.setLayout(new QHBoxLayout(&about));
    about.layout()->setSpacing(30);

    QSvgWidget uc_logo(":res/UCBLACK.svg");
    uc_logo.renderer()->setAspectRatioMode(Qt::KeepAspectRatio);
    about.layout()->addWidget(&uc_logo);

    QLabel label(&about);

    const auto license =
        "Copyright (C) 2021 Theodore Chang, Minghao Li\n\n"
        "This program is free software: you can redistribute it and/or modify\n"
        "it under the terms of the GNU General Public License as published by\n"
        "the Free Software Foundation, either version 3 of the License, or\n"
        "(at your option) any later version.\n\n"
        "This program is distributed in the hope that it will be useful,\n"
        "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
        "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n"
        "GNU General Public License for more details.\n\n"
        "You should have received a copy of the GNU General Public License\n"
        "along with this program. If not, see <http://www.gnu.org/licenses/>.\n\n"
        "PB3D is a simple application to visualise frame models\n"
        "and generate model input files for finite element analysis.\n\n"
        "The program is written and compiled with " +
        tr("Qt version %1.\n").arg(qVersion()) +
        "The program is licensed under GNU GPLv3.\n"
        "The source code can be acquired via GitHub repository.\n\n"
        "The software is written and maintained by\n"
        "\tDr. Minghao Li\t\tminghao.li@canterbury.ac.nz\n"
        "\tDr. Theodore Chang\t\ttlcfem@gmail.com";

    label.setText(license);

    auto font = label.font();
    font.setPixelSize(13);
    label.setFont(font);

    about.layout()->addWidget(&label);

    about.adjustSize();
    about.setFixedSize(about.size());
    about.exec();
}

void ModelBuilder::openFile() {
    on_reset_model_clicked();

    QFileDialog dialog(this);
    dialog.setFileMode(QFileDialog::AnyFile);
    if(dialog.exec()) {
        const auto filename = dialog.selectedFiles();
        if(1 == filename.size() && !model.loadModel(filename.at(0))) {
            QMessageBox msg(this);
            msg.setText(tr("Fail to read file %1.\n").arg(filename.at(0)) + "Please make sure the input file is correct.\nOtherwise contact the authors.\n");
            msg.exec();
        }
        else
            saved_file_name = filename.at(0);
    }

    ui->input_node_tag->setValue(model.getNextNodeTag());
    ui->input_element_tag->setValue(model.getNextElementTag());

    updateNodeList();
    updateFrameSectionList();
    updateWallSectionList();
    updateElementList();
    updateAnalysisSetting();

    ui->canvas->repaint();

    saved = true;
}

void ModelBuilder::updateNodeList() const {
    ui->box_node->clear();
    ui->box_node_i->clear();
    ui->box_node_j->clear();
    ui->box_node_load->clear();
    ui->box_node->addItem("");
    ui->box_node_i->addItem("");
    ui->box_node_j->addItem("");
    ui->box_node_load->addItem("");
    for(const auto& I : model.getNodeTag()) {
        ui->box_node->addItem(QString::number(I));
        ui->box_node_i->addItem(QString::number(I));
        ui->box_node_j->addItem(QString::number(I));
        ui->box_node_load->addItem(QString::number(I));
    }

    ui->canvas->repaint();
}

void ModelBuilder::updateElementList() const {
    ui->box_element->clear();
    ui->box_element->addItem("");

    for(const auto& I : model.getElementTag()) ui->box_element->addItem(QString::number(I));

    ui->canvas->repaint();
}

void ModelBuilder::updateAnalysisSetting() const {
    ui->box_unit->setCurrentIndex(model.unit_system - 1);

    if(1 != ui->box_unit->currentIndex()) {
        ui->slider_bc_size->setValue(5000);
        ui->slider_axis_size->setValue(100000);
        ui->slider_grid_size->setValue(17);
    }

    ui->input_damping->setText(QString::number(model.damping_ratio));
    ui->input_period->setText(QString::number(model.natural_period));
    ui->input_scale->setText(QString::number(model.scale_factor));
    ui->input_qfx->setValue(model.quadrature_frame.at(0));
    ui->input_qfy->setValue(model.quadrature_frame.at(1));
    ui->input_qfz->setValue(model.quadrature_frame.at(2));
    ui->input_qwx->setValue(model.quadrature_wall.at(0));
    ui->input_qwy->setValue(model.quadrature_wall.at(1));

    if(!model.acc_record.at(0).isEmpty()) {
        ui->check_accx->setChecked(true);
        ui->input_accx->setText(model.acc_record.at(0));
    }
    if(!model.acc_record.at(1).isEmpty()) {
        ui->check_accy->setChecked(true);
        ui->input_accy->setText(model.acc_record.at(1));
    }

    ui->input_relf->setText(QString::number(model.tolerance.at(0)));
    ui->input_relx->setText(QString::number(model.tolerance.at(1)));
    ui->input_relu->setText(QString::number(model.tolerance.at(2)));
    ui->input_absf->setText(QString::number(model.tolerance.at(3)));
    ui->input_absx->setText(QString::number(model.tolerance.at(4)));
    ui->input_absu->setText(QString::number(model.tolerance.at(5)));
}

void ModelBuilder::updateWallSectionList() {
    ui->input_wall_section_tag->setValue(model.getNextWallSectionTag());

    ui->box_wall_section->clear();
    for(const auto& I : model.getWallSectionTag()) { ui->box_wall_section->addItem(QString::number(I)); }

    if(ui->box_element_type->currentText() == "Wall") {
        ui->box_section->clear();
        for(const auto& I : model.getWallSectionTag()) { ui->box_section->addItem(QString::number(I)); }
    }

    ui->canvas->repaint();
}

void ModelBuilder::updateFrameSectionList() {
    ui->input_frame_section_tag->setValue(model.getNextFrameSectionTag());

    ui->box_frame_section->clear();
    for(const auto& I : model.getFrameSectionTag()) { ui->box_frame_section->addItem(QString::number(I)); }

    if(ui->box_element_type->currentText() != "Wall") {
        ui->box_section->clear();
        for(const auto& I : model.getFrameSectionTag()) { ui->box_section->addItem(QString::number(I)); }
    }

    ui->canvas->repaint();
}

void ModelBuilder::on_box_unit_currentIndexChanged(int index) {
    model.changeUnit(index + 1);

    saved = false;
}

void ModelBuilder::on_input_scale_textChanged(const QString& F) {
    model.changeScale(F);

    saved = false;
}

void ModelBuilder::on_input_damping_textChanged(const QString& F) {
    model.changeDamping(F);

    saved = false;
}

void ModelBuilder::on_input_period_textChanged(const QString& F) {
    model.changePeriod(F);

    saved = false;
}

void ModelBuilder::on_input_qfx_textChanged(const QString& qfx) {
    if(qfx.toInt() > 5) ui->input_qfx->setValue(5);
    model.quadrature_frame[0] = std::min(5, qfx.toInt());

    saved = false;
}

void ModelBuilder::on_input_qfy_textChanged(const QString& qfy) {
    model.quadrature_frame[1] = qfy.toInt();

    saved = false;
}

void ModelBuilder::on_input_qfz_textChanged(const QString& qfz) {
    model.quadrature_frame[2] = qfz.toInt();

    saved = false;
}

void ModelBuilder::on_input_qwx_textChanged(const QString& qwx) {
    model.quadrature_wall[0] = qwx.toInt();

    saved = false;
}

void ModelBuilder::on_input_qwy_textChanged(const QString& qwy) {
    model.quadrature_wall[1] = qwy.toInt();

    saved = false;
}

void ModelBuilder::on_box_modify_type_currentIndexChanged(const int type) const {
    ui->box_node->setCurrentIndex(0);
    ui->input_modify_node_a->setText("");
    ui->input_modify_node_b->setText("");
    ui->input_modify_node_c->setText("");
    if(0 == type) {
        ui->label_modify_node_a->setText("");
        ui->label_modify_node_b->setText("");
        ui->label_modify_node_c->setText("");
        ui->label_modify_node_a->setDisabled(true);
        ui->label_modify_node_b->setDisabled(true);
        ui->label_modify_node_c->setDisabled(true);
        ui->input_modify_node_a->setDisabled(true);
        ui->input_modify_node_b->setDisabled(true);
        ui->input_modify_node_c->setDisabled(true);
        ui->label_valid_node->setDisabled(false);
        ui->box_node->setDisabled(false);
    }
    else {
        ui->label_modify_node_a->setDisabled(false);
        ui->label_modify_node_b->setDisabled(false);
        ui->label_modify_node_c->setDisabled(false);
        ui->input_modify_node_a->setDisabled(false);
        ui->input_modify_node_b->setDisabled(false);
        ui->input_modify_node_c->setDisabled(false);
        if(1 == type) {
            ui->label_modify_node_a->setText("Start");
            ui->label_modify_node_b->setText("Interval");
            ui->label_modify_node_c->setText("End");
            ui->label_valid_node->setDisabled(true);
            ui->box_node->setDisabled(true);
        }
        else if(2 == type) {
            ui->label_modify_node_a->setText("A");
            ui->label_modify_node_b->setText("B");
            ui->label_modify_node_c->setText("C");
            ui->label_valid_node->setDisabled(true);
            ui->box_node->setDisabled(true);
        }
        else if(3 == type) {
            ui->label_modify_node_a->setText("New X");
            ui->label_modify_node_b->setText("New Y");
            ui->label_modify_node_c->setText("New Z");
            ui->label_valid_node->setDisabled(false);
            ui->box_node->setDisabled(false);
        }
    }
}

void ModelBuilder::on_button_split_element_clicked() {
    const auto tag = ui->box_element->currentText().toInt();
    const auto segment = ui->input_split->text().toInt();

    if(segment < 2) return;

    model.splitElement(tag, segment);

    ui->input_node_tag->setValue(model.getNextNodeTag());
    updateNodeList();

    ui->input_element_tag->setValue(model.getNextElementTag());
    updateElementList();

    saved = false;
}

void ModelBuilder::on_button_remove_wall_section_clicked() {
    const auto tag = ui->box_wall_section->currentText().toInt();

    model.removeWallSection(tag);

    updateWallSectionList();

    saved = false;
}

void ModelBuilder::on_button_remove_frame_section_clicked() {
    const auto tag = ui->box_frame_section->currentText().toInt();

    model.removeFrameSection(tag);

    updateFrameSectionList();

    saved = false;
}

void ModelBuilder::on_button_remove_element_clicked() {
    const auto tag = ui->box_element->currentText().toInt();

    model.removeElement(tag);

    updateElementList();

    saved = false;
}

void ModelBuilder::on_button_remove_all_element_clicked() {
    model.removeElement();

    ui->input_element_tag->setValue(1);
    updateElementList();

    saved = false;
}

void ModelBuilder::on_button_modify_node_clicked() {
    const auto index = ui->box_modify_type->currentIndex();
    const auto tag = ui->box_node->currentText().toInt();

    if(0 == index) { model.removeNode(tag); }
    else if(1 == index || 2 == index) {
        for(const auto& I : highlighted_group) model.removeNode(I);
    }
    else if(3 == index) {
        const auto a = ui->input_modify_node_a->text().toFloat();
        const auto b = ui->input_modify_node_b->text().toFloat();
        const auto c = ui->input_modify_node_c->text().toFloat();
        model.get<Database::Node>(tag).position = QVector3D{a, b, c};
    }

    ui->input_modify_node_a->setText("");
    ui->input_modify_node_b->setText("");
    ui->input_modify_node_c->setText("");

    updateNodeList();

    saved = false;
}

void ModelBuilder::on_button_add_node_clicked() {
    const auto tag = ui->input_node_tag->value();
    const auto x = ui->input_x->text().toFloat();
    const auto y = ui->input_y->text().toFloat();
    const auto z = ui->input_z->text().toFloat();
    const auto dx = ui->input_dx->text().toFloat();
    const auto dy = ui->input_dy->text().toFloat();
    const auto dz = ui->input_dz->text().toFloat();
    const auto nx = ui->input_nx->value();
    const auto ny = ui->input_ny->value();
    const auto nz = ui->input_nz->value();

    if(0 == tag) return;

    for(auto I = 0; I < nx; ++I)
        for(auto J = 0; J < ny; ++J)
            for(auto K = 0; K < nz; ++K) model.add(model.getNextNodeTag(), Database::Node{QVector3D{x + dx * static_cast<float>(I), y + dy * static_cast<float>(J), z + dz * static_cast<float>(K)}});

    ui->input_node_tag->setValue(model.getNextNodeTag());

    updateNodeList();

    saved = false;
}

void ModelBuilder::on_button_change_section_clicked() {
    const auto tag = ui->box_element->currentText().toInt();
    const auto sec_tag = ui->box_section_2->currentText().toInt();

    model.changeSection(tag, sec_tag);

    ui->canvas->repaint();

    saved = false;
}

void ModelBuilder::on_button_clear_bc_clicked() {
    for(auto& [fst, snd] : model.getNodePool()) model.changeFixity(fst, QVector<bool>(6, false));

    ui->box_node_load->setCurrentIndex(0);

    ui->canvas->repaint();

    saved = false;
}

void ModelBuilder::on_button_clear_load_clicked() {
    const auto type = ui->box_load_type->currentText();

    if(type == "Mass")
        for(auto& [fst, snd] : model.getNodePool()) model.changeMass(fst, 0.);
    else if(type == "Displacement")
        for(auto& [fst, snd] : model.getNodePool()) model.changeDisplacement(fst, QVector<double>(6, 0.));
    else
        for(auto& [fst, snd] : model.getNodePool()) model.changeLoad(fst, QVector<double>(6, 0.));

    ui->box_node_load->setCurrentIndex(0);

    ui->canvas->repaint();

    saved = false;
}

void ModelBuilder::on_button_add_bc_clicked() {
    const auto tag = ui->box_node_load->currentText().toInt();
    const auto increx = ui->input_bc_increx->text().toInt();
    const auto increy = ui->input_bc_increy->text().toInt();
    const auto increz = ui->input_bc_increz->text().toInt();
    const auto repeatx = ui->input_bc_repeatx->text().toInt();
    const auto repeaty = ui->input_bc_repeaty->text().toInt();
    const auto repeatz = ui->input_bc_repeatz->text().toInt();

    const auto x = ui->box_x->checkState() == Qt::Checked;
    const auto y = ui->box_y->checkState() == Qt::Checked;
    const auto z = ui->box_z->checkState() == Qt::Checked;
    const auto rx = ui->box_rx->checkState() == Qt::Checked;
    const auto ry = ui->box_ry->checkState() == Qt::Checked;
    const auto rz = ui->box_rz->checkState() == Qt::Checked;

    for(auto I = 0; I < repeatx; ++I)
        for(auto J = 0; J < repeaty; ++J)
            for(auto K = 0; K < repeatz; ++K) model.changeFixity(tag + I * increx + J * increy + K * increz, QVector<bool>{x, y, z, rx, ry, rz});

    ui->box_node_load->setCurrentIndex(0);

    ui->canvas->repaint();

    saved = false;
}

void ModelBuilder::on_button_add_load_clicked() {
    const auto tag = ui->box_node_load->currentText().toInt();
    const auto increx = ui->input_bc_increx->text().toInt();
    const auto increy = ui->input_bc_increy->text().toInt();
    const auto increz = ui->input_bc_increz->text().toInt();
    const auto repeatx = ui->input_bc_repeatx->text().toInt();
    const auto repeaty = ui->input_bc_repeaty->text().toInt();
    const auto repeatz = ui->input_bc_repeatz->text().toInt();
    const auto type = ui->box_load_type->currentText();

    const auto x = ui->input_loadx->text().toDouble();
    const auto y = ui->input_loady->text().toDouble();
    const auto z = ui->input_loadz->text().toDouble();
    const auto rx = ui->input_loadrx->text().toDouble();
    const auto ry = ui->input_loadry->text().toDouble();
    const auto rz = ui->input_loadrz->text().toDouble();

    if(type == "Mass")
        for(auto I = 0; I < repeatx; ++I)
            for(auto J = 0; J < repeaty; ++J)
                for(auto K = 0; K < repeatz; ++K) model.changeMass(tag + I * increx + J * increy + K * increz, x);
    else if(type == "Force")
        for(auto I = 0; I < repeatx; ++I)
            for(auto J = 0; J < repeaty; ++J)
                for(auto K = 0; K < repeatz; ++K) model.changeLoad(tag + I * increx + J * increy + K * increz, QVector<double>{x, y, z, rx, ry, rz});
    else if(type == "Displacement")
        for(auto I = 0; I < repeatx; ++I)
            for(auto J = 0; J < repeaty; ++J)
                for(auto K = 0; K < repeatz; ++K) model.changeDisplacement(tag + I * increx + J * increy + K * increz, QVector<double>{x, y, z, rx, ry, rz});

    ui->box_node_load->setCurrentIndex(0);

    ui->canvas->repaint();

    saved = false;
}

void ModelBuilder::on_button_add_element_clicked() {
    const auto tag = ui->input_element_tag->text().toInt();
    const auto sec_tag = ui->box_section->currentText().toInt();
    const auto nodei_tag = ui->box_node_i->currentText().toInt();
    const auto nodej_tag = ui->box_node_j->currentText().toInt();
    const auto increix = ui->input_increix->text().toInt();
    const auto increjx = ui->input_increjx->text().toInt();
    const auto increiy = ui->input_increiy->text().toInt();
    const auto increjy = ui->input_increjy->text().toInt();
    const auto increiz = ui->input_increiz->text().toInt();
    const auto increjz = ui->input_increjz->text().toInt();
    const auto repeati = ui->input_repeati->value();
    const auto repeatj = ui->input_repeatj->value();
    const auto repeatk = ui->input_repeatk->value();
    const auto type = ui->box_element_type->currentText();
    const auto orient = type == "Wall" ? ui->box_orient->currentIndex() + 1 : 0;

    if(tag == 0 || nodei_tag == 0 || nodej_tag == 0 || sec_tag == 0) return;

    for(auto I = 0; I < repeati; ++I)
        for(auto J = 0; J < repeatj; ++J)
            for(auto K = 0; K < repeatk; ++K) {
                const auto new_i = nodei_tag + I * increix + J * increiy + K * increiz;
                const auto new_j = nodej_tag + I * increjx + J * increjy + K * increjz;

                model.add(model.getNextElementTag(), Database::Element(sec_tag, QVector<int>{new_i, new_j}, type, orient));
            }

    ui->input_element_tag->setValue(model.getNextElementTag());

    updateElementList();

    ui->box_node_i->setCurrentIndex(0);
    ui->box_node_j->setCurrentIndex(0);

    saved = false;
}

void ModelBuilder::on_button_add_wall_section_clicked() {
    const auto tag = ui->input_wall_section_tag->text().toInt();

    if(tag == 0) return;

    const auto l = ui->input_ls->text().toDouble();
    const auto d = ui->input_ds->text().toDouble();
    const auto e = ui->input_es->text().toDouble();
    const auto ys = ui->input_ys->text().toDouble();
    const auto th = ui->input_th->text().toDouble();
    const auto dl = ui->input_dl->text().toDouble();

    const auto q0t = ui->input_q0t->text().toDouble();
    const auto q1t = ui->input_q1t->text().toDouble();
    const auto q2t = ui->input_q2t->text().toDouble();
    const auto xkt = ui->input_xkt->text().toDouble();
    const auto dmaxt = ui->input_dmaxt->text().toDouble();
    const auto sdft = ui->input_sdft->text().toDouble();

    const auto q0p = ui->input_q0p->text().toDouble();
    const auto q1p = ui->input_q1p->text().toDouble();
    const auto q2p = ui->input_q2p->text().toDouble();
    const auto xkp = ui->input_xkp->text().toDouble();
    const auto dmaxp = ui->input_dmaxp->text().toDouble();
    const auto sdfp = ui->input_sdfp->text().toDouble();

    model.add<Database::WallSection>(tag, Database::WallSection{QVector<double>{l, d, e, ys, th, dl, q0t, q1t, q2t, xkt, dmaxt, sdft, q0p, q1p, q2p, xkp, dmaxp, sdfp}});

    ui->input_wall_section_tag->setValue(model.getNextWallSectionTag());

    updateWallSectionList();

    saved = false;
}

void ModelBuilder::on_button_add_frame_section_clicked() {
    const auto tag = ui->input_frame_section_tag->text().toInt();

    if(tag == 0) return;

    const auto elastic_modulus = ui->input_e->text().toDouble();
    const auto shear_modulus = ui->input_g->text().toDouble();
    const auto w = ui->input_w->text().toDouble();
    const auto h = ui->input_h->text().toDouble();
    const auto type = ui->box_material_type->currentText();

    model.add<Database::FrameSection>(tag, Database::FrameSection{type, QVector<double>{elastic_modulus, shear_modulus, w, h}});

    ui->input_frame_section_tag->setValue(model.getNextFrameSectionTag());

    updateFrameSectionList();

    saved = false;
}

void ModelBuilder::on_reset_model_clicked() {
    model = Database();

    updateNodeList();
    updateElementList();
    updateFrameSectionList();
    updateWallSectionList();

    ui->input_node_tag->setValue(1);
    ui->input_frame_section_tag->setValue(1);
    ui->input_wall_section_tag->setValue(1);
    ui->input_element_tag->setValue(1);
    ui->label_section_info->clear();
    ui->label_section_info_2->clear();

    saved = false;
}

void ModelBuilder::on_box_element_type_currentTextChanged(const QString& F) {
    if(F == "Wall") {
        updateWallSectionList();
        ui->box_orient->setEnabled(true);
    }
    else {
        updateFrameSectionList();
        ui->box_orient->setEnabled(false);
    }
}

void ModelBuilder::on_box_element_currentTextChanged(const QString& F) {
    ui->box_section_2->clear();

    const auto ele_tag = F.toInt();

    if(ele_tag == 0) return;

    const auto ele_type = model.get<Database::Element>(F.toInt()).type;

    for(const auto& I : (ele_type == Database::Element::Type::Wall ? model.getWallSectionTag() : model.getFrameSectionTag())) { ui->box_section_2->addItem(QString::number(I)); }
}

void ModelBuilder::on_box_load_type_currentTextChanged(const QString& F) const {
    if(F == "Mass") {
        ui->label_loadx->setText("Mass");
        ui->label_loady->setText("");
        ui->label_loadz->setText("");
        ui->label_loadrx->setText("");
        ui->label_loadry->setText("");
        ui->label_loadrz->setText("");
        ui->input_loady->setDisabled(true);
        ui->input_loadz->setDisabled(true);
        ui->input_loadrx->setDisabled(true);
        ui->input_loadry->setDisabled(true);
        ui->input_loadrz->setDisabled(true);
        ui->label_loady->setDisabled(true);
        ui->label_loadz->setDisabled(true);
        ui->label_loadrx->setDisabled(true);
        ui->label_loadry->setDisabled(true);
        ui->label_loadrz->setDisabled(true);
    }
    else {
        ui->label_loadx->setText("X");
        ui->label_loady->setText("Y");
        ui->label_loadz->setText("Z");
        ui->label_loadrx->setText("RX");
        ui->label_loadry->setText("RY");
        ui->label_loadrz->setText("RZ");
        ui->input_loady->setDisabled(false);
        ui->input_loadz->setDisabled(false);
        ui->input_loadrx->setDisabled(false);
        ui->input_loadry->setDisabled(false);
        ui->input_loadrz->setDisabled(false);
        ui->label_loady->setDisabled(false);
        ui->label_loadz->setDisabled(false);
        ui->label_loadrx->setDisabled(false);
        ui->label_loadry->setDisabled(false);
        ui->label_loadrz->setDisabled(false);
    }
}

void ModelBuilder::on_input_relf_textChanged(const QString& t) {
    model.tolerance[0] = t.toDouble();

    saved = false;
}

void ModelBuilder::on_input_relx_textChanged(const QString& t) {
    model.tolerance[1] = t.toDouble();

    saved = false;
}

void ModelBuilder::on_input_relu_textChanged(const QString& t) {
    model.tolerance[2] = t.toDouble();

    saved = false;
}

void ModelBuilder::on_input_absf_textChanged(const QString& t) {
    model.tolerance[3] = t.toDouble();

    saved = false;
}

void ModelBuilder::on_input_absx_textChanged(const QString& t) {
    model.tolerance[4] = t.toDouble();

    saved = false;
}

void ModelBuilder::on_input_absu_textChanged(const QString& t) {
    model.tolerance[5] = t.toDouble();

    saved = false;
}

void ModelBuilder::on_box_section_textHighlighted(const QString& F) {
    ui->label_section_info->clear();

    const auto sec_tag = F.toInt();

    if(sec_tag == 0) return;

    QString text;

    if(ui->box_element_type->currentText() != "Wall") {
        const auto& t_para = model.get<Database::FrameSection>(sec_tag).parameter;
        text = tr("Selected Frame Section Info:\nTag:\t%1\nParameters:\n\t%2\t%3\t%4\t%5").arg(sec_tag).arg(t_para.at(0)).arg(t_para.at(1)).arg(t_para.at(2)).arg(t_para.at(3));
    }
    else {
        const auto& t_para = model.get<Database::WallSection>(sec_tag).parameter;
        text = tr("Selected Wall Section Info:\nTag:\t%1\nParameters:\n\t%2\t%3\t%4\t%5\t%6\t%7").arg(sec_tag).arg(t_para.at(0)).arg(t_para.at(1)).arg(t_para.at(2)).arg(t_para.at(3)).arg(t_para.at(4)).arg(t_para.at(5));
        text += tr("\n\t%1\t%2\t%3\t%4\t%5\t%6").arg(t_para.at(6)).arg(t_para.at(7)).arg(t_para.at(8)).arg(t_para.at(9)).arg(t_para.at(10)).arg(t_para.at(11));
        text += tr("\n\t%1\t%2\t%3\t%4\t%5\t%6").arg(t_para.at(12)).arg(t_para.at(13)).arg(t_para.at(14)).arg(t_para.at(15)).arg(t_para.at(16)).arg(t_para.at(17));
    }

    ui->label_section_info->setText(text);
}

void ModelBuilder::on_box_section_2_textHighlighted(const QString& F) {
    ui->label_section_info->clear();

    const auto sec_tag = F.toInt();

    if(sec_tag == 0) return;

    const auto ele_type = model.get<Database::Element>(ui->box_element->currentText().toInt()).type;

    QString text;

    if(ele_type != Database::Element::Type::Wall) {
        const auto& t_para = model.get<Database::FrameSection>(sec_tag).parameter;
        text = tr("Selected Frame Section Info:\nTag:\t%1\nParameters:\n\t%2\t%3\t%4\t%5").arg(sec_tag).arg(t_para.at(0)).arg(t_para.at(1)).arg(t_para.at(2)).arg(t_para.at(3));
    }
    else {
        const auto& t_para = model.get<Database::WallSection>(sec_tag).parameter;
        text = tr("Selected Wall Section Info:\nTag:\t%1\nParameters:\n\t%2\t%3\t%4\t%5\t%6\t%7").arg(sec_tag).arg(t_para.at(0)).arg(t_para.at(1)).arg(t_para.at(2)).arg(t_para.at(3)).arg(t_para.at(4)).arg(t_para.at(5));
        text += tr("\n\t%1\t%2\t%3\t%4\t%5\t%6").arg(t_para.at(6)).arg(t_para.at(7)).arg(t_para.at(8)).arg(t_para.at(9)).arg(t_para.at(10)).arg(t_para.at(11));
        text += tr("\n\t%1\t%2\t%3\t%4\t%5\t%6").arg(t_para.at(12)).arg(t_para.at(13)).arg(t_para.at(14)).arg(t_para.at(15)).arg(t_para.at(16)).arg(t_para.at(17));
    }

    ui->label_section_info->setText(text);
}

void ModelBuilder::on_input_wall_section_tag_textChanged(const QString&) const {
    ui->input_ls->clear();
    ui->input_ds->clear();
    ui->input_es->clear();
    ui->input_ys->clear();
    ui->input_th->clear();
    ui->input_dl->clear();

    ui->input_q0t->clear();
    ui->input_q1t->clear();
    ui->input_q2t->clear();
    ui->input_xkt->clear();
    ui->input_dmaxt->clear();
    ui->input_sdft->clear();

    ui->input_q0p->clear();
    ui->input_q1p->clear();
    ui->input_q2p->clear();
    ui->input_xkp->clear();
    ui->input_dmaxp->clear();
    ui->input_sdfp->clear();
}

void ModelBuilder::on_input_frame_section_tag_textChanged(const QString&) const {
    ui->input_e->clear();
    ui->input_g->clear();
    ui->input_w->clear();
    ui->input_h->clear();
}

void ModelBuilder::on_box_frame_section_textHighlighted(const QString& F) {
    const auto& t_para = model.get<Database::FrameSection>(F.toInt()).parameter;
    const QString text = tr("Selected Frame Section Info:\nTag:\t%1\nParameters:\n\t%2\t%3\t%4\t%5").arg(F.toInt()).arg(t_para.at(0)).arg(t_para.at(1)).arg(t_para.at(2)).arg(t_para.at(3));

    ui->label_section_info_2->setText(text);
}

void ModelBuilder::on_box_wall_section_textHighlighted(const QString& F) {

    const auto& t_para = model.get<Database::WallSection>(F.toInt()).parameter;
    auto text = tr("Selected Wall Section Info:\nTag:\t%1\nParameters:\n\t%2\t%3\t%4\t%5\t%6\t%7").arg(F.toInt()).arg(t_para.at(0)).arg(t_para.at(1)).arg(t_para.at(2)).arg(t_para.at(3)).arg(t_para.at(4)).arg(t_para.at(5));
    text += tr("\n\t%1\t%2\t%3\t%4\t%5\t%6").arg(t_para.at(6)).arg(t_para.at(7)).arg(t_para.at(8)).arg(t_para.at(9)).arg(t_para.at(10)).arg(t_para.at(11));
    text += tr("\n\t%1\t%2\t%3\t%4\t%5\t%6").arg(t_para.at(12)).arg(t_para.at(13)).arg(t_para.at(14)).arg(t_para.at(15)).arg(t_para.at(16)).arg(t_para.at(17));

    ui->label_section_info_2->setText(text);
}

void ModelBuilder::on_box_element_textHighlighted(const QString& F) {
    const auto ele_tag = F.toInt();

    if(ele_tag == 0) return;

    const auto& t_ele = model.get<Database::Element>(ele_tag);

    auto text = tr("Selected Element Info:\nTag:\t%1").arg(ele_tag);
    text += tr("\nConnects nodes:\t%1\t%2").arg(t_ele.encoding.at(0)).arg(t_ele.encoding.at(1));

    if(t_ele.type == Database::Element::Type::Wall) {
        const auto& t_para = model.get<Database::WallSection>(t_ele.section_tag).parameter;
        text += tr("\nWall Section Info:\nTag:\t%1\nParameters:\n\t%2\t%3\t%4\t%5\t%6\t%7").arg(t_ele.section_tag).arg(t_para.at(0)).arg(t_para.at(1)).arg(t_para.at(2)).arg(t_para.at(3)).arg(t_para.at(4)).arg(t_para.at(5));
        text += tr("\n\t%1\t%2\t%3\t%4\t%5\t%6").arg(t_para.at(6)).arg(t_para.at(7)).arg(t_para.at(8)).arg(t_para.at(9)).arg(t_para.at(10)).arg(t_para.at(11));
        text += tr("\n\t%1\t%2\t%3\t%4\t%5\t%6").arg(t_para.at(12)).arg(t_para.at(13)).arg(t_para.at(14)).arg(t_para.at(15)).arg(t_para.at(16)).arg(t_para.at(17));
        text += tr("\nSection Orientation:\t%1").arg(t_ele.orient == 1 ? "X" : "Y");
    }
    else {
        const auto& t_para = model.get<Database::FrameSection>(t_ele.section_tag).parameter;
        text += tr("\nFrame Section Info:\nTag:\t%1\nParameters:\n\t%2\t%3\t%4\t%5").arg(t_ele.section_tag).arg(t_para.at(0)).arg(t_para.at(1)).arg(t_para.at(2)).arg(t_para.at(3));
    }

    ui->label_section_info->setText(text);
}

void ModelBuilder::on_check_accx_clicked(const bool checked) {
    if(!checked) return;

    QFileDialog dialog(this);
    dialog.setFileMode(QFileDialog::ExistingFile);
    if(dialog.exec()) {
        const auto filename = dialog.selectedFiles();
        if(1 == filename.size()) {
            const auto name = filename.at(0).split("/").last();
            ui->input_accx->setText(name);
            model.changeAccxRecord(name);

            saved = false;
        }
    }
}

void ModelBuilder::on_check_accy_clicked(const bool checked) {
    if(!checked) return;

    QFileDialog dialog(this);
    dialog.setFileMode(QFileDialog::ExistingFile);
    if(dialog.exec()) {
        const auto filename = dialog.selectedFiles();
        if(1 == filename.size()) {
            const auto name = filename.at(0).split("/").last();
            ui->input_accy->setText(name);
            model.changeAccyRecord(name);

            saved = false;
        }
    }
}

void ModelBuilder::on_button_select_clicked() {
    QFileDialog dialog(this);
    dialog.setFileMode(QFileDialog::ExistingFile);
    if(dialog.exec()) {
        const auto filename = dialog.selectedFiles();
        if(1 == filename.size()) ui->label_exe->setText(filename.at(0));
    }
}

void ModelBuilder::on_button_run_clicked() {
    auto path = ui->label_exe->text();

    if(path.isEmpty()) {
        path = QDir::currentPath() + "/PB3D.CORE.exe";
        QFile exe_file(path);
        if(!exe_file.exists())
            path.clear();
        else {
            ui->label_exe->setText(path);
            ui->label_exe->update();
        }
    }

    if(path.isEmpty() || (!path.endsWith(".exe") && !path.endsWith(".EXE"))) {
        QMessageBox msg(QMessageBox::Critical, tr("Error"), tr("Please selectc the correct executable first."), QMessageBox::Ok, this);
        msg.exec();
    }
    else {
        if(saved_file_name.isEmpty()) {
            QMessageBox msg(QMessageBox::Information, tr("Note"), tr("Please save the model first."), QMessageBox::Ok, this);
            msg.exec();
            writeOutput();
        }
        if(saved_file_name.isEmpty()) {
            QMessageBox msg(QMessageBox::Critical, tr("Error"), tr("Fail to load the file."), QMessageBox::Ok, this);
            msg.exec();
        }
        else if(saved) {
            QDir pwd(saved_file_name);
            pwd.cdUp();
            path.prepend(" && ");
            path.prepend(pwd.absolutePath());
            path.prepend("/C cd ");
            path.append(' ');
            path.append(saved_file_name);
            QProcess process;
            process.startDetached("cmd.exe", QStringList() << path);
        }
    }
}
