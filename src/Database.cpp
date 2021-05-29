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

#include "Database.h"
#include <QFile>
#include <QRegularExpression>

const std::unordered_map<int, Database::Node>& Database::getNodePool() const { return node_pool; }

const std::unordered_map<int, Database::WallSection>& Database::getWallSectionPool() const { return wall_section_pool; }

const std::unordered_map<int, Database::FrameSection>& Database::getFrameSectionPool() const { return frame_section_pool; }

const std::unordered_map<int, Database::Element>& Database::getElementPool() const { return element_pool; }

QVector<int> Database::getNodeTag() const {
	QVector<int> pool;

	pool.reserve(static_cast<int>(node_pool.size()));
	for(auto& [fst, snd] : node_pool) pool.append(fst);
	std::sort(pool.begin(), pool.end());

	return pool;
}

QVector<int> Database::getWallSectionTag() const {
	QVector<int> pool;

	pool.reserve(static_cast<int>(wall_section_pool.size()));
	for(auto& [fst, snd] : wall_section_pool) pool.append(fst);
	std::sort(pool.begin(), pool.end());

	return pool;
}

QVector<int> Database::getFrameSectionTag() const {
	QVector<int> pool;

	pool.reserve(static_cast<int>(frame_section_pool.size()));
	for(auto& [fst, snd] : frame_section_pool) pool.append(fst);
	std::sort(pool.begin(), pool.end());

	return pool;
}

QVector<int> Database::getElementTag() const {
	QVector<int> pool;

	pool.reserve(static_cast<int>(element_pool.size()));
	for(auto& [fst, snd] : element_pool) pool.append(fst);
	std::sort(pool.begin(), pool.end());

	return pool;
}

int Database::getNextNodeTag() {
	auto tag = 0;

	for(auto& [fst, snd] : node_pool) if(fst > tag) tag = fst;

	return tag + 1;
}

int Database::getNextWallSectionTag() {
	auto tag = 0;

	for(auto& [fst, snd] : wall_section_pool) if(fst > tag) tag = fst;

	return tag + 1;
}

int Database::getNextFrameSectionTag() {
	auto tag = 0;

	for(auto& [fst, snd] : frame_section_pool) if(fst > tag) tag = fst;

	return tag + 1;
}

int Database::getNextElementTag() {
	auto tag = 0;

	for(auto& [fst, snd] : element_pool) if(fst > tag) tag = fst;

	return tag + 1;
}

bool Database::removeNode(const int T) {
	auto I = element_pool.begin();
	while(I != element_pool.end()) {
		if(I->second.encoding.count(T) > 0) I = element_pool.erase(I);
		else ++I;
	}

	return 1 == node_pool.erase(T);
}

bool Database::removeWallSection(const int T) {
	auto I = element_pool.begin();
	while(I != element_pool.end()) {
		if(I->second.type == Element::Type::Wall && I->second.section_tag == T) I = element_pool.erase(I);
		else ++I;
	}

	return 1 == wall_section_pool.erase(T);
}

bool Database::removeFrameSection(const int T) {
	auto I = element_pool.begin();
	while(I != element_pool.end()) {
		if(I->second.type != Element::Type::Wall && I->second.section_tag == T) I = element_pool.erase(I);
		else ++I;
	}

	return 1 == frame_section_pool.erase(T);
}

bool Database::removeElement(const int T) { return 1 == element_pool.erase(T); }

void Database::changePosition(const int tag, QVector3D&& position) { if(node_pool.find(tag) != node_pool.end()) node_pool[tag].position = position; }

void Database::changeFixity(const int tag, QVector<bool>&& fixity) { if(node_pool.find(tag) != node_pool.end()) node_pool[tag].fixity = fixity; }

void Database::changeLoad(const int tag, QVector<double>&& load) { if(node_pool.find(tag) != node_pool.end()) node_pool[tag].load = load; }

void Database::changeMass(const int tag, const double mass) { if(node_pool.find(tag) != node_pool.end()) node_pool[tag].mass = mass; }

void Database::changeDisplacement(const int tag, QVector<double>&& displacement) { if(node_pool.find(tag) != node_pool.end()) node_pool[tag].displacement = displacement; }

void Database::changeSection(const int ele, const int sec) {
	if(element_pool.find(ele) == element_pool.end()) return;

	if(element_pool.at(ele).type == Element::Type::Wall) { if(wall_section_pool.find(sec) == wall_section_pool.end()) return; }
	else if(frame_section_pool.find(sec) == frame_section_pool.end()) return;

	element_pool[ele].section_tag = sec;
}

void Database::changeUnit(const int F) { unit_system = F; }

void Database::changeAnalysisType(const int F) {
	if(F != 0 && F != 1) return;

	analysis_type = F;
}

void Database::changeDamping(const QString& F) { damping_ratio = std::max(0., F.toDouble()); }

void Database::changePeriod(const QString& F) { natural_period = std::max(0., F.toDouble()); }

void Database::changeScale(const QString& F) { scale_factor = std::max(0., F.toDouble()); }

void Database::changeAccxRecord(const QString& F) { acc_record[0] = F; }

void Database::changeAccyRecord(const QString& F) { acc_record[1] = F; }

void Database::splitElement(const int tag, const int segment) {
	if(element_pool.find(tag) == element_pool.end()) return;

	auto& t_element = element_pool[tag];

	const auto& coor_i = node_pool[t_element.encoding.at(0)].position;
	const auto& coor_j = node_pool[t_element.encoding.at(1)].position;

	auto node_tag = getNextNodeTag();

	for(auto I = 1; I < segment; ++I) {
		const auto new_coor = static_cast<float>(I) / static_cast<float>(segment) * (coor_j - coor_i) + coor_i;
		add(node_tag + I - 1, Node{QVector3D{new_coor.x(), new_coor.y(), new_coor.z()}});
	}

	auto element_tag = getNextElementTag();

	t_element.highlighted = false;

	auto element_copy = t_element;

	element_copy.encoding[1] = node_tag;

	element_pool.try_emplace(element_tag++, element_copy);

	element_copy = t_element;

	element_copy.encoding[0] = node_tag + segment - 2;

	element_pool.try_emplace(element_tag++, element_copy);

	for(auto I = 0; I < segment - 2; ++I) {
		element_copy = t_element;

		element_copy.encoding[0] = node_tag++;
		element_copy.encoding[1] = node_tag;

		element_pool.try_emplace(element_tag++, element_copy);
	}

	removeElement(tag);
}

void Database::removeElement() { element_pool.clear(); }

bool Database::loadModel(const QString& file_name) {
	QFile file(file_name);
	file.open(QIODevice::ReadOnly);
	QTextStream script(&file);

	QString command;
	QStringList pool;

	auto skip_blank = [&](const int size = 0) {
		while(true) {
			command = remove_comment(script.readLine());
			if(!command.isEmpty()) {
				pool = command.split(" ");

				auto I = pool.begin();

				while(I != pool.end()) {
					if(I->isEmpty()) I = pool.erase(I);
					else ++I;
				}

				if(!pool.isEmpty()) break;
			}
		}

		return pool.size() >= size;
	};

	// title
	skip_blank();

	// frame quadrature
	if(!skip_blank(3)) return false;
	quadrature_frame = {pool.at(0).toInt(), pool.at(1).toInt(), pool.at(2).toInt()};

	// wall quadrature
	if(!skip_blank(2)) return false;
	quadrature_wall = {pool.at(0).toInt(), pool.at(1).toInt()};

	// unit
	if(!skip_blank(1)) return false;
	changeUnit(pool.at(0).toInt());

	if(!skip_blank(1)) return false;
	changeAnalysisType(pool.at(0).toInt());

	if(!skip_blank(1)) return false;
	const auto accx = pool.at(0).toInt() == 1;
	if(!skip_blank(1)) return false;
	const auto accy = pool.at(0).toInt() == 1;
	if(!skip_blank(2)) return false;
	changeDamping(pool.at(0));
	changePeriod(pool.at(1));
	if(!skip_blank(1)) return false;
	changeScale(pool.at(0));
	if(accx) {
		if(!skip_blank(1)) return false;
		changeAccxRecord(pool.at(0));
	}
	if(accy) {
		if(!skip_blank(1)) return false;
		changeAccyRecord(pool.at(0));
	}

	if(!skip_blank(6)) return false;
	const auto node_num = pool.at(0).toInt();
	const auto beam_num = pool.at(1).toInt();
	const auto brace_num = pool.at(2).toInt();
	const auto wall_num = pool.at(3).toInt();
	const auto frame_type_num = pool.at(4).toInt();
	const auto wall_type_num = pool.at(5).toInt();

	for(auto I = 0; I < node_num; ++I) {
		if(!skip_blank(4)) return false;
		add<Node>(pool.at(0).toInt(), Node{QVector3D{pool.at(1).toFloat(), pool.at(2).toFloat(), pool.at(3).toFloat()}});
	}

	for(auto I = 0; I < frame_type_num; ++I) {
		if(!skip_blank(5)) return false;
		add<FrameSection>(I + 1, FrameSection{pool.at(0), QVector<double>{pool.at(1).toDouble(), pool.at(2).toDouble(), pool.at(3).toDouble(), pool.at(4).toDouble()}});
	}

	for(auto I = 0; I < wall_type_num; ++I) {
		add<WallSection>(I + 1, WallSection{QVector<double>{}});
		auto& t_section = wall_section_pool[I + 1];
		for(auto J = 0; J < 3; ++J) {
			if(!skip_blank(6)) return false;
			for(auto K = 0; K < 6; ++K) t_section.parameter.append(pool.at(K).toDouble());
		}
	}

	for(auto I = 0; I < beam_num; ++I) {
		if(!skip_blank(5)) return false;
		add<Element>(getNextElementTag(), Element{pool.at(4).toInt(), QVector<int>{pool.at(2).toInt(), pool.at(3).toInt()}, "Frame", 0});
	}

	for(auto I = 0; I < brace_num; ++I) {
		if(!skip_blank(5)) return false;
		add<Element>(getNextElementTag(), Element{pool.at(4).toInt(), QVector<int>{pool.at(2).toInt(), pool.at(3).toInt()}, "Brace", 0});
	}

	for(auto I = 0; I < wall_num; ++I) {
		if(!skip_blank(5)) return false;
		add<Element>(getNextElementTag(), Element{pool.at(3).toInt(), QVector<int>{pool.at(1).toInt(), pool.at(2).toInt()}, "Wall", pool.at(4).toInt()});
	}

	if(!skip_blank(1)) return false;

	const auto mass_num = pool.at(0).toInt();

	for(auto I = 0; I < mass_num; ++I) {
		if(!skip_blank(3)) return false;
		changeMass(pool.at(1).toInt(), pool.at(2).toDouble());
	}

	skip_blank();

	const auto bc_num = pool.at(0).toInt();

	for(auto I = 0; I < bc_num; ++I) {
		if(!skip_blank(10)) return false;
		changeFixity(pool.at(1).toInt(), QVector<bool>{pool.at(4).toInt() != 0, pool.at(5).toInt() != 0, pool.at(6).toInt() != 0, pool.at(7).toInt() != 0, pool.at(8).toInt() != 0, pool.at(9).toInt() != 0});
	}

	skip_blank();

	skip_blank();

	if(!skip_blank(3)) return false;
	tolerance[0] = pool.at(0).toDouble();
	tolerance[1] = pool.at(1).toDouble();
	tolerance[2] = pool.at(2).toDouble();

	if(!skip_blank(3)) return false;
	tolerance[3] = pool.at(0).toDouble();
	tolerance[4] = pool.at(1).toDouble();
	tolerance[5] = pool.at(2).toDouble();

	return true;
}

bool Database::saveModel(const QString& file_name) {
	QFile file(file_name);
	file.open(QIODevice::WriteOnly);
	QTextStream output(&file);

    output.setRealNumberNotation(QTextStream::FixedNotation);
	output.setRealNumberPrecision(6);
    output.setRealNumberNotation(QTextStream::ScientificNotation);
    output.setRealNumberPrecision(3);

	output << "MODEL GENERATED BY PB3D\n";
	output << quadrature_frame.at(0) << ' ' << quadrature_frame.at(1) << ' ' << quadrature_frame.at(2) << '\n';
	output << quadrature_wall.at(0) << ' ' << quadrature_wall.at(1) << '\n';
	output << unit_system << "\n\n";
	output << analysis_type << '\n';
	output << (acc_record.at(0).isEmpty() ? 0 : 1) << '\n';
	output << (acc_record.at(1).isEmpty() ? 0 : 1) << '\n';
	output << damping_ratio << ' ' << natural_period << '\n';
	output << scale_factor << '\n';
	if(!acc_record.at(0).isEmpty()) output << acc_record.at(0) << '\n';
	if(!acc_record.at(1).isEmpty()) output << acc_record.at(1) << '\n';

	output << '\n';

	compress();

	output << node_pool.size() << ' ';

	auto frame_size = 0;
	auto brace_size = 0;
	auto wall_size = 0;
	for(auto& [fst, snd] : element_pool) {
		if(snd.type == Element::Type::Frame) frame_size++;
		else if(snd.type == Element::Type::Brace) brace_size++;
		else if(snd.type == Element::Type::Wall) wall_size++;
	}
	output << frame_size << ' ';
	output << brace_size << ' ';
	output << wall_size << ' ';

	output << frame_section_pool.size() << ' ';
	output << wall_section_pool.size() << '\n';

	output << '\n';

	serialize<Node>(output);
	serialize<FrameSection>(output);
	serialize<WallSection>(output);
	serialize<Element>(output);

	serializeMass(output);
	serializeBC(output);

	output << tolerance.at(0) << ' ';
	output << tolerance.at(1) << ' ';
	output << tolerance.at(2) << '\n';
	output << tolerance.at(3) << ' ';
	output << tolerance.at(4) << ' ';
	output << tolerance.at(5) << '\n';

	return true;
}

void Database::serializeMass(QTextStream& output) {
	auto counter = 0;
	for(auto& [fst, snd] : node_pool) if(snd.mass > 0.) ++counter;

	output << counter << " ! TOTAL NUMBER OF NODES APPILED WITH MASS\n";

	auto t_node = node_pool.cbegin();
	for(auto I = 1; I <= counter; ++I) {
		while(t_node->second.mass <= 0.) ++t_node;
		output << I << ' ' << t_node->first << ' ' << t_node->second.mass << '\n';
		++t_node;
	}

	output << "\n";
}

void Database::serializeBC(QTextStream& output) {
	auto counter = 0;
	for(auto& [fst, snd] : node_pool)
		for(auto& I : snd.fixity)
			if(I) {
				++counter;
				break;
			}

	output << counter << " ! TOTAL NUMBER OF NODES APPLIED WITH BC\n";

	QVector<int> bc_list;

	auto t_node = node_pool.cbegin();
	for(auto I = 1; I <= counter; ++I) {
		long long t_num;
		while(true) {
			t_num = t_node->second.fixity.count(true);
			if(t_num != 0) break;
			++t_node;
		}

		bc_list.append(t_node->first);

		output << I << ' ' << t_node->first << ' ' << t_num << " 0";

		auto idx = 1;
		for(auto& J : t_node->second.fixity) output << ' ' << (J ? idx++ : 0);

		output << "\n";

		++t_node;
	}

	output << "\n";

	output << bc_list.size() << " ! NUMBER OF NODES USED FOR CALCULATING BASE SHEAR\n";
	for(auto I = bc_list.cbegin(); I != bc_list.cend(); ++I) output << *I << ((I != bc_list.cend() - 1) ? ' ' : '\n');

	output << '\n';
}

void Database::compress() {
    auto node_tag = getNodeTag();
	for(auto I = 0, J = 1; I < node_tag.size(); ++I, ++J) compress_node(node_tag.at(I), J);

    node_tag = getNodeTag();
    std::vector<bool> label(node_tag.size() + 1, false);

    for(auto& [fst, snd] : element_pool) {
        if(!label.at(snd.encoding.at(0))) label.at(snd.encoding.at(0)) = true;
        if(!label.at(snd.encoding.at(1))) label.at(snd.encoding.at(1)) = true;
    }

    for(auto I = 1; I < label.size(); ++I) if(!label.at(I)) node_pool.erase(I);

	const auto wall_section_tag = getWallSectionTag();
	for(auto I = 0, J = 1; I < wall_section_tag.size(); ++I, ++J) compress_wall_section(wall_section_tag.at(I), J);

	const auto frame_section_tag = getFrameSectionTag();
	for(auto I = 0, J = 1; I < frame_section_tag.size(); ++I, ++J) compress_frame_section(frame_section_tag.at(I), J);
}

void Database::compress_node(const int old_tag, const int new_tag) {
	if(old_tag == new_tag) return;

	auto t_node = node_pool.extract(old_tag);
	t_node.key() = new_tag;
	node_pool.insert(std::move(t_node));

	for(auto& [fst, snd] : element_pool) {
		if(snd.encoding.at(0) == old_tag) snd.encoding[0] = new_tag;
		if(snd.encoding.at(1) == old_tag) snd.encoding[1] = new_tag;
	}
}

void Database::compress_wall_section(const int old_tag, const int new_tag) {
	if(old_tag == new_tag) return;

	auto t_section = wall_section_pool.extract(old_tag);
	t_section.key() = new_tag;
	wall_section_pool.insert(std::move(t_section));

	for(auto& [fst, snd] : element_pool) {
		if(snd.type != Element::Type::Wall) continue;
		if(snd.section_tag == old_tag) snd.section_tag = new_tag;
	}
}

void Database::compress_frame_section(const int old_tag, const int new_tag) {
	if(old_tag == new_tag) return;

	auto t_section = frame_section_pool.extract(old_tag);
	t_section.key() = new_tag;
	frame_section_pool.insert(std::move(t_section));

	for(auto& [fst, snd] : element_pool) {
		if(snd.type == Element::Type::Wall) continue;
		if(snd.section_tag == old_tag) snd.section_tag = new_tag;
	}
}

QString Database::remove_comment(QString in) {
	in.replace(QString("\t"), QString(" "));
	in.replace(QString(","), QString(" "));
	in.remove(QRegularExpression("!.*$"));
	return in;
}

template<typename T> void Database::highlight(int, bool) { throw; }

template<> void Database::highlight<Database::Node>(const int tag, const bool highlighted) {
	const auto t_node = node_pool.find(tag);
	if(t_node == node_pool.end()) return;
	t_node->second.highlighted = highlighted;
}

template<> void Database::highlight<Database::Element>(const int tag, const bool highlighted) {
	const auto t_element = element_pool.find(tag);
	if(t_element == element_pool.end()) return;
	t_element->second.highlighted = highlighted;
}

template<typename T> bool Database::add(int, T&&) { throw; }

template<> bool Database::add<Database::Node>(const int tag, Node&& obj) { return node_pool.try_emplace(tag, std::forward<Node>(obj)).second; }

template<> bool Database::add<Database::WallSection>(const int tag, WallSection&& obj) { return wall_section_pool.try_emplace(tag, std::forward<WallSection>(obj)).second; }

template<> bool Database::add<Database::FrameSection>(const int tag, FrameSection&& obj) { return frame_section_pool.try_emplace(tag, std::forward<FrameSection>(obj)).second; }

template<> bool Database::add<Database::Element>(const int tag, Element&& obj) {
	for(auto& I : obj.encoding) if(node_pool.find(I) == node_pool.end()) return false;

	return element_pool.try_emplace(tag, std::forward<Element>(obj)).second;
}

template<typename T> T& Database::get(int) { throw; }

template<> Database::Node& Database::get<Database::Node>(const int tag) { return node_pool.at(tag); }

template<> Database::WallSection& Database::get<Database::WallSection>(const int tag) { return wall_section_pool.at(tag); }

template<> Database::FrameSection& Database::get<Database::FrameSection>(const int tag) { return frame_section_pool.at(tag); }

template<> Database::Element& Database::get<Database::Element>(const int tag) { return element_pool.at(tag); }

template<typename T> void Database::serialize(QTextStream&) { throw; }

template<> void Database::serialize<Database::Node>(QTextStream& output) {
	for(auto& I : getNodeTag()) {
		output << I << ' ';
		auto& snd = node_pool.at(I);
		output << snd.x() << ' ';
		output << snd.y() << ' ';
		output << snd.z() << '\n';
	}
    output << "\n";
}

template<> void Database::serialize<Database::WallSection>(QTextStream& output) {
	for(auto& I : getWallSectionTag()) {
		auto& snd = wall_section_pool.at(I);
		output << snd.parameter.at(0) << ' ';
		output << snd.parameter.at(1) << ' ';
		output << snd.parameter.at(2) << ' ';
		output << snd.parameter.at(3) << ' ';
		output << snd.parameter.at(4) << ' ';
		output << snd.parameter.at(5) << '\n';
		output << snd.parameter.at(6) << ' ';
		output << snd.parameter.at(7) << ' ';
		output << snd.parameter.at(8) << ' ';
		output << snd.parameter.at(9) << ' ';
		output << snd.parameter.at(10) << ' ';
		output << snd.parameter.at(11) << '\n';
		output << snd.parameter.at(12) << ' ';
		output << snd.parameter.at(13) << ' ';
		output << snd.parameter.at(14) << ' ';
		output << snd.parameter.at(15) << ' ';
		output << snd.parameter.at(16) << ' ';
		output << snd.parameter.at(17) << "\n\n";
    }
}

template<> void Database::serialize<Database::FrameSection>(QTextStream& output) {
	for(auto& I : getFrameSectionTag()) {
		auto& snd = frame_section_pool.at(I);
		output << static_cast<int>(snd.type) + 1 << ' ';
		output << snd.parameter.at(0) << ' ';
		output << snd.parameter.at(1) << ' ';
		output << snd.parameter.at(2) << ' ';
		output << snd.parameter.at(3) << '\n';
	}

    output << "\n";
}

template<> void Database::serialize<Database::Element>(QTextStream& output) {
	for(auto& I : getElementTag()) {
		auto& snd = element_pool.at(I);
		if(snd.type != Element::Type::Frame) continue;
		output << I << " 1 ";
		output << snd.encoding.at(0) << ' ';
		output << snd.encoding.at(1) << ' ';
		output << snd.section_tag << '\n';
	}
    output << "\n";
	for(auto& I : getElementTag()) {
		auto& snd = element_pool.at(I);
		if(snd.type != Element::Type::Brace) continue;
		output << I << " 2 ";
		output << snd.encoding.at(0) << ' ';
		output << snd.encoding.at(1) << ' ';
		output << snd.section_tag << '\n';
	}
    output << "\n";
	for(auto& I : getElementTag()) {
		auto& snd = element_pool.at(I);
		if(snd.type != Element::Type::Wall) continue;
		output << I << ' ';
		output << snd.encoding.at(0) << ' ';
		output << snd.encoding.at(1) << ' ';
		output << snd.section_tag << ' ' << snd.orient << '\n';
	}
    output << "\n";
}

Database::Element::Element(const int st, QVector<int> e, const QString& t, const int o)
	: section_tag(st)
	, encoding(std::forward<QVector<int>>(e))
	, orient(o) {
	if(t == "Wall") type = Type::Wall;
	else if(t == "Brace") type = Type::Brace;
	else if(t == "Frame") type = Type::Frame;
}

Database::FrameSection::FrameSection(const QString& T, QVector<double>&& P)
	: parameter(std::forward<QVector<double>>(P)) {
	if(T == "Steel" || T == "2") type = Type::Steel;
	else if(T == "Wood" || T == "1") type = Type::Wood;
}
