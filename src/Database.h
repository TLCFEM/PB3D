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

#ifndef DATABASE_H
#define DATABASE_H

#include <QTextStream>
#include <QVector3D>
#include <QVector>

class Database {
public:
	struct Node {
		QVector3D position = QVector3D(0, 0, 0);
		QVector<bool> fixity = QVector<bool>(6, false);
		QVector<double> load = QVector<double>(6, 0.);
		QVector<double> displacement = QVector<double>(6, 0.);
		double mass = 0.;
		bool highlighted = false;

		[[nodiscard]] float x() const { return position.x(); }

		[[nodiscard]] float y() const { return position.y(); }

		[[nodiscard]] float z() const { return position.z(); }
	};

	struct WallSection {
		QVector<double> parameter = QVector<double>();
	};

	struct FrameSection {
		enum class Type : int {
			Wood,
			Steel
		};

		FrameSection(const QString&, QVector<double>&&);

		Type type = Type::Steel;
		QVector<double> parameter = QVector<double>();
	};

	struct Element {
		enum class Type : int {
			Wall,
			Brace,
			Frame
		};

		explicit Element(int = 0, QVector<int>  = QVector<int>(), const QString& = "Frame", int = 1);

		int section_tag = 0;
		QVector<int> encoding = QVector<int>();
		Type type = Type::Frame;
		int orient = 1;
		bool highlighted = false;
	};

	[[nodiscard]] const std::unordered_map<int, Node>& getNodePool() const;
	[[nodiscard]] const std::unordered_map<int, WallSection>& getWallSectionPool() const;
	[[nodiscard]] const std::unordered_map<int, FrameSection>& getFrameSectionPool() const;
	[[nodiscard]] const std::unordered_map<int, Element>& getElementPool() const;

	[[nodiscard]] QVector<int> getNodeTag() const;
	[[nodiscard]] QVector<int> getWallSectionTag() const;
	[[nodiscard]] QVector<int> getFrameSectionTag() const;
	[[nodiscard]] QVector<int> getElementTag() const;

	int getNextNodeTag();
	int getNextWallSectionTag();
	int getNextFrameSectionTag();
	int getNextElementTag();

	template<typename T> bool add(int, T&&);
	template<typename T> void highlight(int, bool);
	template<typename T> T& get(int);

	bool removeNode(int);
	bool removeWallSection(int);
	bool removeFrameSection(int);
	bool removeElement(int);

	void changePosition(int, QVector3D&&);
	void changeFixity(int, QVector<bool>&&);
	void changeLoad(int, QVector<double>&&);
	void changeMass(int, double);
	void changeDisplacement(int, QVector<double>&&);
	void changeSection(int, int);
	void splitElement(int, int);
	void removeElement();

	void changeUnit(int);
	void changeAnalysisType(int);
	void changeDamping(const QString&);
	void changePeriod(const QString&);
	void changeScale(const QString&);
	void changeAccxRecord(const QString&);
	void changeAccyRecord(const QString&);

	bool loadModel(const QString&);
	bool saveModel(const QString&);

	template<typename T> void serialize(QTextStream&);

	void serializeMass(QTextStream&);
	void serializeBC(QTextStream&);

    QVector<int> quadrature_frame = QVector<int>(3, 5);
    QVector<int> quadrature_wall = QVector<int>(2, 5);

	int unit_system = 1;
	int analysis_type = 1;
	double damping_ratio = 5.;
	double natural_period = 1.;
	double scale_factor = 1.;

	QVector<QString> acc_record = QVector<QString>(2, "");

	QVector<double> tolerance = QVector<double>(6, 1E-3);

protected:
	std::unordered_map<int, Node> node_pool;
	std::unordered_map<int, WallSection> wall_section_pool;
	std::unordered_map<int, FrameSection> frame_section_pool;
	std::unordered_map<int, Element> element_pool;

	void compress();
	void compress_node(int, int);
	void compress_wall_section(int, int);
	void compress_frame_section(int, int);

	static QString remove_comment(QString);
};

template<> bool Database::add<Database::Node>(int, Node&&);
template<> bool Database::add<Database::WallSection>(int, WallSection&&);
template<> bool Database::add<Database::FrameSection>(int, FrameSection&&);
template<> bool Database::add<Database::Element>(int, Element&&);
template<> void Database::serialize<Database::Node>(QTextStream&);
template<> void Database::serialize<Database::WallSection>(QTextStream&);
template<> void Database::serialize<Database::FrameSection>(QTextStream&);
template<> void Database::serialize<Database::Element>(QTextStream&);

#endif // DATABASE_H
