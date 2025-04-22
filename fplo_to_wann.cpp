#define BOHR_TO_ANG 0.529177249
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <eigen3/Eigen/Dense>
#include "fplo_to_wann.h"

using namespace std;
using namespace Eigen;

bool Read_block(ifstream &file, map<array<int, 5>, array<double, 2>> &structure, array<double, 3> &max, array<double, 3> &min, Matrix3d &A){
	string line;
	getline(file, line);
	stringstream ss_1(line);
	string dummy;
	ss_1 >> dummy;
	map<array<int, 3>, int> count;
	if (dummy != "Tij,"){
		return false;
	}
	array<int ,5> key;
	getline(file, line);
	stringstream ss_2(line);
	ss_2 >> dummy;
	key[3] = stoi(dummy);
	ss_2 >> dummy;
	key[4] = stoi(dummy);
	while(getline(file, line)){
		stringstream ss(line);
		bool switch_away = false;
		array<double, 2> values;
		Vector3d b;
		array<int ,3 > position;
		for (int ii = 0; ii < 5; ++ii){
			ss >> dummy;
			if (dummy == "end"){
				switch_away = true;
				break;
			}
			if (ii == 4 || ii == 3){
				values[ii - 3] = stod(dummy);
			}
			else{
				b(ii) = stod(dummy);
			}
		}
			Vector3d solution = A.fullPivLu().solve(b);
			for (int ii = 0; ii < 3; ++ii){
				int signm = 0;
				if (abs(solution[ii]) > 0.01){
					signm = (solution[ii] > 0) ? 1:-1;
				}
				if (signm == 1){
				key[ii] = ceil(abs(solution[ii]) - 0.01) * signm;
				}
				else if (signm == -1){
				key[ii] = floor(abs(solution[ii]) + 0.01) * signm;
				}
				else if (signm == 0){
				key[ii] = 0;
				}
				position[ii] = key[ii];
				if (key[ii] > max[ii]){
					max[ii] = key[ii];
				}
				if (key[ii] < min[ii]){
					min[ii] = key[ii];
				}
			}
		if (switch_away){
			break;
		}
		structure[key] = values;
	}
	return true;
}

void Write_ham(ofstream &file, map<array<int, 5>, array<double, 2>> &structure, int num_wann, array<double, 3> &max, array<double, 3> &min){
	file << "Some text\n";
	int num_lines = (abs(max[0]) + abs(min[0]) + 1) * (abs(max[1]) + abs(min[1]) + 1) * (abs(max[2]) + abs(min[2]) + 1);
	file << num_wann << "\n";
	file << num_lines << "\n";
	for (int ii = 0; ii < num_lines; ++ii){
		if (ii % 15 == 0 && ii != 0){
			file << "\n";
		}
		file << setw(5) << 1;
	}
	file << "\n";
	array<int, 5> indices;
	array<double, 2> to_print;
	for (int ii = min[0]; ii <= max[0]; ++ii){
		for (int jj = min[1]; jj <= max[1]; ++jj){
			for (int kk = min[2]; kk <= max[2]; ++kk){
				for (int nn = 1; nn <= num_wann; ++nn){
					for (int mm = 1; mm <= num_wann; ++mm){
						indices = {ii, jj, kk, nn, mm};
						for (int hh = 0; hh < 5; ++hh){
							file << setw(5) << indices[hh];
						}
						if(structure.contains(indices)){
							to_print = structure[indices];
						}
						else{
						to_print = {0, 0};
						}
						file << setw(20) << fixed << setprecision(6) << to_print[0];
						file << setw(20) << fixed << setprecision(6) << to_print[1];
						file << endl;
						}
					}
				}
			}
		}
	}
void To_file_ham(map<array<int, 5>, array<double, 2>> &structure_1, map<array<int, 5>, array<double, 2>> &structure_2, ofstream &file, vector<string> &orbitals, array<int, 3> &cell, int num_wann){
	array<int, 5> indices;
	indices[0] = cell[0];
	indices[1] = cell[1];
	indices[2] = cell[2];
	file << "x" << "\t";
	for (auto &str : orbitals){
		file << str << "\t";
	}
	file << "x" << "\t";
	for (auto &str : orbitals){
		file << str << "\t";
	}
	file << "\n";
	for (int ii = 1; ii <= num_wann; ++ii){
		file << orbitals[ii - 1] << "\t";
		indices[3] = ii;
		for (int jj = 1; jj <= num_wann; ++jj){
			indices[4] = jj; 
			file << structure_1[indices][0] << "\t";
		}

		file << orbitals[ii -1] << "\t";
		indices[3] = ii;
		for (int jj = 1; jj <= num_wann; ++jj){
			indices[4] = jj; 
			file << structure_2[indices][0] << "\t";
		}
	file << "\n";
	}
}
void Read_orbs(ifstream &file, vector<string> &orbitals, int num_wann){
	for (int ii = 0; ii < num_wann; ++ii){
		string line;
		getline(file, line);
		stringstream ss(line);
		string dummy, orb;
		ss >> dummy;
		ss >> orb;
		orbitals.push_back(orb);
	}
}

Matrix3d Get_to_centres(ifstream &data, ofstream &POSCAR, vector<string> &orbs, int &num_wann, int &num_spin){
	Matrix3d basis_2;
	string line;
	POSCAR << "SOme text" << "\n" << 1 << "\n";
	while(getline(data, line)){
		stringstream ss(line);
		string first;
		ss >> first;
		if (first == "wannames:"){
			Read_orbs(data, orbs, num_wann);
		}
		if (first == "nspin:"){
			getline(data, line);
			stringstream ss_1(line);
			ss_1 >> first;
			num_spin = stoi(first);
			
		}
		if (first == "nwan:"){
			getline(data, line);
			stringstream ss_1(line);
			ss_1 >> first;
			num_wann = stoi(first);
		}
		if (first == "lattice_vectors:"){
			for (int ii = 0; ii < 3; ++ii){
				getline(data, line);
				stringstream ss_1(line);
				for (int jj = 0; jj < 3; ++jj){
					ss_1 >> first;
					basis_2(jj, ii) = stod(first);
					POSCAR << setw(17) << fixed <<setprecision(5) << stod(first) * BOHR_TO_ANG;
				}
				POSCAR << "\n";
			}
		}
		if (first == "wancenters:"){
			break;
		}
	}
	return basis_2;
}

map<string, int> Get_elements(string filename){
	int garbage;
	string file_name = "SAVE_COUNT_AND_NAMES_OF_FPLO_ELEMENTS.txt";
	string command1 = "touch "+filename;
	string command2 = "grep -E '^([A-Za-z])+[1-9]+ ' \"" + filename + "\"| sed -E 's/^([A-Za-z]+[1-9]+).+/\\1/'| sort -u |sed -E 's/^([A-Za-z]+)[1-9]+/\\1/'| uniq -c > "+file_name;
	string command3 = "rm "+file_name;
	garbage = system(command1.c_str());
	garbage = system(command2.c_str());
	ifstream reader(file_name);
	string line;
	map<string, int> to_return;
	while(getline(reader, line)){
		stringstream ss(line);
		string count_str, element;
		ss >> count_str;
		ss >> element;
		int count = stoi(count_str);
		to_return[element] = count;
	}
	garbage = system(command3.c_str());
	return to_return;
}
vector<string> Get_elements_order(string filename){
	int garbage;
	string file_name = "SAVE_ORDER_OF_FPLO_ELEMENTS.txt";
	string command1 = "touch "+file_name;
	string command2 = "grep -E '^([A-Za-z])+[1-9]+ ' \""+ filename +"\"| sed -E 's/^([A-Za-z]+[1-9]+).+/\\1/'| sort -u| sed -E 's/^([A-Za-z]+)[1-9]+/\\1/' > "+file_name;
	string command3 = "rm "+file_name;
	garbage = system(command1.c_str());
	garbage = system(command2.c_str());
	ifstream reader(file_name);
	string line;
	vector<string> to_return;
	while(getline(reader, line)){
		stringstream ss(line);
		string element;
		ss >> element;
		to_return.push_back(element);
	}
	garbage = system(command3.c_str());
	return to_return;
}

vector<array<double,3>> Write_centres(ifstream &data, ofstream &Output, vector<string> &ordered_elements, int num_wann){
	Output << setw(6) << num_wann + ordered_elements.size()  << endl;
	Output << "Some text" << endl;
	string line;
	vector<array<double, 3>> tosafe;
	while(getline(data, line)){
		array<double, 3> center;
		stringstream ss(line);
		string number;
		bool switch_away = false;
		for (int ii = 0; ii < 3; ++ii){
			ss >> number;
			if (number == "symmetry:"){
				switch_away = true;
				break;
			}
			if (ii == 0){
				Output <<setw(4) << left << "X" << right;
			}
			center[ii] = stod(number) * BOHR_TO_ANG;
			Output << setw(17) << setprecision(8) << fixed <<center[ii];
		}
		if (switch_away){
			break;
		}
		Output << endl;
		if (tosafe.empty()){
			tosafe.push_back(center);
		}
		else{
			bool switch_last = false;
			for (int ii = 0; ii < 3; ++ii){
				if(center[ii] != tosafe.back()[ii]){
					switch_last = true;
				}
			}
			if (switch_last){

				tosafe.push_back(center);
			}
		}



	}
	for (unsigned int ii = 0; ii < tosafe.size(); ++ii){
		for (int jj = 0; jj < 1; ++jj){
			Output << setw(4) << left << ordered_elements[ii] << right;
			for (int kk = 0; kk < 3; ++kk){
				Output << setw(17) << setprecision(8) << fixed << tosafe[ii + jj][kk];
			}
			Output << endl;
		}
	}
	return tosafe;
}

void Write_to_POSCAR(ofstream &POSCAR, map<string, int> &elements, vector<string> ordered_elements, vector<array<double, 3>> centres){
	for (const auto& el : elements){
		POSCAR << setw(6) << el.first;
	}
	POSCAR << "\n";
	for (const auto& el : elements){
		POSCAR << setw(6) << el.second;
	}
	POSCAR << "\n";
	POSCAR << "Cartesian\n";
	for (const auto& el : elements){
		for (unsigned int ii = 0; ii < ordered_elements.size(); ++ii){
			if (el.first == ordered_elements[ii]){
				for (int jj = 0; jj < 3; ++jj){
					POSCAR << setw(17) << fixed <<setprecision(5) << centres[ii][jj];
				}
				POSCAR << "\n";
			}
		}
	}
}

void Find_hamiltonians(ifstream &data){
	string line;
	while(getline(data, line)){
		string first;
		stringstream ss(line);
		ss >> first;
		if (first == "spin:"){
			break;
		}
	}
	getline(data, line);
}
