#include <map>
#include <array>
#include <fstream>
#include<Eigen/Dense>
#include <string>
#include <vector>
bool Read_block(std::ifstream &file, std::map<std::array<int, 5>, std::array<double, 2>> &structure, double lattice, std::array<double, 3> &max, std::array<double, 3> &min, Eigen::Matrix3d &A);
void Write_ham(std::ofstream &file, std::map<std::array<int, 5>, std::array<double, 2>> &structure, int num_wann, std::array<double, 3> &max, std::array<double, 3> &min);
void To_file_ham(std::map<std::array<int, 5>, std::array<double, 2>> &structure_1, std::map<std::array<int, 5>, std::array<double, 2>> &structure_2, std::ofstream &file, std::vector<std::string> &orbitals, std::array<int, 3> &cell, int num_wann);
void Read_orbs(std::ifstream &file, std::vector<std::string> &orbitals, int num_wann);
