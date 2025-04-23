#if __has_include(<Eigen/Dense>)
    #include <Eigen/Dense>
#elif __has_include(<eigen3/Eigen/Dense>)
    #include <eigen3/Eigen/Dense>
#else
    #error "Could not find Eigen library. Please install Eigen or set the include path correctly."
#endif

#include <map>
#include <array>
#include <fstream>
#include <string>
#include <vector>
bool Read_block(std::ifstream &file, std::map<std::array<int, 5>, std::array<double, 2>> &structure, std::array<double, 3> &max, std::array<double, 3> &min, Eigen::Matrix3d &A);
void Write_ham(std::ofstream &file, std::map<std::array<int, 5>, std::array<double, 2>> &structure, int num_wann, std::array<double, 3> &max, std::array<double, 3> &min);
void To_file_ham(std::map<std::array<int, 5>, std::array<double, 2>> &structure_1, std::map<std::array<int, 5>, std::array<double, 2>> &structure_2, std::ofstream &file, std::vector<std::string> &orbitals, std::array<int, 3> &cell, int num_wann);
void Read_orbs(std::ifstream &file, std::vector<std::string> &orbitals, int num_wann);
Eigen::Matrix3d Get_to_centres(std::ifstream &data, std::ofstream &POSCAR, std::vector<std::string> &orbs, int &num_wann, int &num_spin);
std::vector<std::string> Get_elements_order(std::string filename);
std::map<std::string, int> Get_elements(std::string filename);
std::vector<std::array<double,3>> Write_centres(std::ifstream &data, std::ofstream &Output, std::vector<std::string> &ordered_elements, int num_wann);
void Write_to_POSCAR(std::ofstream &POSCAR, std::map<std::string, int> &elements, std::vector<std::string> ordered_elements, std::vector<std::array<double, 3>> centres);
void Find_hamiltonians(std::ifstream &data);
