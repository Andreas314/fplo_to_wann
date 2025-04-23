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
#include <unistd.h>
#include <eigen3/Eigen/Dense>
#include "fplo_to_wann.h"
using namespace std;
using namespace Eigen;

int main(int argc, char* argv[])
{
	cout << "fplo_to_wann: Translation of fplo +hamdata to Wannier90 _hr.dat file.\n";
	//Open +hamdata and check if it exists
	ifstream data("+hamdata");
	if (!data.is_open()){
		throw runtime_error("Error: No +hamdata in current directory!");
	}

	//A block processing input command
	//So far, only accepts -H and three integers, which tells the programme number of cell, whose hamiltonian
	//is to be printed to a file H_X_Y_Z.dat
	bool write_ham_to_file = false;
	string inp_dumm;
	array<int, 3> print_cell; //cell number
	string dummy_1 = "";
	if (argc > 1){
		dummy_1 = argv[1];
	}
	if (dummy_1 == "-H"){
		write_ham_to_file = true;
		//Read the cell number
		for (int ii = 0; ii < 3; ++ii){
			dummy_1 = argv[ii + 2];
			try{
				print_cell[ii] = stoi(dummy_1);
			}
			catch(...){
				invalid_argument("Invalid cell index: " + dummy_1);
			}
		}
	}
	else if (argc > 1){
		throw invalid_argument("Invalid flag: " + dummy_1);
	}

	//Two functions reading number and order of elements using bash functions and regex
	//Important to correctly write POSCAR and centres files
	map<string, int> elements = Get_elements("+hamdata"); //Stores and element names and counts
	vector<string> ordered_elements = Get_elements_order("+hamdata"); //Store element order in +hamdata
	vector<string> orbs; //orbital names from +hamdata
	int num_wann, num_spin; //number of wannier functions and spin
	int garbage; //garbage value to store output of system command
	
	garbage = system("test -d fplo_to_wann_results_files && rm -r fplo_to_wann_results_files");
	garbage = system("mkdir fplo_to_wann_results_files");
	garbage = chdir("fplo_to_wann_results_files");
	garbage = system("touch POSCAR");
	ofstream POSCAR("POSCAR");
	//Read beginning of hamdata -> number of spins and wannier functions, name of orbitals 
	//and lattice vectors in real space which are also written to POSCAR and outputted
	//Read until you get data to wanncentres
	Matrix3d basis = Get_to_centres(data, POSCAR, orbs, num_wann, num_spin);
	
	cout << "Started reading +hamdata with:\n npin: " << num_spin << "\n nwann: " << num_wann <<"\n";
	
	//Block creating output file for wannier centres
	string centres_name;
	if (num_spin == 2){
		centres_name = "fplo_to_wann.up_centres.xyz";
	}
	else{
		centres_name = "fplo_to_wann_centres.xyz";
	}
	string command = "touch " + centres_name;
	garbage = system(command.c_str());
	ofstream Output(centres_name.c_str());
	
	//Write wannier centres to OUTPUT file, then take the uniqe ones and output them as 
	//position of ionts in the cell
	vector<array<double,3>> centres = Write_centres(data, Output, ordered_elements, num_wann);
	
	//Write number and name of atoms of each element and position of atoms to POSCAR
	Write_to_POSCAR(POSCAR, elements, ordered_elements, centres);
	
	if (num_spin == 2){
		garbage = system("cp fplo_to_wann.up_centres.xyz fplo_to_wann.down_centres.xyz");
	}
	Output.close();

	cout << "Wannier centres written, start processing hopping terms.\n";
 
	//Find block with Hamiltonian energies -> skip unimportant lines with symmetries
	Find_hamiltonians(data);
	//Variables containing hopping terms keys are in order [3*cell indices, 2*orbital indices],
	//values are real and imagianry part of the hopping term
	map<array<int, 5>, array<double, 2>> Hamiltonian_spin1, Hamiltonian_spin2;
	//Variables to store minimal and maximal values of cell indices to set boundaries for 
	//the loop writting the result to a file
	array<double, 3> min_1 = {0, 0, 0}, min_2 = {0, 0, 0}, max_1 = {0, 0, 0}, max_2 = {0, 0, 0};
	bool switch_away = true;
	
	//Read both spins
	while(switch_away){
		//Read block for one pair of orbitals and keep track of min and max
		switch_away = Read_block(data, Hamiltonian_spin1, max_1, min_1, basis);
	}
	string hr_name;
	if (num_spin == 2){
		hr_name = "fplo_to_wann.up_hr.dat";
	}
	else{
		hr_name = "fplo_to_wann_hr.dat";
	}
	command = "touch "+hr_name;
	garbage = system(command.c_str());
	ofstream Output_2(hr_name);
	Write_ham(Output_2, Hamiltonian_spin1, num_wann, max_1, min_1);
	Output_2.close();
	if (num_spin == 2){
		string line;
		getline(data, line);
		getline(data, line);
		switch_away = true;
		while(switch_away){
			switch_away = Read_block(data, Hamiltonian_spin2, max_2, min_2, basis);
		}

		garbage = system("touch fplo_to_wann.down_hr.dat");
		Output_2.open("fplo_to_wann.down_hr.dat");
		Write_ham(Output_2, Hamiltonian_spin2, num_wann, max_2, min_2);
		Output_2.close();
	}
	if (write_ham_to_file){
		string file_name = "HAMILTONIAN_" + to_string(print_cell[0]) + "_"+ to_string(print_cell[1]) + "_"+ to_string(print_cell[2]) + ".dat";
		string command = "touch " + file_name;
		garbage = system(command.c_str());
		Output_2.open(file_name);
		To_file_ham(Hamiltonian_spin1, Hamiltonian_spin2, Output_2, orbs, print_cell, num_wann);	
	}
	cout << "All processing finished.\n";
}
