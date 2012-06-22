#include <string>
#include <vector>
#include <iostream>
#include <stdio.h>
#include "EmuPlusPlus.h" // load the emuplusplus header

using namespace std;

int main (int argc, char** argv){

	if(argc < 2){
		fprintf(stderr, "run with path to statefile as first argument\n");
		return(EXIT_FAILURE);
	}
		
	string filename(argv[1]);

	cout << "# loading emulator from: " << filename << endl;

	emulator my_emu(filename); // construct the emulator
	
	vector<double> the_point;
	vector<double> the_mean;
	vector<double> the_err;
	double dtemp;
	int expected_params = my_emu.number_params;
	int r;
	//int expected_r = 1;
	
	// now read locations in sample space from stdio and loop until eof
	while(!feof(stdin)){
		//for(int i = 0; (i <  expected_params) && ( r == expected_r);  i++){
		for(int i = 0; (i <  expected_params);  i++){
			r = fscanf(stdin, "%lf%*c", &dtemp);
			the_point.push_back(dtemp);
		}

		// if( r < expected_r) {
		// 	break;
		// }

		cout << "# the_point: ";
		for(int i = 0; i < expected_params; i++)
			cout << the_point[i] << " "; 
		cout << endl;

		my_emu.QueryEmulator(the_point, the_mean, the_err);

		cout << "# mean: ";
		for(int i = 0; i < the_mean.size(); i++)
			cout << the_mean[i] << " "; 
		cout << endl;

		cout << "# err: ";
		for(int i = 0; i < the_err.size(); i++)
			cout << the_err[i] << " "; 
		cout << endl;
		
		// now clear up for the next point
		the_point.clear();
		the_mean.clear();
		the_err.clear();
	}			
		
	
	return EXIT_SUCCESS;
}
