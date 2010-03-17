#include "Riostream.h"
#include <math.h>

/** 
 * read cov-errors.txt files
 * 
 * formatted like x,y,z,.... mean(position), var(position), truevalue(position), error(position)
 */

int emu_ascii_read(char* in_dir=0x0, const int nparams, char* filename=0x0, char* out_dir=0x0){
	const int CHAR_BUF_SIZE = 128;
	const int READ_BUF_SIZE = 512;
	char charbuf[CHAR_BUF_SIZE];	
	char buf[READ_BUF_SIZE];
	char *split_string;
	int linecount = 0;
	int maxlines = 500; // hack
	int i,j;
	ifstream file;

	sprintf(charbuf, "%s/%s", in_dir, filename);
	file.open(charbuf);

	sprintf(charbuf, "%s/%s.root", out_dir, filename);
	//TFile* rootfile = new TFile("rootfile.root", "RECREATE") //if you want to overwrite, "NEW" if you want to create but not overwrite
	TFile *rfile = TFile::Open(charbuf, "RECREATE"); // if it's not recreate you cant write to the tree for some reason?
	//TFile rfile(charbuf, "CREATE", "testing lpm data analysis");
	TTree* tree = new TTree("tree", "emu-data-tree");
	tree->SetDirectory(rfile);

	Float_t mean, var, truevalue, error;
	Float_t params[nparams];

	for(i = 0; i < nparams; i++) params[i] =0.0;
	mean = 0.0; 
	var = 0.0;
	truevalue = 0.0;
	error = 0.0;
	
	// assign the param branches
	for(i = 0;i < nparams; i++){
		sprintf(charbuf, "x%d", i);
		sprintf(buf, "%s/F", charbuf); // need two buffers here
		tree->Branch(charbuf, &(params[i]), buf);
	}
	tree->Branch("mean", &mean, "mean/F");
	tree->Branch("var", &var, "var/F");
	tree->Branch("truevalue", &truevalue, "truevalue/F");
	tree->Branch("error", &error, "error/F");
	tree->Print();

	for(i = 0; i < maxlines; i++){
		file.getline(buf,READ_BUF_SIZE); 
		printf("%s\n", buf);
		split_string = strtok(buf, "\t ");
		for(j = 0; j < nparams; j++){ // read all the params in 
			sscanf(split_string, "%f", &params[j]);
			split_string = strtok(NULL, "\t ");
			printf("%f ", params[j]);
		}
		sscanf(split_string, "%f", &mean);
		split_string = strtok(NULL, "\t ");
		sscanf(split_string, "%f", &var);
		split_string = strtok(NULL, "\t ");
		sscanf(split_string, "%f", &truevalue);
		split_string = strtok(NULL, "\t ");
		sscanf(split_string, "%f", &error);
		printf("%f %f %f %f \n", mean, var, truevalue, error);
		tree->Fill();
	}
	
	tree->Print();
	tree->Write();
	file.close();
	delete rfile;
	return(0);
}
