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


int stacked_hist(char* in_dir=0x0, char* filename=0x0, char* out_dir=0x0){
	// open a root file and make a contour plot
	const int CHAR_BUF_SIZE = 128;
	char charbuf[CHAR_BUF_SIZE];
	
	TCanvas *c2h = new TCanvas("c2h","testing",10,10,800,600);
	c2h->Divide(1,1); // dont need to divide here
	c2h->SetFillColor(17);
	c2h->cd(1);
	sprintf(charbuf, "%s/%s", in_dir, filename);
	TFile *rfile = TFile::Open(charbuf); // if it's not recreate you cant write to the tree for some reason?
	int i;
	Float_t x,y,z,w;
	TBranch *branch1 = tree->GetBranch("x0");
	TBranch *branch2 = tree->GetBranch("x1");
	TBranch *branch3 = tree->GetBranch("mean");
	TBranch *branch4 = tree->GetBranch("var");
	branch1->SetAddress(&x);
	branch2->SetAddress(&y);
	branch3->SetAddress(&z);
	branch4->SetAddress(&w);
	
	// histogram stack
	THStack *stack = new THStack("stack", "stacked 2d hists");
	
	// the mean value
	TH2F *h2sta = new TH2F("h2sta", "h2sta", 20, 0,1,20,0,1);
	h2sta->SetFillColor(38);
	for(i = 0; i < tree->GetEntries(); i++){
		branch1->GetEntry(i);
		branch2->GetEntry(i);		
		branch3->GetEntry(i);
		h2sta->Fill(x,y,z);
	}

	// var plus value
	TH2F *h2stb = new TH2F("h2stb", "h2stb", 20, 0,1,20,0,1);
	h2stb->SetFillColor(39);
	for(i = 0; i < tree->GetEntries(); i++){
		branch1->GetEntry(i);
		branch2->GetEntry(i);		
		branch3->GetEntry(i);
		branch4->GetEntry(i);
		h2stb->Fill(x,y,z+(w/2));
	}

	// var minus
	TH2F *h2stc = new TH2F("h2stc", "h2stc", 20, 0,1,20,0,1);
	h2stc->SetFillColor(40);
	for(i = 0; i < tree->GetEntries(); i++){
		branch1->GetEntry(i);
		branch2->GetEntry(i);		
		branch3->GetEntry(i);
		branch4->GetEntry(i);
		h2stc->Fill(x,y,z-(w/2));
	}

	stack->Add(h2sta);
	stack->Add(h2stb);
	stack->Add(h2stc);
	stack->Draw("colz");
	return c2h;


}

int contour_maker(char* in_dir=0x0, char* filename=0x0, char* out_dir=0x0){
	// open a root file and make a contour plot
	const int CHAR_BUF_SIZE = 128;
	char charbuf[CHAR_BUF_SIZE];
	
	TCanvas *c2h = new TCanvas("c2h","testing",10,10,800,600);
	c2h->Divide(1,2);
	c2h->SetFillColor(17);
	c2h->cd(1);

	sprintf(charbuf, "%s/%s", in_dir, filename);
	TFile *rfile = TFile::Open(charbuf); // if it's not recreate you cant write to the tree for some reason?
	int i;
	Float_t x,y,z,w;
	TBranch *branch1 = tree->GetBranch("x0");
	TBranch *branch2 = tree->GetBranch("x1");
	TBranch *branch3 = tree->GetBranch("mean");
	TBranch *branch4 = tree->GetBranch("truevalue");
	branch1->SetAddress(&x);
	branch2->SetAddress(&y);
	branch3->SetAddress(&z);
	branch4->SetAddress(&w);

	TGraph2D *graph = new TGraph2D();

	for(i = 0; i < tree->GetEntries(); i++){
		branch1->GetEntry(i);
		branch2->GetEntry(i);		
		branch3->GetEntry(i);
		branch4->GetEntry(i);
		//cout << x << y << z << endl;
		graph->SetPoint(i,x,y,z);
		//graph->SetPointError(i,0,0,w);
	}
	gStyle->SetPalette(1);
	graph->Draw("SURF1");

	c2h->cd(2);
	TGraph2D *graph2 = new TGraph2D();

	for(i = 0; i < tree->GetEntries(); i++){
		branch1->GetEntry(i);
		branch2->GetEntry(i);		
		branch3->GetEntry(i);
		branch4->GetEntry(i);
		//cout << x << y << z << endl;
		graph2->SetPoint(i,x,y,w);
		//graph->SetPointError(i,0,0,w);
	}
	gStyle->SetPalette(1);
	graph2->Draw("SURF1");
	c2h->Update();
	

}
	

