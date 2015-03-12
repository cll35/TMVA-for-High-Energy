#include <TSystem.h>
#include <TChain.h>
#include <TROOT.h>
#include <TRint.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TEventList.h>
#include "TClonesArray.h"
#include <TH1D.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TF1.h>
#include "TLorentzVector.h"
#include "TCanvas.h"
#include <TMath.h>
#include <TProfile.h>
#include "TStyle.h"
#include <time.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include "TLatex.h"
#include "TLegend.h"
#include "Utilities.h"
#include "TVectorD.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "DelphesClasses.h"

using namespace std;

double MH2=100.0*100.0;

double lepton_invariant_mass(Muon *ptlep, MissingET *miss)
{
	double plx=(ptlep->PT)*cos(ptlep->Phi);
	double ply=(ptlep->PT)*sin(ptlep->Phi);
	double misspx=(miss->MET)*cos(miss->Phi);
	double misspy=(miss->MET)*sin(miss->Phi);
	
	double result=sqrt(2*(ptlep->PT)*(miss->MET)-2*(plx*misspx+ply*misspy ));
	return result;
}

double lepton_invariant_mass(Electron *ptlep, MissingET *miss)
{
	double plx=(ptlep->PT)*cos(ptlep->Phi);
	double ply=(ptlep->PT)*sin(ptlep->Phi);
	double misspx=(miss->MET)*cos(miss->Phi);
	double misspy=(miss->MET)*sin(miss->Phi);
	
	double result=sqrt(2*(ptlep->PT)*(miss->MET)-2*(plx*misspx+ply*misspy ));
	return result;
}

//Calculates the top invariant mass using ELECTRON
//and missing energy and higgs mass
double top_invariant_mass(Electron *ptlep, MissingET *miss)
{
	double plx=(ptlep->PT)*cos(ptlep->Phi);
	double ply=(ptlep->PT)*sin(ptlep->Phi);
	double plz=(ptlep->PT)*sinh(ptlep->Eta);
	double pl  =(ptlep->PT)*cosh(ptlep->Eta); //massless lepton energy
	
	double mh2=MH2;
	double mt=173.7;
	double misspx=(miss->MET)*cos(miss->Phi);
	double misspy=(miss->MET)*sin(miss->Phi);
	double a=pow(((double)plz/pl),2)-1;
	double b=2*((plx*misspx+ply*misspy)/(pl)+mh2/(2*pl))*(plz/pl);
	double c=((plx*misspx+ply*misspy)/(pl)+mh2/(2*pl))-pow((miss->MET),2);
	double res_sqrt=(b*b-4*a*c);
	
	//cout << a << "\t" << b << "\t" << c << "\t"  << res_sqrt << endl;
	if (res_sqrt<0) return 10.0;
	else
	{
		double result1=(-b- sqrt(res_sqrt))/(2*a);
		double result2=(-b+sqrt(res_sqrt))/(2*a);
		if (result1>0&&result2<0) return result1;
		else if (result1<0&&result2>0) return result2;
		else if (result1>0&&result2>0) return min(result1-mt,result2-mt) +mt;
		return 0;
	}
}

double top_invariant_mass(Muon *ptlep, MissingET *miss)
{
	double plx=(ptlep->PT)*cos(ptlep->Phi);
	double ply=(ptlep->PT)*sin(ptlep->Phi);
	double plz=(ptlep->PT)*sinh(ptlep->Eta);
	double pl=(ptlep->PT)*cosh(ptlep->Eta);
	
	double mh2=MH2;
	double mt=173.7;
	double misspx=(miss->MET)*cos(miss->Phi);
	double misspy=(miss->MET)*sin(miss->Phi);
	double a=pow(((double)plz/pl),2)-1;
	double b=2*((plx*misspx+ply*misspy)/(pl)+mh2/(2*pl))*(plz/pl);
	double c=((plx*misspx+ply*misspy)/(pl)+mh2/(2*pl))-pow((miss->MET),2);
	double res_sqrt=(b*b-4*a*c);
	
	//cout << a << "\t" << b << "\t" << c << "\t"  << res_sqrt << endl;
	if (res_sqrt<0) return 10.0;
	else
	{
		double result1=(-b- sqrt(res_sqrt))/(2*a);
		double result2=(-b+sqrt(res_sqrt))/(2*a);
		if (result1>0&&result2<0) return result1;
		else if (result1<0&&result2>0) return result2;
		else if (result1>0&&result2>0) return min(result1-mt,result2-mt) +mt;
		return 0;
	}
}


//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

int main(int argc, char*argv[])
{
	//----------------------------------------------------------------------------------------
	/////////////////////////////////////////////  Read Input Parameters  ////////////////////
	//----------------------------------------------------------------------------------------
	ofstream myfile;
	myfile.open ("Signal.csv", ios::out | ios::app);
	
	myfile << "Signal File Index,"<< "M_lv,"<<"Missing_Energy,"<<"P_btagged,"<<"P_Tl,"<<"M_lvb,"<<"N_jets,"<<endl;
	
	if (argc < 2) {printf("******Error in input parameters\n");return 1;}
	
	CommandLine c1;
	c1.parse(argc,argv);
	gROOT->ProcessLine("#include >vector<");
	
	string InputFileName	= c1.getValue<string>	("InputFileName");
	string OutputFileName	= c1.getValue<string>	("OutputFileName");
	//string OutputFileTag	= c1.getValue<string>	("OutputFileTag");
	//string JetAlgo		= c1.getValue<string>	("JetAlgo");
	//vector<string> HLTbit	= c1.getVector<string> 	("HLTbit","");
	int NENTRIES			= c1.getValue<int>	("NEntries");
	
	if (!c1.check()) return 0;
	c1.print(); // Printing the options
	
	//string outputfile = OutputFileName + "_" + OutputFileTag;
	//TFile *outf = new TFile(outputfile.c_str(),"RECREATE");
	TFile *outf = new TFile(OutputFileName.c_str(),"RECREATE");
	
	cout << "________________________________________________________________\n";
	cout << "\n";
	cout << "time start  " << endl;
	gSystem->Exec("date '+%H:%M:%S'");
	
	//----------------------------------------------------------------------------------------
	///////////////////////////////////////////////  Histogram Output  ///////////////////////
	//----------------------------------------------------------------------------------------
	// Book histograms
	TFile hfile("tmva_class_exampleS.root","RECREATE");
	
	TTree *tree = new TTree("TreeS","Example");
	Int_t split = 0;
	Int_t bsize = 64000;
	
	Float_t Pt_lepton;
	Float_t E_t;
	Float_t Pt_btag;
	Float_t Mt_lv;
	Float_t Mi_lvb;
	Float_t N_jets;
	TBranch *var1 = tree->Branch("var1",&Pt_lepton,bsize,split);//Pt_lepton
	TBranch *var2 = tree->Branch("var2",&E_t,bsize,split);// missing et
	TBranch *var3 = tree->Branch("var3",&Pt_btag,bsize,split);//Ptbtag
	TBranch *var4 = tree->Branch("var4",&Mt_lv,bsize,split);//lepton+neutrino enerjisi
	TBranch *var5 = tree->Branch("var5",&Mi_lvb,bsize,split);//lepton+neutrino+btag enerjisi
	TBranch *var6 = tree->Branch("var6",&N_jets,bsize,split);//lepton+neutrino+btag enerjisi
	
	TH1 *histJet_pt[10];
	//TH1 *histJet_eta[10];
	//TH1 *histJet_phi[10];
	char hist_name[100];
	
	//----------------------------------------------------------------------------------------
	TDirectory *jetdir= outf->mkdir("jets");
	jetdir->cd();
	
	for (int i=0; i<10; i++)
	{
		sprintf(hist_name,"jet_pt%i",i);
		histJet_pt[i] = new TH1F(hist_name, "jet P_{T}", 1000, 0.0, 1000.0);
		
		sprintf(hist_name,"jet_eta%i",i);
		//histJet_eta[i] = new TH1F(hist_name, "jet eta", 100, -5.0, 5.0);
		
		sprintf(hist_name,"jet_phi%i",i);
		//histJet_phi[i] = new TH1F(hist_name, "jet phi", 100, -5.0, 5.0);
	}
	
	TH1 *jet_size 		 = new TH1F("jet_size",     "Number of Jets", 10, 0, 10.0);
	TH1 *histJet_btag 	 = new TH1F("histJet_btag", "Number of B-tagged jets", 10, 0, 10.0);
	TH1 *histJet_btag_pt = new TH1F("histJet_btag_pt", "PT of B-tagged jets", 100, 0.0, 500.0);
	TH1 *jet_size_cut8   = new TH1F("jet_size_cut8", "Number of Jets with 30<PT GeV", 10, 0, 10.0);
	TH1 *jet_size_cut9   = new TH1F("jet_size_cut9", "Number of Jets with 15<PT<30 GeV", 10, 0, 10.0);
	TH1 *hist_before_jet_eta = new TH1F("hist_before_jet_eta", "Jet Eta Before cut 9", 100, -5.0, 5.0);
	
	TH1 *jet_alpha_t 	= new TH1F("jet_alpha_t",     "Alpha_T for leading Jets", 100, 0, 2);

	//----------------------------------------------------------------------------------------
	TDirectory *lepdir= outf->mkdir("lepton");
	lepdir->cd();
	
	TH1 *hist_gen_lepton  	     = new TH1F("gen_lepton"		     , "Number of GEN Leptons ", 10, 0, 10.0);
	TH1 *hist_lepton_before_trig = new TH1F("numb_lepton_before_trig", "Number of SIM Leptons |eta|<2.5 and PT>20/30(elec/muon)", 10, 0, 10.0);
	TH1 *hist_lepton_pass_trig   = new TH1F("numb_lepton_pass_trig"	 , "Number of SIM Leptons |eta|<2.5 and PT>20/30(elec/muon)", 10, 0, 10.0);
	TH1 *hist_lepton_pass_10gev  = new TH1F("numb_lepton_pass_10gev" , "Number of SIM Leptons and PT>10(elec/muon)", 10, 0, 10.0);
	
	TH1 *histElec_pt 	= new TH1F("elec_pt1"	 , "1st elec P_{T}", 100, 0.0, 500.0);
	TH1 *histElec_phi	= new TH1F("elec_pt1_phi", "1st elec Phi  ", 100, -5.0, 5.0);
	TH1 *histElec_eta	= new TH1F("elec_pt1_eta", "1st elec Eta  ", 100, -5.0, 5.0);
	
	TH1 *histMuon_pt  = new TH1F("muon_pt1"     , "1st mu P_{T}", 100, 0.0, 500.0);
	TH1 *histMuon_phi = new TH1F("muon_pt1_phi" , "1st mu Phi  ", 100, -5.0, 5.0);
	TH1 *histMuon_eta = new TH1F("muon_pt1_eta" , "1st mu Eta  ", 100, -5.0, 5.0);
	
	TH1 *histLepton_pt 	= new TH1F("lepton_pt"  , "lepton P_{T} ", 100, 0.0, 500.0);
	TH1 *histLepton_eta	= new TH1F("lepton_eta" , "lepton Eta   ", 100, -5.0, 5.0);
	TH1 *histLepton_phi	= new TH1F("lepton_phi" , "lepton  Phi  ", 100, -5.0, 5.0);
	
	TH1 *lepton_invmass = new TH1F("lepton_invmass", "lepton inv mass", 100, 0, 500);
	TH1 *top_invmass 	= new TH1F("top_invmass"   , "top inv mass   ", 100, 0, 500);
	
	TH2 *bjet_lepton_delta_eta 	= new TH2F("bjet_lepton_delta_eta"   , "Delta Eta between lepton and btagjet ", 100, 0, 250, 100, 0,   +5);
	TH2 *bjet_lepton_delta_phi 	= new TH2F("bjet_lepton_delta_phi"   , "Delta Phi between lepton and btagjet ", 100, 0, 250, 100, 0, +6.3);
	TH2 *bjet_lepton_deltaR     = new TH2F("bjet_lepton_deltaR"      , "Delta R between lepton and btagjet   ", 100, 0, 250, 100, 0,  +10);
	
	//----------------------------------------------------------------------------------------
	TDirectory *metdir= outf->mkdir("met");
	metdir->cd();
	
	TH1 *histMET_et  = new TH1F("histMET_et" , "MET",  500,  0.0, 500.0);
	TH1 *histMET_eta = new TH1F("histMET_eta", "MET Eta", 100, -5.0, 5.0);
	TH1 *histMET_phi = new TH1F("histMET_phi", "histMET_phi Phi", 100, -5.0, 5.0);
	
	//----------------------------------------------------------------------------------------
	////////////////////////////////////////  My Data STructure  /////////////////////////////
	//----------------------------------------------------------------------------------------
	
	TH1F::SetDefaultSumw2(true);
	
	gSystem->Load("libExRootAnalysis");
	gSystem->Load("libDelphes");
	
	//string OutputFileName="skimming_delphes.root";
	//TFile *outf = new TFile(OutputFileName.c_str(),"RECREATE");
	
	int loop1[5]={0, 10, 20, 30, 40};
	int loop11[4]={4114, 120428, 104556, 7458};
	int iiii=0;
	while (iiii!=4){
	
	//ofstream myfile;
	//myfile.open ("skimming_mass70_1.txt");
	
	// Create chain of root trees
	TChain chain("Delphes");
	//chain.Add(InputFileName.c_str());
	char filename[1000];
	FILE *input;
	input = fopen(InputFileName.c_str(),"r");
	if (input != NULL)
	{
	  
	  int iii=0;
		// lets read each line and get the filename from it
		while ((fscanf(input,"%s\n",filename) != EOF) && iii !=loop1[iiii+1])
		{
		  iii++;
		  if (iii>loop1[iiii]){
			printf("%s\n",filename);
			chain.Add(filename);
		  }
		}
	}
	
	// Create object of class ExRootTreeReader
	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
	Long64_t numberOfEntries = loop11[iiii]; //treeReader->GetEntries();
	
	// Get pointers to branches used in this analysis
	TClonesArray *branchJet 	= treeReader->UseBranch("Jet");
	TClonesArray *branchParticle = treeReader->UseBranch("Particle");
	TClonesArray *branchElectron = treeReader->UseBranch("Electron");
	TClonesArray *branchMuon 	= treeReader->UseBranch("Muon");
	//TClonesArray *branchPhoton 	= treeReader->UseBranch("Photon");
	TClonesArray *branchMET 	= treeReader->UseBranch("MissingET");
	
	//----------------------------------------------------------------------------------------
	/////////////////////////////////////  LOOP  Over the EVENTS  ////////////////////////////
	//----------------------------------------------------------------------------------------
	
	
	cout << "Set Branch Addresses" << endl;
	
	int decade = 0;
	unsigned entries = 0;
	//unsigned int counter =0;
	
	if (NENTRIES==-1) entries=numberOfEntries;
	else entries=NENTRIES;
	
	outf->cd();
	TVectorD v(entries);
	v[0] = entries;
	v.Write("nevent");
	
	cout << "Reading TREE: " << numberOfEntries << " events available, \n";
	cout << "\t\t" << NENTRIES << " of them will be analyzed." << endl;
	
	GenParticle *part, *mother;//, *status;
	Jet *jet[10];
	Electron *elec[4]; //, *elec2, *elec3;
	Muon *mu[4]; //, *mu2, *mu3;
	MissingET *met;
	
	int event_counter1_after_trigger=0;
	int event_counter2_after_10gev=0;
	int event_counter3_after_leptonpt55=0;
	int event_counter4_after_met=0;
	int event_counter5_after_btag=0;
	int event_counter6_after_topinvmass=0;
	int event_counter7_after_leptinvmass=0;
	int event_counter8_after_onejet=0;
	int event_counter9_after_onejet=0;
	int event_counter10_after_jeteta=0;
	Float_t pt_bt=0;
	// Loop over all events
	for(Long64_t entry = 0; entry < entries; ++entry)
	{
		double progress = 10.0*entry/(1.0*entries);
		int k = TMath::FloorNint(progress);
		if (k > decade) {   cout << 10*k << "%\t"; gSystem->Exec("date '+%H:%M:%S'"); cout << endl;	}
		decade = k;
		
		// Load selected branches with data from specified event
		treeReader->ReadEntry(entry);
		
		//////////////////////////////////////////
		// TRIGGER EMULATION
		//////////////////////////////////////////
		
		int numb_jet=0;
		int counter_btag=0;
		Electron * this_elec=0;
		Muon     * this_muon=0;
		//int numb_muon=0;
		
		//filter lepton=1 events
		int numb_elec_pass_cuts=0;
		int numb_muon_pass_cuts=0;
		
		if( branchElectron->GetEntries() > 0)
		{
			for(int i=0; i<branchElectron->GetEntries(); i++)
			{
				elec[i] = (Electron *) branchElectron->At(i);
				if(  (elec[i]->PT)>30 && abs(elec[i]->Eta)<2.5 )
				{
					numb_elec_pass_cuts++;
					this_elec=elec[i];
				}
			}
		}
		
		if( branchMuon->GetEntries() > 0)
		{
			for(int i=0; i<branchMuon->GetEntries(); i++)
			{
				mu[i] = (Muon *) branchMuon->At(i);
				if(  (mu[i]->PT)>20 && abs(mu[i]->Eta)<2.5 )
				{
					numb_muon_pass_cuts++;
					this_muon=mu[i];
				}
			}
		}
		
		// number of lepton before the trigger
		hist_lepton_before_trig->Fill((numb_muon_pass_cuts + numb_elec_pass_cuts));
		
		// if the total number of muons is higher than 0 fire for this event
		if ((numb_muon_pass_cuts + numb_elec_pass_cuts)==0) continue;
		event_counter1_after_trigger++;
		
		// number of leptons after the trigger cut
		hist_lepton_pass_trig->Fill(numb_muon_pass_cuts + numb_elec_pass_cuts);
		
		//////////////////////////////////////////
		//////////////////////////////////////////
		//	ELECTRON AND MUON BRANCH
		//////////////////////////////////////////
		
		// filter the event if there is lepton PT>55
		
		if ( (numb_elec_pass_cuts+numb_muon_pass_cuts) != 1 ) continue;
		event_counter3_after_leptonpt55++;
		
		if(numb_elec_pass_cuts==1)
		{
			//cout << this_elec->PT << endl;
			histElec_pt->Fill(this_elec->PT);
			histElec_eta->Fill(this_elec->Eta);
			histElec_phi->Fill(this_elec->Phi);
			histLepton_pt->Fill(this_elec->PT);
			histLepton_eta->Fill(this_elec->Eta);
			histLepton_phi->Fill(this_elec->Phi);
		}
		else
		if(numb_muon_pass_cuts ==1 )
		{
				//mu1 = (Muon *) branchMuon->At(this_muon);
				//cout << mu1->PT << endl;
				histMuon_pt->Fill(this_muon->PT);
				histMuon_eta->Fill(this_muon->Eta);
				histMuon_phi->Fill(this_muon->Phi);
				histLepton_pt->Fill(this_muon->PT);
				histLepton_eta->Fill(this_muon->Eta);
				histLepton_phi->Fill(this_muon->Phi);
		}
		
		
		
		//////////////////////////////////////////
		
		// TODO
		// make a comparison between the signal MET and background MET
		// make a comparison and employ a cut on the background
		
		
		double this_met=0;
		if(branchMET->GetEntries() > 0)
		{
			met = (MissingET *) branchMET->At(0);
			histMET_et->Fill(met->MET);
			histMET_eta->Fill(met->Eta);
			histMET_phi->Fill(met->Phi);
			this_met=met->MET;
			
		}
		
		// TODO filter over this met ???
		//		if (this_met < ) continue;
		// FILTER events having MET <50
		//if (this_met < 50. ) continue;
		event_counter4_after_met++;
		
		
		////////////////////////////////////////////////////////////////////
		//	JET BRANCH
		////////////////////////////////////////////////////////////////////
		Jet *this_jet;

		double  alpha_t;
		if(branchJet->GetEntries() > 1 )
		{
		  
		    jet[0] = (Jet*) branchJet->At(0);
		    jet[1] = (Jet*) branchJet->At(1);
		  
		    //double jet_inv_mass12=sqrt(pow(jet[1]->ET+jet[2]->ET,2)-pow(jet[1]->PT+jet[2]->PT,2) );
		    double jet_inv_mass12=sqrt(2*jet[0]->PT*jet[1]->PT*(cosh(jet[0]->Eta-jet[1]->Eta)-cos(jet[0]->Phi-jet[1]->Phi)));
		    alpha_t=(double)(jet[1]->PT)/jet_inv_mass12; 
		 }
		 
		 // go to the end of the loop for alpha_t>0.5, 
		 // because we assume these are qcd events
		 if(alpha_t <0.5) continue;
		 
		 
		 
		 
		
		
		
		int counter_btagfilter=0;
		// filter bjet tagged events
		if(branchJet->GetEntries() > 0)
		{
			for (int i=0; i< branchJet->GetEntries(); i++)
			{
				jet[i] = (Jet*) branchJet->At(i);
				if (jet[i]->BTag==1)
				{
					histJet_btag_pt->Fill(jet[i]->PT);
					counter_btag++;
					if (jet[i]->PT < 100)
					{
						this_jet=jet[i];
						counter_btagfilter++;
					}
				}
			}
			histJet_btag->Fill(counter_btag);
		}
		// or the number of b-tagged jets is not 1
		if ( counter_btagfilter!=1 ) continue;
		event_counter5_after_btag++;
		
		
		
		
		
	///////////////////////////////////////////////////////////////////////////
	///////////////////////// end of PRESELCTION
	///////////////////////////////////////////////////////////////////////////
	
	
		if(numb_elec_pass_cuts==1)
		{
			//cout << this_elec->PT << endl;
			Pt_lepton=this_elec->PT;
		}
		else
		if(numb_muon_pass_cuts ==1 )
		{
				//mu1 = (Muon *) branchMuon->At(this_muon);
				//cout << mu1->PT << endl;
				Pt_lepton=this_muon->PT;
		}
	
		


		if (numb_elec_pass_cuts==1)
		{
			double deltaR2=pow(abs( (this_jet->Eta)-(this_elec->Eta) ),2)+pow(abs( (this_jet->Phi)-(this_elec->Phi) ),2);
			bjet_lepton_delta_eta->Fill( this_elec->PT, abs( (this_jet->Eta)-(this_elec->Eta) ) );
			bjet_lepton_delta_phi->Fill( this_elec->PT, abs( (this_jet->Phi)-(this_elec->Phi) ));
			bjet_lepton_deltaR->Fill( this_elec->PT, sqrt(deltaR2) );
		}
		else if (numb_muon_pass_cuts==1)
		{
			double deltaR2=pow( (this_jet->Eta)-(this_muon->Eta) ,2)+pow( (this_jet->Phi)-(this_muon->Phi) ,2);
			bjet_lepton_delta_eta->Fill( this_muon->PT, abs( (this_jet->Eta)-(this_muon->Eta) ) );
			bjet_lepton_delta_phi->Fill( this_muon->PT, abs( (this_jet->Phi)-(this_muon->Phi) ));
			bjet_lepton_deltaR->Fill( this_muon->PT, sqrt(deltaR2) );
		}
		
		// TODO make a comparison plot for batg-jet
		// multiplicity and pt distribution before and after
		
		E_t = this_met;
		
		//////////////////////////////////////////////////////////////////////////////////////
		// TOP INVARIANT MASS and LEPTONIC TRANSVERSE MASS
		//////////////////////////////////////////////////////////////////////////////////////
		double topinvariantmass=0;
		double leptoninvariantmass=0;
		
		if ( numb_muon_pass_cuts == 1 )	// there is a elec in the event
		{
			//cout << "e  "<< top_invariant_mass(elec1, this_met);
			topinvariantmass=top_invariant_mass(this_muon, met);
			top_invmass->Fill(topinvariantmass);
			Mi_lvb=topinvariantmass;
		}
		else
			if ( numb_elec_pass_cuts == 1 )	// there is a muon in the event
			{
				//cout << "m  "<< top_invariant_mass(mu1, this_met);
				topinvariantmass=top_invariant_mass(this_elec, met);
				top_invmass->Fill(topinvariantmass);
				Mi_lvb=topinvariantmass;
			}
		//cout<< endl;
		
		if ( numb_muon_pass_cuts == 1 )	// there is a elec in the event
		{
			//cout << "e  "<< lepton_invariant_mass(elec1, this_met);
			leptoninvariantmass=lepton_invariant_mass(this_muon, met);
			lepton_invmass->Fill(leptoninvariantmass);
			Mt_lv=leptoninvariantmass;
		}
		else
			if ( numb_elec_pass_cuts == 1 )	// there is a muon in the event
			{
				//cout << "m  "<< lepton_invariant_mass(mu1, this_met);
				leptoninvariantmass=lepton_invariant_mass(this_elec, met);
				lepton_invmass->Fill(leptoninvariantmass);
				Mt_lv=leptoninvariantmass;
			}
		
		
		/*if ( topinvariantmass < 280  ) continue;
		event_counter6_after_topinvmass++;
		
		if ( topinvariantmass < 85  ) continue;
		event_counter7_after_leptinvmass++;*/
		
		//cout<< endl;
		
		////////////////////////////////////////////////////////////////////
		//	JET BRANCH
		////////////////////////////////////////////////////////////////////
		
		numb_jet=0;
		if(branchJet->GetEntries() > 0)
		{
			for (int i=0; i< branchJet->GetEntries(); i++)
			{
				jet[i] = (Jet*) branchJet->At(i);
				if ( (jet[i]->PT)>20 && abs(jet[i]->Eta)< 4.9 ) numb_jet++;
			}
		}
		
		Pt_btag=this_jet->PT;
		
		N_jets=numb_jet;
		
		tree->Fill();
		myfile << iiii << "," << Mt_lv;
		myfile << "," << E_t;
		myfile << "," << Pt_btag;
		myfile << "," << Pt_lepton;
		myfile << "," << Mi_lvb;
		myfile << "," << N_jets << "\n";
	
		
	}	//end of event loop
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	
	//---------------------------------------------------------------------------------------------------------
	////////////////////////////////  Saving the histograms into output  //////////////////////////////////////////
	//---------------------------------------------------------------------------------------------------------
	// end of event loop
	/*
	cout << "nevents : " << entries << endl;
	cout << "event after trigger     : "    << event_counter1_after_trigger     << "\t" << (double)event_counter1_after_trigger/entries <<  endl;
	cout << "event after 10gev       : "    << event_counter2_after_10gev       << "\t" << (double)event_counter2_after_10gev/entries << endl;
	cout << "event after Lep(PT)>55  : "    << event_counter3_after_leptonpt55  << "\t" << (double)event_counter3_after_leptonpt55/entries << endl;
	cout << "event after met > 50    : "    << event_counter4_after_met         << "\t" << (double)event_counter4_after_met/entries << endl;
	cout << "event after btag = 1    : "    << event_counter5_after_btag        << "\t" << (double)event_counter5_after_btag/entries << endl;
	cout << "event after topinvmass  : "    << event_counter6_after_topinvmass  << "\t" << (double)event_counter6_after_topinvmass/entries << endl;
	cout << "event after leptinvmass : "    << event_counter7_after_leptinvmass << "\t" << (double)event_counter7_after_leptinvmass/entries << endl;
	cout << "event after numb_jet =1 : "    << event_counter8_after_onejet      << "\t" << (double)event_counter8_after_onejet/entries << endl;
	cout << "event after numb_jet =1 : "    << event_counter9_after_onejet      << "\t" << (double)event_counter9_after_onejet/entries << endl;
	cout << "event after jet_eta<4.9 : "    << event_counter10_after_jeteta     << "\t" << (double)event_counter10_after_jeteta/entries << endl;
	*/
	cout <<"CELAL" << "," << (double)event_counter1_after_trigger/entries << "," << (double)event_counter3_after_leptonpt55/entries<< "," << (double)event_counter4_after_met/entries<< "," << (double)event_counter5_after_btag/entries<< "\n";
	
	//cout<<pt_bt<<endl;
	iiii++;
	}
	
	tree->Write("",TObject::kOverwrite);
	hfile.Write("",TObject::kOverwrite);
	//myfile.close();
	outf->Write();
	outf->Close();
	myfile.close();
	
	//end of main loop
}


# TMVA-for-High-Energy
# TMVA-for-High-Energy
