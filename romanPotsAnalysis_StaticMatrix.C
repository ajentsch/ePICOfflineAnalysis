//-------------------------
//
// Simple analysis code to analyze EPIC simulation output for Roman Pots
//
// Author: Alex Jentsch
//
// Date of last author update: 2/16/2022
//
//------------------------


using namespace std;


void romanPotsAnalysis_StaticMatrix(){

	TString fileList = "./inputFileList_DVCS_18x275_updated.list";
	
	TString outputName = "EPIC_DVCS_Simulation_Output_";	

	TString date = "2_16_2023_";
	
	TString run  = "run_0";

	cout << "Input FileList: " << fileList << endl;
	TString fileType_ROOT = ".root";
	TString outputFileName = outputName + date + run + fileType_ROOT;
	string fileName;
	TFile * inputRootFile;
	TTree * rootTree;
	cout << "Output file: " << outputFileName << endl;


	ifstream fileListStream;
	fileListStream.open(fileList);
	if(!fileListStream) { cout << "NO_LIST_FILE " << fileList << endl; return;}

	
    //---------------------Roman Pots reconstruction constants------------------
	
	//N.B. this is all bullshit for right now while we solve the online reco problem
	
    double local_x_offset_station_1 = -833.3878326;
    double local_x_offset_station_2 = -924.342804;
    double local_x_slope_offset = -0.00622147;
    double local_y_slope_offset = -0.0451035;
    double crossingAngle; // -0.025
    double nomMomentum = 275.0;
	
	//simplified transfer matrix for X_L = 1.0 (beam orbit)

    const double aXRP[2][2] = {{2.102403743, 29.11067626},
                               {0.186640381, 0.192604619}};
    const double aYRP[2][2] = {{0.0000159900, 3.94082098},
                               {0.0000079946, -0.1402995}};

    double aXRPinv[2][2] = {{0.0, 0.0},
                            {0.0, 0.0}};
    double aYRPinv[2][2] = {{0.0, 0.0},
                            {0.0, 0.0}};
	
	double det = aXRP[0][0] * aXRP[1][1] - aXRP[0][1] * aXRP[1][0];
	
	if (det == 0) {
		cout << "ERROR: Reco matrix determinant = 0! Matrix cannot be inverted! Double-check matrix!" << endl;
		return;
	}

    aXRPinv[0][0] = aXRP[1][1] / det;
    aXRPinv[0][1] = -aXRP[0][1] / det;
    aXRPinv[1][0] = -aXRP[1][0] / det;
    aXRPinv[1][1] = aXRP[0][0] / det;

    det = aYRP[0][0] * aYRP[1][1] - aYRP[0][1] * aYRP[1][0];
    aYRPinv[0][0] = aYRP[1][1] / det;
    aYRPinv[0][1] = -aYRP[0][1] / det;
    aYRPinv[1][0] = -aYRP[1][0] / det;
    aYRPinv[1][1] = aYRP[0][0] / det;

	//--------------------------------------------------------------------------

	//histograms
	
    TH1D* eta_MC = new TH1D("h_eta",";Pseudorapidity, #eta",100,0.0,15.0);
    TH1D* energy_MC = new TH1D("h_energy_MC",";E_{MC} [GeV]",100,0,300);
	TH1D* px_MC = new TH1D("px_MC", ";p_{x} [GeV/c]", 100, -10.0, 10.0);
	TH1D* py_MC = new TH1D("py_MC", ";p_{y} [GeV/c]", 100, -10.0, 10.0);
	TH1D* pt_MC = new TH1D("pt_MC", ";p_{t} [GeV/c]", 100, 0.0, 2.0);
	TH1D* pz_MC = new TH1D("pz_MC", ";p_{z} [GeV/c]", 100, 0.0, 320.0);
	
	TH2D* dvcs_photon_thetaPhi_MC = new TH2D("dvcs_photon_theta_phi_MC", ";Pseudorapidity, #eta; Azimuthal Angle, #phi [rad.]", 100, -4.0, 4.0, 100, -TMath::Pi()-0.1, TMath::Pi()+0.1);
	TH2D* dvcs_proton_thetaPhi_MC = new TH2D("dvcs_proton_theta_phi_MC", ";Pseudorapidity, #eta; Azimuthal Angle, #phi [rad.]", 100, -10.0, 10.0, 100, -TMath::Pi()-0.1, TMath::Pi()+0.1);
	
	TH2D* emcal_hit_map_theta_phi_reco = new TH2D("emcal_hit_map_theta_phi", ";Pseudorapidity, #eta; Azimuthal Angle, #phi [rad.]", 100, -4.0, 4.0, 100, -TMath::Pi()-0.1, TMath::Pi()+0.1);
	
	TH1D* px_smearedFF = new TH1D("px_smearedFF", ";p_{x} [GeV/c]", 100, -10.0, 10.0);
	TH1D* py_smearedFF = new TH1D("py_smearedFF", ";p_{y} [GeV/c]", 100, -10.0, 10.0);
	TH1D* pt_smearedFF = new TH1D("pt_smearedFF", ";p_{t} [GeV/c]", 100, 0.0, 2.0);
	TH1D* pz_smearedFF = new TH1D("pz_smearedFF", ";p_{z} [GeV/c]", 100, 0.0, 320.0);
	
	TH1D* px_RomanPots = new TH1D("px_RomanPots", ";p_{x} [GeV/c]", 100, -10.0, 10.0);
	TH1D* py_RomanPots = new TH1D("py_RomanPots", ";p_{y} [GeV/c]", 100, -10.0, 10.0);
	TH1D* pt_RomanPots = new TH1D("pt_RomanPots", ";p_{t} [GeV/c]", 100, 0.0, 2.0);
	TH1D* pz_RomanPots = new TH1D("pz_RomanPots", ";p_{z} [GeV/c]", 100, 0.0, 320.0);
	
	TH2D* rp_occupancy_map = new TH2D("Roman_pots_occupancy_map", "Roman_pots_occupancy_map", 100, -1100, -700, 100, -70, -70);

	int fileCounter = 0;
	int iEvent = 0;

	while(getline(fileListStream, fileName)){

	    TString tmp = fileName;

	    cout << "Input file " << fileCounter << ": " << fileName << endl;

	    inputRootFile = new TFile(tmp);
	    if(!inputRootFile){ cout << "MISSING_ROOT_FILE"<< fileName << endl; continue;}
		
		fileCounter++;

		TTree * evtTree = (TTree*)inputRootFile->Get("events");


		int numEvents = evtTree->GetEntries();

    	TTreeReader tree_reader(evtTree);       // !the tree reader

		//MC particles
    
    	TTreeReaderArray<float> mc_px_array = {tree_reader, "MCParticles.momentum.x"};
    	TTreeReaderArray<float> mc_py_array = {tree_reader, "MCParticles.momentum.y"};
    	TTreeReaderArray<float> mc_pz_array = {tree_reader, "MCParticles.momentum.z"};
    	TTreeReaderArray<double> mc_mass_array = {tree_reader, "MCParticles.mass"};
    	TTreeReaderArray<int> mc_pdg_array = {tree_reader, "MCParticles.PDG"};
		
		//generated particles
    
    	TTreeReaderArray<float> gen_px_array = {tree_reader, "GeneratedParticles.momentum.x"};
    	TTreeReaderArray<float> gen_py_array = {tree_reader, "GeneratedParticles.momentum.y"};
    	TTreeReaderArray<float> gen_pz_array = {tree_reader, "GeneratedParticles.momentum.z"};
    	TTreeReaderArray<double> gen_mass_array = {tree_reader, "GeneratedParticles.mass"};
    	TTreeReaderArray<int> gen_pdg_array = {tree_reader, "GeneratedParticles.PDG"};


		// Roman Pots
	
   	 	TTreeReaderArray<float> reco_hits_RP_x = {tree_reader, "ForwardRomanPotRecHits.position.x"};
    	TTreeReaderArray<float> reco_hits_RP_y = {tree_reader, "ForwardRomanPotRecHits.position.y"};
    	TTreeReaderArray<float> reco_hits_RP_z = {tree_reader, "ForwardRomanPotRecHits.position.z"};


		// FastSmearing plugin
	
    	TTreeReaderArray<float> reco_fast_FarForward_px = {tree_reader, "SmearedFarForwardParticles.momentum.x"};
    	TTreeReaderArray<float> reco_fast_FarForward_py = {tree_reader, "SmearedFarForwardParticles.momentum.y"};
    	TTreeReaderArray<float> reco_fast_FarForward_pz = {tree_reader, "SmearedFarForwardParticles.momentum.z"};
		TTreeReaderArray<int> reco_fast_FarForward_PDG = {tree_reader, "SmearedFarForwardParticles.PDG"};
		
		//EMCAL branches
		
		TTreeReaderArray<float> barrelEMCAL_phi = {tree_reader, "EcalBarrelSciGlassTruthClusters.intrinsicPhi"};
		TTreeReaderArray<float> barrelEMCAL_theta = {tree_reader, "EcalBarrelSciGlassTruthClusters.intrinsicTheta"};
		
		TTreeReaderArray<float> negEndCapEMCAL_phi = {tree_reader, "EcalEndcapNTruthClusters.intrinsicPhi"};
		TTreeReaderArray<float> negEndCapEMCAL_theta = {tree_reader, "EcalEndcapNTruthClusters.intrinsicTheta"};
		
		TTreeReaderArray<float> posEndCapEMCAL_phi = {tree_reader, "EcalEndcapPTruthClusters.intrinsicPhi"};
		TTreeReaderArray<float> posEndCapEMCAL_theta = {tree_reader, "EcalEndcapPTruthClusters.intrinsicTheta"};
	

		cout << "file has " << evtTree->GetEntries() <<  " events..." << endl;

		
		
		tree_reader.SetEntriesRange(0, evtTree->GetEntries());
		while (tree_reader.Next()) {

			cout << "Reading event: " << iEvent << endl;

	    	//MCParticles
	        //finding the far-forward proton;
	    	//TLorentzVector scatMC(0,0,0,0);
			TVector3 mctrk;
			bool goodHit1 = false;
			bool goodHit2 = false;
			
			//skip events entirely with no RP hits to remove acceptance loss
			
			if(reco_hits_RP_x.GetSize() < 1){iEvent++; continue; }
			
			if(reco_hits_RP_x.GetSize() > 3){
				
				double goodHitX[2] = {0.0, 0.0};
				double goodHitY[2] = {0.0, 0.0};
				double goodHitZ[2] = {0.0, 0.0};
		
				for(int iHit = 0; iHit < reco_hits_RP_x.GetSize(); iHit++){
					
					cout << "hit "<< iHit <<" - (z, x) = (" << reco_hits_RP_z[iHit] << ", " << reco_hits_RP_x[iHit] <<")"<< endl;
					
					if(!goodHit2 && reco_hits_RP_z[iHit] > 27099.0 && reco_hits_RP_z[iHit] < 28022.0){ 
					
						goodHitX[1] = reco_hits_RP_x[iHit];
						goodHitY[1] = reco_hits_RP_y[iHit];
						goodHitZ[1] = reco_hits_RP_z[iHit];
						goodHit2 = true;
					
					}
					if(!goodHit1 && reco_hits_RP_z[iHit] > 25099.0 && reco_hits_RP_z[iHit] < 26022.0){ 
					
						goodHitX[0] = reco_hits_RP_x[iHit];
						goodHitY[0] = reco_hits_RP_y[iHit];
						goodHitZ[0] = reco_hits_RP_z[iHit];
						goodHit1 = true;
					
					}
				}
			}
			
			if(!goodHit1 || !goodHit2){ iEvent++; continue; }
			
	    	double maxPt=-99.;
	    	for(int imc=0;imc<mc_px_array.GetSize();imc++){
	    		mctrk.SetXYZ(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);	
				
	    		if(mc_pdg_array[imc] == 2212 && mctrk.Perp() > 0.1){
	    			
					mctrk.RotateY(0.025);
					
					eta_MC->Fill(mctrk.Eta());
					
					px_MC->Fill(mctrk.Px());
					py_MC->Fill(mctrk.Py());
					pt_MC->Fill(mctrk.Perp());
					pz_MC->Fill(mctrk.Pz());
				}
	    		if(mc_pdg_array[imc] == 22){
	    			
					dvcs_photon_thetaPhi_MC->Fill(mctrk.Eta(), mctrk.Phi());
				}
				
	    	}
			
	    	for(int igen=0;igen<gen_px_array.GetSize();igen++){
	    		mctrk.SetXYZ(gen_px_array[igen], gen_py_array[igen], gen_pz_array[igen]);	
				
	    		if(gen_pdg_array[igen] == 22){
	    			
					dvcs_photon_thetaPhi_MC->Fill(mctrk.Eta(), mctrk.Phi());
				}
	    		if(gen_pdg_array[igen] == 2212){
	    			
					dvcs_proton_thetaPhi_MC->Fill(mctrk.Eta(), mctrk.Phi());
				}
				
	    	}
			
			for(int ical=0;ical<barrelEMCAL_phi.GetSize();ical++){
				
				double radians = barrelEMCAL_theta[ical];
				double eta = -TMath::Log(TMath::Tan(radians/2.0));
		
				emcal_hit_map_theta_phi_reco->Fill(eta, barrelEMCAL_phi[ical]);
			}
			for(int ical=0;ical<negEndCapEMCAL_phi.GetSize();ical++){
		
				double radians = negEndCapEMCAL_theta[ical];
				double eta = -TMath::Log(TMath::Tan(radians/2.0));
				emcal_hit_map_theta_phi_reco->Fill(eta, negEndCapEMCAL_phi[ical]);
			}
			for(int ical=0;ical<posEndCapEMCAL_phi.GetSize();ical++){
		
				double radians = posEndCapEMCAL_theta[ical];
				double eta = -TMath::Log(TMath::Tan(radians/2.0));
				emcal_hit_map_theta_phi_reco->Fill(eta, posEndCapEMCAL_phi[ical]);
			}
		
			
			if(mctrk.Perp() > 0.1){
			
		    	for(int iSFFPart = 0; iSFFPart < reco_fast_FarForward_px.GetSize(); iSFFPart++){
	    		
					TVector3 smearedFFtrk(reco_fast_FarForward_px[iSFFPart], reco_fast_FarForward_py[iSFFPart], reco_fast_FarForward_pz[iSFFPart]);	
	    		
					if(reco_fast_FarForward_PDG[iSFFPart]==2212){
	    			
						smearedFFtrk.RotateY(0.025);
					
						px_smearedFF->Fill(smearedFFtrk.Px());
						py_smearedFF->Fill(smearedFFtrk.Py());
						pt_smearedFF->Fill(smearedFFtrk.Perp());
						pz_smearedFF->Fill(smearedFFtrk.Pz());
					
		    		}
		    	}
			
				goodHit1 = false;
				goodHit2 = false;
			
				//for(int iRpHit = 0; iRpHit < reco_hits_RP_x.GetSize(); iRpHit++){
				if(reco_hits_RP_x.GetSize() > 3){ 
			
					double goodHitX[2] = {0.0, 0.0};
					double goodHitY[2] = {0.0, 0.0};
					double goodHitZ[2] = {0.0, 0.0};
					
					//bool goodHit1 = false;
					//bool goodHit2 = false;
			
					for(int iHit = 0; iHit < reco_hits_RP_x.GetSize(); iHit++){
						
						cout << "hit "<< iHit <<" - (z, x) = (" << reco_hits_RP_z[iHit] << ", " << reco_hits_RP_x[iHit] <<")"<< endl;
						
						if(!goodHit2 && reco_hits_RP_z[iHit] > 27099.0 && reco_hits_RP_z[iHit] < 28022.0){ 
						
							goodHitX[1] = reco_hits_RP_x[iHit];
							goodHitY[1] = reco_hits_RP_y[iHit];
							goodHitZ[1] = reco_hits_RP_z[iHit];
							goodHit2 = true;
						
						}
						if(!goodHit1 && reco_hits_RP_z[iHit] > 25099.0 && reco_hits_RP_z[iHit] < 26022.0){ 
						
							goodHitX[0] = reco_hits_RP_x[iHit];
							goodHitY[0] = reco_hits_RP_y[iHit];
							goodHitZ[0] = reco_hits_RP_z[iHit];
							goodHit1 = true;
						
						}
					}
					
					if(!goodHit1 || !goodHit2){ continue; }
				
					rp_occupancy_map->Fill(goodHitX[0], goodHitY[0]); 
					
					cout << "good hit 1st layer: " << " - (z, x) = (" << goodHitZ[0] << ", " << goodHitX[0] <<")"<< endl;
					cout << "good hit 2nd layer: " << " - (z, x) = (" << goodHitZ[1] << ", " << goodHitX[1] <<")"<< endl;
					
			        // extract hit, subtract orbit offset – this is to get the hits in the coordinate system of the orbit
		            // trajectory
		            double XL[2] = {goodHitX[0] - local_x_offset_station_1, goodHitX[1] - local_x_offset_station_2};
		            double YL[2] = {goodHitY[0], goodHitY[1]};

		            double base = goodHitZ[1] - goodHitZ[0];

					cout << "Z-distance = " << base << endl;

		            if (base == 0) {
		                cout << "ERROR: Detector separation = 0! Cannot calculate slope!" << endl;
		                continue;
		            }

		            double Xip[2] = {0.0, 0.0};
		            double Xrp[2] = {XL[1], (1000 * (XL[1] - XL[0]) / (base)) - local_x_slope_offset}; //- _SX0RP_;
		            double Yip[2] = {0.0, 0.0};
		            double Yrp[2] = {YL[1], (1000 * (YL[1] - YL[0]) / (base)) - local_y_slope_offset}; //- _SY0RP_;

		            // use the hit information and calculated slope at the RP + the transfer matrix inverse to calculate the
		            // Polar Angle and deltaP at the IP

		            for (unsigned i0 = 0; i0 < 2; i0++) {
		                for (unsigned i1 = 0; i1 < 2; i1++) {
		                    Xip[i0] += aXRPinv[i0][i1] * Xrp[i1];
		                    Yip[i0] += aYRPinv[i0][i1] * Yrp[i1];
		                }
		            }

		            // convert polar angles to radians
		            double rsx = Xip[1] / 1000.;
		            double rsy = Yip[1] / 1000.;

		            // calculate momentum magnitude from measured deltaP – using thin lens optics.
		            double p = nomMomentum * (1 + 0.01 * Xip[0]);
		            double norm = std::sqrt(1.0 + rsx * rsx + rsy * rsy);

		            //double prec[3] = {(p * rsx / norm), (p * rsy / norm), (p / norm)};
				
					TVector3 prec_romanpots((p * rsx / norm), (p * rsy / norm), (p / norm)); 
				
					px_RomanPots->Fill(prec_romanpots.Px());
					py_RomanPots->Fill(prec_romanpots.Py());
					pt_RomanPots->Fill(prec_romanpots.Perp());
					pz_RomanPots->Fill(prec_romanpots.Pz());
				
				}
			}
			iEvent++;
		}// event loop
		
		inputRootFile->Close();
		
	}// input file loop
		
	cout << "Check integrals: " << endl;
	cout << "pt_mc integral = " << pt_MC->Integral() << endl;
	cout << "pt_RP_reco integral = " << pt_RomanPots->Integral() << endl;

	TFile * outputFile = new TFile(outputFileName, "RECREATE");

	eta_MC->Write();
	px_MC->Write();
	py_MC->Write();
	pt_MC->Write();
	pz_MC->Write();
	px_smearedFF->Write();
	py_smearedFF->Write();
	pt_smearedFF->Write();
	pz_smearedFF->Write();
	px_RomanPots->Write();
	py_RomanPots->Write();
	pt_RomanPots->Write();
	pz_RomanPots->Write();
	
	energy_MC->Write();
	
	rp_occupancy_map->Write();
	
	dvcs_photon_thetaPhi_MC->Write();
	dvcs_proton_thetaPhi_MC->Write();
	emcal_hit_map_theta_phi_reco->Write();

	outputFile->Close();

	

    return;

}

