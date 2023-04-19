/*************************************************************************
    > File Name: _Read_file.cpp
    > Author: Ma PengXiong
    > Mail: 1182905413@qq.com 
    > Created Time: Tue 04 Apr 2023 12:57:07 PM CST
 ************************************************************************/

#include <fstream>
#include <iostream>
#include "time.h"
#include "string.h"
#include "TFile.h"
#include "TSystem.h"
#include "TString.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMinuit.h"
#include "TAxis.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TVector3.h"
#include "TStyle.h"
#include <vector>
#include "TVirtualFFT.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
using namespace std;
void _Readfile(const char* , const char* ,int ,int );

int main(int argc,const char* argv[]){
 
	int plotif = atoi(argv[3]);
	int duid = atoi(argv[4]);
	_Readfile(argv[1],argv[2],plotif,duid);
	return 0;

}

void _Readfile(const char* input, const char* output, int pic, int du)
{
	TCanvas *TC = new TCanvas("cc","tc",900,300);
	TC->SetFillColor(45);
	TString prefixname = "/home/mapx/mapx/Reco/ROOTFile/";
	TString prefixname1;
	TC->Print(Form(prefixname+"%s"+".pdf[",output ));
	TFile *f1 = new TFile(Form(prefixname1+"%s",input));
	TTree *trevoltage = (TTree*)f1->Get("teventvoltage");
	TTree *treadc = (TTree*)f1->Get("teventadc");
	TFile *fcre = new TFile(Form(prefixname+"%s",output ),"recreate");
	TTree *tree = new TTree("_tree","Comparison with old");
	TTree *tree_eff = new TTree("_trigger_eff","trigger effective");
	int Nentries = treadc->GetEntries();
	float RecoZenith;
	float RecoAzimuth;
	float Delta_time;
	float Delta_T0[13]={-9999};
	float Delta_P2P[13]={-999999};
	float Delta_XP2P[13]={-999999};
	float Delta_YP2P[13]={-999999};
	float Delta_ZP2P[13]={-999999};
	float P2P_Value[13]={-999999};
	float TZero[13]={-9999999};
	float P2P_T0[13]={-9999};
	float Ant_Dis[13]={-9999};
	float RealZenith;
	float RealAzimuth;
	float FootSize;
	float PrimaryE;
	float P2P =0.0;
	float covariance2 =0.0;
	float CPUTime =0.0;
	int TriggerN = 0;
	float D2Xmax=-1e4;
	float Eff_trigger = -99999;
	float N_tr = 0;
	float Pri_E = -9999999;
	float Pri_Theta = -9999999;
	float Pri_Phi = -9999999;
	int TotalDU[5]={0};
	tree_eff->Branch("Eff_trigger",&Eff_trigger,"Eff_trigger/F");
	tree_eff->Branch("Pri_E",&Pri_E,"Pri_E/F");
	tree_eff->Branch("Pri_Theta",&Pri_Theta,"Pri_Theta/F");
	tree_eff->Branch("Pri_Phi",&Pri_Phi,"Pri_Phi/F");

	tree->Branch("RecoZenith",&RecoZenith,"RecoZenith/F");
	tree->Branch("RecoAzimuth",&RecoAzimuth,"RecoAzimuth/F");
	tree->Branch("RealZenith",&RealZenith,"RealZenith/F");
	tree->Branch("RealAzimuth",&RealAzimuth,"RealAzimuth/F");
	tree->Branch("Delta_time",&Delta_time,"Delta_time/F");
	tree->Branch("Delta_T0[13]",Delta_T0,"Delta_T0[13]/F");
	tree->Branch("TZero[13]",TZero,"TZero[13]/F");
	tree->Branch("Delta_P2P[13]",Delta_P2P,"Delta_P2P[13]/F");
	tree->Branch("Delta_XP2P[13]",Delta_XP2P,"Delta_XP2P[13]/F");
	tree->Branch("Delta_YP2P[13]",Delta_YP2P,"Delta_YP2P[13]/F");
	tree->Branch("Delta_ZP2P[13]",Delta_ZP2P,"Delta_ZP2P[13]/F");
	tree->Branch("P2P_Value[13]",P2P_Value,"P2P_Value[13]/F");
	tree->Branch("P2P_T0[13]",P2P_T0,"P2P_T0[13]/F");
	tree->Branch("Ant_Dis[13]",Ant_Dis,"Ant_Dis[13]/F");
	tree->Branch("TotalDU[5]",TotalDU,"TotalDU[5]/I");
	tree->Branch("FootSize",&FootSize,"FootSize/F");
	tree->Branch("PrimaryE",&PrimaryE,"PrimaryE/F");
	tree->Branch("P2P",&P2P,"P2P/F");
	tree->Branch("CPUTime",&CPUTime,"CPUTime/F");
	tree->Branch("TriggerN",&TriggerN,"TriggerN/I");
	tree->Branch("D2Xmax",&D2Xmax,"D2Xmax/F");
	tree->Branch("covariance2",&covariance2,"covariance2/F");
	TH1D *Histo = new TH1D("_Signal",";Signal;",500,0,5000);
	char FileName[63];
	vector<vector<short>> trace_x;
	vector<vector<short>> trace_y;
	vector<vector<short>> trace_z;
	vector<vector<short>> trace_a;
	vector<float> pos_x;
	vector<float> pos_y;
	vector<float> pos_z;
	vector<unsigned short> du_id;
	vector<unsigned short> gps_temp;
	vector<unsigned short> atm_temperature;
	vector<unsigned short> atm_pressure;
	vector<unsigned short> acceleration_x;
	vector<unsigned short> acceleration_y;
	vector<unsigned short> acceleration_z;
	vector<unsigned short> battery_level;
	vector<unsigned short> atm_humidity;
	vector<unsigned int> gps_time;
	vector<unsigned int> du_nanoseconds;
	double xmax_alt= -9999;
	double xmax_pos_shc[3] = {0};
	float shower_core_pos[4]={0};
	float ground_alt;
	float shower_energy;
	float shower_zenith;
	float shower_azimuth;
	vector <unsigned short> Detectors_det_id;
	float magnetic_field[3];
	vector<float> pos ;
	int Nstar=0;
	TF1* fx;
	TF1* fy;
	TF1* fz;

	vector<vector<short>> *branch_trace_x = &trace_x;
	vector<vector<short>> *branch_trace_y = &trace_y;
	vector<vector<short>> *branch_trace_z = &trace_z;
	vector<vector<short>> *branch_trace_a = &trace_a;
	vector<float> *branch_pos_x = &pos_x;
	vector<float> *branch_pos_y = &pos_y;
	vector<float> *branch_pos_z = &pos_z;
	vector<unsigned short> *branch_du_id = &du_id;
	vector<unsigned int> *du_gps = &gps_time;
	vector<unsigned short> *du_temp = &gps_temp;
	vector<unsigned short> *at_hum = &atm_humidity;
	vector<unsigned short> *at_temp = &atm_temperature;
	vector<unsigned short> *at_press = &atm_pressure;
	vector<unsigned short> *du_ax = &acceleration_x;
	vector<unsigned short> *du_ay = &acceleration_y;
	vector<unsigned short> *du_az = &acceleration_z;
	vector<unsigned short> *at_bat = &battery_level;
	vector<unsigned int> *du_gps_nano = &du_nanoseconds;
//	trevoltage->SetBranchAddress("pos_x",&branch_pos_x);//empty is not avaible to read from file!!!!! aborted
//	trevoltage->SetBranchAddress("pos_y",&branch_pos_y);
//	trevoltage->SetBranchAddress("pos_z",&branch_pos_z);
	treadc->SetBranchAddress("trace_1",&branch_trace_x);
	treadc->SetBranchAddress("trace_2",&branch_trace_y);
	treadc->SetBranchAddress("trace_3",&branch_trace_z);
	treadc->SetBranchAddress("trace_0",&branch_trace_a);
	treadc->SetBranchAddress("du_id",&branch_du_id);
	treadc->SetBranchAddress("gps_time",&du_gps);
	treadc->SetBranchAddress("gps_temp",&du_temp);
	treadc->SetBranchAddress("du_nanoseconds",&du_gps_nano);
	treadc->SetBranchAddress("atm_humidity",&at_hum);
	treadc->SetBranchAddress("atm_temperature",&at_temp);
	treadc->SetBranchAddress("acceleration_x",&du_ax);
	treadc->SetBranchAddress("acceleration_y",&du_ay);
	treadc->SetBranchAddress("acceleration_z",&du_az);
	treadc->SetBranchAddress("battery_level",&at_bat);
	treadc->SetBranchAddress("atm_pressure",&at_press);
	int totot_id=0;
	for(int i=Nstar ;i< Nentries ;i+=1)
	{
	    
		vector<vector<short>> Detectors_trace_V1;
		vector<vector<short>> Detectors_trace_V2;
		vector<vector<short>> Detectors_trace_V3;
		vector<vector<short>> Detectors_trace_V0;
		vector<vector<short>> Detectors_t_0 ;

//		treSimu->GetEntry(i);
		treadc->GetEntry(i);
		double xmaxDis[3];
		
		xmaxDis[0] = xmax_pos_shc[0];
		xmaxDis[1] = xmax_pos_shc[1];
		xmaxDis[2] = xmax_pos_shc[2];
		float distanceXmax = sqrt( pow(xmaxDis[0] - shower_core_pos[0],2.0) + pow(xmaxDis[1] - shower_core_pos[1],2.0) + pow(xmaxDis[2] - ground_alt,2.0) )/1e3;

		//cout<<"distanceXmax "<<distanceXmax<<endl;
		float XmaxAlt = xmax_alt/1e3;

		Detectors_trace_V1 =trace_x;
		Detectors_trace_V2 =trace_y;
		Detectors_trace_V3 =trace_z;
		Detectors_trace_V0 =trace_a;
		int ite=0;
		int ntrace = Detectors_trace_V1.size();
		if(du_id[0]==1078) TotalDU[0]+=1;
		if(du_id[0]==1080) TotalDU[1]+=1;
		if(du_id[0]==1081) TotalDU[2]+=1;
		if(du_id[0]==1082) TotalDU[3]+=1;
		if(du_id[0]==1094) TotalDU[4]+=1;
		if(i != Nentries-1 )
			continue;
		cout<<TotalDU[0]<<"\t"<<TotalDU[1]<<"\t"<<TotalDU[2]<<"\t"<<TotalDU[3]<<"\t"<<TotalDU[4]<<endl;
	    TH1D* noiseX[ntrace];// = new TH1D(Form("NoiseX_%d_%d",i,it),"_XNoise",200,-1000,1000);
	    TH1D* noiseY[ntrace];// = new TH1D(Form("NoiseY_%d_%d",i,it),"_YNoise",200,-1000,1000);
	    TH1D* noiseZ[ntrace];// = new TH1D(Form("NoiseZ_%d_%d",i,it),"_ZNoise",200,-1000,1000);
	    TH1D* noiseA[ntrace];// = new TH1D(Form("NoiseZ_%d_%d",i,it),"_ZNoise",200,-1000,1000);
	    TH1D* noiseX_2[ntrace];// = new TH1D(Form("NoiseX_%d_%d",i,it),"_XNoise",200,-1000,1000);
	    TH1D* noiseY_2[ntrace];// = new TH1D(Form("NoiseY_%d_%d",i,it),"_YNoise",200,-1000,1000);
	    TH1D* noiseZ_2[ntrace];// = new TH1D(Form("NoiseZ_%d_%d",i,it),"_ZNoise",200,-1000,1000);
	    TH1D* noiseA_2[ntrace];// = new TH1D(Form("NoiseZ_%d_%d",i,it),"_ZNoise",200,-1000,1000);
		TH1F *adcx[ntrace];
		TH1F *adcy[ntrace];
		TH1F *adcz[ntrace];
		TH1F *adc0[ntrace];
		for(int it= 0;it<ntrace;it++)
		{
		
			float T0 = 0;
			float Tpre = -250;
			float MaxS=-10;
			int Index=-10;
			
            int npoints = Detectors_trace_V1.at(it).size();
            double V0[npoints]={0};
            double VX[npoints]={0};
            double VY[npoints]={0};
            double VZ[npoints]={0};
            float FFTX[npoints];//={0};
            float FFTY[npoints];//={0};
            float FFTZ[npoints];//={0};
			float tFFt[513]={0};
			double V3D[npoints]={0};
            float EX[npoints]={0};
            float EY[npoints]={0};
            float EZ[npoints]={0};
            double Tt[npoints]={0};            

			double Maxsignal_X = 0;
			double Maxsignal_Y = 0;
			double Maxsignal_Z = 0;
			int Signal_T_X = -999;
			int Signal_T_Y = -999;
			int Signal_T_Z = -999;

			noiseX[it] = new TH1D(Form("NoiseX_%dDu_%d",i,it),"_XNoise;ADC;",200,-600,600);
			noiseY[it] = new TH1D(Form("NoiseY_%dDu_%d",i,it),"_YNoise;ADC;",200,-600,600);
			noiseZ[it] = new TH1D(Form("NoiseZ_%dDu_%d",i,it),"_ZNoise;ADC;",200,-600,600);
			noiseA[it] = new TH1D(Form("NoiseADC_%dDu_%d",i,it),"Pure_ADCNoise;ADC;",200,-600,600);
			noiseX[it]->AddDirectory(kFALSE);
			noiseY[it]->AddDirectory(kFALSE);
			noiseZ[it]->AddDirectory(kFALSE);
			noiseA[it]->AddDirectory(kFALSE);
			noiseX_2[it] = new TH1D(Form("NoiseX2_%dDu_%d",i,it),"_XNoise;ADC;",200,-600,600);
			noiseY_2[it] = new TH1D(Form("NoiseY2_%dDu_%d",i,it),"_YNoise;ADC;",200,-600,600);
			noiseZ_2[it] = new TH1D(Form("NoiseZ2_%dDu_%d",i,it),"_ZNoise;ADC;",200,-600,600);
			noiseA_2[it] = new TH1D(Form("NoiseADC2_%dDu_%d",i,it),"Pure_ADCNoise;ADC;",200,-600,600);
			noiseX_2[it]->AddDirectory(kFALSE);
			noiseY_2[it]->AddDirectory(kFALSE);
			noiseZ_2[it]->AddDirectory(kFALSE);
			noiseA_2[it]->AddDirectory(kFALSE);
		
			adcx[it] = new TH1F(Form("ChannelX_%d",it),"_X;ADC;",npoints+1,0,npoints);
			adcy[it] = new TH1F(Form("ChannelY_%d",it),"_Y;ADC;",npoints+1,0,npoints);
			adcz[it] = new TH1F(Form("ChannelZ_%d",it),"_Z;ADC;",npoints+1,0,npoints);
			adc0[it] = new TH1F(Form("ChannelF_%d",it),"Float;ADC;",npoints+1,0,npoints);
			if(npoints < 512 ) continue;
			//loop the trace for time
			for (int j=0; j<Detectors_trace_V1.at(it).size(); j++)
			{
               	float X=Detectors_trace_V1.at(it).at(j);
               	float Y=Detectors_trace_V2.at(it).at(j);
               	float Z=Detectors_trace_V3.at(it).at(j);
               	float F=Detectors_trace_V0.at(it).at(j);
				//temporary option for very large ADC 2**13
				if(TMath::Abs(X) > 8192) 
					X = 0.0;
				if(TMath::Abs(Y) > 8192) 
					Y = 0.0;
				if(TMath::Abs(Z) > 8192) 
					Z = 0.0;
				if(TMath::Abs(F) > 8192) 
					F = 0.0;
//				if( sqrt(X*X + Y*Y+ Z*Z) > PeakV)
//					PeakV = sqrt(X*X + Y*Y+ Z*Z) ;
               	Tt[j]=T0 + Tpre + 2.0*j;
          		VX[j]=X;
          		VY[j]=Y;
          		VZ[j]=Z;
          		V0[j]=F;
               	TVector3 V(X,Y,Z);
          		V3D[j]=V.Mag();
				if(abs(X) > Maxsignal_X )
				{
					Maxsignal_X = abs(X);
					Signal_T_X = j;
				}
				if(abs(Y) > Maxsignal_Y )
				{
					Maxsignal_Y = abs(Y);
					Signal_T_Y = j;
				}
				if(abs(Z) > Maxsignal_Z )
				{
					Maxsignal_Z = abs(Z);
					Signal_T_Z = j;
				}
               	if(abs(V.Mag()) > MaxS)
               	{			
               		Index = j;
               		MaxS = abs (V.Mag() );
               	}
				noiseX[it]->Fill(X);
				noiseY[it]->Fill(Y);
				noiseZ[it]->Fill(Z);
				noiseA[it]->Fill(F);
				noiseA_2[it]->Fill(F);
				adcx[it]->SetBinContent(j+1,X);
				adcy[it]->SetBinContent(j+1,Y);
				adcz[it]->SetBinContent(j+1,Z);
				adc0[it]->SetBinContent(j+1,F);

			}
		    Histo->Fill(Maxsignal_X) ;
			

			//ChannelX
			if(Signal_T_X - 512 > 0 )
			{
				for(int s=0;s< Signal_T_X - 150;s++)
				{
					noiseX_2[it]->Fill(Detectors_trace_V1.at(it).at(s));
				}
			}
			else
			{
				for(int s=Signal_T_X+150;s<npoints;s++)
					noiseX_2[it]->Fill(Detectors_trace_V1.at(it).at(s));
			}
			//ChannelY

			if(Signal_T_Y - 512 > 0 )
			{
				for(int s=0;s< Signal_T_Y - 150;s++)
				{
					noiseY_2[it]->Fill(Detectors_trace_V2.at(it).at(s));
				}
			}
			else
			{
				for(int s=Signal_T_Y+150;s<npoints;s++)
					noiseY_2[it]->Fill(Detectors_trace_V2.at(it).at(s));
			}

			if(abs(Signal_T_Y - Signal_T_X) > 50 ) continue;
			//ChannelZ

			if(Signal_T_Z - 512 > 0 )
			{
				for(int s=0;s< Signal_T_Z - 150;s++)
				{
					noiseZ_2[it]->Fill(Detectors_trace_V3.at(it).at(s));
				}
			}
			else
			{
				for(int s=Signal_T_Z+150;s<npoints;s++)
					noiseZ_2[it]->Fill(Detectors_trace_V3.at(it).at(s));
			}

			double SigmaX = noiseX[it]->GetRMS();
			double SigmaY = noiseY[it]->GetRMS();
			double SigmaZ = noiseZ[it]->GetRMS();
			double SigmaA = noiseA[it]->GetRMS();
			double Sigma_2X = noiseX_2[it]->GetRMS();
			double Sigma_2Y = noiseY_2[it]->GetRMS();
			double Sigma_2Z = noiseZ_2[it]->GetRMS();
			double Sigma_2A = noiseA_2[it]->GetRMS();
			double X_max = TMath::Max(SigmaX,Sigma_2X);
			double X_min = TMath::Min(SigmaX,Sigma_2X);
			double Y_max = TMath::Max(SigmaY,Sigma_2Y);
			double Y_min = TMath::Min(SigmaY,Sigma_2Y);
			double Z_max = TMath::Max(SigmaZ,Sigma_2Z);
			double Z_min = TMath::Min(SigmaZ,Sigma_2Z);
			//1678766540 is to Mar 14 4:00 UTC
			if(gps_time[2] < 1678766540 - 10*86400 || gps_time[2] > 1678766540+5*86400) continue;
			//if(SigmaX > 30 && SigmaY > 30  ) 
			if(Maxsignal_X/Sigma_2X > 5.0 && Maxsignal_Y/Sigma_2Y > 5.0  ) 
			// T   X  Y   Z  A  gps_temp accx accy accz amb_temp press humidity battery
			// 0   1  2   3  4  5         5    7    8     9        10    11      12
			{
				//cout<<"Signal-To-Noise "<<Maxsignal_X/Sigma_2X<<"\t"<<Maxsignal_X<<"\t"<<Sigma_2X<<endl;
			    //cout<<gps_time[2]<<"\t"<<SigmaX<<"\t"<<SigmaY<<"\t"<<SigmaZ<<"\t"<<SigmaA<<"\t"<<gps_temp[0]<<"\t";
			    //cout<<atm_humidity[0]*0.61035156/460.375 - 1650/460.375<<"\t"<<atm_pressure[0]*0.61035156/460.375 - 1650/460.375<<"\t"<<atm_temperature[0]*0.61035156/460.375 - 1650/460.375<<"\t";
			    //cout<<acceleration_x[0]*0.61035156/19.5 - 400/19.5<<"\t"<<acceleration_y[0]*0.61035156/3300./0.009+0.095/0.009<<"\t"<<acceleration_z[0]*0.61035156*0.048 - 500*0.048<<"\t"<<battery_level[0]*0.61035156*0.109/18<<endl;
			}
			else 
				continue;
		
			//*****************************fit&*****************************
			//cout<<"SNR "<< Maxsignal_X/SigmaX << "\t"<< Maxsignal_Y/SigmaY <<"\t"<< Maxsignal_Z/SigmaZ <<endl;
//			if(i==0 && it==0)noiseX[it]->Write();
//			if(i==0 && it==0)noiseY[it]->Write();
//			if(i==0 && it==0)noiseZ[it]->Write();
			noiseX[it]->Delete();
			noiseY[it]->Delete();
			noiseZ[it]->Delete();
			noiseA[it]->Delete();
			noiseX_2[it]->Delete();
			noiseY_2[it]->Delete();
			noiseZ_2[it]->Delete();
			noiseA_2[it]->Delete();


		}
	}
}
