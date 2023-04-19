/*************************************************************************
    > File Name: _Readfile.cpp
    > Author: Ma PengXiong
    > Mail: mapx@pmo.ac.cnn 
    > Created Time: Thu Jul  8 21:26:40 2021
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
#include "TH2F.h"
#include "TF1.h"
#include "TVector3.h"
#include "TStyle.h"
#include "TTimeStamp.h"
#include <vector>
#include "TVirtualFFT.h"
#include <fftw3.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex>
using namespace std;
#define PI 3.14159265358979323
void graph_var_markers(TGraph*, float*,float*,float*);
void _Readfile(const char* , const char* ,int ,int,const char* );
void _bandpass(double*, double*,double,double);
void _filter(double*, double*,double,double);
void FFT(std::complex<double> *, int , int );
void bandpassFilter(double *, int , double, double, int );
//void fft_real_f(int, const double*, float*);
int64_t getTimeStamp();
time_t StringToDatetime(const char*); 
vector<double> TP1; 
vector<double> VPos;
vector<float> t0;
vector<float> t1;
vector<float> detT0;
vector<float> position;
vector<float> position2;
vector<float> p2p;
vector<float> ppv;
vector<vector <float>> Vpos;
//TVector3 position_middle;
//TVector3 positionFirst;
float neff;
double par[2];
float sol = 3e8;
float d2r = TMath::DegToRad();
float r2d = TMath::RadToDeg();
void _bandpass(double *data, double *outdata, double f1, double f2)
{
	const int n = 1024;
	double fs = 500.0;
	double w1 = 2 * PI * f1 / fs;
	double w2 = 2 * PI * f2 / fs;
	double A = 1;//gain bandpass
	double b[3], a[3];
    b[0] = (w2 - w1) / PI;
    b[1] = 0;
    b[2] = -(w2 - w1) / PI;
    a[0] = 1 + A * (w2 - w1) / PI;
    a[1] = -2 * cos((w1 + w2) / 2) / a[0];
    a[2] = (1 - A * (w2 - w1) / PI) / a[0];
	// apply filter to input signal
    double y[n];
    y[0] = b[0] * data[0];
    y[1] = b[0] * data[1] + b[1] * data[0] - a[1] * y[0];
    for (int i = 2; i < n; i++) {
        y[i] = b[0] * data[i] + b[1] * data[i-1] + b[2] * data[i-2] - a[1] * y[i-1] - a[2] * y[i-2];
    }
	for (int i = 0; i < n; i++)
	{
		outdata[i] = y[i];
	}
	return ;
}
//  FFT 
void FFT(std::complex<double> *data, int length, int isign)
{
    int n, m, j, i;
    std::complex<double> w, wx, temp;
    double theta, wr, wi;

    n = length << 1;
    j = 1;
    for (i = 1; i < n; i += 2)
    {
        if (j > i)
        {
            std::swap(data[j - 1], data[i - 1]);
            std::swap(data[j], data[i]);
        }
        m = n >> 1;
        while (m >= 2 && j > m)
        {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    for (m = 2; m <= n; m <<= 1)
    {
        theta = isign * (2 * M_PI / m);
        wr = 1.0;
        wi = 0.0;
        wx = std::complex<double>(cos(theta), sin(theta));
        for (j = 1; j <= m >> 1; j++)
        {
            for (i = j; i <= n; i += m)
            {
                int k = i + (m >> 1);
                temp = data[k - 1] * std::complex<double>(wr, wi);
                data[k - 1] = data[i - 1] - temp;
                data[i - 1] += temp;
            }
            wr = (temp = wr * wx).real();
            wi = temp.imag();
        }
    }
}
void bandpassFilter(double *signal, int length, double f_lower, double f_upper, int sample_rate)
{
    std::complex<double> *data = new std::complex<double>[length];
    for (int i = 0; i < length; i++)
    {
        data[i] = std::complex<double>(signal[i], 0);
    }
    FFT(data, length, 1);

    double df = sample_rate / static_cast<double>(length);

    for (int i = 0; i < length / 2; i++)
    {
        double freq = i * df;
        if (freq > f_lower && freq < f_upper)
        {
            data[i] = std::complex<double>(0, 0);
            data[length - i - 1] = std::complex<double>(0, 0);
        }
    }

    FFT(data, length, -1);
    for (int i = 0; i < length; i++)
    {
        signal[i] = data[i].real() / static_cast<double>(length);
    }

    delete[] data;
}


void _filter(double *data, double *outdata, double f1, double f2)
{
	const int n = 1024;
	double fs = 500.0;
	double w1 = 2 * PI * f1 / fs;
	double w2 = 2 * PI * f2 / fs;
	double A = 1;//gain bandpass
	double b[3], a[3];
    b[0] = (w2 - w1) / PI;
    b[1] = 0;
    b[2] = -(w2 - w1) / PI;
    a[0] = 1 + A * (w2 - w1) / PI;
    a[1] = -2 * cos((w1 + w2) / 2) / a[0];
    a[2] = (1 - A * (w2 - w1) / PI) / a[0];
	// apply filter to input signal
    double y[n];
    y[0] = b[0] * data[0];
    y[1] = b[0] * data[1] + b[1] * data[0] - a[1] * y[0];
    for (int i = 2; i < n; i++) {
        y[i] = b[0] * data[i] + b[1] * data[i-1] + b[2] * data[i-2] - a[1] * y[i-1] - a[2] * y[i-2];
    }
	for (int i = 0; i < n; i++)
	{
		outdata[i] = y[i];
	}
	return ;
}

time_t StringToDatetime(const char* str)
{
    tm tm_;
    int year, month, day, hour, minute, second;
    sscanf(str, "%d/%d/%d %d:%d:%d", &day, &month, &year, &hour, &minute, &second);
    tm_.tm_year = year - 1900;
    tm_.tm_mon = month - 1;
    tm_.tm_mday = day;
    tm_.tm_hour = hour;
    tm_.tm_min = minute;
    tm_.tm_sec = second;
    tm_.tm_isdst = 0;

    time_t t_ = mktime(&tm_); 
    return t_; //second  
}




int main(int argc,const char* argv[]){
 
	time_t now;
    int unixTime = (int)time(&now);
    time_t tick = (time_t)unixTime;

    struct tm tm;
    char s[100];
    tm = *localtime(&tick);

    strftime(s, sizeof(s), "%Y-%m-%d %H:%M:%S", &tm);
    //printf("%d: %s\n", (int)unixTime, s);

//    return 0;

//	int64_t timeStamp = getTimeStamp();


//    std::chrono::time_point<std::chrono::system_clock, std::chrono::seconds> tpMicro = std::chrono::time_point_cast<std::chrono::seconds>(std::chrono::system_clock::now());
//    time_t timeStamp2 = tpMicro.time_since_epoch().count();


//    std::cout << timeStamp << std::endl;
//    std::cout << timeStamp2 << std::endl;
   
    std::string a = "01/4/2023 00:00:20";
    std::string b = "01/4/2023 00:00:01";
    time_t t1 = StringToDatetime(a.c_str());
    time_t t2 = StringToDatetime(b.c_str());

    //std::cout << t2<<"\t"<<t1<<"\t"<<t2-t1<<endl;


//	gSystem->Load("libMinuit");
	int plotif = atoi(argv[3]);
	int duid = atoi(argv[4]);
	_Readfile(argv[1],argv[2],plotif,duid,argv[5]);
	return 0;
}
void _Readfile(const char* input, const char* output, int pic, int du, const char* date)
{
	TCanvas *TC = new TCanvas("cc","tc",900,300);
	TC->SetFillColor(45);
	TC->SetTickx();
	TC->SetTicky();
	TString prefixname = "/home/mapx/mapx/Chao/ROOTFile/";
	TString prefixname1;
	TC->Print(Form(prefixname+"%s"+".pdf[",output ));
	TFile *f1 = new TFile(Form(prefixname1+"%s",input));
//	TTree *treSimu = (TTree*)f1->Get("teventshower");
	TTree *trevoltage = (TTree*)f1->Get("teventvoltage");
	TTree *treadc = (TTree*)f1->Get("teventadc");
//	cout<<"Entries = "<<treadc->GetEntries()<<endl;
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
	float FootprintSize;
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
	tree->Branch("FootprintSize",&FootprintSize,"FootprintSize/F");
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
	vector<unsigned short> gps_long;
	vector<unsigned short> gps_lat;
	vector<unsigned short> gps_alt;
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
	vector<unsigned short> *du_long = &gps_long;
	vector<unsigned short> *du_lat = &gps_lat;
	vector<unsigned short> *du_alt = &gps_alt;
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
	treadc->SetBranchAddress("gps_long",&du_long);
	treadc->SetBranchAddress("gps_lat",&du_lat);
	treadc->SetBranchAddress("gps_alt",&du_alt);
	treadc->SetBranchAddress("du_nanoseconds",&du_gps_nano);
	treadc->SetBranchAddress("atm_humidity",&at_hum);
	treadc->SetBranchAddress("atm_temperature",&at_temp);
	treadc->SetBranchAddress("acceleration_x",&du_ax);
	treadc->SetBranchAddress("acceleration_y",&du_ay);
	treadc->SetBranchAddress("acceleration_z",&du_az);
	treadc->SetBranchAddress("battery_level",&at_bat);
	treadc->SetBranchAddress("atm_pressure",&at_press);
	TH2F *Time_FreX;
	TH2F *Time_FreY;
	TH2F *Time_FreZ;
	TH1F *his_FX;
	TH1F *his_FY;
	TH1F *his_FZ;
	int Skip=1;//if you want much quickly go through all file. set bigger number.
	int tag=0;
	unsigned int start_time=0;
	int Ndu=2;//total working DUs in the root file.
	int tag_i = 0;
	int Stop_pos = Nentries ;//default last event.
	int Nstart=0;//default first event 
	for(int i=Nstart ;i< Stop_pos ;i+=Skip)
	{
		tag_i++;
	    p2p.clear();
	    ppv.clear();
	    Vpos.clear();
	    t0.clear();
	    t1.clear();
	    position.clear();
	    position2.clear();
	    
		vector<vector<short>> Detectors_trace_V1;
		vector<vector<short>> Detectors_trace_V2;
		vector<vector<short>> Detectors_trace_V3;
		vector<vector<short>> Detectors_trace_V0;
		vector<vector<short>> Detectors_t_0 ;
		if(tag_i==1) 
		{
			//(Stop_pos - Nstart)/Ndu*10, where 10 is to 10 seconds periodical DAQ mode, which is necessnary to convert the fromat of label to time format, 
			//Note if DAQ is on the other modes, please insure the correct conversion from gps timestamp to lable of timeformat. bin-to-bin timeBin-to-timeBin
			Time_FreX = new TH2F("Time_freq_X",";Time;",(Stop_pos - Nstart)/Ndu,0,(Stop_pos - Nstart)/Ndu*10,512,-0.5,250.5);
			Time_FreY = new TH2F("Time_freq_Y",";Time;",(Stop_pos - Nstart)/Ndu,0,(Stop_pos - Nstart)/Ndu*10,512,-0.5,250.5);
			Time_FreZ = new TH2F("Time_freq_Z",";Time;",(Stop_pos - Nstart)/Ndu,0,(Stop_pos - Nstart)/Ndu*10,512,-0.5,250.5);
		}
		if(tag_i==1) 
		{
			his_FX = new TH1F("_freq_X","Channel X;Freq.[MHz];Magnitude",512,-0.5,250.5);
			his_FY = new TH1F("_freq_Y","Channel Y;Freq.[MHz];Magnitude",512,-0.5,250.5);
			his_FZ = new TH1F("_freq_Z","Channel Z;Freq.[MHz];Magnitude",512,-0.5,250.5);
		}

//		treSimu->GetEntry(i);
		treadc->GetEntry(i);
		double xmaxDis[3];
		
		xmaxDis[0] = xmax_pos_shc[0];
		xmaxDis[1] = xmax_pos_shc[1];
		xmaxDis[2] = xmax_pos_shc[2];
		float distanceXmax = sqrt( pow(xmaxDis[0] - shower_core_pos[0],2.0) + pow(xmaxDis[1] - shower_core_pos[1],2.0) + pow(xmaxDis[2] - ground_alt,2.0) )/1e3;

		//cout<<"distanceXmax "<<distanceXmax<<endl;
		float XmaxAlt = xmax_alt/1e3;
		neff = 1+1E-6*325*TMath::Exp(-0.1218 * XmaxAlt);

		Detectors_trace_V1 =trace_x;
		Detectors_trace_V2 =trace_y;
		Detectors_trace_V3 =trace_z;
		Detectors_trace_V0 =trace_a;
		int ntrace = Detectors_trace_V1.size();
		if(ntrace > 13) continue;
		if(du_id[0] != du ) continue;
		TGraph *vx[ntrace];
		TGraph *v_fx[ntrace];
		TGraph *vy[ntrace];
		TGraph *vz[ntrace];
		TGraph *v0[ntrace];
		TH1F *adcx[ntrace];
		TH1F *adcy[ntrace];
		TH1F *adcz[ntrace];
		TH1F *adc0[ntrace];
		TH1D *histox[ntrace];
		TH1D *histoy[ntrace];
		TH1D *histoz[ntrace];
		TH1D *hbpx[ntrace];
		TH1D *hbpy[ntrace];
		TH1D *hbpz[ntrace];
		TGraph *ex[ntrace];
		TGraph *ey[ntrace];
		TGraph *ez[ntrace];
		TGraph *v3D[ntrace];
		vector<float> PosX;
		vector<float> PosY;
		vector<float> PosZ;
		double posX[ntrace]={-999999};
		double posY[ntrace]={-999999};
		double posZ[ntrace]={-999999};

		float omega[ntrace];
		float vpeak[ntrace];
		double ant_pos_evB[ntrace];
		double ant_pos_evvB[ntrace];
		double ant_pos_ev[ntrace];
		TVector3 ant_pos_proj[ntrace];

		float Magni[ntrace]={0};
		float Time[ntrace] ={0};
		float Type[ntrace] ={0};
		float Tmin=1e10;
		float Tmax= -1e10;

		int IndexP=0;
		float residual=0;
		float xmin=-999999;
		float ymin=-999999;
		float zmin=-999999;
		float xmax=-999999;
		float ymax=-999999;
		float zmax=-999999;
		
		P2P = 0;
		CPUTime = 0.0;
		PosX = pos_x;
		PosY = pos_y;
		PosZ = pos_z;
		Detectors_det_id = du_id;
		float covariance=0;
			
		float T0 ;
		double PeakV =0;
	    TH1D* noiseX[ntrace];// = new TH1D(Form("NoiseX_%d_%d",i,it),"_XNoise",200,-1000,1000);
	    TH1D* noiseY[ntrace];// = new TH1D(Form("NoiseY_%d_%d",i,it),"_YNoise",200,-1000,1000);
	    TH1D* noiseZ[ntrace];// = new TH1D(Form("NoiseZ_%d_%d",i,it),"_ZNoise",200,-1000,1000);
	    TH1D* noiseA[ntrace];// = new TH1D(Form("NoiseZ_%d_%d",i,it),"_ZNoise",200,-1000,1000);
	    TH1D* noiseX_2[ntrace];// = new TH1D(Form("NoiseX_%d_%d",i,it),"_XNoise",200,-1000,1000);
	    TH1D* noiseY_2[ntrace];// = new TH1D(Form("NoiseY_%d_%d",i,it),"_YNoise",200,-1000,1000);
	    TH1D* noiseZ_2[ntrace];// = new TH1D(Form("NoiseZ_%d_%d",i,it),"_ZNoise",200,-1000,1000);
	    TH1D* noiseA_2[ntrace];// = new TH1D(Form("NoiseZ_%d_%d",i,it),"_ZNoise",200,-1000,1000);
		time_t t_start ;
		if(tag_i<Ndu*Skip + Ndu*2 )
		{
			//20230417123123 [0,15]
			std::string specific_date;
			specific_date += date[6];
			specific_date += date[7];
			specific_date += "/";
			specific_date += date[4];
			specific_date += date[5];
			specific_date += "/";
			specific_date += date[0];
			specific_date += date[1];
			specific_date += date[2];
			specific_date += date[3];
			specific_date += " ";
			specific_date += date[8];
			specific_date += date[9];
			specific_date += ":";
			specific_date += date[10];
			specific_date += date[11];
			specific_date += ":";
			specific_date += date[12];
			specific_date += date[13];
			// 17/04/2023 12:31:23
			t_start = StringToDatetime(specific_date.c_str());

			//abandaon those event with significant strange wrong timestamp
			if(gps_time[2] < t_start || gps_time[2] > t_start + 10*3600) continue;
			if(tag > 0) continue;
			tag++;
			start_time = gps_time[2];		
			//UTC+8 Beijing Time Zone. To build the fast relations of temperature, solar strength and local time.
		    TTimeStamp aa = TTimeStamp(gps_time[2]+3600*8);
		    unsigned int Date = aa.GetDate();
		    unsigned int Moment = aa.GetTime();
		    unsigned int month = aa.GetMonth();
		    int year = Date/10000;//2023
		    int date = Date%year;//April 8 ;408
		    int day = date%100;// 8
		    int hour= Moment/10000;
		    int minute =0;
		    int seconds = 0;
			int day2hour=10000;
			int day2min=100;
		    if (hour != 0) {
		        minute=Moment%(hour*day2hour)/day2min;
				if(minute ==0 )
					seconds=Moment%(hour*day2hour);
				else 	
					seconds=Moment%(hour*day2hour)%(day2min*minute);
		    }
		    else
		    {
		        minute=Moment/day2min;
		        seconds = Moment% day2min;
		    }
		    TDatime da(year,month,day,hour,minute,seconds);
		    //printf("%d%02d%02d-%02d%02d%02d\n",year,month,day,hour,minute,seconds);
		    gStyle->SetTimeOffset(da.Convert());
		    Time_FreX->SetTitle(Form("ChannelX DU%d UTC+8h %d%02d%02d-%02d:%02d:%02d",du,year,month,day,hour,minute,seconds));
		    Time_FreX->SetTitleFont(20);
		    Time_FreX->GetXaxis()->SetTimeDisplay(1);
		    Time_FreX->GetXaxis()->SetLabelFont(20);
		    Time_FreX->GetXaxis()->SetLabelSize(0.05);
		    Time_FreX->GetYaxis()->SetLabelFont(20);
		    Time_FreX->GetYaxis()->SetLabelSize(0.05);
		    Time_FreX->GetYaxis()->SetLabelColor(kRed);
		    Time_FreX->GetXaxis()->SetLabelColor(kRed);
		    Time_FreX->GetYaxis()->SetTitleColor(kRed);
		    Time_FreX->GetXaxis()->SetTitleColor(kRed);
		    Time_FreX->GetYaxis()->SetTitleSize(0.05);
		    Time_FreX->GetXaxis()->SetTitleSize(0.04);
		    Time_FreX->GetXaxis()->SetTimeFormat("%m\/%d-%H:%M");
		    Time_FreY->SetTitle(Form("ChannelY DU%d UTC+8h %d%02d%02d-%02d:%02d:%02d",du,year,month,day,hour,minute,seconds));
		    Time_FreY->SetTitleFont(20);
		    Time_FreY->GetXaxis()->SetTimeDisplay(1);
		    Time_FreY->GetXaxis()->SetLabelFont(20);
		    Time_FreY->GetXaxis()->SetLabelSize(0.05);
		    Time_FreY->GetYaxis()->SetLabelFont(20);
		    Time_FreY->GetYaxis()->SetLabelSize(0.05);
		    Time_FreY->GetYaxis()->SetLabelColor(kRed);
		    Time_FreY->GetXaxis()->SetLabelColor(kRed);
		    Time_FreY->GetYaxis()->SetTitleColor(kRed);
		    Time_FreY->GetXaxis()->SetTitleColor(kRed);
		    Time_FreY->GetYaxis()->SetTitleSize(0.05);
		    Time_FreY->GetXaxis()->SetTitleSize(0.04);
		    Time_FreY->GetXaxis()->SetTimeFormat("%m\/%d-%H:%M");
		    Time_FreZ->SetTitle(Form("ChannelZ DU%d UTC+8h %d%02d%02d-%02d:%02d:%02d",du,year,month,day,hour,minute,seconds));
		    Time_FreZ->SetTitleFont(20);
		    Time_FreZ->GetXaxis()->SetTimeDisplay(1);
		    Time_FreZ->GetXaxis()->SetLabelFont(20);
		    Time_FreZ->GetXaxis()->SetLabelSize(0.05);
		    Time_FreZ->GetYaxis()->SetLabelFont(20);
		    Time_FreZ->GetYaxis()->SetLabelSize(0.05);
		    Time_FreZ->GetYaxis()->SetLabelColor(kRed);
		    Time_FreZ->GetXaxis()->SetLabelColor(kRed);
		    Time_FreZ->GetYaxis()->SetTitleColor(kRed);
		    Time_FreZ->GetXaxis()->SetTitleColor(kRed);
		    Time_FreZ->GetYaxis()->SetTitleSize(0.05);
		    Time_FreZ->GetXaxis()->SetTitleSize(0.04);
		    Time_FreZ->GetXaxis()->SetTimeFormat("%m\/%d-%H:%M");
			//his_FX->SetTitle(Form("ChannelX DU%d UTC+8h %d%02d%02d-%02d:%02d:%02d",du,year,month,day,hour,minute,seconds));
		    his_FX->SetTitleFont(20);
		    his_FX->GetXaxis()->SetLabelFont(20);
		    his_FX->GetXaxis()->SetLabelSize(0.05);
		    his_FX->GetYaxis()->SetLabelFont(20);
		    his_FX->GetYaxis()->SetLabelSize(0.05);
		    his_FX->GetYaxis()->SetTitleSize(0.05);
		    his_FX->GetXaxis()->SetTitleOffset(1.5);
		    his_FX->GetXaxis()->SetTitleSize(0.05);
			his_FX->GetXaxis()->SetTitleColor(kRed);
			his_FX->GetYaxis()->SetTitleColor(kRed);
			his_FX->GetXaxis()->SetLabelColor(kRed);
			his_FX->GetYaxis()->SetLabelColor(kRed);
			//his_FY->SetTitle(Form("ChannelY DU%d UTC+8h %d%02d%02d-%02d:%02d:%02d",du,year,month,day,hour,minute,seconds));
		    his_FY->SetTitleFont(20);
		    his_FY->GetXaxis()->SetLabelFont(20);
		    his_FY->GetXaxis()->SetLabelSize(0.05);
		    his_FY->GetYaxis()->SetLabelFont(20);
		    his_FY->GetYaxis()->SetLabelSize(0.05);
		    his_FY->GetYaxis()->SetTitleSize(0.05);
		    his_FY->GetXaxis()->SetTitleOffset(1.5);
		    his_FY->GetXaxis()->SetTitleSize(0.05);
			his_FY->GetXaxis()->SetTitleColor(kRed);
			his_FY->GetYaxis()->SetTitleColor(kRed);
			his_FY->GetXaxis()->SetLabelColor(kRed);
			his_FY->GetYaxis()->SetLabelColor(kRed);
			//his_FZ->SetTitle(Form("ChannelZ DU%d UTC+8h %d%02d%02d-%02d:%02d:%02d",du,year,month,day,hour,minute,seconds));
		    his_FZ->SetTitleFont(20);
		    his_FZ->GetXaxis()->SetLabelFont(20);
		    his_FZ->GetXaxis()->SetLabelSize(0.05);
		    his_FZ->GetYaxis()->SetLabelFont(20);
		    his_FZ->GetYaxis()->SetLabelSize(0.05);
		    his_FZ->GetYaxis()->SetTitleSize(0.05);
		    his_FZ->GetXaxis()->SetTitleOffset(1.5);
		    his_FZ->GetXaxis()->SetTitleSize(0.05);
			his_FZ->GetXaxis()->SetTitleColor(kRed);
			his_FZ->GetYaxis()->SetTitleColor(kRed);
			his_FZ->GetXaxis()->SetLabelColor(kRed);
			his_FZ->GetYaxis()->SetLabelColor(kRed);
		}
		for(int it= 0;it<ntrace;it++)
		{
			T0 = 0;
			float Tpre = 0;
			float MaxS=-10;
			int Index=-10;
			
            int npoints = Detectors_trace_V1.at(it).size();
			//currently, we set the length of trace to 1024. 
			if(npoints != 1024) continue;
		    if(Detectors_trace_V2.at(it).size() != 1024) continue;
		    if(Detectors_trace_V3.at(it).size() != 1024) continue;
		    if(Detectors_trace_V0.at(it).size() != 1024) continue;
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
			int W_X = 0;
			int W_Y = 0;
			int W_Z = 0;
			int W_F = 0;
			//loop the trace for time
			for (int j=0; j<Detectors_trace_V1.at(it).size(); j++)
			{
               	float X=Detectors_trace_V1.at(it).at(j);
               	float Y=Detectors_trace_V2.at(it).at(j);
               	float Z=Detectors_trace_V3.at(it).at(j);
               	float F=Detectors_trace_V0.at(it).at(j);
				//temporary option for very large ADC 2**13
				//if(TMath::Abs(X) > 8192) 
				//	X = 0.0;
				//if(TMath::Abs(Y) > 8192) 
				//	Y = 0.0;
				//if(TMath::Abs(Z) > 8192) 
				//	Z = 0.0;
				//if(TMath::Abs(F) > 8192) 
				//	F = 0.0;
				//	counting the number of wrong ADC!!!  
				if(TMath::Abs(X) > 8192) 
					W_X += 1 ;
				if(TMath::Abs(Y) > 8192) 
					W_Y += 1 ;
				if(TMath::Abs(Z) > 8192) 
					W_Z += 1 ;
				if(TMath::Abs(F) > 8192) 
					W_F += 1 ;
				ppv.push_back( sqrt(X*X + Y*Y+ Z*Z) );
				if( sqrt(X*X + Y*Y+ Z*Z) > PeakV)
					PeakV = sqrt(X*X + Y*Y+ Z*Z) ;
               	Tt[j]=T0 + Tpre + 1.0*j;
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

			// coincidence in time domain
			//if(abs(Signal_T_Y - Signal_T_X) > 50 ) continue;
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

			//ChannelF temporarally not used here

			//if(Signal_T_X - 512 > 0 )
			//{
			//	for(int s=0;s< Signal_T_X - 150;s++)
			//	{
			//		noiseX_2[it]->Fill(Detectors_trace_V1.at(it).at(s));
			//	}
			//}
			//else
			//{
			//	for(int s=Signal_T_X+150;s<npoints;s++)
			//		noiseX_2[it]->Fill(Detectors_trace_V1.at(it).at(s));
			//}
            /*
				for(int s=0;s<513;s++)
				{
					tFFt[s] = 0.5*s;
					noiseA_2[it]->Fill(Detectors_trace_V0.at(it).at(s));
			//		cout<<Detectors_trace_V0.at(it).at(s)<<"\t";
				}
			//	cout<<endl;
			//	*/
            noiseX[it]->Fit("gaus","Q0");
            noiseY[it]->Fit("gaus","Q0");
            noiseZ[it]->Fit("gaus","Q0");
            noiseA[it]->Fit("gaus","Q0");
            noiseX_2[it]->Fit("gaus","Q0");
            noiseY_2[it]->Fit("gaus","Q0");
            noiseZ_2[it]->Fit("gaus","Q0");
            noiseA_2[it]->Fit("gaus","Q0");
            TF1 *funcx = (TF1*)noiseX[it]->GetFunction("gaus");
            TF1 *funcy = (TF1*)noiseY[it]->GetFunction("gaus");
            TF1 *funcz = (TF1*)noiseZ[it]->GetFunction("gaus");
            TF1 *funca = (TF1*)noiseA[it]->GetFunction("gaus");
            TF1 *func2_x = (TF1*)noiseX_2[it]->GetFunction("gaus");
            TF1 *func2_y = (TF1*)noiseY_2[it]->GetFunction("gaus");
            TF1 *func2_z = (TF1*)noiseZ_2[it]->GetFunction("gaus");
            TF1 *func2_a = (TF1*)noiseA_2[it]->GetFunction("gaus");

			double SigmaX = funcx->GetParameter(2);
			double SigmaY = funcy->GetParameter(2);
			double SigmaZ = funcz->GetParameter(2);
			double SigmaA = funca->GetParameter(2);
			double MeanX = funcx->GetParameter(1);
			double MeanY = funcy->GetParameter(1);
			double MeanZ = funcz->GetParameter(1);
			double MeanA = funca->GetParameter(1);
			//double SigmaX = noiseX[it]->GetRMS(); 
			//double SigmaY = noiseY[it]->GetRMS(); 
			//double SigmaZ = noiseZ[it]->GetRMS(); 
			//double SigmaA = noiseA[it]->GetRMS(); 
			double Sigma_2X = func2_x->GetParameter(2);
			double Sigma_2Y = func2_y->GetParameter(2);
			double Sigma_2Z = func2_z->GetParameter(2);
			double Sigma_2A = func2_a->GetParameter(2);
			//double Sigma_2X = noiseX_2[it]->GetRMS(); 
			//double Sigma_2Y = noiseY_2[it]->GetRMS(); 
			//double Sigma_2Z = noiseZ_2[it]->GetRMS(); 
			//double Sigma_2A = noiseA_2[it]->GetRMS(); 
			double X_max = TMath::Max(SigmaX,Sigma_2X);
			double X_min = TMath::Min(SigmaX,Sigma_2X);
			double Y_max = TMath::Max(SigmaY,Sigma_2Y);
			double Y_min = TMath::Min(SigmaY,Sigma_2Y);
			double Z_max = TMath::Max(SigmaZ,Sigma_2Z);
			double Z_min = TMath::Min(SigmaZ,Sigma_2Z);
			//if(TMath::Abs( SigmaA) < 30  ) continue;
			//cout<<gps_time[2]<<"\t"<<start_time<<endl;
			if( gps_time[2] - start_time  > 3600*9 ) continue; 
			if( gps_time[2] - start_time  < 0 ) continue;
			//if(SigmaX < 300 && SigmaY < 300  ) 
			////if(Maxsignal_X/Sigma_2X > 5.0 && Maxsignal_Y/Sigma_2Y > 5.0  ) 
			//// T   X  Y   Z  A  gps_temp accx accy accz amb_temp press humidity battery
			//// 0   1  2   3  4  5         5    7    8     9        10    11      12
			//{
			//  //cout<<"Signal-To-Noise "<<Maxsignal_X/Sigma_2X<<"\t"<<Maxsignal_X<<"\t"<<Sigma_2X<<endl;
			//  cout<<gps_time[2]<<"\t"<<SigmaX<<"\t"<<SigmaY<<"\t"<<SigmaZ<<"\t"<<SigmaA<<"\t"<<gps_temp[0]<<"\t";
			//  cout<<atm_humidity[0]*0.61035156/460.375 - 1650/460.375<<"\t"<<atm_pressure[0]*0.61035156/460.375 - 1650/460.375<<"\t"<<atm_temperature[0]*0.61035156/460.375 - 1650/460.375<<"\t";
			//  cout<<acceleration_x[0]*0.61035156/19.5 - 400/19.5<<"\t"<<acceleration_y[0]*0.61035156/3300./0.009+0.095/0.009<<"\t"<<acceleration_z[0]*0.61035156*0.048 - 500*0.048<<"\t"<<battery_level[0]*0.61035156*0.109/18<<endl;
			//}
			//else 
			//	continue;
		
			//*****************************fit&*****************************
			//cout<<"SNR "<< Maxsignal_X/SigmaX << "\t"<< Maxsignal_Y/SigmaY <<"\t"<< Maxsignal_Z/SigmaZ <<endl;
			if(i==0 && it==0)noiseX[it]->Write();
			if(i==0 && it==0)noiseY[it]->Write();
			if(i==0 && it==0)noiseZ[it]->Write();
			noiseX[it]->Delete();
			noiseY[it]->Delete();
			noiseZ[it]->Delete();
			noiseA[it]->Delete();
			noiseX_2[it]->Delete();
			noiseY_2[it]->Delete();
			noiseZ_2[it]->Delete();
			noiseA_2[it]->Delete();
			if(pic == 1)
			{
				//**************FFT*****************
			    /*fft_real_f(npoints,VX, FFTX);	
			    fft_real_f(npoints,VY, FFTY);	
			    fft_real_f(npoints,VZ, FFTZ);	
				*/
				//**************FFT*****************
	            TPad *pad1 = new TPad("fpad1","my pad1",0.0, 0.51,0.32,0.99);
	            TPad *pad2 = new TPad("fpad2","my pad2",0.34,0.51,0.66,0.99);
	            TPad *pad3 = new TPad("fpad3","my pad3",0.68,0.51,0.99,0.99);
	            TPad *pad4 = new TPad("fpad4","my pad4",0.0, 0.0,0.32,0.5);
	            TPad *pad5 = new TPad("fpad5","my pad5",0.34,0.0,0.66,0.5);
	            TPad *pad6 = new TPad("fpad6","my pad6",0.68,0.0,0.99,0.5);
	            pad1->Draw();
	            pad2->Draw();
	            pad3->Draw();
	            pad4->Draw();
	            pad5->Draw();
	            pad6->Draw();
				gStyle->SetOptStat(0);
				gStyle->SetOptFit(0);
		       
		        TTimeStamp t_tmp = TTimeStamp(gps_time[2]+3600*8);//UTC+8 Beijing Time Zone
		        unsigned int g_Moment = t_tmp.GetTime();
				
				vx[it] = new TGraph(npoints,Tt,VX);
             	//vx[it]->SetMarkerColor(kRed-3);
             	//vx[it]->SetLineColor(kRed-3);
             	vx[it]->SetMarkerSize(0.2);
             	vx[it]->SetMarkerStyle(kFullCircle);
             	//vx[it]->SetName(Form("_Event%d_Antenna%d_x",i,du_id[it]));
             	vx[it]->SetTitle(Form("_Evt%d_Ant%d_X_%06i;Time [ns];ADC ",i,du_id[it],g_Moment));
				vx[it]->GetXaxis()->SetLabelSize(0.05);
				vx[it]->GetXaxis()->SetTitleSize(0.06);
				vx[it]->GetYaxis()->SetTitleSize(0.07);
				vx[it]->GetYaxis()->SetLabelSize(0.05);
				vx[it]->GetYaxis()->SetTitleOffset(1);
             	//vx[it]->Write();
    //			vx[it]->Draw("AP");
             	vy[it] = new TGraph(npoints,Tt,VY);
             	//vy[it]->SetMarkerColor(kSpring-7);
             	//vy[it]->SetLineColor(kSpring-7);
             	vy[it]->SetMarkerSize(0.2);
             	vy[it]->SetMarkerStyle(kFullCircle);
             	vy[it]->SetName(Form("_Event%d_Antenna%d_y",i,du_id[it]));
             	//vy[it]->SetTitle(Form("_Event%d_Antenna%d_y;Time [ns];ADC ",i,du_id[it]));
             	vy[it]->SetTitle(Form("_Evt%d_Ant%d_Y_%06d;Time [ns];ADC ",i,du_id[it],g_Moment));
				vy[it]->GetXaxis()->SetLabelSize(0.05);
				vy[it]->GetXaxis()->SetTitleSize(0.06);
				vy[it]->GetYaxis()->SetLabelSize(0.05);
    			//vy[it]->Write();
				//         	vy[it]->Draw("Psame");
             	vz[it] = new TGraph(npoints,Tt,VZ);
             	//vz[it]->SetMarkerColor(kGreen);
             	//vz[it]->SetLineColor(kGreen);
             	vz[it]->SetMarkerSize(0.2);
             	vz[it]->SetMarkerStyle(kFullCircle);
             	vz[it]->SetName(Form("_Event%d_Antenna%d_z",i,du_id[it]));
             	//vz[it]->SetTitle(Form("_Event%d_Antenna%d_z;Time [ns];ADC ",i,du_id[it]));
             	vz[it]->SetTitle(Form("_Evt%d_Ant%d_Z_%06d;Time [ns];ADC ",i,du_id[it],g_Moment));
				vz[it]->GetXaxis()->SetLabelSize(0.05);
				vz[it]->GetXaxis()->SetTitleSize(0.06);
				vz[it]->GetYaxis()->SetLabelSize(0.05);
    			//vz[it]->Write();
    			v0[it] = new TGraph(npoints,Tt,V0);
             	v0[it]->SetMarkerColor(kCyan-3);
             	v0[it]->SetMarkerSize(0.2);
             	v0[it]->SetLineColor(kCyan-3);
             	v0[it]->SetMarkerStyle(kFullCircle);
             	v0[it]->SetName(Form("_Event%d_Antenna%d_f",i,du_id[it]));
             	v0[it]->SetTitle(Form("_Event%d_Antenna%d_f;Time [ns];ADC ",i,du_id[it]));
             	v3D[it] = new TGraph(npoints,Tt,V3D);
             	v3D[it]->SetMarkerColor(kMagenta);
             	v3D[it]->SetMarkerSize(0.2);
             	v3D[it]->SetLineColor(kMagenta);
             	v3D[it]->SetMarkerStyle(kFullCircle);
             	v3D[it]->SetName(Form("_Event%d_Antenna%d_3D",i,du_id[it]));
             	v3D[it]->SetTitle(Form("_Event%d_Antenna%d_3D;Time [ns];ADC ",i,du_id[it]));
             	//v3D[it]->Write();
				TVirtualFFT::SetTransform(0);
				TH1* fftx=0;
				TH1D fft_x("fx",";Freq[MHz];Magnitude",npoints/2,-0.5,250.5);		
				fftx = adcx[it]->FFT(fftx, "MAG");
                fftx->SetTitle("X Magnitude of the 1st transform");
				fft_x.GetXaxis()->SetLabelSize(0.05);
				fft_x.GetYaxis()->SetLabelSize(0.07);
				fft_x.GetYaxis()->SetTitleSize(0.07);
				fft_x.GetYaxis()->SetTitleOffset(1);
				fft_x.SetLineWidth(1);
				fft_x.SetMarkerSize(0.2);
				for(int x=0;x<npoints/2;x++)
				{
					fft_x.SetBinContent(x+1,fftx->GetBinContent(x+1)/sqrt(npoints));
					Time_FreX->SetBinContent( (tag_i)/Ndu,x+1,fftx->GetBinContent(x+1)/sqrt(npoints));
				}
				his_FX->Add(&fft_x);
				cout<<i<<"\tMaxBin "<<fft_x.GetMaximumBin()<<endl;
				//if(fftx->GetBinContent(245)/fftx->GetBinContent(250) > 5) 
				//{
				//	fftx->SetBinContent(245,fftx->GetBinContent(250));
				//}
				
				//-----------------------
				/*TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
				Double_t *re_full = new Double_t[npoints];
				Double_t *im_full = new Double_t[npoints];
				fft->GetPointsComplex(re_full,im_full);
				TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &npoints, "C2R M K");
				fft_back->SetPointsComplex(re_full,im_full);
				fft_back->Transform();
				TH1 *hb = 0;
				hb = TH1::TransformHisto(fft_back,hb,"Re");
				hb->SetTitle("The backward transform result");
				hb->SetLineColor(kMagenta-2);
				hb->SetMarkerColor(kMagenta-2);
				*/
				//-----------------------
				//hb->SetTitle("The bandpass filter result");
				//hb->SetLineColor(kMagenta-2);
				//hb->SetMarkerColor(kMagenta-2);
				

				TH1* ffty=0;
				TH1D fft_y("fy",";Freq[MHz];",npoints/2,-0.5,250.5);		
				ffty = adcy[it]->FFT(ffty, "MAG");
                ffty->SetTitle("Y Magnitude of the 1st transform");
				fft_y.GetXaxis()->SetLabelSize(0.05);
				fft_y.GetYaxis()->SetLabelSize(0.05);
				fft_y.SetLineWidth(1);
				fft_y.SetMarkerSize(0.2);
				for(int x=0;x<npoints/2+1;x++)
				{
					fft_y.SetBinContent(x+1,ffty->GetBinContent(x+1)/sqrt(npoints));
					Time_FreY->SetBinContent( (tag_i )/Ndu,x+1,ffty->GetBinContent(x+1)/sqrt(npoints));
				}
				his_FY->Add(&fft_y);

				TH1* fftz=0;
				TH1D fft_z("fz",";Freq[MHz];",npoints/2,-0.5,250.5);		
				fftz = adcz[it]->FFT(fftz, "MAG");
                fftz->SetTitle("Z Magnitude of the 1st transform");
				fft_z.GetXaxis()->SetLabelSize(0.05);
				fft_z.GetYaxis()->SetLabelSize(0.05);
				fft_z.SetLineWidth(1);
				fft_z.SetMarkerSize(0.2);
				for(int x=0;x<npoints/2+1;x++)
				{
					fft_z.SetBinContent(x+1,fftz->GetBinContent(x+1)/sqrt(npoints));
					Time_FreZ->SetBinContent( (tag_i )/Ndu,x+1,fftz->GetBinContent(x+1)/sqrt(npoints));
				}
				his_FZ->Add(&fft_z);
				
				//----------------------------------------------------------------------------------------
				hbpx[it] = new TH1D(Form("_filter_ADCX_Du%d_%d",du_id[it],i),Form("_ADCX_%dDu_Eve%d;time;ADC",du_id[it],i),npoints+1,0,npoints);
				hbpy[it] = new TH1D(Form("_filter_ADCY_Du%d_%d",du_id[it],i),Form("_ADCY_%dDu_Eve%d;time;ADC",du_id[it],i),npoints+1,0,npoints);
				hbpz[it] = new TH1D(Form("_filter_ADCZ_Du%d_%d",du_id[it],i),Form("_ADCZ_%dDu_Eve%d;time;ADC",du_id[it],i),npoints+1,0,npoints);
				histox[it] = new TH1D(Form("NoiseX_Du%d_%d",du_id[it],i),";ADCX;",200,-600,600);
				histoy[it] = new TH1D(Form("NoiseY_Du%d_%d",du_id[it],i),";ADCY;",200,-600,600);
				histoz[it] = new TH1D(Form("NoiseZ_Du%d_%d",du_id[it],i),";ADCZ;",200,-600,600);
			    histox[it]->AddDirectory(kFALSE);
			    histoy[it]->AddDirectory(kFALSE);
			    histoz[it]->AddDirectory(kFALSE);
			    hbpx[it]->AddDirectory(kFALSE);
			    hbpy[it]->AddDirectory(kFALSE);
			    hbpz[it]->AddDirectory(kFALSE);

				double F_X1[npoints];
				double F_X2[npoints];
				double F_X[npoints]={0};
				double F_Y[npoints]={0};
				double F_Z[npoints]={0};
				
				if( (fft_x.GetMaximumBin()*0.48828125 < 45 && fft_x.GetMaximumBin()*0.48828125 > 25 ) || ( fft_y.GetMaximumBin()*0.48828125 < 45 && fft_y.GetMaximumBin()*0.48828125 > 25 ) 
				||  (fft_x.GetMaximumBin()*0.48828125 < 145 && fft_x.GetMaximumBin()*0.48828125 > 115 ) || ( fft_y.GetMaximumBin()*0.48828125 < 145 && fft_y.GetMaximumBin()*0.48828125 > 115 ) )				
				{      
						_bandpass(VX,F_X, (fft_x.GetMaximumBin()-25)*0.48828125, (fft_x.GetMaximumBin()+25)*0.48828125);
						_bandpass(VY,F_Y, (fft_y.GetMaximumBin()-25)*0.48828125, (fft_y.GetMaximumBin()+25)*0.48828125);
						_bandpass(VZ,F_Z, (fft_x.GetMaximumBin()-25)*0.48828125, (fft_x.GetMaximumBin()+25)*0.48828125);
				    for(int x=0;x<npoints;x++)
				    {
				    	F_X[x] = VX[x] - F_X[x];
				    	F_Y[x] = VY[x] - F_Y[x];
				    	F_Z[x] = VZ[x] - F_Z[x];
				    	//F_X[x] = F_X1[x]+F_X2[x];
				    	hbpx[it]->SetBinContent(x+1,F_X[x]);
				    	hbpy[it]->SetBinContent(x+1,F_Y[x]);
				    	hbpz[it]->SetBinContent(x+1,F_Z[x]);
						//histox[it]->Fill(F_X[x]);
						//histoy[it]->Fill(F_Y[x]);
						//histoz[it]->Fill(F_Z[x]);
				    //	hbpx[it]->SetBinContent(x+1,F_X1[x]+F_X2[x]);
				    }
				}
				else 
				{ 
				    for(int x=0;x<npoints;x++)
				    {
				    	F_X[x] = VX[x];
				    	F_Y[x] = VY[x];
				    	F_Z[x] = VZ[x];
				    	hbpx[it]->SetBinContent(x+1,F_X[x]);
				    	hbpy[it]->SetBinContent(x+1,F_Y[x]);
				    	hbpz[it]->SetBinContent(x+1,F_Z[x]);
						//histox[it]->Fill(F_X[x]);
						//histoy[it]->Fill(F_Y[x]);
						//histoz[it]->Fill(F_Z[x]);
				    }
				}
				//-----------------------New Bankground nose ----------------------------

			    if(Signal_T_X - 512 > 0 )
			    {
			    	for(int s=0;s< Signal_T_X - 150;s++)
			    	{
			    		histox[it]->Fill(F_X[s]);
			    	}
			    }
			    else
			    {
			    	for(int s=Signal_T_X+150;s<npoints;s++)
			    		histox[it]->Fill(F_X[s]);
			    }
			    //ChannelY

			    if(Signal_T_Y - 512 > 0 )
			    {
			    	for(int s=0;s< Signal_T_Y - 150;s++)
			    	{
			    		histoy[it]->Fill(F_Y[s]);
			    	}
			    }
			    else
			    {
			    	for(int s=Signal_T_Y+150;s<npoints;s++)
			    		histoy[it]->Fill(F_Y[s]);
			    }

			    //if(abs(Signal_T_Y - Signal_T_X) > 50 ) continue;
			    //ChannelZ

			    if(Signal_T_Z - 512 > 0 )
			    {
			    	for(int s=0;s< Signal_T_Z - 150;s++)
			    	{
						histoz[it]->Fill(F_Z[s]);
			    	}
			    }
			    else
			    {
			    	for(int s=Signal_T_Z+150;s<npoints;s++)
						histoz[it]->Fill(F_Z[s]);
			    }
                histox[it]->Fit("gaus","Q0","");
                histoy[it]->Fit("gaus","Q0","");
                histoz[it]->Fit("gaus","Q0","");
				//Test!!!!!!!!!//GausDistri
                TF1 *F_hisx = (TF1*)histox[it]->GetFunction("gaus");
                TF1 *F_hisy = (TF1*)histoy[it]->GetFunction("gaus");
                TF1 *F_hisz = (TF1*)histoz[it]->GetFunction("gaus");
                double NewSigmaX = F_hisx->GetParameter(2);
                double NewSigmaY = F_hisy->GetParameter(2);
                double NewSigmaZ = F_hisz->GetParameter(2);
                double NewMeanX = F_hisx->GetParameter(1);
                double NewMeanY = F_hisy->GetParameter(1);
                double NewMeanZ = F_hisz->GetParameter(1);
				//5 times Signal-To-Noise ratios
				//if(F_X[Signal_T_X]/NewSigmaX < 5.0 || F_Y[Signal_T_Y]/NewSigmaY < 5.0) continue;
				//if (NewSigmaX < 25) continue;
				pad1->cd();
		        //pad1->SetRightMargin(0.05);
		        //pad1->SetGrid();
		        pad1->SetTickx();
		        pad1->SetTicky();
             	hbpx[it]->SetMarkerColor(kRed-3);
             	hbpx[it]->SetMarkerSize(0.1);
             	hbpx[it]->SetLineColor(kRed-3);
             	hbpx[it]->SetMarkerStyle(kFullCircle);
				hbpx[it]->GetXaxis()->SetLabelSize(0.05);
				hbpx[it]->GetYaxis()->SetLabelSize(0.05);
			    vx[it]->Draw();
			    hbpx[it]->Draw("same");
		        pad2->cd();
		        //pad2->SetLeftMargin(0.15);
		        //pad2->SetRightMargin(0.02);
		        //pad2->SetGrid();
		        pad2->SetTickx();
		        pad2->SetTicky();
             	hbpy[it]->SetMarkerColor(kCyan);
             	hbpy[it]->SetMarkerSize(0.1);
             	hbpy[it]->SetLineColor(kCyan);
             	hbpy[it]->SetMarkerStyle(kFullCircle);
				hbpy[it]->GetXaxis()->SetLabelSize(0.05);
				hbpy[it]->GetYaxis()->SetLabelSize(0.05);
			    vy[it]->Draw();
			    hbpy[it]->Draw("same");
		        pad3->cd();
		        //pad3->SetLeftMargin(0.15);
		        //pad3->SetRightMargin(0.02);
		        //pad3->SetRightMargin(0.05);
		        //pad3->SetGrid();
		        pad3->SetTickx();
		        pad3->SetTicky();
             	hbpz[it]->SetMarkerColor(kGreen);
             	hbpz[it]->SetMarkerSize(0.1);
             	hbpz[it]->SetLineColor(kGreen);
             	hbpz[it]->SetMarkerStyle(kFullCircle);
				hbpz[it]->GetXaxis()->SetLabelSize(0.05);
				hbpz[it]->GetYaxis()->SetLabelSize(0.05);
			    vz[it]->Draw();
			    hbpz[it]->Draw("same");
				pad4->cd();
		        pad4->SetLogy();
		        pad4->SetTickx();
		        pad4->SetTicky();
				histox[it]->SetLineColor(kMagenta-3);
				histox[it]->Draw();
				histox[it]->GetXaxis()->SetLabelSize(0.05);
				histox[it]->GetYaxis()->SetLabelSize(0.05);
				histoy[it]->Draw("same");
				pad4->BuildLegend(0.13,0.1,0.43,0.35);
				//fft_x.Draw();
				pad5->cd();
		        pad5->SetLogy();
		        pad5->SetTickx();
		        pad5->SetTicky();
				fft_x.SetLineColor(kMagenta-3);
				fft_x.SetMarkerColor(kMagenta-3);
				fft_x.Draw("L");
				fft_y.Draw("Lsame");
				pad5->BuildLegend(0.25,0.1,0.35,0.25);
				pad6->cd();
				pad6->SetTickx();
				pad6->SetTicky();
		        pad6->SetLogy();
				fft_z.Draw("L");
				//TC->Print(Form(prefixname+"%s"+".pdf",output ));
                //double NewSigmaX = histox[it]->GetRMS();
                //double NewSigmaY = histoy[it]->GetRMS();
                //double NewSigmaZ = histoz[it]->GetRMS();
				//----------------------------------------------------------------------------------------
				//TC->SaveAs(Form(prefixname+"%s_"+"Event%d_Du%d"+".pdf",output,i,du_id[it] ));
				//TC->Print(Form(prefixname+"%s"+".pdf",output ));
	            hbpx[it]->Delete();
	            hbpy[it]->Delete();
	            hbpz[it]->Delete();
	            histox[it]->Delete();
	            histoy[it]->Delete();
	            histoz[it]->Delete();
	            pad1->Delete();
	            pad2->Delete();
	            pad3->Delete();
	            pad4->Delete();
	            pad5->Delete();
	            pad6->Delete();
			    cout<<gps_time[2]<<"\t"<<NewSigmaX<<"\t"<<NewSigmaY<<"\t"<<NewSigmaZ<<"\t"<<SigmaA<<"\t"<<gps_temp[0]<<"\t";
			    cout<<atm_humidity[0]*0.61035156/460.375 - 1650/460.375<<"\t"<<atm_pressure[0]*0.61035156/460.375 - 1650/460.375<<"\t"<<atm_temperature[0]*0.61035156/460.375 - 1650/460.375<<"\t";
			    cout<<acceleration_x[0]*0.61035156/19.5 - 400/19.5<<"\t"<<acceleration_y[0]*0.61035156/3300./0.009+0.095/0.009<<"\t"<<acceleration_z[0]*0.61035156*0.048 - 500*0.048<<"\t"<<battery_level[0]*0.61035156*0.109/18<<"\t"<<NewMeanX<<"\t"<<NewMeanY<<"\t"<<NewMeanZ<<"\t"<<MeanA<<"\t"<<W_X<<"\t"<<W_Y<<"\t"<<W_Z<<"\t"<<W_F<<endl;
			}
		}
	}
	TPad *Pad1 = new TPad("frepad1","my Pad1",0.01, 0.0,0.32,0.99);
	TPad *Pad2 = new TPad("frepad2","my Pad2",0.32,0.0,0.99,0.99);
	Pad1->SetLeftMargin(0.15);
	Pad1->SetRightMargin(0.05);
	Pad2->SetLeftMargin(0.01);
	Pad1->SetTickx();
	Pad1->SetTicky();
	Pad1->SetGridy();
	Pad2->SetTickx();
	Pad2->SetTicky();
	Pad1->Draw();
	Pad2->Draw();
	Pad2->SetLogz();
	Pad1->cd();
    his_FX->Draw("hbar");
	Pad1->RedrawAxis();
	Pad1->SetLogx();
	Pad2->cd();
    Time_FreX->GetYaxis()->SetLabelColor(kWhite);	
    Time_FreX->Draw("colz");	
	TC->Print(Form(prefixname+"%s"+".pdf",output ));
	TC->SaveAs(Form(prefixname+"%s"+"_2D_ChannelX.png",output ));
	Pad1->cd();
    his_FY->Draw("hbar");	
	Pad1->RedrawAxis();
	Pad1->SetLogx();
	Pad2->cd();
	Time_FreY->GetYaxis()->SetLabelColor(kWhite);
    Time_FreY->Draw("colz");	
	TC->Print(Form(prefixname+"%s"+".pdf",output ));
	TC->SaveAs(Form(prefixname+"%s"+"_2D_ChannelY.png",output ));
	Pad1->cd();
    his_FZ->Draw("hbar");	
	Pad1->RedrawAxis();
	Pad1->SetLogx();
	Pad2->cd();
	Time_FreZ->GetYaxis()->SetLabelColor(kWhite);
    Time_FreZ->Draw("colz");	
	TC->Print(Form(prefixname+"%s"+".pdf",output ));
	TC->SaveAs(Form(prefixname+"%s"+"_2D_ChannelZ.png",output ));
	tree_eff->Fill();
	TC->Print(Form(prefixname+"%s"+".pdf]",output ));
	fcre->Write();
	fcre->Close();
}
