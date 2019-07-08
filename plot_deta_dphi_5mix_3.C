//Author: zhangliuyao
//Date: 2017/12/2
//compute the projecion from 2D figures. for fig.6
void  plot_deta_dphi_5mix_3()
{
	gROOT->Reset();
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetOptTitle(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetFillColor(0);
	gStyle->SetTitleColor(1);
	gStyle->SetStatColor(0);
	gStyle->SetTextFont(132);
	gStyle->SetLabelFont(132,"X");
	gStyle->SetTitleFont(132,"X");
	gStyle->SetLabelFont(132,"Y");
	gStyle->SetTitleFont(132,"Y");
	gStyle->SetStatFont(132);

	char caseC12[3][32] = {"nofluc","Chain","Triangle"};
	Double_t PI=TMath::Pi();
	TH2D *h2S_1[3];
	TH2D *h2S_2[3];
	TH2D *h2S_3[3];
	TH2D *h2S_4[3];
	TH2D *h2S_5[3];
	TH2D *h2B_1[3];
	TH2D *h2B_2[3];
	TH2D *h2B_3[3];
	TH2D *h2B_4[3];
	TH2D *h2B_5[3];

	char h2SBName1[64];
	char h2SBName2[64];
	char h2SBName3[64];
	char h2SBName4[64];
	char h2SBName5[64];
	for(int k1=0; k1<3; k1++)
	{
		sprintf(h2SBName1,"h2S_1_%d",k1);
		h2S_1[k1] = new TH2D(h2SBName1,h2SBName1,100,-1.6,1.6,200,-1.5,2.*PI-1.5);
		sprintf(h2SBName1,"h2B_1_%d",k1);
		h2B_1[k1] = new TH2D(h2SBName1,h2SBName1,100,-1.6,1.6,200,-1.5,2.*PI-1.5);
	}

	for(int k2=0; k2<3; k2++)
	{
		sprintf(h2SBName2,"h2S_2_%d",k2);
		h2S_2[k2] = new TH2D(h2SBName2,h2SBName2,100,-1.6,1.6,200,-1.5,2.*PI-1.5);
		sprintf(h2SBName2,"h2B_2_%d",k2);
		h2B_2[k2] = new TH2D(h2SBName2,h2SBName2,100,-1.6,1.6,200,-1.5,2.*PI-1.5);
	}

	for(int k3=0; k3<3; k3++)
	{
		sprintf(h2SBName3,"h2S_3_%d",k3);
		h2S_3[k3] = new TH2D(h2SBName3,h2SBName3,100,-1.6,1.6,200,-1.5,2.*PI-1.5);
		sprintf(h2SBName3,"h2B_3_%d",k3);
		h2B_3[k3] = new TH2D(h2SBName3,h2SBName3,100,-1.6,1.6,200,-1.5,2.*PI-1.5);
	}

	for(int k4=0; k4<3; k4++)
	{
		sprintf(h2SBName4,"h2S_4_%d",k4);
		h2S_4[k4] = new TH2D(h2SBName4,h2SBName4,100,-1.6,1.6,200,-1.5,2.*PI-1.5);
		sprintf(h2SBName4,"h2B_4_%d",k4);
		h2B_4[k4] = new TH2D(h2SBName4,h2SBName4,100,-1.6,1.6,200,-1.5,2.*PI-1.5);
	}
	for(int k5=0; k5<3; k5++)
	{
		sprintf(h2SBName5,"h2S_5_%d",k5);
		h2S_5[k5] = new TH2D(h2SBName5,h2SBName5,100,-1.6,1.6,200,-1.5,2.*PI-1.5);
		sprintf(h2SBName5,"h2B_5_%d",k5);
		h2B_5[k5] = new TH2D(h2SBName5,h2SBName5,100,-1.6,1.6,200,-1.5,2.*PI-1.5);
	}

	TFile *fin1[3];
	TFile *fin2[3];
	TFile *fin3[3];
	TFile *fin4[3];
	TFile *fin5[3];
	TH2D *htmp1;
	TH2D *htmp2;
	TH2D *htmp3;
	TH2D *htmp4;
	TH2D *htmp5;
	char finName1[512];
	char finName2[512];
	char finName3[512];
	char finName4[512];
	char finName5[512];
	char histName1[128];
	char histName2[128];
	char histName3[128];
	char histName4[128];
	char histName5[128];
	Double_t theEta=0;
	Double_t theEta1=0;
	Double_t theEta2=0;
	Double_t theEta3=0;
	Double_t theEta4=0;
	Double_t theEta5=0;
	Double_t sumS_1[3]={0};
	Double_t sumS_2[3]={0};
	Double_t sumS_3[3]={0};
	Double_t sumS_4[3]={0};
	Double_t sumS_5[3]={0};
	Double_t sumB_1[3]={0};
	Double_t sumB_2[3]={0};
	Double_t sumB_3[3]={0};
	Double_t sumB_4[3]={0};
	Double_t sumB_5[3]={0};
	Double_t sumS_flow_1[3] = {0};
	Double_t sumB_flow_1[3] = {0};
	Double_t sumS_flow_2[3] = {0};
	Double_t sumB_flow_2[3] = {0};
	Double_t sumS_flow_3[3] = {0};
	Double_t sumB_flow_3[3] = {0};
	Double_t sumS_flow_4[3] = {0};
	Double_t sumB_flow_4[3] = {0};
	Double_t sumS_flow_5[3] = {0};
	Double_t sumB_flow_5[3] = {0};

	TF1 *constant = new TF1("constant","1.",-10,10);
	TH1D *h1S_1_1[3];
	TH1D *h1S_2_1[3];
	TH1D *h1B_1_1[3];
	TH1D *h1B_2_1[3];
	TH1D *h1S_1_2[3];
	TH1D *h1S_2_2[3];
	TH1D *h1B_1_2[3];
	TH1D *h1B_2_2[3];
	TH1D *h1S_1_3[3];
	TH1D *h1S_2_3[3];
	TH1D *h1B_1_3[3];
	TH1D *h1B_2_3[3];
	TH1D *h1S_1_4[3];
	TH1D *h1S_2_4[3];
	TH1D *h1B_1_4[3];
	TH1D *h1B_2_4[3];
	TH1D *h1S_1_5[3];
	TH1D *h1S_2_5[3];
	TH1D *h1B_1_5[3];
	TH1D *h1B_2_5[3];
	Double_t etaCut[3] = {-1.3,0.,1.3};//-1,1 -> -2,2
	Int_t etaBin1[3];
	Int_t etaBin2[3];
	Int_t etaBin3[3];
	Int_t etaBin4[3];
	Int_t etaBin5[3];
	Double_t vnn[3][6]={0};
	Double_t Dphi=0, CDphi=0, vnnC=0;
	TF1 *fCDPhi[3];
	TF1 *fCDPhi_nn[3][3];//c12Case, nn=1,2,3
	char f1Name[64];
	//======================================
	//difPt:signal pi:0.2,kaon:0.3,proton:0.5,lambda:0.6,sigma:0.7
	//===================================== 
	char ntmax2_fitPt[100][100] = { "2_pionplus_pionminus_Pt_0.2_2.5", "2_pionplus_pionplus_Pt_0.2_2.5", "2_pionminus_pionminus_Pt_0.2_2.5",
                                       "2_kaonplus_kaonminus_Pt_0.3_2.5", "2_kaonplus_kaonplus_Pt_0.3_2.5", "2_kaonminus_kaonminus_Pt_0.3_2.5",
                                       "2_proton_protonbar_Pt_0.5_2.5", "2_proton_proton_Pt_0.5_2.5", "2_protonbar_protonbar_Pt_0.5_2.5",
                                       "2_lambda_lambdabar_Pt_0.6_2.5", "2_lambda_lambda_Pt_0.6_2.5", "2_lambdabar_lambdabar_Pt_0.6_2.5"};
     char ntmax2_default_fitPt[100][100] = { "2_pionplus_pionminus_Pt_0.2_2.5_default", "2_pionplus_pionplus_Pt_0.2_2.5_default", "2_pionminus_pionminus_Pt_0.2_2.5_default",
                                         "2_kaonplus_kaonminus_Pt_0.3_2.5_default", "2_kaonplus_kaonplus_Pt_0.3_2.5_default", "2_kaonminus_kaonminus_Pt_0.3_2.5_default",
                                         "2_proton_protonbar_Pt_0.5_2.5_default", "2_proton_proton_Pt_0.5_2.5_default", "2_protonbar_protonbar_Pt_0.5_2.5_default",
                                         "2_lambda_lambdabar_Pt_0.6_2.5_default", "2_lambda_lambda_Pt_0.6_2.5_default", "2_lambdabar_lambdabar_Pt_0.6_2.5_default"};

        char ntmax2_tau0_difPt[100][100] = {"2_pionplus_pionminus_PtMore0.2_tau0","2_pionplus_pionplus_PtMore0.2_tau0","2_pionminus_pionminus_PtMore0.2_tau0",
		"2_kaonplus_kaonminus_PtMore0.3_tau0","2_kaonplus_kaonplus_PtMore0.3_tau0","2_kaonminus_kaonminus_PtMore0.3_tau0",
		"2_proton_protonbar_PtMore0.5_tau0","2_proton_proton_PtMore0.5_tau0","2_protonbar_protonbar_PtMore0.5_tau0",
		"2_lambda_lambdabar_PtMore0.6_tau0","2_lambda_lambda_PtMore0.6_tau0_80mi","2_lambdabar_lambdabar_PtMore0.6_tau0_80mi",
		"2_proton_lambda_PtMore0.5_tau0_80mi","2_protonbar_lambdabar_PtMore0.5_tau0_80mi",
		"2_proton_lambdabar_PtMore0.5_tau0_80mi","2_protonbar_lambda_PtMore0.5_tau0_80mi"};
	char ntmax2_tau0_amongPt[100][100] = {"2_proton_proton_Pt_0.05_1.0_tau0",     "2_protonbar_protonbar_Pt_0.05_1.0_tau0",
		"2_lambda_lambda_Pt_0.05_1.0_tau0_80mi","2_lambdabar_lambdabar_Pt_0.05_1.0_tau0_80mi",
		"2_proton_lambda_Pt_0.05_1.0_tau0_80mi","2_protonbar_lambdabar_Pt_0.05_1.0_tau0_80mi",
		"2_proton_lambdabar_Pt_0.05_1.0_tau0_80mi","2_protonbar_lambda_Pt_0.05_1.0_tau0_80mi"};
	char ntmax2_tau0_fitPt[100][100] = { "2_pionplus_pionminus_Pt_0.2_2.5_tau0_40mi", "2_pionplus_pionplus_Pt_0.2_2.5_tau0_40mi", "2_pionminus_pionminus_Pt_0.2_2.5_tau0_40mi", 
		"2_kaonplus_kaonminus_Pt_0.3_2.5_tau0_40mi", "2_kaonplus_kaonplus_Pt_0.3_2.5_tau0_40mi", "2_kaonminus_kaonminus_Pt_0.3_2.5_tau0_40mi", 
		"2_proton_protonbar_Pt_0.5_2.5_tau0_80mi_again", "2_proton_proton_Pt_0.5_2.5_tau0_80mi_again", "2_protonbar_protonbar_Pt_0.5_2.5_tau0_80mi", 
		"2_lambda_lambdabar_Pt_0.6_2.5_tau0_80mi", "2_lambda_lambda_Pt_0.6_2.5_tau0_80mi", "2_lambdabar_lambdabar_Pt_0.6_2.5_tau0_80mi",
		"2_proton_lambda_Pt_0.5_2.5_tau0_80mi", "2_protonbar_lambdabar_Pt_0.5_2.5_tau0_80mi"};                                      
	char ntmax2_tau0_default_difPt[100][100] = {"2_pionplus_pionminus_PtMore0.2_tau0_default","2_pionplus_pionplus_PtMore0.2_tau0_default","2_pionminus_pionminus_PtMore0.2_tau0_default",
		"2_kaonplus_kaonminus_PtMore0.3_tau0_default","2_kaonplus_kaonplus_PtMore0.3_tau0_default","2_kaonminus_kaonminus_PtMore0.3_tau0_default",
		"2_proton_protonbar_PtMore0.5_tau0_default","2_proton_proton_PtMore0.5_tau0_default","2_protonbar_protonbar_PtMore0.5_tau0_default",
		"2_lambda_lambdabar_PtMore0.6_tau0_default","2_lambda_lambda_PtMore0.6_tau0_default","2_lambdabar_lambdabar_PtMore0.6_tau0_default",
		"2_proton_lambda_PtMore0.5_tau0_default","2_protonbar_lambdabar_PtMore0.5_tau0_default",
		"2_proton_lambdabar_PtMore0.5_tau0_default","2_protonbar_lambda_PtMore0.5_tau0_default"};

	char ntmax2_tau0_default_amongPt[100][100] = {"2_proton_proton_Pt_0.05_1.0_tau0_default","2_protonbar_protonbar_Pt_0.05_1.0_tau0_default",
		"2_lambda_lambda_Pt_0.05_1.0_tau0_default","2_lambdabar_lambdabar_Pt_0.05_1.0_tau0_default",
		"2_proton_lambda_Pt_0.05_1.0_tau0_default","2_protonbar_lambdabar_Pt_0.05_1.0_tau0_default",
		"2_proton_lambdabar_Pt_0.05_1.0_tau0_default","2_protonbar_lambda_Pt_0.05_1.0_tau0_default"};
	char ntmax2_tau0_default_fitPt[100][100] = { "2_pionplus_pionminus_Pt_0.2_2.5_tau0_default", "2_pionplus_pionplus_Pt_0.2_2.5_tau0_default", "2_pionminus_pionminus_Pt_0.2_2.5_tau0_default", 
		"2_kaonplus_kaonminus_Pt_0.3_2.5_tau0_default", "2_kaonplus_kaonplus_Pt_0.3_2.5_tau0_default", "2_kaonminus_kaonminus_Pt_0.3_2.5_tau0_default", 
		"2_proton_protonbar_Pt_0.5_2.5_tau0_default_again", "2_proton_proton_Pt_0.5_2.5_tau0_default_again", "2_protonbar_protonbar_Pt_0.5_2.5_tau0_default", 
		"2_lambda_lambdabar_Pt_0.6_2.5_tau0_default", "2_lambda_lambda_Pt_0.6_2.5_tau0_default", "2_lambdabar_lambdabar_Pt_0.6_2.5_tau0_default",
		"2_proton_lambda_Pt_0.5_2.5_tau0_default", "2_protonbar_lambdabar_Pt_0.5_2.5_tau0_default"};                                      

	char ntmax2_fitPt_newCoal[100][100] = {
                "2_pionplus_pionplus_Pt_0.2_2.5_newCoal_0mb","2_pionplus_pionplus_Pt_0.2_2.5_newCoal_1.5mb","2_pionplus_pionplus_Pt_0.2_2.5_newCoal_3mb","2_pionplus_pionplus_Pt_0.2_2.5_newCoal_6mb",
		"2_proton_proton_Pt_0.5_2.5_newCoal_0mb","2_proton_proton_Pt_0.5_2.5_newCoal_1.5mb","2_proton_proton_Pt_0.5_2.5_newCoal_3mb","2_proton_proton_Pt_0.5_2.5_newCoal_6mb", 
                "2_protonbar_protonbar_Pt_0.5_2.5_newCoal_0mb","2_protonbar_protonbar_Pt_0.5_2.5_newCoal_1.5mb","2_protonbar_protonbar_Pt_0.5_2.5_newCoal_3mb","2_protonbar_protonbar_Pt_0.5_2.5_newCoal_6mb",
                "2_protonbar_protonbar_Pt_0.01_2.5_newCoal_0mb","2_protonbar_protonbar_Pt_0.01_2.5_newCoal_1.5mb","2_protonbar_protonbar_Pt_0.01_2.5_newCoal_3mb",
                "2_protonbar_protonbar_Pt_0.01_1.0_newCoal_0mb","2_protonbar_protonbar_Pt_0.01_1.0_newCoal_1.5mb","2_protonbar_protonbar_Pt_0.01_1.0_newCoal_3mb"}; 

	char ntmax100_difPt[100][100] = {"100_pionplus_pionminus_PtMore0.2","100_pionplus_pionplus_PtMore0.2","100_pionminus_pionminus_PtMore0.2",
		"100_kaonplus_kaonminus_PtMore0.3","100_kaonplus_kaonplus_PtMore0.3","100_kaonminus_kaonminus_PtMore0.3",
		"100_proton_protonbar_PtMore0.5","100_proton_proton_PtMore0.5","100_protonbar_protonbar_PtMore0.5",
		"100_lambda_lambdabar_PtMore0.6","100_lambda_lambda_PtMore0.6","100_lambdabar_lambdabar_PtMore0.6",
		"100_proton_labmda_PtMore0.5", "100_protonbar_lambdabar_PtMore0.5",
		"100_proton_lambdabar_PtMore0.5","100_protonbar_lambda_PtMore0.5",
		"100_sigmaplus_sigmaminus_PtMore0.6","100_sigmaplus_sigmaplus_PtMore0.6","100_sigmaminus_sigmaminus_PtMore0.6"}; 
	char ntmax100_amongPt[100][100] = {"100_proton_proton_Pt_0.05_1.0","100_protonbar_protonbar_Pt_0.05_1.0",
		"100_lambda_lambda_Pt_0.05_1.0","100_lambdabar_lambdabar_Pt_0.05_1.0",
		"100_proton_lambda_Pt_0.05_1.0","100_protonbar_lambdabar_Pt_0.05_1.0",
		"100_proton_lambdabar_Pt_0.05_1.0","100_protonbar_lambda_Pt_0.05_1.0"
	};
	char ntmax100_fitPt[100][100] = { "100_pionplus_pionminus_Pt_0.2_2.5", "100_pionplus_pionplus_Pt_0.2_2.5", "100_pionminus_pionminus_Pt_0.2_2.5", 
		"100_kaonplus_kaonminus_Pt_0.3_2.5", "100_kaonplus_kaonplus_Pt_0.3_2.5", "100_kaonminus_kaonminus_Pt_0.3_2.5", 
		"100_proton_protonbar_Pt_0.5_2.5_again", "100_proton_proton_Pt_0.5_2.5_again", "100_protonbar_protonbar_Pt_0.5_2.5", 
		"100_lambda_lambdabar_Pt_0.5_2.5_again", "100_lambda_lambda_Pt_0.6_2.5", "100_lambdabar_lambdabar_Pt_0.6_2.5",
		"100_lambda_lambda_Pt_0.6_2.5_fixed_background",
		"100_proton_lambda_Pt_0.5_2.5", "100_protonbar_lambdabar_Pt_0.5_2.5"};
        char ntmax100_fitPt_newCoal[100][100] = { "100_pionplus_pionminus_Pt_0.2_2.5_newCoal", "100_pionplus_pionplus_Pt_0.2_2.5_newCoal", "100_pionminus_pionminus_Pt_0.2_2.5_newCoal",
                                                 "100_kaonplus_kaonminus_Pt_0.3_2.5_newCoal", "100_kaonplus_kaonplus_Pt_0.3_2.5_newCoal", "100_kaonminus_kaonminus_Pt_0.3_2.5_newCoal",
                                                 "100_proton_protonbar_Pt_0.5_2.5_newCoal", "100_proton_proton_Pt_0.5_2.5_newCoal", "100_protonbar_protonbar_Pt_0.5_2.5_newCoal",
                                                 "100_lambda_lambdabar_Pt_0.6_2.5_newCoal", "100_lambda_lambda_Pt_0.6_2.5_newCoal", "100_lambdabar_lambdabar_Pt_0.6_2.5_newCoal",
                                                 "100_proton_lambda_Pt_0.5_2.5_newCoal", "100_protonbar_lambdabar_Pt_0.5_2.5_newCoal"};
 
        char ntmax100_amongPt_newCoal[100][100] = {"100_proton_proton_Pt_0.01_1.0_newCoal","100_protonbar_protonbar_Pt_0.01_1.0_newCoal",
                                          "100_lambda_lambda_Pt_0.01_1.0_newCoal","100_lambdabar_lambdabar_Pt_0.01_1.0_newCoal"};

	char ntmax100_default_difPt[100][100] = {"100_pionplus_pionminus_PtMore0.2_default","100_pionplus_pionplus_PtMore0.2_default","100_pionminus_pionminus_PtMore0.2_default",
		"100_kaonplus_kaonminus_PtMore0.3_default","100_kaonplus_kaonplus_PtMore0.3_default","100_kaonminus_kaonminus_PtMore0.3_default",
		"100_proton_protonbar_PtMore0.5_default","100_proton_proton_PtMore0.5_default","100_protonbar_protonbar_PtMore0.5_default",
		"100_lambda_lambdabar_PtMore0.6_default","100_lambda_lambda_PtMore0.6_default","100_lambdabar_lambdabar_PtMore0.6_default",
		"100_proton_lambda_PtMore0.5_default", "100_protonbar_lambdabar_PtMore0.5_default",
		"100_proton_lambdabar_PtMore0.5_default","100_protonbar_lambda_PtMore0.5_default"}; 

	char ntmax100_default_amongPt[100][100] = {"100_proton_proton_Pt_0.05_1.0_default","100_protonbar_protonbar_Pt_0.05_1.0_default",
		"100_lambda_lambda_Pt_0.05_1.0_default","100_lambdabar_lambdabar_Pt_0.05_1.0_default",
		"100_proton_lambda_Pt_0.05_1.0_default","100_protonbar_lambdabar_Pt_0.05_1.0_default",
		"100_proton_lambdabar_Pt_0.05_1.0_default","100_protonbar_lambda_Pt_0.05_1.0_default"};
	char ntmax100_default_fitPt[100][100] = { "100_pionplus_pionminus_Pt_0.2_2.5_default", "100_pionplus_pionplus_Pt_0.2_2.5_default", "100_pionminus_pionminus_Pt_0.2_2.5_default", 
		"100_kaonplus_kaonminus_Pt_0.3_2.5_default", "100_kaonplus_kaonplus_Pt_0.3_2.5_default", "100_kaonminus_kaonminus_Pt_0.3_2.5_default", 
		"100_proton_protonbar_Pt_0.5_2.5_default_again", "100_proton_proton_Pt_0.5_2.5_default_again", "100_protonbar_protonbar_Pt_0.5_2.5_default", 
		"100_lambda_lambdabar_Pt_0.6_2.5_default", "100_lambda_lambda_Pt_0.6_2.5_default", "100_lambdabar_lambdabar_Pt_0.6_2.5_default",
		"100_proton_lambda_Pt_0.5_2.5_default", "100_protonbar_lambdabar_Pt_0.5_2.5_default"};                                      
	char ntmax300_amongPt[100][100] = {"300_proton_proton_Pt_0.05_1.0","300_protonbar_protonbar_Pt_0.05_1.0",
		"300_lambda_lambda_Pt_0.05_1.0","300_lambdabar_lambdabar_Pt_0.05_1.0"};

	char ntmax300_default_difPt[100][100] = {"300_pionplus_pionminus_PtMore0.2_default","300_pionplus_pionplus_PtMore0.2_default","300_pionminus_pionminus_PtMore0.2_default",
		"300_kaonplus_kaonminus_PtMore0.3_default","300_kaonplus_kaonplus_PtMore0.3_default","300_kaonminus_kaonminus_PtMore0.3_default",
		"300_proton_protonbar_PtMore0.5_default","300_proton_proton_PtMore0.5_default","300_protonbar_protonbar_PtMore0.5_default",
		"300_lambda_lambdabar_PtMore0.6_default","300_lambda_lambda_PtMore0.6_default","300_lambdabar_lambdabar_PtMore0.6_default",
		"300_sigmaplus_sigmaminus_PtMore0.6_default","300_sigmaplus_sigmaplus_PtMore0.6_default","300_sigmaminus_sigmaminus_PtMore0.6_default"};

	char title_name[100][100] = {" #pi^{+}-#pi^{-} pairs"," #pi^{+}-#pi^{+} pairs"," #pi^{-}-#pi^{-} pairs",
		" K^{+}-K^{-} pairs"," K^{+}-K^{+} pairs"," K^{-}-K^{-} pairs",
		" p-#bar{p} pairs"," p-p pairs", " #bar{p}-#bar{p} pairs",
		" #Lambda-#bar{#Lambda} paris", " #Lambda-#Lambda pairs"," #bar{#Lambda}-#bar{#Lambda} pairs",
		" p-#Lambda pairs ", " #bar{p}-#bar{#Lambda} pairs", " p-#bar{#Lambda} pairs", "#bar{p}-#Lambda", 
		"(e1) #Sigma^{+}#Sigma^{-}","(e2) #Sigma^{+}#Sigma^{+}","(e3) #Sigma^{-}#Sigma^{-}"
			"p#Lambda","#bar{p}#bar{#Lambda}", "p#bar{#Lambda}",
		"#bar{p}#Lambda"}; 	
	char file_name[100][100] = {"pionplus_pionminus_fit","pionplus_pionplus_newCoal_last","pionminus_pionminus_newCoal_last",
		"kaonplus_kaonminus_fit","kaonplus_kaonplus_newCoal_last","kaonminus_kaonminus_newCoal_last",
		"proton_protonbar_fit_again","proton_proton_lowPt_newCoal_last","protonbar_protonbar_lowPt_newCoal_last",
		"lambda_lambdabar_fit","lambda_lambda_lowPt_newCoal_last","lambdabar_lambdabar_lowPt_newCoal_last",
		"proton_lambda_fit", "protonbar_lambdabar_fit", "proton_lambdabar", "protonbar_lambda", 
		"sigmaplus_sigmaminus","sigmaplus_sigmaplus","sigmaminus_sigmaminus"};
	int ip_number=11; 
	int ip_number1=4; 
	int ip_number2=3; 
	int ip_number3=11; 
	int ip_number4=3; 
	int ip_number5=10; 
	int title_number=10; 
	//==========================================================================
	//data declaration; 
	TGraphAsymmErrors *graph[3];
	TGraphAsymmErrors *grae[3];
	TDirectory *f2;
	TDirectory *dir;

	char mul_y[3][20] = {"y1", "y2", "y3"};
	char tit[3][10] = {"y1", "y2", "y3"};
	char finName[50];
	char title_y[100];
	Double_t graph_exl[3][29] = {};
	Double_t graph_exh[3][29] = {};
	Double_t graph_eyl[3][29] = {};
	Double_t graph_eyh[3][29] = {};
	//-----------------------------------------
	for(int ic=0; ic<1; ic++)//  only one file,so ic=3->1 
	{ 
		//================================================================================================================================================	
		 //for ntmax2_fitPt
	sprintf(finName1,"/home/zhangliuyao/zhangliuyao_scratch/analysis-flow-2PC-etaGap4/results/result-pp_ntmax2_PseudoRapidity/result-pp_ntmax2_difPt_rootfile/result-pp_ntmax%s.root", ntmax2_fitPt[ip_number1]);
        sprintf(finName2,"/home/zhangliuyao/zhangliuyao_scratch/analysis-flow-2PC-etaGap2/results/result-pp_ntmax100_PseudoRapidity/result-pp_ntmax100_difPt_rootfile/result-pp_ntmax%s.root", ntmax100_amongPt[ip_number2]);

		//for ntmamx2_default_fitPt
       sprintf(finName3,"/home/zhangliuyao/zhangliuyao_scratch/analysis-flow-2PC-etaGap4/results/result-pp_ntmax2_PseudoRapidity/result-pp_ntmax2_default_difPt_rootfile/result-pp_ntmax%s.root", ntmax2_default_fitPt[ip_number3]);
                // for ntmax2_tau0_difPt
		//sprintf(finName1,"/home/zhangliuyao/zhangliuyao_scratch/analysis-flow-2PC-etaGap4/results/result-pp_ntmax2_PseudoRapidity/result-pp_ntmax2_tau0_difPt_rootfile/result-pp_ntmax%s.root", ntmax2_tau0_difPt[ip_number1]);
		//================================================================================================================================================	
		// for ntmax2_tau0_amongPt
		//sprintf(finName1,"/home/zhangliuyao/zhangliuyao_scratch/analysis-flow-2PC-etaGap4/results/result-pp_ntmax2_PseudoRapidity/result-pp_ntmax2_tau0_difPt_rootfile/result-pp_ntmax%s.root", ntmax2_tau0_amongPt[ip_number1]);
		//sprintf(finName1,"/home/zhangliuyao/zhangliuyao_scratch/analysis-flow-2PC-etaGap4/results/result-pp_ntmax2_PseudoRapidity/result-pp_ntmax2_newCoal_rootfile/result-pp_ntmax%s.root", ntmax2_fitPt_newCoal[ip_number1]);
		//sprintf(finName2,"/home/zhangliuyao/zhangliuyao_scratch/analysis-flow-2PC-etaGap4/results/result-pp_ntmax2_PseudoRapidity/result-pp_ntmax2_newCoal_rootfile/result-pp_ntmax%s.root", ntmax2_fitPt_newCoal[ip_number2]);
		////================================================================================================================================================	
		//// for ntmax2_tau0_default_difPt
		//sprintf(finName3,"/home/zhangliuyao/zhangliuyao_scratch/analysis-flow-2PC-etaGap4/results/result-pp_ntmax2_PseudoRapidity/result-pp_ntmax2_newCoal_rootfile/result-pp_ntmax%s.root", ntmax2_fitPt_newCoal[ip_number3]);
		//===============================================================================================================================================
		//for ntmax100_ntmax100_amongPt
		//sprintf(finName2,"/home/zhangliuyao/zhangliuyao_scratch/analysis-flow-2PC-etaGap2/results/result-pp_ntmax100_PseudoRapidity/result-pp_ntmax100_difPt_rootfile/result-pp_ntmax%s.root", ntmax100_fitPt[ip_number2]);
		//===============================================================================================================================================
		//for ntmax100_difPt
		//sprintf(finName2,"/home/zhangliuyao/zhangliuyao_scratch/analysis-flow-2PC-etaGap2/results/result-pp_ntmax100_PseudoRapidity/result-pp_ntmax100_difPt_rootfile/result-pp_ntmax%s.root", ntmax100_difPt[ip_number2]);
		//===============================================================================================================================================
		//for ntmax100_difPt
		sprintf(finName5,"/home/zhangliuyao/zhangliuyao_scratch/analysis-flow-2PC-etaGap2/results/result-pp_ntmax100_PseudoRapidity/result-pp_ntmax100_difPt_rootfile/result-pp_ntmax%s.root", ntmax100_fitPt[ip_number5]);
		//===============================================================================================================================================
		sprintf(finName4,"/home/zhangliuyao/zhangliuyao_scratch/analysis-flow-2PC-etaGap2/results/result-pp_ntmax100_PseudoRapidity/result-pp_ntmax100_difPt_rootfile_newCoal/result-pp_ntmax%s.root", ntmax100_amongPt_newCoal[ip_number4]);
		//===================================================================================================================================================
		//for ntmax300_difPt
		//sprintf(finName3,"/home/zhangliuyao/zhangliuyao_scratch/analysis-flow-2PC-etaGap/results/result-pp_ntmax300_PseudoRapidity/result-pp_ntmax300_difPt_rootfile/result-pp_ntmax%s.root",ntmax300_amongPt[ip_number3]);
		//===================================================================================================================================================
		//for ntmax300_amongPt
		//sprintf(finName3,"/home/zhangliuyao/zhangliuyao_scratch/analysis-flow-2PC-etaGap/results/result-pp_ntmax300_PseudoRapidity/result-pp_ntmax300_difPt_rootfile/result-pp_ntmax%s.root",ntmax300_amongPt[ip_number3]);
		//===================================================================================================================================================
		//for ntmax300_default_difPt
		//sprintf(finName3,"/home/zhangliuyao/zhangliuyao_scratch/analysis-flow-2PC-etaGap/results/result-pp_ntmax300_PseudoRapidity/result-pp_ntmax300_difPt_default_rootfile/result-pp_ntmax%s.root",ntmax300_default_difPt[ip_number3]);
		cout << finName1 << endl; 
		cout << finName2 << endl; 
		cout << finName3 << endl; 
		cout << finName4 << endl; 
		cout << finName5 << endl; 
		fin1[ic] = new TFile(finName1,"READ");
		fin2[ic] = new TFile(finName2,"READ");
		fin3[ic] = new TFile(finName3,"READ");
		fin4[ic] = new TFile(finName4,"READ");
		fin5[ic] = new TFile(finName5,"READ");

		for(int ih1=0; ih1<100; ih1++)
		{
			sprintf(histName1,"Sh2_Deta_Dphi_%d",ih1);
			htmp1 = (TH2D *) fin1[ic]->Get(histName1);
			h2S_1[ic]->Add(htmp1);

			sprintf(histName1,"Bh2_Deta_Dphi_%d",ih1);
			htmp1 = (TH2D *) fin1[ic]->Get(histName1);
			h2B_1[ic]->Add(htmp1);
		}
		for(int ih2=0; ih2<1; ih2++)
		{
			sprintf(histName2,"Sh2_Deta_Dphi_%d",ih2);
			htmp2 = (TH2D *) fin2[ic]->Get(histName2);
			h2S_2[ic]->Add(htmp2);

			sprintf(histName2,"Bh2_Deta_Dphi_%d",ih2);
			htmp2 = (TH2D *) fin2[ic]->Get(histName2);
			h2B_2[ic]->Add(htmp2);
		}
		for(int ih3=0; ih3<100; ih3++)
		{
			sprintf(histName3,"Sh2_Deta_Dphi_%d",ih3);
			htmp3 = (TH2D *) fin3[ic]->Get(histName3);
			h2S_3[ic]->Add(htmp3);

			sprintf(histName3,"Bh2_Deta_Dphi_%d",ih3);
			htmp3 = (TH2D *) fin3[ic]->Get(histName3);
			h2B_3[ic]->Add(htmp3);
		}
		for(int ih4=0; ih4<100; ih4++)
		{
			sprintf(histName4,"Sh2_Deta_Dphi_%d",ih4);
			htmp4 = (TH2D *) fin4[ic]->Get(histName4);
			h2S_4[ic]->Add(htmp4);

			sprintf(histName4,"Bh2_Deta_Dphi_%d",ih4);
			htmp4 = (TH2D *) fin4[ic]->Get(histName4);
			h2B_4[ic]->Add(htmp4);
		}
		for(int ih5=0; ih5<1; ih5++)
		{
			sprintf(histName5,"Sh2_Deta_Dphi_%d",ih5);
			htmp5 = (TH2D *) fin5[ic]->Get(histName5);
			h2S_5[ic]->Add(htmp5);

			sprintf(histName5,"Bh2_Deta_Dphi_%d",ih5);
			htmp5 = (TH2D *) fin5[ic]->Get(histName5);
			h2B_5[ic]->Add(htmp5);
		}
		//fin.Close();
		// his is for one dimensional delta phi distribution. 
		for(int i1=1; i1<=h2S_1[ic]->GetNbinsX(); i1++)
		{
			theEta1 = h2S_1[ic]->GetXaxis()->GetBinCenter(i1);
			if(theEta1>1.3 || theEta1<-1.3) continue;
			for(int j1=1; j1<=h2S_1[ic]->GetNbinsY(); j1++)
			{
				sumS_flow_1[ic] += h2S_1[ic]->GetBinContent(i1,j1);
				sumB_flow_1[ic] += h2B_1[ic]->GetBinContent(i1,j1);
			}
		}
		for(int i2=1; i2<=h2S_2[ic]->GetNbinsX(); i2++)
		{
			theEta2 = h2S_2[ic]->GetXaxis()->GetBinCenter(i2);
			if(theEta2>1.3 || theEta2<-1.3) continue;
			for(int j2=1; j2<=h2S_2[ic]->GetNbinsY(); j2++)
			{
				sumS_flow_2[ic] += h2S_2[ic]->GetBinContent(i2,j2);
				sumB_flow_2[ic] += h2B_2[ic]->GetBinContent(i2,j2);
			}
		}
		for(int i3=1; i3<=h2S_3[ic]->GetNbinsX(); i3++)
		{
			theEta3 = h2S_3[ic]->GetXaxis()->GetBinCenter(i3);
			if(theEta3>1.3 || theEta3<-1.3) continue;
			for(int j3=1; j3<=h2S_3[ic]->GetNbinsY(); j3++)
			{
				sumS_flow_3[ic] += h2S_3[ic]->GetBinContent(i3,j3);
				sumB_flow_3[ic] += h2B_3[ic]->GetBinContent(i3,j3);
			}
		}
		for(int i4=1; i4<=h2S_4[ic]->GetNbinsX(); i4++)
		{
			theEta4 = h2S_4[ic]->GetXaxis()->GetBinCenter(i4);
			if(theEta4>1.3 || theEta4<-1.3) continue;
			for(int j4=1; j4<=h2S_4[ic]->GetNbinsY(); j4++)
			{
				sumS_flow_4[ic] += h2S_4[ic]->GetBinContent(i4,j4);
				sumB_flow_4[ic] += h2B_4[ic]->GetBinContent(i4,j4);
			}
		}
		for(int i5=1; i5<=h2S_5[ic]->GetNbinsX(); i5++)
		{
			theEta5 = h2S_5[ic]->GetXaxis()->GetBinCenter(i5);
			if(theEta5>1.3 || theEta5<-1.3) continue;
			for(int j5=1; j5<=h2S_5[ic]->GetNbinsY(); j5++)
			{
				sumS_flow_5[ic] += h2S_5[ic]->GetBinContent(i5,j5);
				sumB_flow_5[ic] += h2B_5[ic]->GetBinContent(i5,j5);
			}
		}

		cout << "number of axis:" << h2S_1[ic]->GetNbinsX()<< endl;  
		cout << "number of axis:" << h2S_2[ic]->GetNbinsX()<< endl;  
		cout << "number of axis:" << h2S_3[ic]->GetNbinsX()<< endl;  
		cout << "number of axis:" << h2S_4[ic]->GetNbinsX()<< endl;  
		cout << "number of axis:" << h2S_5[ic]->GetNbinsX()<< endl;  
		for(int ie1=0; ie1<3; ie1++) etaBin1[ie1] = h2S_1[0]->GetXaxis()->FindBin(etaCut[ie1]);
		for(int ie2=0; ie2<3; ie2++) etaBin2[ie2] = h2S_2[0]->GetXaxis()->FindBin(etaCut[ie2]);
		for(int ie3=0; ie3<3; ie3++) etaBin3[ie3] = h2S_3[0]->GetXaxis()->FindBin(etaCut[ie3]);
		for(int ie4=0; ie4<3; ie4++) etaBin4[ie4] = h2S_4[0]->GetXaxis()->FindBin(etaCut[ie4]);
		for(int ie5=0; ie5<3; ie5++) etaBin5[ie5] = h2S_5[0]->GetXaxis()->FindBin(etaCut[ie5]);
		h1S_1_1[ic] = (TH1D *) h2S_1[ic]->ProjectionY("h1S_1_1[ic]",etaBin1[0],etaBin1[1]-1,"");
		h1S_2_1[ic] = (TH1D *) h2S_1[ic]->ProjectionY("h1S_2_1[ic]",etaBin1[1]-1,etaBin1[2]-1,"");
		h1S_1_1[ic]->Add(h1S_2_1[ic]);//prc86,014907(2012) :formula 15 
		h1B_1_1[ic] = (TH1D *) h2B_1[ic]->ProjectionY("h1B_1_1[ic]",etaBin1[0],etaBin1[1]-1,"");
		h1B_2_1[ic] = (TH1D *) h2B_1[ic]->ProjectionY("h1B_2_1[ic]",etaBin1[1]-1,etaBin1[2]-1,"");
		h1B_1_1[ic]->Add(h1B_2_1[ic]);


		h1S_1_2[ic] = (TH1D *) h2S_2[ic]->ProjectionY("h1S_1_2[ic]",etaBin2[0],etaBin2[1]-1,"");
		h1S_2_2[ic] = (TH1D *) h2S_2[ic]->ProjectionY("h1S_2_2[ic]",etaBin2[1]-1,etaBin2[2]-1,"");
		h1S_1_2[ic]->Add(h1S_2_2[ic]);//prc86,014907(2012) :formula 15 
		h1B_1_2[ic] = (TH1D *) h2B_2[ic]->ProjectionY("h1B_1_2[ic]",etaBin2[0],etaBin2[1]-1,"");
		h1B_2_2[ic] = (TH1D *) h2B_2[ic]->ProjectionY("h1B_2_2[ic]",etaBin2[1]-1,etaBin2[2]-1,"");
		h1B_1_2[ic]->Add(h1B_2_2[ic]);


		h1S_1_3[ic] = (TH1D *) h2S_3[ic]->ProjectionY("h1S_1_3[ic]",etaBin3[0],etaBin3[1]-1,"");
		h1S_2_3[ic] = (TH1D *) h2S_3[ic]->ProjectionY("h1S_2_3[ic]",etaBin3[1]-1,etaBin3[2]-1,"");
		h1S_1_3[ic]->Add(h1S_2_3[ic]);//prc86,014907(2012) :formula 15 
		h1B_1_3[ic] = (TH1D *) h2B_3[ic]->ProjectionY("h1B_1_3[ic]",etaBin3[0],etaBin3[1]-1,"");
		h1B_2_3[ic] = (TH1D *) h2B_3[ic]->ProjectionY("h1B_2_3[ic]",etaBin3[1]-1,etaBin3[2]-1,"");
		h1B_1_3[ic]->Add(h1B_2_3[ic]);

		h1S_1_4[ic] = (TH1D *) h2S_4[ic]->ProjectionY("h1S_1_4[ic]",etaBin4[0],etaBin4[1]-1,"");
		h1S_2_4[ic] = (TH1D *) h2S_4[ic]->ProjectionY("h1S_2_4[ic]",etaBin4[1]-1,etaBin4[2]-1,"");
		h1S_1_4[ic]->Add(h1S_2_4[ic]);//prc86,014907(2012) :formula 15 
		h1B_1_4[ic] = (TH1D *) h2B_4[ic]->ProjectionY("h1B_1_4[ic]",etaBin4[0],etaBin4[1]-1,"");
		h1B_2_4[ic] = (TH1D *) h2B_4[ic]->ProjectionY("h1B_2_4[ic]",etaBin4[1]-1,etaBin4[2]-1,"");
		h1B_1_4[ic]->Add(h1B_2_4[ic]);

		h1S_1_5[ic] = (TH1D *) h2S_5[ic]->ProjectionY("h1S_1_5[ic]",etaBin5[0],etaBin5[1]-1,"");
		h1S_2_5[ic] = (TH1D *) h2S_5[ic]->ProjectionY("h1S_2_5[ic]",etaBin5[1]-1,etaBin5[2]-1,"");
		h1S_1_5[ic]->Add(h1S_2_5[ic]);//prc86,015907(2012) :formula 15 
		h1B_1_5[ic] = (TH1D *) h2B_5[ic]->ProjectionY("h1B_1_5[ic]",etaBin5[0],etaBin5[1]-1,"");
		h1B_2_5[ic] = (TH1D *) h2B_5[ic]->ProjectionY("h1B_2_5[ic]",etaBin5[1]-1,etaBin5[2]-1,"");
		h1B_1_5[ic]->Add(h1B_2_5[ic]);

		h1S_1_1[ic]->Sumw2();
		h1B_1_1[ic]->Sumw2();
		h1S_1_1[ic]->Divide(constant,sumS_flow_1[ic]/sumB_flow_1[ic]);//TH1::Divide(TF1*f1,Double_t c1), perform the operation : this = this/c1*f1. 
		h1S_1_1[ic]->Divide(h1B_1_1[ic]);

		h1S_1_2[ic]->Sumw2();
		h1B_1_2[ic]->Sumw2();
		h1S_1_2[ic]->Divide(constant,sumS_flow_2[ic]/sumB_flow_2[ic]);//TH1::Divide(TF1*f1,Double_t c1), perform the operation : this = this/c1*f1. 
		h1S_1_2[ic]->Divide(h1B_1_2[ic]);

		h1S_1_3[ic]->Sumw2();
		h1B_1_3[ic]->Sumw2();
		h1S_1_3[ic]->Divide(constant,sumS_flow_3[ic]/sumB_flow_3[ic]);//TH1::Divide(TF1*f1,Double_t c1), perform the operation : this = this/c1*f1. 
		h1S_1_3[ic]->Divide(h1B_1_3[ic]);

		h1S_1_4[ic]->Sumw2();
		h1B_1_4[ic]->Sumw2();
		h1S_1_4[ic]->Divide(constant,sumS_flow_4[ic]/sumB_flow_4[ic]);//TH1::Divide(TF1*f1,Double_t c1), perform the operation : this = this/c1*f1. 
		h1S_1_4[ic]->Divide(h1B_1_4[ic]);

		h1S_1_5[ic]->Sumw2();
		h1B_1_5[ic]->Sumw2();
		h1S_1_5[ic]->Divide(constant,sumS_flow_5[ic]/sumB_flow_5[ic]);//TH1::Divide(TF1*f1,Double_t c1), perform the operation : this = this/c1*f1. 
		h1S_1_5[ic]->Divide(h1B_1_5[ic]);
		//for plot
		h1S_1_1[ic]->Rebin(20);   //after Rebin(8),then divide 8 . 
		h1S_1_1[ic]->Divide(constant,20);
		h1S_1_1[ic]->SetMarkerStyle(25);
		h1S_1_1[ic]->SetMarkerSize(1.3);
		h1S_1_1[ic]->SetMarkerColor(6);
		h1S_1_1[ic]->SetLineColor(6);
		h1S_1_1[ic]->SetLineWidth(2);

		h1S_1_2[ic]->Rebin(20);   //after Rebin(8),then divide 8 . 
		h1S_1_2[ic]->Divide(constant,20);
		h1S_1_2[ic]->SetMarkerStyle(26);
		h1S_1_2[ic]->SetMarkerSize(1.3); 
		h1S_1_2[ic]->SetMarkerColor(2); 
		h1S_1_2[ic]->SetLineColor(2); 
		h1S_1_2[ic]->SetLineWidth(2); 

		h1S_1_3[ic]->Rebin(10);   //after Rebin(5),then divide 8 . 
		h1S_1_3[ic]->Divide(constant,10);
		h1S_1_3[ic]->SetMarkerStyle(24);
		h1S_1_3[ic]->SetMarkerSize(1.3); 
		h1S_1_3[ic]->SetMarkerColor(1); 
		h1S_1_3[ic]->SetLineColor(1); 
		h1S_1_3[ic]->SetLineWidth(2); 

		h1S_1_4[ic]->Rebin(20);   //after Rebin(8),then divide 8 . 
		h1S_1_4[ic]->Divide(constant,20);
		h1S_1_4[ic]->SetMarkerStyle(28);
		h1S_1_4[ic]->SetMarkerSize(1.3); 
		h1S_1_4[ic]->SetMarkerColor(4); 
		h1S_1_4[ic]->SetLineColor(4); 
		h1S_1_4[ic]->SetLineWidth(2); 

		h1S_1_5[ic]->Rebin(5);   //after Rebin(8),then divide 8 . 
		h1S_1_5[ic]->Divide(constant,5);
		h1S_1_5[ic]->SetMarkerStyle(30);
		h1S_1_5[ic]->SetMarkerSize(1.3); 
		h1S_1_5[ic]->SetMarkerColor(4); 
		h1S_1_5[ic]->SetLineColor(4); 
		h1S_1_5[ic]->SetLineWidth(2); 
		//this is two distribution of delta phi and delta eta 
		h2S_1[ic]->Rebin2D(4,4);
		h2S_2[ic]->Rebin2D(4,4);
		h2S_3[ic]->Rebin2D(4,4);
		h2S_4[ic]->Rebin2D(4,4);
		h2S_5[ic]->Rebin2D(4,4);
		h2B_1[ic]->Rebin2D(4,4);
		h2B_2[ic]->Rebin2D(4,4);
		h2B_3[ic]->Rebin2D(4,4);
		h2B_4[ic]->Rebin2D(4,4);
		h2B_5[ic]->Rebin2D(4,4);

		for(int ie1=1; ie1<=h2S_1[ic]->GetNbinsX(); ie1++)
		{
			theEta = h2S_1[ic]->GetXaxis()->GetBinCenter(ie1);
			if(theEta>1.3 || theEta<-1.3) continue;
			for(int jp1=1; jp1<=h2S_1[ic]->GetNbinsY(); jp1++)
			{
				sumS_1[ic] += h2S_1[ic]->GetBinContent(ie1,jp1);
				sumB_1[ic] += h2B_1[ic]->GetBinContent(ie1,jp1);
			}
		}
		h2S_1[ic]->Divide(constant,sumS_1[ic]/sumB_1[ic]);
		h2S_1[ic]->Divide(h2B_1[ic]);

		for(int ie2=1; ie2<=h2S_2[ic]->GetNbinsX(); ie2++)
		{
			theEta = h2S_2[ic]->GetXaxis()->GetBinCenter(ie2);
			if(theEta>1.3 || theEta<-1.3) continue;
			for(int jp2=1; jp2<=h2S_2[ic]->GetNbinsY(); jp2++)
			{
				sumS_2[ic] += h2S_2[ic]->GetBinContent(ie2,jp2);
				sumB_2[ic] += h2B_2[ic]->GetBinContent(ie2,jp2);
			}
		}
		h2S_2[ic]->Divide(constant,sumS_2[ic]/sumB_2[ic]);
		h2S_2[ic]->Divide(h2B_2[ic]);

		for(int ie3=1; ie3<=h2S_3[ic]->GetNbinsX(); ie3++)
		{
			theEta = h2S_3[ic]->GetXaxis()->GetBinCenter(ie3);
			if(theEta>1.3 || theEta<-1.3) continue;
			for(int jp3=1; jp3<=h2S_3[ic]->GetNbinsY(); jp3++)
			{
				sumS_3[ic] += h2S_3[ic]->GetBinContent(ie3,jp3);
				sumB_3[ic] += h2B_3[ic]->GetBinContent(ie3,jp3);
			}
		}
		h2S_3[ic]->Divide(constant,sumS_3[ic]/sumB_3[ic]);
		h2S_3[ic]->Divide(h2B_3[ic]);

		for(int ie4=1; ie4<=h2S_4[ic]->GetNbinsX(); ie4++)
		{
			theEta = h2S_4[ic]->GetXaxis()->GetBinCenter(ie4);
			if(theEta>1.3 || theEta<-1.3) continue;
			for(int jp4=1; jp4<=h2S_4[ic]->GetNbinsY(); jp4++)
			{
				sumS_4[ic] += h2S_4[ic]->GetBinContent(ie4,jp4);
				sumB_4[ic] += h2B_4[ic]->GetBinContent(ie4,jp4);
			}
		}
		h2S_4[ic]->Divide(constant,sumS_4[ic]/sumB_4[ic]);
		h2S_4[ic]->Divide(h2B_4[ic]);

		for(int ie5=1; ie5<=h2S_5[ic]->GetNbinsX(); ie5++)
		{
			theEta = h2S_5[ic]->GetXaxis()->GetBinCenter(ie5);
			if(theEta>1.3 || theEta<-1.3) continue;
			for(int jp5=1; jp5<=h2S_5[ic]->GetNbinsY(); jp5++)
			{
				sumS_5[ic] += h2S_5[ic]->GetBinContent(ie5,jp5);
				sumB_5[ic] += h2B_5[ic]->GetBinContent(ie5,jp5);
			}
		}
		h2S_5[ic]->Divide(constant,sumS_5[ic]/sumB_5[ic]);
		h2S_5[ic]->Divide(h2B_5[ic]);

	}// deadline for ic loop 

	//=====================================================
	//processing for data 
	// experiment data from Alice sqrt{7} Tev |delta-eta|<1.3 .From arxiv: 1612.08975v1 
	//N1,fig.5 
	for(int ex=0;ex<2;ex++)
	{
		TFile *f = new TFile("./experiment_data/HEPData-ins1507157-v1-root.root", "Read");
		TDirectory *dir;
		dir= (TDirectory *)f->Get("Table 6");
		sprintf(finName, "Graph1D_%s", mul_y[ex]);
		graph[ex] = (TGraphAsymmErrors*)dir->Get(finName);
		Double_t *graph_x = graph[ex]->GetX();
		Double_t *graph_y = graph[ex]->GetY();
		grae[ex] = new TGraphAsymmErrors(29);
		for(Int_t i=0; i<graph[ex]->GetN();i++)
		{
			graph_exl[ex][i] = graph[ex]->GetErrorXlow(i);
			graph_exh[ex][i] = graph[ex]->GetErrorXhigh(i);
			graph_eyl[ex][i] = graph[ex]->GetErrorYlow(i);
			graph_eyh[ex][i] = graph[ex]->GetErrorYhigh(i);
			sprintf(title_y,"Graph1D_%s",tit[ex]);
			grae[ex]->SetPoint(i,graph_x[i],graph_y[i]);
			grae[ex]->SetPointError(i,graph_exl[ex][i],graph_exh[ex][i], graph_eyl[ex][i], graph_eyh[ex][i]);
			grae[ex]->SetMarkerStyle(20); 
			grae[ex]->SetMarkerSize(1.2); 
			grae[ex]->SetMarkerColor(1);
			grae[ex]->SetLineWidth(2);
		}
	}
	//end of processing for data
	//==========================================================
	TH3F *back = new TH3F("back","back",100,-1.6,1.6,200,-1.5,2.*PI-1.5,10,0,1.5);
	back->GetXaxis()->SetNdivisions(505);
	back->GetYaxis()->SetNdivisions(506);
	back->GetZaxis()->SetNdivisions(506);
	back->GetXaxis()->CenterTitle(true);
	back->GetYaxis()->CenterTitle(true);
	back->GetZaxis()->CenterTitle(true);
	back->GetXaxis()->SetTitleOffset(1.3);
	back->GetYaxis()->SetTitleOffset(1.3);
	back->GetZaxis()->SetTitleOffset(1.3);
	back->GetXaxis()->SetLabelOffset(0.015);
	back->GetYaxis()->SetLabelOffset(0.015);
	back->GetZaxis()->SetLabelOffset(0.015);
	back->GetXaxis()->SetTitle("#Delta#eta");
	back->GetYaxis()->SetTitle("#Delta#phi");
	back->GetZaxis()->SetTitle("C(#Delta#eta,#Delta#phi)");
	back->GetXaxis()->SetTitleFont(132);
	back->GetYaxis()->SetTitleFont(132);
	back->GetZaxis()->SetTitleFont(132);
	back->GetXaxis()->SetLabelSize(0.07);
	back->GetYaxis()->SetLabelSize(0.07);
	back->GetZaxis()->SetLabelSize(0.07);
	back->GetXaxis()->SetTitleSize(0.07);
	back->GetYaxis()->SetTitleSize(0.07);
	back->GetZaxis()->SetTitleSize(0.07);

	TH2F *back_flow = new TH2F("back_flow","back_flow",200,-1.5,2.*PI-1.5,11,0.70,1.48);
	back_flow->GetXaxis()->SetNdivisions(506);
	back_flow->GetYaxis()->SetNdivisions(506);
	back_flow->GetXaxis()->CenterTitle(true);
	back_flow->GetYaxis()->CenterTitle(true);
	back_flow->GetXaxis()->SetTitleOffset(0.95);
	back_flow->GetYaxis()->SetTitleOffset(0.95);// this represents the division between title and label. 
	back_flow->GetXaxis()->SetLabelOffset(0.01);
	back_flow->GetYaxis()->SetLabelOffset(0.01);
	back_flow->GetXaxis()->SetTitle("#Delta#phi");
	back_flow->GetYaxis()->SetTitle("C(#Delta#phi)");
	back_flow->GetXaxis()->SetTitleFont(132);
	back_flow->GetYaxis()->SetTitleFont(132);
	back_flow->GetXaxis()->SetLabelSize(0.08);
	back_flow->GetYaxis()->SetLabelSize(0.08);
	back_flow->GetXaxis()->SetTitleSize(0.08);
	back_flow->GetYaxis()->SetTitleSize(0.08);

	TLegend *leg[3];
	for(int iL=0; iL<3; iL++)
	{
		leg[iL] = new TLegend(0.4, 0.9, 0.7, 0.97,NULL,"brNDC");
		leg[iL]->SetName("leg[iL]");
		leg[iL]->SetBorderSize(0);
		leg[iL]->SetTextSize(0.09);
		leg[iL]->SetLineColor(0);
		leg[iL]->SetLineStyle(0);
		leg[iL]->SetLineWidth(0);
		leg[iL]->SetFillColor(0);
		leg[iL]->SetFillStyle(0);
		//leg[iL]->AddEntry(graph_v2[0][0][0],"v_{2}{PP}","p");
	}
	//"nofluc","Chain","Triangle"
	leg[0]->SetHeader(title_name[title_number]);
	leg[1]->SetHeader("(b) Chain");
	leg[2]->SetHeader("(c) Triangle");

	TLegend *leg_flow[3];
	for(int iL=0; iL<3; iL++)
	{
		if(iL==0) leg_flow[iL] = new TLegend(0.15, 0.65, 0.6, 0.87,NULL,"brNDC");
		else leg_flow[iL] = new TLegend(0.2, 0.93, 0.5, 0.95,NULL,"brNDC");
		leg_flow[iL]->SetName("leg_flow[iL]");
		leg_flow[iL]->SetBorderSize(0);
		leg_flow[iL]->SetTextSize(0.05);//from 0.07
		leg_flow[iL]->SetLineColor(0);
		leg_flow[iL]->SetLineStyle(0);
		leg_flow[iL]->SetLineWidth(0);
		leg_flow[iL]->SetFillColor(0);
		leg_flow[iL]->SetFillStyle(0);
		//leg_flow[iL]->AddEntry(graph_v2[0][0][0],"v_{2}{PP}","p");
	}
	//"nofluc","Chain","Triangle"
	leg_flow[0]->AddEntry(grae[0],"Alice: pp #sqrt{s}=7TeV,|#Delta#eta|<1.3","pE"); 
	leg_flow[0]->AddEntry(h1S_1_2[0],"AMPT-Melting, t_{H}=20(fm/c)","pE"); 
	leg_flow[0]->AddEntry(h1S_1_4[0],"AMPT-Melting, newCoal. t_{H}=20(fm/c)","pE"); 
	//leg_flow[0]->AddEntry(h1S_1_4[0],"AMPT-Default, t(hadron)=20(fm/c)","pE"); 
	//leg_flow[0]->AddEntry(h1S_1_5[0],"Including #Lambda, #Sigma^{0} weak decay","pE"); 

	TLatex *tex = new TLatex( -0.5,0.75," p_{T} < 1.0 GeV/c(AMPT)"); 
	tex->SetTextSize(0.05); 
	tex->SetLineWidth(0); 
	tex->SetTextFont(42); 
	tex->SetTextAngle(0); 
	//TLatex *textit = new TLatex( 2.0,0.75,"(a1) #pi^{+}-#pi^{+}"); 
	//TLatex *textit = new TLatex( 2.0,0.75,"(a2) #pi^{-}-#pi^{-}"); 
	//TLatex *textit = new TLatex( 2.0,0.75,"(b1) K^{+}-K^{+}"); 
	//TLatex *textit = new TLatex( 2.0,0.75,"(b2) K^{-}-K^{-}"); 
	//TLatex *textit = new TLatex( 3.0,0.85,"(a) p-p"); 
	//TLatex *textit = new TLatex( 3.0,0.85,"(b) #bar{p}-#bar{p}"); 
	//TLatex *textit = new TLatex( 3.0,0.85,"(c) #Lambda-#Lambda"); 
	TLatex *textit = new TLatex( 3.0,0.85,"(d) #bar{#Lambda}-#bar{#Lambda}");
	textit->SetTextSize(0.06);
	textit->SetLineWidth(0);
	textit->SetTextFont(42);
	textit->SetTextAngle(0);


	TCanvas *c1 = new TCanvas("c1", "c1",50,0,550,550);
	//c1->SetLogx(1);
	//c1->SetLogy(1);
	//  c1->SetGrid(1,1);
	//  c1->SetFillStyle(4000);
	//  c1->SetFrameFillStyle(4000);
	c1->SetBorderMode(0);
	c1->SetFillColor(10);
	c1->SetFrameFillColor(0);
	c1->SetFrameBorderMode(0);
	c1->SetLeftMargin(0.15);//0.15
	c1->SetTopMargin(0.05);
	c1->SetRightMargin(0.02);
	c1->SetBottomMargin(0.1);
	//c1->Divide(2,2,0.01,0.01);
	//c1->Divide(2,1,0.0,0.0);

	//	//==============================================================
	Double_t padX[2] = {0.0,1.0};
	Double_t padY[3] = {0.0,1.0};

	TPad *pad[1][2];
	char padName[32];
	for(int i=0; i<1; i++)//the pad 3*2   i=3->1 ,j=2->2
	{
		for(int j=0; j<2; j++)
		{
			sprintf(padName,"pad%d_%d",i,j);

			pad[i][j] = new TPad(padName,padName,padX[i],padY[j],padX[i+1],padY[j+1]);

			pad[i][j]->Draw();
			pad[i][j]->cd();
			pad[i][j]->Range(0,0,1.0,1.0);
			//    pad[i]->SetTickx(1);
			//    pad[i]->SetTicky(1);
			//pad[i]->SetLogx(1);

			if(j==1) pad[i][j]->SetBottomMargin(0.15);//this is the division of the two pads. 
			else pad[i][j]->SetBottomMargin(0.30);//this is the bottom margin of second pad. 
			pad[i][j]->SetLeftMargin(0.16);//0.18
			pad[i][j]->SetRightMargin(0.02);
			pad[i][j]->SetTopMargin(0.1);

			if(j==1)
			{
				//back->Draw("c");
				//h2S[i]->Draw("SURF1SAME");
				back_flow->Draw();
				//leg[0]->Draw("same");
				//leg_flow[i]->Draw("same");
				//tex->Draw("same"); 
				textit->Draw("same"); 
				//h1S_1_1[i]->Draw("same");
				h1S_1_2[i]->Draw("same");
				//h1S_1_3[i]->Draw("same");
				h1S_1_4[i]->Draw("same");
				//h1S_1_5[i]->Draw("same");
				grae[0]->Draw("Epsame");
				//fCDPhi[i]->Draw("same");
				//for(int nn=0; nn<3; nn++) fCDPhi_nn[i][nn]->Draw("same");
			}
			pad[i][j]->Modified();
			c1->cd();
		}
	}
	char fin_pro[1000]; 
	sprintf(fin_pro,"./figs_PseudoRapidity_paper/figs_Mul_projection/figs6/figs_%s.eps",file_name[ip_number]); 
	c1->SaveAs(fin_pro);
}
