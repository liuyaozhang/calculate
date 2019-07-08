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

  gStyle->SetPalette(1);

  char caseC12[3][32] = {"nofluc","Chain","Triangle"};

  Double_t PI=TMath::Pi();
  TH2D *h2S[3];
  TH2D *h2B[3];

  char h2SBName[64];
  for(int i=0; i<3; i++)
  {
    sprintf(h2SBName,"h2S_%d",i);
    h2S[i] = new TH2D(h2SBName,h2SBName,200,-5,5,200,-1.5,2.*PI-1.5);
    sprintf(h2SBName,"h2B_%d",i);
    h2B[i] = new TH2D(h2SBName,h2SBName,200,-5,5,200,-1.5,2.*PI-1.5);
  }

  TFile *fin[3];
  TH2D *htmp;

  char finName[512];
  char histName[128];
  Double_t theEta=0;
  Double_t sumS[3]={0};
  Double_t sumB[3]={0};
  Double_t sumS_flow[3] = {0};
  Double_t sumB_flow[3] = {0};

  TF1 *constant = new TF1("constant","1.",-10,10);
  TH1D *h1S_1[3];
  TH1D *h1S_2[3];
  TH1D *h1B_1[3];
  TH1D *h1B_2[3];
  Double_t etaCut[4] = {-5.,-2.,2.,5.};//-1,1 -> -2,2
  Int_t etaBin[4];

  Double_t vnn[3][6]={0};
  //Double_t vnnErr[6]={0};
  Double_t Dphi=0, CDphi=0, vnnC=0;

  TF1 *fCDPhi[3];
  TF1 *fCDPhi_nn[3][3];//c12Case, nn=1,2,3
  char f1Name[64];

  for(int ic=0; ic<3; ic++)
  {
    sprintf(finName,"../results/result-minibias/pT1-3GeV/result-2PC-etaGap-c12-Au197-200GeV-%s/result-all.root",caseC12[ic]);
    fin[ic] = new TFile(finName,"READ");

    for(int ih=0; ih<10; ih++)
    {
	sprintf(histName,"Sh2_Deta_Dphi_%d",20+ih);
	htmp = (TH2D *) fin[ic]->Get(histName);
	h2S[ic]->Add(htmp);

	sprintf(histName,"Bh2_Deta_Dphi_%d",20+ih);
	htmp = (TH2D *) fin[ic]->Get(histName);
	h2B[ic]->Add(htmp);
    }
    //fin.Close();

    for(int i=1; i<=h2S[ic]->GetNbinsX(); i++)
    {
	theEta = h2S[ic]->GetXaxis()->GetBinCenter(i);
	if(theEta>-2.&&theEta<2.) continue;
	for(int j=1; j<=h2S[ic]->GetNbinsY(); j++)
	{
	  sumS_flow[ic] += h2S[ic]->GetBinContent(i,j);
	  sumB_flow[ic] += h2B[ic]->GetBinContent(i,j);
	}
    }
    for(int ie=0; ie<4; ie++) etaBin[ie] = h2S[0]->GetXaxis()->FindBin(etaCut[ie]);
    h1S_1[ic] = (TH1D *) h2S[ic]->ProjectionY("h1S_1[ic]",etaBin[0],etaBin[1]-1,"");
    h1S_2[ic] = (TH1D *) h2S[ic]->ProjectionY("h1S_2[ic]",etaBin[2]-1,etaBin[3]-1,"");
    h1S_1[ic]->Add(h1S_2[ic]);
    h1B_1[ic] = (TH1D *) h2B[ic]->ProjectionY("h1B_1[ic]",etaBin[0],etaBin[1]-1,"");
    h1B_2[ic] = (TH1D *) h2B[ic]->ProjectionY("h1B_2[ic]",etaBin[2]-1,etaBin[3]-1,"");
    h1B_1[ic]->Add(h1B_2[ic]);

    h1S_1[ic]->Sumw2();
    h1B_1[ic]->Sumw2();
    h1S_1[ic]->Divide(constant,sumS_flow[ic]/sumB_flow[ic]);
    h1S_1[ic]->Divide(h1B_1[ic]);

    vnnC=0;
    for(int ip=1; ip<=h1S_1[ic]->GetNbinsX(); ip++)
    {
	Dphi = h1S_1[ic]->GetBinCenter(ip);
	CDphi = h1S_1[ic]->GetBinContent(ip);
	for(int nn=1; nn<=6; nn++)
	{
	  vnn[ic][nn-1] += cos(Double_t(nn)*Dphi)*CDphi;
	}
	vnnC   += CDphi;
    }

    for(int nn=1; nn<=6; nn++) vnn[ic][nn-1] /= vnnC;

    sprintf(f1Name,"fCDphi_%d",ic);
    fCDPhi[ic] = new TF1(f1Name,"[0]*(1.+2.*([1]*cos(x)+[2]*cos(2.*x)+[3]*cos(3.*x)+[4]*cos(4.*x)+[5]*cos(5.*x)+[6]*cos(6.*x)))",-10,10);
    fCDPhi[ic]->SetLineStyle(4);
    fCDPhi[ic]->SetLineColor(kRed);
    fCDPhi[ic]->SetParameter(0,1.);
    for(int nn=1; nn<=3; nn++) fCDPhi[ic]->SetParameter(nn,vnn[ic][nn-1]);

    sprintf(f1Name,"fCDphi_%d_0",ic);
    fCDPhi_nn[ic][0] = new TF1(f1Name,"1.+2.*[0]*cos(x)",-10,10);
    fCDPhi_nn[ic][0]->SetParameter(0,vnn[ic][0]);
    fCDPhi_nn[ic][1] = new TF1(f1Name,"1.+2.*[0]*cos(2.*x)",-10,10);
    fCDPhi_nn[ic][1]->SetParameter(0,vnn[ic][1]);
    fCDPhi_nn[ic][2] = new TF1(f1Name,"1.+2.*[0]*cos(3.*x)",-10,10);
    fCDPhi_nn[ic][2]->SetParameter(0,vnn[ic][2]);
    fCDPhi_nn[ic][0]->SetParameter(0,vnn[ic][0]);

    fCDPhi_nn[ic][0]->SetLineStyle(1);
    fCDPhi_nn[ic][1]->SetLineStyle(2);
    fCDPhi_nn[ic][2]->SetLineStyle(3);

    fCDPhi_nn[ic][0]->SetLineColor(kBlack);
    fCDPhi_nn[ic][1]->SetLineColor(kBlue);
    fCDPhi_nn[ic][2]->SetLineColor(kGreen+4);

    //for plot
    h1S_1[ic]->Rebin(8);
    h1S_1[ic]->Divide(constant,8);
    h1S_1[ic]->SetMarkerStyle(24);
    //

    //
    h2S[ic]->Rebin2D(8,8);
    h2B[ic]->Rebin2D(8,8);

    for(int ie=1; ie<=h2S[ic]->GetNbinsX(); ie++)
    {
	theEta = h2S[ic]->GetXaxis()->GetBinCenter(ie);
	if(theEta>-2.&&theEta<2.) continue;
	for(int jp=1; jp<=h2S[ic]->GetNbinsY(); jp++)
	{
	  sumS[ic] += h2S[ic]->GetBinContent(ie,jp);
	  sumB[ic] += h2B[ic]->GetBinContent(ie,jp);
	}
    }

    h2S[ic]->Divide(constant,sumS[ic]/sumB[ic]);
    h2S[ic]->Divide(h2B[ic]);
  } // for (ic;;)

  TH3F *back = new TH3F("back","back",200,-5,5,200,-1.5,2.*PI-1.5,10,0,1.5);
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

  TH2F *back_flow = new TH2F("back_flow","back_flow",200,-1.5,2.*PI-1.5,10,0.9,1.09);
  back_flow->GetXaxis()->SetNdivisions(505);
  back_flow->GetYaxis()->SetNdivisions(506);
  back_flow->GetXaxis()->CenterTitle(true);
  back_flow->GetYaxis()->CenterTitle(true);
  back_flow->GetXaxis()->SetTitleOffset(1.1);
  back_flow->GetYaxis()->SetTitleOffset(1.2);
  back_flow->GetXaxis()->SetLabelOffset(0.015);
  back_flow->GetYaxis()->SetLabelOffset(0.015);
  back_flow->GetXaxis()->SetTitle("#Delta#phi");
  back_flow->GetYaxis()->SetTitle("C(#Delta#phi)");
  back_flow->GetXaxis()->SetTitleFont(132);
  back_flow->GetYaxis()->SetTitleFont(132);
  back_flow->GetXaxis()->SetLabelSize(0.07);
  back_flow->GetYaxis()->SetLabelSize(0.07);
  back_flow->GetXaxis()->SetTitleSize(0.07);
  back_flow->GetYaxis()->SetTitleSize(0.07);

  TLegend *leg[3];
  for(int iL=0; iL<3; iL++)
  {
    leg[iL] = new TLegend(0., 0.93, 0.3, 0.99,NULL,"brNDC");
    leg[iL]->SetName("leg[iL]");
    leg[iL]->SetBorderSize(0);
    leg[iL]->SetTextSize(0.07);
    leg[iL]->SetLineColor(0);
    leg[iL]->SetLineStyle(0);
    leg[iL]->SetLineWidth(0);
    leg[iL]->SetFillColor(0);
    leg[iL]->SetFillStyle(0);
    //leg[iL]->AddEntry(graph_v2[0][0][0],"v_{2}{PP}","p");
  }
  //"nofluc","Chain","Triangle"
  leg[0]->SetHeader("(a) Woods-Saxon");
  leg[1]->SetHeader("(b) Chain");
  leg[2]->SetHeader("(c) Triangle");

  TLegend *leg_flow[3];
  for(int iL=0; iL<3; iL++)
  {
    if(iL==0) leg_flow[iL] = new TLegend(0.2, 0.7, 0.5, 0.95,NULL,"brNDC");
    else leg_flow[iL] = new TLegend(0.2, 0.93, 0.5, 0.95,NULL,"brNDC");
    leg_flow[iL]->SetName("leg_flow[iL]");
    leg_flow[iL]->SetBorderSize(0);
    leg_flow[iL]->SetTextSize(0.07);
    leg_flow[iL]->SetLineColor(0);
    leg_flow[iL]->SetLineStyle(0);
    leg_flow[iL]->SetLineWidth(0);
    leg_flow[iL]->SetFillColor(0);
    leg_flow[iL]->SetFillStyle(0);
    //leg_flow[iL]->AddEntry(graph_v2[0][0][0],"v_{2}{PP}","p");
  }
  //"nofluc","Chain","Triangle"
  leg_flow[0]->SetHeader("(d) Woods-Saxon");
  leg_flow[0]->AddEntry(fCDPhi_nn[0][0],"n=1","l");
  leg_flow[0]->AddEntry(fCDPhi_nn[0][1],"n=2","l");
  leg_flow[0]->AddEntry(fCDPhi_nn[0][2],"n=3","l");
  leg_flow[0]->AddEntry(fCDPhi[0],"sum","l");
  leg_flow[1]->SetHeader("(e) Chain");
  leg_flow[2]->SetHeader("(f) Triangle");


  TCanvas *c1 = new TCanvas("c1", "c1",50,0,700*3,500*2);
  //c1->SetLogx(1);
  //c1->SetLogy(1);
  //  c1->SetGrid(1,1);
  //  c1->SetFillStyle(4000);
  //  c1->SetFrameFillStyle(4000);
  c1->SetBorderMode(0);
  c1->SetFillColor(10);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);
  c1->SetLeftMargin(0.15);
  c1->SetTopMargin(0.05);
  c1->SetRightMargin(0.02);
  c1->SetBottomMargin(0.1);
  //c1->Divide(2,2,0.01,0.01);
  //c1->Divide(2,1,0.0,0.0);


  Double_t padX[4] = {0.0,1./3.,1./3.*2,1.};
  Double_t padY[3] = {0.0,1./2.,1.};

  TPad *pad[3][2];
  char padName[32];
  int xN=0;
  int yN=0;
  for(int i=0; i<3; i++)
  {
    for(int j=0; j<2; j++)
    {
	sprintf(padName,"pad%d_%d",i,j);

	pad[i][j] = new TPad(padName,padName,padX[i],padY[j],padX[i+1],padY[j+1]);

	pad[i][j]->Draw();
	pad[i][j]->cd();
	pad[i][j]->Range(0,0,1,1);
	//    pad[i]->SetTickx(1);
	//    pad[i]->SetTicky(1);
	//pad[i]->SetLogx(1);

	if(j==1) pad[i][j]->SetBottomMargin(0.1);
	else pad[i][j]->SetBottomMargin(0.16);
	pad[i][j]->SetLeftMargin(0.18);
	pad[i][j]->SetRightMargin(0.02);
	pad[i][j]->SetTopMargin(0.01);

	if(j==1)
	{
	  back->Draw("c");
	  leg[i]->Draw();
	  h2S[i]->Draw("SURF1SAME");// 2D correlation distribution 
	}
	else
	{
	  back_flow->Draw("c");
	  leg_flow[i]->Draw();
	  h1S_1[i]->Draw("EPsame");// Delta phi correlation funcion : C(Delta_phi):formual 15 , Physical Review C86,014907(2012)
	  fCDPhi[i]->Draw("same");
	  for(int nn=0; nn<3; nn++) fCDPhi_nn[i][nn]->Draw("same");
	}
	pad[i][j]->Modified();
	c1->cd();
    }
  }

}
