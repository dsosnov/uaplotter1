{
  double x[11] = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5};
  double dx[11];
  for (short unsigned int i=0; i<11; i++){
    dx[11] = .5;
  };
  double ym1[11] = {720,175,93,85,64,61,85,39,165,296,2317};
  double yp1[11] = {744,257,172,249,236,261,326,187,707,692,1516};
  double ym2[11] = {720,175,94,85,64,61,85,44,164,297,2311};
  double yp2[11] = {744,257,168,247,235,261,325,192,705,691, 1522};
  double ym3[11] = {834,298,76,78,60,60,82,42,349,605,1920};
  double yp3[11] = {908,376,126,223,220,245,308,182,845,625,1320};

  TGraphErrors * gm1 = new TGraphErrors(11, x,ym1, dx,0);
  gm1->SetLineColor(1); gm1->SetMarkerStyle(20);
  TGraphErrors * gm2 = new TGraphErrors(11, x,ym2, dx,0);
  gm2->SetLineColor(4);
  TGraphErrors * gm3 = new TGraphErrors(11, x,ym3, dx,0);
  gm3->SetLineColor(2);
  
  TGraphErrors * gp1 = new TGraphErrors(11, x,yp1, dx,0);
  gp1->SetLineColor(1);gp1->SetMarkerStyle(20);
  TGraphErrors * gp2 = new TGraphErrors(11, x,yp2, dx,0);
  gp2->SetLineColor(4);
  TGraphErrors * gp3 = new TGraphErrors(11, x,yp3, dx,0);
  gp3->SetLineColor(2);

  TCanvas * c1 = new TCanvas("c1","c1", 1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  gm1->Draw("APL"); gm2->Draw("PLsames"); gm3->Draw("PLsames"); 
  c1->cd(2);
  gp1->Draw("APL"); gp2->Draw("PLsames"); gp3->Draw("PLsames"); 
}


