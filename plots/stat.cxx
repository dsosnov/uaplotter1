{

  int PbpHsel =  157185;
  int PbpHstat = 0;
  int PbpHplus[11] = {3074, 560, 420, 502, 377, 389, 269, 64, 73, 29, 0};
  int PbpEsel = 107484;
  int PbpEstat = 0;
  int PbpEplus[11] = {2454, 448, 319, 355, 204, 158, 126, 71, 189, 238, 0};
  int PbpDsel= 467874;
  int PbpDstat = 0;
  int PbpDplus[11] = {4374, 1180, 710, 784, 575, 525, 526, 304, 1681, 3119};

  for(int i=0; i<11; i++){
    PbpHstat+=PbpHplus[i];
    PbpEstat+=PbpEplus[i];
    PbpDstat+=PbpDplus[i];
  };
  
  TH1F * Hijing = new TH1F("Hijing", "Hijing", 12, -1, 11);
  TH1F * EPOS = new TH1F("EPOS", "EPOS", 12, -1, 11);
  TH1F * data = new TH1F("data", "data", 12, -1, 11);
  for(int i=0; i<11; i++){
	Hijing->SetBinContent(i+2, float(PbpHplus[i])/PbpHstat);
	Hijing->SetBinError  (i+2, sqrt(float(PbpHplus[i]))/PbpHstat);
	EPOS->SetBinContent  (i+2, float(PbpEplus[i])/PbpEstat);
	EPOS->SetBinError    (i+2, sqrt(float(PbpEplus[i]))/PbpEstat);
	data->SetBinContent  (i+2, float(PbpDplus[i])/PbpDstat);
	data->SetBinError    (i+2, sqrt(float(PbpDplus[i]))/PbpDstat);
  };
  std::cout << "Hjing " << float(PbpHstat)/PbpHsel << "\t" << sqrt(PbpHstat)/PbpHsel << std::endl;
  std::cout << "EPOS  " << float(PbpEstat)/PbpEsel << "\t" << sqrt(PbpEstat)/PbpEsel << std::endl;
  std::cout << "data  " << float(PbpDstat)/PbpDsel << "\t" << sqrt(PbpDstat)/PbpDsel << std::endl;
}
