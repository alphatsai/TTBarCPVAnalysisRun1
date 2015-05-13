////* Caculate ACP
template<class th1>
double caculateACPerror( th1* Oh)
{
   int binOg0=2; 	
   int binOl0=1; 	
  
   double Og0 = Oh->GetBinContent(binOg0);
   double Ol0 = Oh->GetBinContent(binOl0);
  
   double n3  = (Og0+Ol0)*(Og0+Ol0)*(Og0+Ol0);
   double ACPe = 2*sqrt(Og0*Ol0/n3);

   printf("Call caculateACPerror: %f\n", ACPe);	
   return ACPe;
}
template<class th1>
double caculateACP( th1* Oh)
{
   int binOg0=2; 	
   int binOl0=1; 	
  
   double Og0 = Oh->GetBinContent(binOg0);
   double Ol0 = Oh->GetBinContent(binOl0);

   double ACP = (Og0-Ol0)/(Og0+Ol0);
   printf("Call caculateACP: %f\n", ACP);	
   return ACP;	
}
void caculateACP( TFile* f, std::string histName, std::string Oname="O")
{
   int binOg0=2; 	
   int binOl0=1; 	
  
   TH1D* Oh = (TH1D*)f->Get(histName.c_str());

   double Og0 = Oh->GetBinContent(binOg0);
   double Ol0 = Oh->GetBinContent(binOl0);

   double ACP = caculateACP(Oh);
   double ACPe = caculateACPerror(Oh);
  
   printf("%-15s: %s>0=%5.0f, %s<0=%5.0f, ACP=%6.3f +/-%6.3f \n", histName.c_str(), Oname.c_str(), Og0, Oname.c_str(), Ol0, ACP, ACPe);

}
