////* Caculate ACP
template<class th1>
double caculateACPerror( th1* Oh, bool murmur=false)
{
   int binOg0=2; 	
   int binOl0=1; 	
  
   double Og0 = Oh->GetBinContent(binOg0);
   double Ol0 = Oh->GetBinContent(binOl0);
  
   double n3  = (Og0+Ol0)*(Og0+Ol0)*(Og0+Ol0);
   double ACPe = 2*sqrt(Og0*Ol0/n3);

   if ( murmur ) printf("O<0 : %.1f, O>0 : %.1f\n", Ol0, Og0);
   printf("Call caculateACPerror: %f\n", ACPe);	
   return ACPe;
}
template<class th1>
double caculateACPerrorWrt( th1* Oh, bool murmur=false)
{
   int binOg0=2; 	
   int binOl0=1; 	
  
   double Og0 = Oh->GetBinContent(binOg0);
   double Ol0 = Oh->GetBinContent(binOl0);
   double Og0e = Oh->GetBinError(binOg0);
   double Ol0e = Oh->GetBinError(binOl0);
   
   double t1 = 2*Ol0/((Ol0+Og0)*(Ol0+Og0));  
   double t2 = 2*Og0/((Ol0+Og0)*(Ol0+Og0)); 
   double ACPe = sqrt( t1*t1*Og0e*Og0e + t2*t2*Ol0e*Ol0e ); 

   if ( murmur ) printf("O<0 : %.1f, O>0 : %.1f\n", Ol0, Og0);
   printf("Call caculateACPerror: %f\n", ACPe);	
   return ACPe;
}
template<class th1>
double caculateACP( th1* Oh, bool murmur=false)
{
   int binOg0=2; 	
   int binOl0=1; 	
  
   double Og0 = Oh->GetBinContent(binOg0);
   double Ol0 = Oh->GetBinContent(binOl0);

   double ACP = (Og0-Ol0)/(Og0+Ol0);
   if ( murmur ) printf("O<0 : %.1f, O>0 : %.1f\n", Ol0, Og0);
   printf("Call caculateACP: %f\n", ACP);	
   return ACP;	
}

template<class th1>
void caculateACPDetail( th1* h, std::string Oname="O" )
{
   int binOg0=2; 	
   int binOl0=1; 	
  
   double Og0 = h->GetBinContent(binOg0);
   double Ol0 = h->GetBinContent(binOl0);

   double ACP = caculateACP(h);
   double ACPe = caculateACPerrorWrt(h);
   //double ACPe = caculateACPerror(h);
  
   printf("%s>0=%5.0f, %s<0=%5.0f, ACP=%6.3f +/-%6.3f \n", Oname.c_str(), Og0, Oname.c_str(), Ol0, ACP, ACPe);
}
void caculateACPDetail( TFile* f, std::string histName, std::string Oname="O")
 {
   TH1D* Oh = (TH1D*)f->Get(histName.c_str());
   printf("%-15s: ", histName.c_str());	
   caculateACPDetail(Oh, Oname);
}
