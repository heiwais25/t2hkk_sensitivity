#include "BargerPropagator.h"

class ProbWrapper {

    public:

        ProbWrapper(double Delta_M2, double theta23,  double delta_m2, double theta13, double theta12, double Delta);
        ProbWrapper(double Delta_M2  = 2.4e-3, double theta23 = 0.5);
        ProbWrapper(bool isEarthProfile, const char * f);
        virtual ~ProbWrapper();

        // Methods to get the oscillation probabililties
        
        
        double GetProbNuMuNuMu(float energy);
        double GetProbNuMuNuE(float energy);
        double GetProbNuMuBarNuMuBar(float energy);
        double GetProbNuENuE(float energy);
        double GetProbNuEBarNuEBar(float energy);
        double GetProbNuMuBarNuEBar(float energy);

        double GetProbNuMuNuMu(float energy, float cosineZ);
        double GetProbNuMuNuE(float energy, float cosineZ);
        double GetProbNuMuBarNuMuBar(float energy, float cosineZ);
        double GetProbNuENuE(float energy, float cosineZ);
        double GetProbNuEBarNuEBar(float energy, float cosineZ);
        double GetProbNuMuBarNuEBar(float energy, float cosineZ);

        double GetProb(int initialSpecies, int finalSpecies);

        // Method to set whether we will use Earth profile
        // void SetEarthProfile(bool EarthProfile = false){isEarthProfile = EarthProfile};

        // Methods to change the oscillation parameters
        void SetAntiNeutrino(bool isAntiNeutrino); 
        void SetInverseMassHierarchy(bool isInverted); 
        void SetDeltaCP(double angle); 
        void SetOscillationParameters(double Delta_M2  = 2.4e-3, double theta23 = 0.5,  double delta_m2 = 7.6e-5, double theta13 = 0.0241, double theta12 = 0.307, double Delta = 0.0);
        void SetBaseLine(double baseline);
        void SetDensity(double density);

        /// Density and path length
        double BasePath;
        double Density;
                                                                                                
        /// Oscillation Parameters
        bool kSquared;   // using sin^2(x) variables?
        // bool isEarthProfile;
        int    kNuBar;
        double DM2;
        double Theta23;
        double Theta13;
        double dm2;
        double Theta12;
        double delta;

        BargerPropagator *bNu;

};

void ProbWrapper::SetOscillationParameters(double Delta_M2, double theta23, double delta_m2, double theta13, double theta12, double Delta)
{
    DM2     =  Delta_M2;
    Theta23 =  theta23;
    dm2     =  delta_m2;
    Theta13 =  theta13;
    Theta12 =  theta12;
    delta   =  Delta;
}

void ProbWrapper::SetAntiNeutrino(bool isAntiNeutrino)
{
    if (isAntiNeutrino) kNuBar = -1;
    else kNuBar = 1;
}

void ProbWrapper::SetInverseMassHierarchy(bool isInverted)
{
    if (isInverted) DM2 = fabs(DM2)*(-1.0);
    else DM2 = fabs(DM2);
}

// Set Delta CP in units of degrees
void ProbWrapper::SetDeltaCP(double angle)
{
    delta   =  angle * (3.1415926/180.0);
}

void ProbWrapper::SetBaseLine(double baseline){
    BasePath = baseline;
}

void ProbWrapper::SetDensity(double density){
    Density = density;
}


ProbWrapper::ProbWrapper(double Delta_M2, double theta23)
{

    /// Density and path length
    BasePath = 1100.0;
    Density = 2.60;

    /// Oscillation Parameters
    kSquared  = true;   // using sin^2(x) variables?

    kNuBar  =  1;
    DM2     =  Delta_M2;
    Theta23 =  theta23;
    Theta13 =  0.0241;
    dm2     =  7.6e-5;
    Theta12 =  0.307;
    delta   =  0.0;

    std::cout << "Using          " << std::endl
        << "      DM2      " <<  DM2      << std::endl
        << "      Theta23  " <<  Theta23  << std::endl
        << "      Theta13  " <<  Theta13  << std::endl
        << "      dm2      " <<  dm2      << std::endl
        << "      Theta12  " <<  Theta12  << std::endl;

    bNu = new BargerPropagator();
    bNu->UseMassEigenstates( false );

}

ProbWrapper::ProbWrapper(double Delta_M2, double theta23, double delta_m2, double theta13, double theta12, double Delta)
{

    /// Density and path length
    BasePath = 1100.0;
    Density = 2.60;

    /// Oscillation Parameters
    kSquared  = true;   // using sin^2(x) variables?
    kNuBar  =  1;

    SetOscillationParameters(Delta_M2, theta23, delta_m2, theta13, theta12, Delta);

    std::cout << "Using          " << std::endl
        << "      DM2      " <<  DM2      << std::endl
        << "      Theta23  " <<  Theta23  << std::endl
        << "      dm2      " <<  dm2      << std::endl
        << "      Theta13  " <<  Theta13  << std::endl
        << "      Theta12  " <<  Theta12  << std::endl;

    bNu = new BargerPropagator();
    bNu->UseMassEigenstates( false );

}

ProbWrapper::ProbWrapper(bool isEarthProfile, const char * f)
{
    /// Density and path length
    double Delta_M2  = 2.4e-3; 
    double theta23 = 0.5;

    BasePath = 1100.0;
    Density = 2.60;
    /// Oscillation Parameters
    kSquared  = true;   // using sin^2(x) variables?

    kNuBar  =  1;
    DM2     =  Delta_M2;
    Theta23 =  theta23;
    Theta13 =  0.0241;
    dm2     =  7.6e-5;
    Theta12 =  0.307;
    delta   =  0.0;

    std::cout << "Using          " << std::endl
    << "      DM2      " <<  DM2      << std::endl
    << "      Theta23  " <<  Theta23  << std::endl
    << "      Theta13  " <<  Theta13  << std::endl
    << "      dm2      " <<  dm2      << std::endl
    << "      Theta12  " <<  Theta12  << std::endl;

    if(isEarthProfile == true)
    {
        bNu = new BargerPropagator(f);
    }
    else
    {
        bNu = new BargerPropagator();
    }
    bNu->UseMassEigenstates( false );
}

ProbWrapper::~ProbWrapper()
{
    if(bNu!=NULL){
        delete bNu;
        bNu = NULL;
    }
}

double ProbWrapper::GetProb(int initialSpecies, int finalSpecies)
{
    double total_prob = 0.0;
    for(int m=1; m<=3; m++) total_prob += bNu->GetProb(initialSpecies, m); // Normalize the Probabilities //
    if ( total_prob >1.00001 || total_prob<0.99998 ) {std::cerr << "ERROR Prob" << std::endl;}// abort();}

    // Get Muon appearance probability
    double prob = bNu->GetProb(initialSpecies,finalSpecies);     

    return prob;
}

double ProbWrapper::GetProbNuMuNuMu(float energy)
{
    SetAntiNeutrino(0);
    bNu->SetMNS( Theta12,  Theta13, Theta23, dm2, DM2, delta , energy, kSquared, kNuBar ); 
    bNu->propagateLinear( 1*kNuBar, BasePath, Density );

    return GetProb(2,2);
}

double ProbWrapper::GetProbNuMuBarNuMuBar(float energy)
{
    SetAntiNeutrino(1);
    bNu->SetMNS( Theta12,  Theta13, Theta23, dm2, DM2, delta , energy, kSquared, kNuBar ); 
    bNu->propagateLinear( 1*kNuBar, BasePath, Density );

    return GetProb(2,2);
}

double ProbWrapper::GetProbNuMuNuE(float energy)
{
    SetAntiNeutrino(0);
    bNu->SetMNS( Theta12,  Theta13, Theta23, dm2, DM2, delta , energy, kSquared, kNuBar ); 
    bNu->propagateLinear( 1*kNuBar, BasePath, Density );

    return GetProb(2,1);
}

double ProbWrapper::GetProbNuMuBarNuEBar(float energy)
{
    SetAntiNeutrino(1);
    bNu->SetMNS( Theta12,  Theta13, Theta23, dm2, DM2, delta , energy, kSquared, kNuBar ); 
    bNu->propagateLinear( 1*kNuBar, BasePath, Density );

    return GetProb(2,1);
}

double ProbWrapper::GetProbNuENuE(float energy)
{
    SetAntiNeutrino(0);
    bNu->SetMNS( Theta12,  Theta13, Theta23, dm2, DM2, delta , energy, kSquared, kNuBar ); 
    bNu->propagateLinear( 1*kNuBar, BasePath, Density );

    return GetProb(1,1);
}

double ProbWrapper::GetProbNuEBarNuEBar(float energy)
{
    SetAntiNeutrino(1);
    bNu->SetMNS( Theta12,  Theta13, Theta23, dm2, DM2, delta , energy, kSquared, kNuBar ); 
    bNu->propagateLinear( 1*kNuBar, BasePath, Density );

    return GetProb(1,1);
}

/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
    Atmospheric neutrino oscillation probability(using propagate routine in the Bargerpropagate)

 -------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
double ProbWrapper::GetProbNuMuNuMu(float energy, float cosineZ)
{
    SetAntiNeutrino(0);
    bNu->SetMNS( Theta12,  Theta13, Theta23, dm2, DM2, delta , energy, kSquared, kNuBar ); 
    bNu->DefinePath(cosineZ, 25.00 );
    bNu->propagate( 1*kNuBar );

    return GetProb(2,2);
}

double ProbWrapper::GetProbNuMuBarNuMuBar(float energy, float cosineZ)
{
    SetAntiNeutrino(1);
    bNu->SetMNS( Theta12,  Theta13, Theta23, dm2, DM2, delta , energy, kSquared, kNuBar ); 
    bNu->DefinePath( cosineZ, 25.00  );
    bNu->propagate( 1*kNuBar );

    return GetProb(2,2);
}

double ProbWrapper::GetProbNuMuNuE(float energy, float cosineZ)
{
    SetAntiNeutrino(0);
    bNu->SetMNS( Theta12,  Theta13, Theta23, dm2, DM2, delta , energy, kSquared, kNuBar ); 
    bNu->DefinePath( cosineZ, 25.00  );
    bNu->propagate( 1*kNuBar );

    return GetProb(2,1);
}

double ProbWrapper::GetProbNuMuBarNuEBar(float energy, float cosineZ)
{
    SetAntiNeutrino(1);
    bNu->SetMNS( Theta12,  Theta13, Theta23, dm2, DM2, delta , energy, kSquared, kNuBar ); 
    bNu->DefinePath( cosineZ, 25.00  );
    bNu->propagate( 1*kNuBar );

    return GetProb(2,1);
}

double ProbWrapper::GetProbNuENuE(float energy, float cosineZ)
{
    SetAntiNeutrino(0);
    bNu->SetMNS( Theta12,  Theta13, Theta23, dm2, DM2, delta , energy, kSquared, kNuBar ); 
    bNu->DefinePath( cosineZ, 25.00  );
    bNu->propagate( 1*kNuBar );

    return GetProb(1,1);
}

double ProbWrapper::GetProbNuEBarNuEBar(float energy, float cosineZ)
{
    SetAntiNeutrino(1);
    bNu->SetMNS( Theta12,  Theta13, Theta23, dm2, DM2, delta , energy, kSquared, kNuBar ); 
    bNu->DefinePath( cosineZ, 25.00  );
    bNu->propagate( 1*kNuBar );

    return GetProb(1,1);
}

















