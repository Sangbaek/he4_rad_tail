//
// newed.cc
// Developer : Jingyi Zhou
// Based on ###
// History:
//   Feb 2023, J. Zhou, ed event generator, lowest order radiative correction.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TFoam.h"
#include "TRandom2.h"

#include "newhe4.hh"

#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <time.h>
#include <string>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char ** argv)
{
    // inputs
    char mychar[64];

    //std::cout << "Full energy of the incident lepton (MeV): " << std::flush;
    //std::cin.getline(mychar, 64);
    Ei_e = 10.604;

    //std::cout << "Minimum polar angle of the electron (degree): " << std::flush;
    //std::cin.getline(mychar, 64);
    //theta_min = atof(mychar) * deg;
    theta_min = atof(argv[1])* deg;

    //std::cout << "Maximum polar angle of the electron (degree): " << std::flush;
    //std::cin.getline(mychar, 64);
    theta_max = atof(argv[2]) * deg;

    //std::cout << "Minimum energy of bremsstrahlung photons (MeV, -1 to use default value): " << std::flush;
    //std::cin.getline(mychar, 64);
    E_g_min = atof(argv[3]) ;
    E_g_max = atof(argv[4]) ;
    int N = atoi(argv[5]);

    // random number generator
    PseRan = new TRandom2();
    PseRan->SetSeed((int)(time(NULL)));

    // solid angle
    omega = 2.0 * pi * (Cos(theta_min) - Cos(theta_max));

    // initial states, fixed target and beam along +z
    vi_e.SetPxPyPzE(0.0, 0.0, Sqrt(Pow2(Ei_e) - m2), Ei_e);
    vi_d.SetPxPyPzE(0.0, 0.0, 0.0, M);


//    // output file
//    FILE *fp = fopen("ehe4_events.dat", "w");

//    for (int i = 0; i < InterpolPoints; i++) {
//        theta_e = theta_min + i * (theta_max - theta_min) / (InterpolPoints - 1);
//        theta[i] = theta_e;
//        double sinth = Sin(theta_e);
//
//        double e_elastic = ElasticEnergy(theta_e);
//        double e_ef = e_elastic - E_g_min;//maximal e_ef
//        vf_e.SetPxPyPzE(Sqrt(Pow2(e_ef) - m2) * Sin(theta_e), 0, Sqrt(Pow2(e_ef) - m2) * Cos(theta_e), e_ef);
//        v_min = (vi_e + vi_d - vf_e) * (vi_e + vi_d - vf_e) - M2;//maximal e_ef corresponds to minimal v
//        v_min = 0.001;//fixed for now
//
//        if (E_g_max < e_elastic - m)
//            e_ef = e_elastic - E_g_max;
//        else
//            e_ef = m;
//
//        vf_e.SetPxPyPzE(Sqrt(Pow2(e_ef) - m2) * Sin(theta_e), 0, Sqrt(Pow2(e_ef) - m2) * Cos(theta_e), e_ef);
//        v_cut = (vi_e + vi_d - vf_e) * (vi_e + vi_d - vf_e) - M2;
//
//        xs_elastic_sin[i] = ElasticXS_Sin_ed(theta_e) * mkb;
//        xs_born_sin[i] = BornXS_Sin_ed(theta_e) * mkb;
//
//        double sigma_born = xs_born_sin[i] / sinth;
//        double sigma_elastic = xs_elastic_sin[i] / sinth;
//
//        e_ef = e_elastic;
//        vf_e.SetPxPyPzE(Sqrt(Pow2(e_ef) - m2) * Sin(theta_e), 0, Sqrt(Pow2(e_ef) - m2) * Cos(theta_e), e_ef);
//        double q2 = -(vi_e - vf_e) * (vi_e - vf_e);
//        //std::cout<<theta_e * 180.0 / pi<<" "<<q2<<" "<<sigma_born<<" "<<sigma_elastic<<" "<<sigma_elastic / sigma_born<<std::endl;
//        //printf("%8.6lf %10.4le %10.4le %10.4le %8.6lf\n", theta_e * 180.0 / pi, q2, sigma_born, sigma_elastic, sigma_elastic / sigma_born);
//        fprintf(fp, "%8.6lf %12.6le %12.6le %12.6le %8.6lf\n", theta_e * 180.0 / pi, q2, sigma_born, sigma_elastic, sigma_elastic / sigma_born);
//    }
//
//    fclose(fp);
//
//    Interpolator_ElasticXS_Sin_ed.SetData(InterpolPoints, theta, xs_elastic_sin);
//    Interpolator_BornXS_Sin_ed.SetData(InterpolPoints, theta, xs_born_sin);
//
//    TFoam *FoamElastic = new TFoam("FoamElastic");
//    TFoamIntegrand *pFoamElastic = new ElasticIntegrand();
//    FoamElastic->SetkDim(1);
//    FoamElastic->SetnCells(1000); // Set number of cells
//    FoamElastic->SetnSampl(200); // Set number of samples
//    FoamElastic->SetOptRej(1); // Unweighted events in MC generation
//    FoamElastic->SetRho(pFoamElastic); // Set distribution function
//    FoamElastic->SetPseRan(PseRan); // Set random number generator
//    //FoamElastic->SetChat(1); // Set "chat level" in the standard output
//    FoamElastic->Initialize();

    TFoam *FoamBrems = new TFoam("FoamBrems");
    TFoamIntegrand *pFoamBrems = new BremsIntegrand();
    FoamBrems->SetkDim(4);
    FoamBrems->SetnCells(30000); // Set number of cells
    FoamBrems->SetnSampl(1500); // Set number of samples
    FoamBrems->SetOptRej(1); // Unweighted events in MC generation
    FoamBrems->SetRho(pFoamBrems); // Set distribution function
    FoamBrems->SetPseRan(PseRan); // Set random number generator
    //FoamBrems->SetChat(1); // Set "chat level" in the standard output
    FoamBrems->Initialize();

    for (int i = 0; i < 10000000; i++) {
        if (i % 100000 == 0 && i != 0) std::cout << " Initializing:" << i << "\r" << std::flush;

        FoamBrems->MakeEvent();
    }

//    double xsint_elastic = Interpolator_ElasticXS_Sin_ed.Integ(theta_min, theta_max);
//    double xsint_born = Interpolator_BornXS_Sin_ed.Integ(theta_min, theta_max);
    double xsint_brems, xsint_brems_error;
    FoamBrems->GetIntegMC(xsint_brems, xsint_brems_error);

    //double xsint = xsint_elastic + xsint_brems;
    double xsint = xsint_brems;
    double xsint_error = xsint_brems_error;

    int n_elastic = 0;//int(N * (xsint_elastic / xsint));




    FILE * fp = fopen(argv[6], "w");

    int count_elastic = 0, count_brems = 0;

    for (int i = 0; i < N; ++i) {
        if (i % 10000 == 0 && i != 0) std::cout << i << std::endl;

//        if ((n_elastic - count_elastic) > 0 && PseRan->Rndm() < 1.0 * (n_elastic - count_elastic) / (N - i)) {
//            FoamElastic->MakeEvent();
//
//            Ef_e = ElasticEnergy(theta_e);
//
//            vf_e.SetPxPyPzE(Sqrt(Pow2(Ef_e) - m2) * Sin(theta_e), 0.0, Sqrt(Pow2(Ef_e) - m2) * Cos(theta_e), Ef_e);
//            vf_d = vi_e + vi_d - vf_e;
//
//            double phi = 2.0 * pi * PseRan->Rndm();
//
//            vf_e.RotateZ(phi);
//            vf_d.RotateZ(phi);
//
//            Ef_e = vf_e.E();
//            theta_e = vf_e.Theta();
//            phi_e = vf_e.Phi();
//
//            Ef_d = vf_d.E();
//            theta_d = vf_d.Theta();
//            phi_d = vf_d.Phi();
//
//            E_g = 0.0;
//            theta_g = 0.0;
//            phi_g = 0.0;
//
//            count_elastic++;
//        } else {
            FoamBrems->MakeEvent();

            double phi = 2.0 * pi * PseRan->Rndm();

            vf_e.RotateZ(phi);
            vf_d.RotateZ(phi);
            v_g.RotateZ(phi);

            Ef_e = vf_e.E();
            theta_e = vf_e.Theta();
            phi_e = vf_e.Phi();

            Ef_d = vf_d.E();
            theta_d = vf_d.Theta();
            phi_d = vf_d.Phi();

            E_g = v_g.E();
            theta_g = v_g.Theta();
            phi_g = v_g.Phi();

            count_brems++;
//        }

        fprintf(fp, "%11.5lf %10.8lf %11.8lf %11.5lf %10.8lf %11.8lf %11.5lf %10.8lf %11.8lf\n", Ef_e , theta_e, phi_e, Ef_d , theta_d, phi_d, E_g , theta_g, phi_g);
        //printf("%11.5lf %10.8lf %11.8lf %11.5lf %10.8lf %11.8lf %11.5lf %10.8lf %11.8lf\n", Ef_e * 1000, theta_e, phi_e, Ef_d * 1000, theta_d, phi_d, E_g * 1000, theta_g, phi_g);
    }

    fclose(fp);

    fp = fopen(argv[7], "w");

    fprintf(fp, "[\n  {\n    \"beam energy\":%lf,\n     ", Ei_e);
    fprintf(fp, "\"polar angle min\":%lf,\n     ", theta_min/deg);
    fprintf(fp, "\"polar angle max\":%lf,\n     ", theta_max/deg);
    fprintf(fp, "\"polar energy min\":%lf,\n     ", E_g_min);
    fprintf(fp, "\"polar energy max\":%lf,\n     ", E_g_max);
    fprintf(fp, "\"angle acceptance\":%lf,\n     ", omega);
    fprintf(fp, "\"total cross section\":%lf\n  }\n]", xsint);

//    std::cout << xsint_born << " " << xsint  << " " << xsint_elastic << " "  << xsint_brems << " " << std::endl;
//    std::cout <<"Ratio:"<<(xsint/xsint_born)-1.0<<std::endl;


    fclose(fp);
}
