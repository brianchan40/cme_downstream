#include <TFile.h>
#include "error_calc.C"

using namespace std;

TFile *outputfile, *inputfile;

float v2_q[4], v2_q_pair[4], v2_q_err[4], v2_q_pair_err[4], v2_pair_q[4], v2_pair_q_pair[4], v2_pair_q_err[4], v2_pair_q_pair_err[4], v2_single[4], v2_pair[4], v2_single_err[4], v2_pair_err[4];
std::vector<TString> v2_x_labels;

float single_q_single_v2[4], single_q_pair_v2[4], pair_q_single_v2[4], pair_q_pair_v2[4], ensemble_avg[4];
float single_q_single_v2_err[4], single_q_pair_v2_err[4], pair_q_single_v2_err[4], pair_q_pair_v2_err[4], ensemble_avg_err[4];

TCanvas *v2_q2_canvas = new TCanvas("v2_q2_canvas", "v2_q2_canvas", 1500, 800);

std::vector<float> plot_TProfile(const char *profile_name, const char *x_axis_name, const char *y_axis_name, const char *profile_title, bool fit, Color_t mcolor = NULL, int option = 0);
std::vector<float> plot_TProfile(TGraphErrors *temp_profile, const char *x_axis_name, const char *y_axis_name, const char *profile_title, bool fit);
std::vector<float> plot_combined(TProfile *profile1_name, TProfile *profile2_name, const char *x_axis_name, const char *y_axis_name, const char *profile_title, bool fit, int rebin);
std::vector<float> plot_combined(TProfile *profile1_name, TGraphErrors *profile2_name, const char *x_axis_name, const char *y_axis_name, const char *profile_title, bool fit, int rebin);

void plotting_TGraph_Helper(const char *title, const char *x_title, const char *y_title, int npts, std::vector<TString> x_labels, float *y1, float *y1_err, float *y2 = NULL, float *y2_err = NULL, float *y3 = NULL, float *y3_err = NULL, float *y4 = NULL, float *y4_err = NULL, float *y5 = NULL, float *y5_err = NULL, float *y6 = NULL, float *y6_err = NULL);
TGraphErrors* subtract_TProfiles(const char *first_profile, const char *second_profile);
TGraphErrors* multiply_TProfiles(const char *first_profile, const char *second_profile);

void plot_v2_pt_profiles();
void plot_v2_q2_profiles();

void cme_downstream()
{
    int rebin = 2;

    outputfile = new TFile("output.root", "recreate");
    inputfile = new TFile("cmeTest2_17_noRho_NewPlots.root");
    // inputfile = new TFile("cmeTest2_17_NewPlots.root");
    
    v2_x_labels.push_back("Primordial");
    v2_x_labels.push_back("Decay");
    v2_x_labels.push_back("All");
    v2_x_labels.push_back("#rho");

    plot_v2_pt_profiles();
    plot_v2_q2_profiles();

    TGraphErrors *primordial_Gamma_os_Q2 = subtract_TProfiles("primordial_Gamma_os_Q2", "primordial_Gamma_ss_Q2");
    TGraphErrors *both_prim_Gamma_os_Q2 = subtract_TProfiles("both_prim_Gamma_os_Q2", "both_prim_Gamma_ss_Q2");
    TGraphErrors *both_glob_Gamma_os_Q2 = subtract_TProfiles("both_glob_Gamma_os_Q2", "both_glob_Gamma_ss_Q2");
    TGraphErrors *one_prim_Gamma_os_Q2 = subtract_TProfiles("one_prim_Gamma_os_Q2", "one_prim_Gamma_ss_Q2");

    both_glob_Gamma_os_Q2->Fit("pol0");
    both_glob_Gamma_os_Q2->Write("Both_Decay_DeltaGamma_Q2");
    
    plot_TProfile(primordial_Gamma_os_Q2, "Single Pion q^{2}", "#Delta#gamma", "Pion #Delta#gamma vs. Single Pion q^{2}", true); 
    plot_TProfile(both_prim_Gamma_os_Q2, "Single Pion q^{2}", "Both Primordial #Delta#gamma", "Both Primordial Pion #Delta#gamma vs. Single Pion q^{2}", true);
    plot_TProfile(both_glob_Gamma_os_Q2, "Single Pion q^{2}", "Both Decay #Delta#gamma", "Both Decay Pion #Delta#gamma vs. Single Pion q^{2}", true);
    plot_TProfile(one_prim_Gamma_os_Q2, "Single Pion q^{2}", "One Primordial #Delta#gamma", "One Primordial One Decay Pion #Delta#gamma vs. Single Pion q^{2}", true);

    v2_q_temp = plot_combined((TProfile *)inputfile->Get("all_v2_Q2"), primordial_Gamma_os_Q2, "All Pions v_{2}", "#Delta#gamma", "All Pion-Pair #Delta#gamma vs. All Pion v_{2} (combined based on Single Pion q^{2})", true, rebin);
    single_q_single_v2[0] = (v2_q_temp.at(0));
    single_q_single_v2_err[0] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    v2_q_temp = plot_combined((TProfile *)inputfile->Get("all_v2_Q2"), both_prim_Gamma_os_Q2, "All Pions v_{2}", "Both Primordial #Delta#gamma", "Both Primordial Pion-Pair #Delta#gamma vs. All Pion v_{2} (combined based on Single Pion q^{2})", true, rebin);
    single_q_single_v2[2] = (v2_q_temp.at(0));
    single_q_single_v2_err[2] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    v2_q_temp = plot_combined((TProfile *)inputfile->Get("all_v2_Q2"), both_glob_Gamma_os_Q2, "All Pions v_{2}", "Both Decay #Delta#gamma", "Both Decay Pion-Pair #Delta#gamma vs. All Pion v_{2} (combined based on Single Pion q^{2})", true, rebin);
    single_q_single_v2[3] = (v2_q_temp.at(0));
    single_q_single_v2_err[3] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    v2_q_temp = plot_combined((TProfile *)inputfile->Get("all_v2_Q2"), one_prim_Gamma_os_Q2, "All Pions v_{2}", "One Primordial #Delta#gamma", "One Primordial One Decay Pion-Pair #Delta#gamma vs. All Pion v_{2} (combined based on Single Pion q^{2})", true, rebin);
    single_q_single_v2[1] = (v2_q_temp.at(0));
    single_q_single_v2_err[1] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    v2_q_temp = plot_combined((TProfile *)inputfile->Get("all_v2pair_Q2"), primordial_Gamma_os_Q2, "Pion-Pair v_{2}", "#Delta#gamma", "All Pion-Pair #Delta#gamma vs. Pion-Pair v_{2} (combined based on Single Pion q^{2})", true, rebin);
    single_q_pair_v2[0] = (v2_q_temp.at(0));
    single_q_pair_v2_err[0] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    v2_q_temp = plot_combined((TProfile *)inputfile->Get("all_v2pair_Q2"), both_prim_Gamma_os_Q2, "Pion-Pair v_{2}", "Both Primordial #Delta#gamma", "Both Primordial Pion-Pair #Delta#gamma vs. Pion-Pair v_{2} (combined based on Single Pion q^{2})", true, rebin);
    single_q_pair_v2[2] = (v2_q_temp.at(0));
    single_q_pair_v2_err[2] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    v2_q_temp = plot_combined((TProfile *)inputfile->Get("all_v2pair_Q2"), both_glob_Gamma_os_Q2, "Pion-Pair v_{2}", "Both Decay #Delta#gamma", "Both Decay Pion-Pair #Delta#gamma vs. Pion-Pair v_{2} (combined based on Single Pion q^{2})", true, rebin);
    single_q_pair_v2[3] = (v2_q_temp.at(0));
    single_q_pair_v2_err[3] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    v2_q_temp = plot_combined((TProfile *)inputfile->Get("all_v2pair_Q2"), one_prim_Gamma_os_Q2, "Pion-Pair v_{2}", "One Primordial #Delta#gamma", "One Primordial One Decay Pion-Pair #Delta#gamma vs. Pion-Pair v_{2} (combined based on Single Pion q^{2})", true, rebin);
    single_q_pair_v2[1] = (v2_q_temp.at(0));
    single_q_pair_v2_err[1] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    plot_combined((TProfile *)inputfile->Get("all_v2_Q2"), (TProfile *)inputfile->Get("rho_v2_Q2"), "All Pions v_{2}", "#rho v2", "#rho v_{2} vs. All Pion v_{2} (combined based on Single Pion q^{2})", true, 2);
    plot_combined((TProfile *)inputfile->Get("all_v2pair_Q2"), (TProfile *)inputfile->Get("rho_v2_Q2"), "Pion-Pion v_{2}", "#rho v2", "#rho v_{2} vs. Pion-Pair v_{2} (combined based on Single Pion q^{2})", true, 2);
    plot_combined((TProfile *)inputfile->Get("all_v2_Q2"), (TProfile *)inputfile->Get("global_v2_Q2"), "All Pions v_{2}", "Decay Pions v2", "Decay Pions v_{2} vs. All Pion v_{2} (combined based on Single Pion q^{2})", true, 2);
    plot_combined((TProfile *)inputfile->Get("all_v2pair_Q2"), (TProfile *)inputfile->Get("global_v2_Q2"), "Pion-Pion v_{2}", "Decay Pions v2", "Decay Pions v_{2} vs. Pion-Pair v_{2} (combined based on Single Pion q^{2})", true, 2);

    
    TGraphErrors *primordial_Gamma_os_pairQ2 = subtract_TProfiles("primordial_Gamma_os_pairQ2", "primordial_Gamma_ss_pairQ2");
    TGraphErrors *both_prim_Gamma_os_pairQ2 = subtract_TProfiles("both_prim_Gamma_os_pairQ2", "both_prim_Gamma_ss_pairQ2");
    TGraphErrors *both_glob_Gamma_os_pairQ2 = subtract_TProfiles("both_glob_Gamma_os_pairQ2", "both_glob_Gamma_ss_pairQ2");
    TGraphErrors *one_prim_Gamma_os_pairQ2 = subtract_TProfiles("one_prim_Gamma_os_pairQ2", "one_prim_Gamma_ss_pairQ2");

    both_glob_Gamma_os_pairQ2->Fit("pol0");
    both_glob_Gamma_os_pairQ2->Write("Both_Decay_DeltaGamma_pairQ2");

    plot_TProfile(primordial_Gamma_os_pairQ2, "Pion-Pair q^{2}", "#Delta#gamma", "Pion #Delta#gamma vs. Pion-Pair q^{2}", true);
    plot_TProfile(both_prim_Gamma_os_pairQ2, "Pion-Pair q^{2}", "Both Primordial #Delta#gamma", "Both Primordial Pion #Delta#gamma vs. Pion-Pair q^{2}", true);
    plot_TProfile(both_glob_Gamma_os_pairQ2, "Pion-Pair q^{2}", "Both Decay #Delta#gamma", "Both Decay Pion #Delta#gamma vs. Pion-Pair q^{2}", true);
    plot_TProfile(one_prim_Gamma_os_pairQ2, "Pion-Pair q^{2}", "One Primordial #Delta#gamma", "One Primordial One Decay Pion #Delta#gamma vs. Pion-Pair q^{2}", true);
    
    v2_q_temp = plot_combined((TProfile *)inputfile->Get("all_v2_pairQ2"), primordial_Gamma_os_pairQ2, "All Pions v_{2}", "#Delta#gamma", "All Pion-Pair #Delta#gamma vs. All Pion v_{2} (combined based on Pion-Pair q^{2})", true, rebin);
    pair_q_single_v2[0] = (v2_q_temp.at(0));
    pair_q_single_v2_err[0] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    v2_q_temp = plot_combined((TProfile *)inputfile->Get("all_v2_pairQ2"), both_prim_Gamma_os_pairQ2, "All Pions v_{2}", "Both Primordial #Delta#gamma", "Both Primordial Pion-Pair #Delta#gamma vs. All Pion v_{2} (combined based on Pion-Pair q^{2})", true, rebin);
    pair_q_single_v2[2] = (v2_q_temp.at(0));
    pair_q_single_v2_err[2] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    v2_q_temp = plot_combined((TProfile *)inputfile->Get("all_v2_pairQ2"), both_glob_Gamma_os_pairQ2, "All Pions v_{2}", "Both Decay #Delta#gamma", "Both Decay Pion-Pair #Delta#gamma vs. All Pion v_{2} (combined based on Pion-Pair q^{2})", true, rebin);
    pair_q_single_v2[3] = (v2_q_temp.at(0));
    pair_q_single_v2_err[3] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    v2_q_temp = plot_combined((TProfile *)inputfile->Get("all_v2_pairQ2"), one_prim_Gamma_os_pairQ2, "All Pions v_{2}", "One Primordial #Delta#gamma", "One Primordial One Decay Pion-Pair #Delta#gamma vs. All Pion v_{2} (combined based on Pion-Pair q^{2})", true, rebin);
    pair_q_single_v2[1] = (v2_q_temp.at(0));
    pair_q_single_v2_err[1] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    v2_q_temp = plot_combined((TProfile *)inputfile->Get("all_v2pair_pairQ2"), primordial_Gamma_os_pairQ2, "Pion-Pairs v_{2}", "#Delta#gamma", "All Pion-Pair #Delta#gamma vs. Pion-Pair v_{2} (combined based on Pion-Pair q^{2})", true, rebin);
    pair_q_pair_v2[0] = (v2_q_temp.at(0));
    pair_q_pair_v2_err[0] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    v2_q_temp = plot_combined((TProfile *)inputfile->Get("all_v2pair_pairQ2"), both_prim_Gamma_os_pairQ2, "Pion-Pair v_{2}", "#Delta#gamma", "Both Primordial Pion #Delta#gamma vs. Pion-Pair v_{2} (combined based on Pion-Pair q^{2})", true, rebin);
    pair_q_pair_v2[2] = (v2_q_temp.at(0));
    pair_q_pair_v2_err[2] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    v2_q_temp = plot_combined((TProfile *)inputfile->Get("all_v2pair_pairQ2"), both_glob_Gamma_os_pairQ2, "Pion-Pair v_{2}", "#Delta#gamma", "Both Decay Pion #Delta#gamma vs. Pion-Pair v_{2} (combined based on Pion-Pair q^{2})", true, rebin);
    pair_q_pair_v2[3] = (v2_q_temp.at(0));
    pair_q_pair_v2_err[3] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    v2_q_temp = plot_combined((TProfile *)inputfile->Get("all_v2pair_pairQ2"), one_prim_Gamma_os_pairQ2, "Pion-Pair v_{2}", "#Delta#gamma", "One Primordial One Decay Pion #Delta#gamma vs. Pion-Pair v_{2} (combined based on Pion-Pair q^{2})", true, rebin);
    pair_q_pair_v2[1] = (v2_q_temp.at(0));
    pair_q_pair_v2_err[1] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    plot_combined((TProfile *)inputfile->Get("all_v2_pairQ2"), (TProfile *)inputfile->Get("rho_v2_pairQ2"), "All Pions v_{2}", "#rho v2", "#rho v_{2} vs. All Pion v_{2} (combined based on Pion-Pair q^{2})", true, 2);
    plot_combined((TProfile *)inputfile->Get("all_v2pair_pairQ2"), (TProfile *)inputfile->Get("rho_v2_pairQ2"), "All Pion-Pair v_{2}", "#rho v2", "#rho v_{2} vs. Pion-Pair v_{2} (combined based on Pion-Pair q^{2})", true, 2);
    plot_combined((TProfile *)inputfile->Get("all_v2_pairQ2"), (TProfile *)inputfile->Get("global_v2_pairQ2"), "All Pions v_{2}", "Decay Pions v2", "Decay Pions v_{2} vs. All Pion v_{2} (combined based on Pion-Pair q^{2})", true, 2);
    plot_combined((TProfile *)inputfile->Get("all_v2pair_pairQ2"), (TProfile *)inputfile->Get("global_v2_pairQ2"), "Pion-Pair v_{2}", "Decay Pions v2", "Decay Pions v_{2} vs. Pion-Pair v_{2} (combined based on Pion-Pair q^{2})", true, 2);

    plot_combined((TProfile *)inputfile->Get("all_v2_Q2"), (TProfile *)inputfile->Get("flow_samep_Q2"), "All Pions v_{2}", "<cos(phi1 + phi2 - 2*phi_rho)>", "<cos(phi1 + phi2 - 2*phi_rho)> vs. All Pions v_{2} (combined based on All Pion q^{2})", false, 2);
    plot_combined((TProfile *)inputfile->Get("all_v2pair_Q2"), (TProfile *)inputfile->Get("flow_samep_Q2"), "Pion-Pair v_{2}", "<cos(phi1 + phi2 - 2*phi_rho)>", "<cos(phi1 + phi2 - 2*phi_rho)> vs. Pion-Pair v_{2} (combined based on All Pion q^{2})", false, 2);
    plot_combined((TProfile *)inputfile->Get("all_v2_pairQ2"), (TProfile *)inputfile->Get("flow_samep_pair_Q2"), "All Pions v_{2}", "<cos(phi1 + phi2 - 2*phi_rho)>", "<cos(phi1 + phi2 - 2*phi_rho)> vs. All Pions v_{2} (combined based on Pion-Pair q^{2})", false, 2);
    plot_combined((TProfile *)inputfile->Get("all_v2pair_pairQ2"), (TProfile *)inputfile->Get("flow_samep_pair_Q2"), "Pion-Pair v_{2}", "<cos(phi1 + phi2 - 2*phi_rho)>", "<cos(phi1 + phi2 - 2*phi_rho)> vs. Pion-Pair v_{2} (combined based on Pion-Pair q^{2})", false, 2);

    plot_combined((TProfile *)inputfile->Get("all_v2_Q2"), (TProfile *)inputfile->Get("delta_samep_Q2"), "All Pions v_{2}", "<cos(phi1 - phi2)>", "<cos(phi1 - phi2)> vs. All Pions v_{2} (combined based on All Pion q^{2})", false, 2);
    plot_combined((TProfile *)inputfile->Get("all_v2pair_Q2"), (TProfile *)inputfile->Get("delta_samep_Q2"), "Pion-Pair v_{2}", "<cos(phi1 - phi2)>>", "<cos(phi1 - phi2)> vs. Pion-Pair v_{2} (combined based on All Pion q^{2})", false, 2);
    plot_combined((TProfile *)inputfile->Get("all_v2_pairQ2"), (TProfile *)inputfile->Get("delta_samep_pair_Q2"), "All Pions v_{2}", "<cos(phi1 - phi2)>>", "<cos(phi1 - phi2)> vs. All Pions v_{2} (combined based on Pion-Pair q^{2})", false, 2);
    plot_combined((TProfile *)inputfile->Get("all_v2pair_pairQ2"), (TProfile *)inputfile->Get("delta_samep_pair_Q2"), "Pion-Pair v_{2}", "<cos(phi1 - phi2)>", "<cos(phi1 - phi2)> vs. Pion-Pair v_{2} (combined based on Pion-Pair q^{2})", false, 2);

    
    TGraphErrors *Gamma_opp = subtract_TProfiles("Gamma_opp", "Gamma_same");
    cout << "Ensemble Average delta gamma is " << Gamma_opp->GetY()[0] << " +/- " << Gamma_opp->GetErrorY(0) << endl;
    Gamma_opp->Write("Ensemble_Avg_Delta_Gamma");
    ensemble_avg[0] = Gamma_opp->GetY()[0];
    ensemble_avg_err[0] = Gamma_opp->GetErrorY(0);

    TGraphErrors *both_primordial_Gamma_opp = subtract_TProfiles("both_primordial_Gamma_opp", "both_primordial_Gamma_same");
    cout << "Ensemble Average both_primordial delta gamma is " << both_primordial_Gamma_opp->GetY()[0] << " +/- " << both_primordial_Gamma_opp->GetErrorY(0) << endl;
    both_primordial_Gamma_opp->Write("Ensemble_Avg_BothPrim_Delta_Gamma");
    ensemble_avg[2] = both_primordial_Gamma_opp->GetY()[0];
    ensemble_avg_err[2] = both_primordial_Gamma_opp->GetErrorY(0);

    TGraphErrors *both_decay_Gamma_opp = subtract_TProfiles("both_decay_Gamma_opp", "both_decay_Gamma_same");
    cout << "Ensemble Average both_decay delta gamma is " << both_decay_Gamma_opp->GetY()[0] << " +/- " << both_decay_Gamma_opp->GetErrorY(0) << endl;
    both_decay_Gamma_opp->Write();
    ensemble_avg[3] = both_decay_Gamma_opp->GetY()[0];
    ensemble_avg_err[3] = both_decay_Gamma_opp->GetErrorY(0);

    TGraphErrors *one_primordial_Gamma_opp = subtract_TProfiles("one_primordial_Gamma_opp", "one_primordial_Gamma_same");
    cout << "Ensemble Average one_primordial delta gamma is " << one_primordial_Gamma_opp->GetY()[0] << " +/- " << one_primordial_Gamma_opp->GetErrorY(0) << endl;
    one_primordial_Gamma_opp->Write();
    ensemble_avg[1] = one_primordial_Gamma_opp->GetY()[0];
    ensemble_avg_err[1] = one_primordial_Gamma_opp->GetErrorY(0);

    cout << "                          All Pions          One Decay One Primordial          Both Primordial            Both Decay" << endl;
    cout << "single q single v2: ";
    for (int fg = 0; fg < 4; fg++){
        cout << single_q_single_v2[fg] << " +/- " << single_q_single_v2_err[fg] << " ";
    }
    cout << endl;
    cout << "single q pair v2: ";
    for (int fg = 0; fg < 4; fg++){
        cout << single_q_pair_v2[fg] << " +/- " << single_q_pair_v2_err[fg] << " ";
    }
    cout << endl;
    cout << "pair q single v2: ";
    for (int fg = 0; fg < 4; fg++){
        cout << pair_q_single_v2[fg] << " +/- " << pair_q_single_v2_err[fg] << " ";
    }
    cout << endl;
    cout << "pair q pair v2: ";
    for (int fg = 0; fg < 4; fg++){
        cout << pair_q_pair_v2[fg] << " +/- " << pair_q_pair_v2_err[fg] << " ";
    }
    cout << endl;
    cout << "ensemble_avg: ";
    for (int fg = 0; fg < 4; fg++){
        cout << ensemble_avg[fg] << " +/- " << ensemble_avg_err[fg] << " ";
    }
    cout << endl;

    plotting_TGraph_Helper("4 ESE Approaches vs. Ensemble Average", "", "#Delta#gamma_{112}", 4, v2_x_labels, single_q_single_v2, single_q_single_v2_err, single_q_pair_v2, single_q_pair_v2_err, pair_q_single_v2, pair_q_single_v2_err, pair_q_pair_v2, pair_q_pair_v2_err, ensemble_avg, ensemble_avg_err);

    TProfile *gamma_ss_pt_profile = (TProfile *)inputfile->Get("gamma_ss_pt_profile");
    TProfile *gamma_os_pt_profile = (TProfile *)inputfile->Get("gamma_os_pt_profile");
    TProfile *delta_ss_pt_profile = (TProfile *)inputfile->Get("delta_ss_pt_profile");
    TProfile *delta_os_pt_profile = (TProfile *)inputfile->Get("delta_os_pt_profile");
    TProfile *pair_v2_pt_profile = (TProfile *)inputfile->Get("pair_v2_pt_profile");

    gamma_ss_pt_profile->Sumw2();
    gamma_os_pt_profile->Sumw2();
    delta_ss_pt_profile->Sumw2();
    delta_os_pt_profile->Sumw2();
    pair_v2_pt_profile->Sumw2();

    gamma_ss_pt_profile->Write();
    gamma_os_pt_profile->Write();
    delta_ss_pt_profile->Write();
    delta_os_pt_profile->Write();
    pair_v2_pt_profile->Write();

    const int n_pt_bins = 200;
    float delta_gamma[n_pt_bins] = {0.}, delta_delta[n_pt_bins] = {0.}, delta_v2[n_pt_bins] = {0.}, kappa[n_pt_bins] = {0.};
    float delta_gamma_err[n_pt_bins] = {0.}, delta_delta_err[n_pt_bins] = {0.}, delta_v2_err[n_pt_bins] = {0.}, kappa_err[n_pt_bins] = {0.};
    float X[n_pt_bins] = {0.}, X_err[n_pt_bins] = {0.};

    int final_pt_bins = 0;

    for (int o = 1; o <= n_pt_bins; o++)
    {
        if (((delta_os_pt_profile->GetBinContent(o) - delta_ss_pt_profile->GetBinContent(o)) * pair_v2_pt_profile->GetBinContent(o)) != 0)
        {
            X[final_pt_bins] = gamma_os_pt_profile->GetBinCenter(o);

            delta_gamma[final_pt_bins] = gamma_os_pt_profile->GetBinContent(o) - gamma_ss_pt_profile->GetBinContent(o);
            delta_gamma_err[final_pt_bins] = error_add(gamma_os_pt_profile->GetBinError(o), gamma_ss_pt_profile->GetBinError(o));

            delta_delta[final_pt_bins] = delta_os_pt_profile->GetBinContent(o) - delta_ss_pt_profile->GetBinContent(o);
            delta_delta_err[final_pt_bins] = error_add(delta_os_pt_profile->GetBinError(o), delta_ss_pt_profile->GetBinError(o));

            delta_v2[final_pt_bins] = delta_delta[final_pt_bins] * pair_v2_pt_profile->GetBinContent(o);
            delta_v2_err[final_pt_bins] = error_mult(delta_delta[final_pt_bins], pair_v2_pt_profile->GetBinContent(o), delta_delta_err[final_pt_bins], pair_v2_pt_profile->GetBinError(o));

            kappa[final_pt_bins] = delta_gamma[final_pt_bins] / delta_v2[final_pt_bins];
            kappa_err[final_pt_bins] = error_divide(delta_gamma[final_pt_bins], delta_v2[final_pt_bins], delta_gamma_err[final_pt_bins], delta_v2_err[final_pt_bins]);

            final_pt_bins++;
        }
    }

    TGraphErrors *delta_gamma_pt = new TGraphErrors(final_pt_bins, X, delta_gamma, X_err, delta_gamma_err);
    delta_gamma_pt->SetTitle("#Delta#gamma vs. p_{T}");
    delta_gamma_pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    delta_gamma_pt->GetYaxis()->SetTitle("#Delta#gamma");
    delta_gamma_pt->SetMarkerStyle(kFullCircle);
    delta_gamma_pt->Draw("AP");
    TGraphErrors *delta_delta_pt = new TGraphErrors(final_pt_bins, X, delta_delta, X_err, delta_delta_err);
    delta_delta_pt->SetTitle("#Delta#delta vs. p_{T}");
    delta_delta_pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    delta_delta_pt->GetYaxis()->SetTitle("#Delta#delta");
    delta_delta_pt->SetMarkerStyle(kFullCircle);
    delta_delta_pt->Draw("AP");
    TGraphErrors *delta_v2_pt = new TGraphErrors(final_pt_bins, X, delta_v2, X_err, delta_v2_err);
    delta_v2_pt->SetTitle("(#Delta#delta * v_{2}) vs. p_{T}");
    delta_v2_pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    delta_v2_pt->GetYaxis()->SetTitle("(#Delta#delta * v_{2})");
    delta_v2_pt->SetMarkerStyle(kFullCircle);
    delta_v2_pt->Draw("AP");
    TGraphErrors *kappa_pt = new TGraphErrors(final_pt_bins, X, kappa, X_err, kappa_err);
    kappa_pt->SetTitle("#kappa vs. p_{T}");
    kappa_pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    kappa_pt->GetYaxis()->SetTitle("#kappa");
    kappa_pt->SetMarkerStyle(kFullCircle);
    delta_v2_pt->Draw("AP");

    delta_gamma_pt->Write("delta_gamma_pt");
    delta_delta_pt->Write("delta_delta_pt");
    delta_v2_pt->Write("delta_v2_pt");
    kappa_pt->Write("kappa_pt");

    TCanvas v2_flow("v2_flow", "Combined v2_flow vs. gamma", 1500, 800);
    TGraphErrors *v2_flow_graph = multiply_TProfiles("rho_v2_pt_profile", "flow_samep_pT");
    v2_flow_graph->SetLineColor(kRed);
    v2_flow_graph->SetMarkerColor(kRed);
    v2_flow_graph->Draw();
    (TProfile*)inputfile->Get("both_glob_Gamma_opp_pT_samep")->Draw("sames");
    v2_flow.Write();

    TCanvas rho00_measured("rho00_measured", "rho00_measured", 1500, 800);
    TH1D *hCosTheta = (TH1D*) inputfile->Get("hCosTheta");
    TF1 *rho_func = new TF1("rho_func", "[1]* ((1-[0]) + (3*[0]-1)*x*x)", 0, 0.7);
    hCosTheta->Fit("rho_func", "R");
    hCosTheta->Draw();
    rho00_measured.Write();
    

    outputfile->Close();
    inputfile->Close();
}

void plot_v2_pt_profiles(){
    plot_TProfile("primordial_v2_pt_profile", "p_{T} (GeV/c)", "v_{2}", "Primordial Pions v_{2}(p_{T})", false);
    plot_TProfile("global_v2_pt_profile", "p_{T} (GeV/c)", "v_{2}", "Decay Pions v_{2}(p_{T})", false);
    plot_TProfile("all_v2_pt_profile", "p_{T} (GeV/c)", "v_{2}", "All Pions v_{2}(p_{T})", false);
    plot_TProfile("rho_v2_pt_profile", "p_{T} (GeV/c)", "v_{2}", "#rho v_{2}(p_{T})", false);
}

void plot_v2_q2_profiles(){
    std::vector<float> v2_q_temp;

    v2_q_temp = plot_TProfile("primordial_v2_Q2", "Single Pion q^{2}", "v_{2}", "Primordial Pions v_{2}(q^{2})", true, kBlue, 1);
    v2_q[1] = (v2_q_temp.at(0));
    v2_q_err[1] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    v2_q_temp = plot_TProfile("global_v2_Q2", "Single Particle q^{2}", "v_{2}", "Decay Pions v_{2}(q^{2})", true, kRed);
    v2_q[2] = (v2_q_temp.at(0));
    v2_q_err[2] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    v2_q_temp = plot_TProfile("all_v2_Q2", "Single Particle q^{2}", "v_{2}", "All Pions v_{2}(q^{2})", true, kBlack);
    v2_q[0] = (v2_q_temp.at(0));
    v2_q_err[0] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    v2_q_temp = plot_TProfile("rho_v2_Q2", "Single Particle q^{2}", "v_{2}", "#rho v_{2}(q^{2})", true);
    v2_q[3] = (v2_q_temp.at(0));
    v2_q_err[3] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    v2_q2_canvas->Write();

    v2_q_temp = plot_TProfile("primordial_v2_pairQ2", "Pair-Pion q^{2}", "v_{2}", "Primordial Pions v_{2}(q^{2}_{pair})", true);
    v2_q_pair[1] = (v2_q_temp.at(0));
    v2_q_pair_err[1] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    v2_q_temp = plot_TProfile("global_v2_pairQ2", "Pair-Pion q^{2}", "v_{2}", "Decay Pions v_{2}(q^{2}_{pair})", true);
    v2_q_pair[2] = (v2_q_temp.at(0));
    v2_q_pair_err[2] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    v2_q_temp = plot_TProfile("all_v2_pairQ2", "Pair-Pion q^{2}", "v_{2}", "All Pions v_{2}(q^{2}_{pair})", true);
    v2_q_pair[0] = (v2_q_temp.at(0));
    v2_q_pair_err[0] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    v2_q_temp = plot_TProfile("rho_v2_pairQ2", "Pair-Pion q^{2}", "v_{2}", "#rho v_{2}(q^{2}_{pair})", true);
    v2_q_pair[3] = (v2_q_temp.at(0));
    v2_q_pair_err[3] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    v2_q_temp = plot_TProfile("all_v2pair_Q2", "Single Particle q^{2}", "Pion-Pair v_{2}", "All Pions v_{2,pair}(q^{2})", true);
    v2_pair_q[0] = (v2_q_temp.at(0));
    v2_pair_q_err[0] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    v2_pair_q[1] = -10000000;
    v2_pair_q[2] = -10000000;
    v2_pair_q[3] = -10000000;

    v2_q_temp = plot_TProfile("all_v2pair_pairQ2", "Pair-Pion q^{2}", "Pion-Pair v_{2}", "All Pions v_{2,pair}(q^{2}_{pair})", true);
    v2_pair_q_pair[0] = (v2_q_temp.at(0));
    v2_pair_q_pair_err[0] = (v2_q_temp.at(1));
    v2_q_temp.clear();

    v2_pair_q_pair[1] = -1000000;
    v2_pair_q_pair[2] = -1000000;
    v2_pair_q_pair[3] = -1000000;

    TProfile *v2_single_profile = (TProfile *)inputfile->Get("v2_single");
    TProfile *v2_pair_profile = (TProfile *)inputfile->Get("v2_pair");
    v2_single[0] = v2_single_profile->GetBinContent(1) * 100.0;
    v2_single_err[0] = v2_single_profile->GetBinError(1) * 100.0;
    v2_pair[0] = v2_pair_profile->GetBinContent(1) * 100.0;
    v2_pair_err[0] = v2_pair_profile->GetBinError(1) * 100.0;

    v2_single[1] = -1000000;
    v2_single[2] = -1000000;
    v2_single[3] = -1000000;
    v2_pair[1] = -1000000;
    v2_pair[2] = -1000000;
    v2_pair[3] = -1000000;

    plotting_TGraph_Helper("y-intercept of v_{2}(q^{2}) and v_{2}(q^{2}_{pair})", "", "v_{2}%", 4, v2_x_labels, v2_q, v2_q_err, v2_q_pair, v2_q_pair_err, v2_pair_q, v2_pair_q_err, v2_pair_q_pair, v2_pair_q_pair_err, v2_single, v2_single_err, v2_pair, v2_pair_err);
}

TGraphErrors* subtract_TProfiles(const char *first_profile, const char *second_profile){
    
    TProfile *first_TProfile = (TProfile *)inputfile->Get(first_profile);
    TProfile *second_TProfile = (TProfile *)inputfile->Get(second_profile);

    TGraphErrors* return_graph;
    float x[1000] = {0.}, x_err[1000] = {0.}, y[1000] = {0.}, y_err[1000] = {0.};
    
    int n_bins = first_TProfile->GetNbinsX();

    for (int k = 1; k <= n_bins ; k++){
        x[k-1] = first_TProfile->GetBinCenter(k);
        y[k-1] = first_TProfile->GetBinContent(k) - second_TProfile->GetBinContent(k);
        y_err[k-1] = error_add(first_TProfile->GetBinError(k), second_TProfile->GetBinError(k));
        
        /*if ((y[k-1] != 0) && (y_err[k-1] != 0)){
            cout << "x = " << x[k-1] << endl;
            cout << "first = " << first_TProfile->GetBinContent(k) << endl;
            cout << "second = " << second_TProfile->GetBinContent(k) << endl;
            cout << "y = " << y[k-1] << " +/- " << y_err[k-1] << endl;
        }*/
    }

    return_graph = new TGraphErrors(n_bins, x, y, x_err, y_err);
    return_graph->SetName(first_profile);

    return return_graph;
}

TGraphErrors* multiply_TProfiles(const char *first_profile, const char *second_profile){
    
    TProfile *first_TProfile = (TProfile *)inputfile->Get(first_profile);
    TProfile *second_TProfile = (TProfile *)inputfile->Get(second_profile);

    TGraphErrors* return_graph;
    float x[1000] = {0.}, x_err[1000] = {0.}, y[1000] = {0.}, y_err[1000] = {0.};
    
    int n_bins = first_TProfile->GetNbinsX();

    for (int k = 1; k <= n_bins ; k++){
        x[k-1] = first_TProfile->GetBinCenter(k);
        y[k-1] = first_TProfile->GetBinContent(k) * second_TProfile->GetBinContent(k);
        y_err[k-1] = error_mult(first_TProfile->GetBinContent(k), second_TProfile->GetBinContent(k), first_TProfile->GetBinError(k), second_TProfile->GetBinError(k));
        
        /*if ((y[k-1] != 0) && (y_err[k-1] != 0)){
            cout << "x = " << x[k-1] << endl;
            cout << "first = " << first_TProfile->GetBinContent(k) << endl;
            cout << "second = " << second_TProfile->GetBinContent(k) << endl;
            cout << "y = " << y[k-1] << " +/- " << y_err[k-1] << endl;
        }*/
    }

    return_graph = new TGraphErrors(n_bins, x, y, x_err, y_err);
    return_graph->SetName(first_profile);

    return return_graph;
}

std::vector<float> plot_TProfile(const char *profile_name, const char *x_axis_name, const char *y_axis_name, const char *profile_title, bool fit, Color_t draw_color, int option)
{
    TProfile *temp_profile = (TProfile *)inputfile->Get(profile_name);

    outputfile->cd();
    TCanvas c1(profile_name, profile_title, 1500, 800);
    TF1 *function_tmp = new TF1("function_tmp", "[0] + [1] * x + [2] * x * x", 0, 1.5);

    float y_intercept_tmp = 0, y_intercept_tmp_err = 0;

    if (fit)
    {
        temp_profile->Fit(function_tmp, "QNR");
        y_intercept_tmp = function_tmp->GetParameter(0) * 100.0;
        y_intercept_tmp_err = function_tmp->GetParError(0) * 100.0;
        cout << setprecision(3);
        cout << profile_title << " y-intercept = " << y_intercept_tmp << " +/- " << y_intercept_tmp_err << endl;
        function_tmp->Draw("same");
    }
    temp_profile->GetXaxis()->SetTitle(x_axis_name);
    temp_profile->GetYaxis()->SetTitle(y_axis_name);
    temp_profile->SetTitle(profile_title);
    temp_profile->Draw();
    if (fit)
    {
        TLegend *leg = new TLegend(0.1, 0.7, 0.48, 0.9);
        leg->AddEntry((TObject *)0, TString::Format("y-intercept = %.3e +/- %.3e %%", y_intercept_tmp, y_intercept_tmp_err).Data(), "");
        leg->Draw();
    }

    TLine l1(0, 0, 15, 0);
    l1.Draw();

    c1.Write();

    temp_profile->Write(TString::Format("%s_profile", profile_name));

    std::vector<float> temp_vector_TProfile;
    temp_vector_TProfile.push_back(y_intercept_tmp);
    temp_vector_TProfile.push_back(y_intercept_tmp_err);

    return temp_vector_TProfile;
}

std::vector<float> plot_combined(TProfile *profile1_name, TProfile *profile2_name, const char *x_axis_name, const char *y_axis_name, const char *profile_title, bool fit, int rebin)
{
    // TProfile *temp_profile1 = (TProfile*) inputfile->Get(profile1_name);
    // TProfile *temp_profile2 = (TProfile*) inputfile->Get(profile2_name);

    float new_x[1000] = {0.}, new_y[1000] = {0.}, new_x_err[1000] = {0.}, new_y_err[1000] = {0.};
    int N_bins = profile1_name->GetNbinsX() / 5;
    float new_x_numerator[1000] = {0.}, new_x_denomenator[1000] = {0.}, new_y_numerator[1000] = {0.}, new_y_denomenator[1000] = {0.};

    for (int i = 0; i < N_bins / rebin; i++)
    {
        for (int j = 0; j < rebin; j++)
        {
            if (profile1_name->GetBinError(rebin * i + 1 + j) == 0 || profile2_name->GetBinError(rebin * i + 1 + j) == 0)
                continue;
            new_x_numerator[i] += profile1_name->GetBinContent(rebin * i + 1 + j) / pow(profile1_name->GetBinError(rebin * i + 1 + j), 2);
            new_x_denomenator[i] += 1.0 / pow(profile1_name->GetBinError(rebin * i + 1 + j), 2);
            new_y_numerator[i] += profile2_name->GetBinContent(rebin * i + 1 + j) / pow(profile2_name->GetBinError(rebin * i + 1 + j), 2);
            new_y_denomenator[i] += 1.0 / pow(profile2_name->GetBinError(rebin * i + 1 + j), 2);
        }

        new_x[i] = new_x_numerator[i] / new_x_denomenator[i];
        new_x_err[i] = sqrt(1.0 / new_x_denomenator[i]);
        new_y[i] = new_y_numerator[i] / new_y_denomenator[i];
        new_y_err[i] = sqrt(1.0 / new_y_denomenator[i]);

        // cout << "new_x[" << i << "] = " << new_x[i] << endl;
        // cout << "new_y[" << i << "] = " << new_y[i] << endl;
    }

    TGraphErrors *temp_graph = new TGraphErrors(N_bins / rebin, new_x, new_y, new_x_err, new_y_err);

    outputfile->cd();
    TCanvas c1(profile_title, profile_title, 1500, 800);
    TF1 *function_tmp = new TF1("function_tmp", "[0] + [1] * x + [2] * x * x", 0, 0.08);
    if (fit)
    {
        temp_graph->Fit(function_tmp, "QNR");
        cout << setprecision(3);
        cout << profile_title << " y-intercept = " << function_tmp->GetParameter(0) << " +/- " << function_tmp->GetParError(0) << endl;
    }
    temp_graph->GetXaxis()->SetTitle(x_axis_name);
    temp_graph->GetXaxis()->SetRangeUser(0., 20.0);
    temp_graph->GetYaxis()->SetTitle(y_axis_name);
    temp_graph->SetTitle(profile_title);
    temp_graph->SetMarkerStyle(21);
    temp_graph->Draw("AP");
    function_tmp->Draw("same");
    if (fit)
    {
        TLegend *leg = new TLegend(0.1, 0.7, 0.48, 0.9);
        leg->AddEntry((TObject *)0, TString::Format("y-intercept = %.3e +/- %.3e", function_tmp->GetParameter(0), function_tmp->GetParError(0)).Data(), "");
        leg->Draw();
    }
    c1.Write();

    std::vector<float> temp_vector_TProfile;
    temp_vector_TProfile.push_back(function_tmp->GetParameter(0));
    temp_vector_TProfile.push_back(function_tmp->GetParError(0));

    temp_graph->Write(TString::Format("%s_%s", profile1_name->GetName(), profile2_name->GetName()));

    return temp_vector_TProfile;
}

std::vector<float> plot_combined(TProfile *profile1_name, TGraphErrors *profile2_name, const char *x_axis_name, const char *y_axis_name, const char *profile_title, bool fit, int rebin)
{
    float new_x[1000] = {0.}, new_y[1000] = {0.}, new_x_err[1000] = {0.}, new_y_err[1000] = {0.};
    int N_bins = profile1_name->GetNbinsX() / 5;
    float new_x_numerator[1000] = {0.}, new_x_denomenator[1000] = {0.}, new_y_numerator[1000] = {0.}, new_y_denomenator[1000] = {0.};

    for (int i = 0; i < N_bins / rebin; i++)
    {

        for (int j = 0; j < rebin; j++)
        {
            if (profile1_name->GetBinError(rebin * i + 1 + j) == 0 || profile2_name->GetErrorY(rebin * i + j) == 0)
                continue;
            new_x_numerator[i] += profile1_name->GetBinContent(rebin * i + 1 + j) / pow(profile1_name->GetBinError(rebin * i + 1 + j), 2);
            new_x_denomenator[i] += 1.0 / pow(profile1_name->GetBinError(rebin * i + 1 + j), 2);
            new_y_numerator[i] += profile2_name->GetY()[(rebin * i  + j)] / pow(profile2_name->GetErrorY(rebin * i + j), 2);
            new_y_denomenator[i] += 1.0 / pow(profile2_name->GetErrorY(rebin * i + j), 2);
        }

        new_x[i] = new_x_numerator[i] / new_x_denomenator[i];
        new_x_err[i] = sqrt(1.0 / new_x_denomenator[i]);
        new_y[i] = new_y_numerator[i] / new_y_denomenator[i];
        new_y_err[i] = sqrt(1.0 / new_y_denomenator[i]);

        // cout << "new_x[" << i << "] = " << new_x[i] << endl;
        // cout << "new_y[" << i << "] = " << new_y[i] << endl;
    }

    TGraphErrors *temp_graph = new TGraphErrors(N_bins / rebin, new_x, new_y, new_x_err, new_y_err);

    outputfile->cd();
    TCanvas c1(profile_title, profile_title, 1500, 800);
    TF1 *function_tmp = new TF1("function_tmp", "[0] + [1] * x + [2] * x * x", 0, 0.08);
    if (fit)
    {
        temp_graph->Fit(function_tmp, "QNR");
        cout << setprecision(3);
        cout << profile_title << " y-intercept = " << function_tmp->GetParameter(0) << " +/- " << function_tmp->GetParError(0) << endl;
    }
    temp_graph->GetXaxis()->SetTitle(x_axis_name);
    temp_graph->GetXaxis()->SetRangeUser(0., 20.0);
    temp_graph->GetYaxis()->SetTitle(y_axis_name);
    temp_graph->SetTitle(profile_title);
    temp_graph->SetMarkerStyle(21);
    temp_graph->Draw("AP");
    function_tmp->Draw("same");
    if (fit)
    {
        TLegend *leg = new TLegend(0.1, 0.7, 0.48, 0.9);
        leg->AddEntry((TObject *)0, TString::Format("y-intercept = %.3e +/- %.3e", function_tmp->GetParameter(0), function_tmp->GetParError(0)).Data(), "");
        leg->Draw();
    }
    c1.Write();

    std::vector<float> temp_vector_TProfile;
    temp_vector_TProfile.push_back(function_tmp->GetParameter(0));
    temp_vector_TProfile.push_back(function_tmp->GetParError(0));

    temp_graph->Write(TString::Format("%s_%s", profile1_name->GetName(), profile2_name->GetName()));

    return temp_vector_TProfile;
}

std::vector<float> plot_TProfile(TGraphErrors *temp_profile, const char *x_axis_name, const char *y_axis_name, const char *profile_title, bool fit)
{

    outputfile->cd();
    TCanvas c1(profile_title, profile_title, 1500, 800);
    TF1 *function_tmp = new TF1("function_tmp", "[0] + [1] * x + [2] * x * x", 0, 5);
    if (fit)
    {
        temp_profile->Fit(function_tmp, "QNR");
        cout << setprecision(3);
        cout << profile_title << " y-intercept = " << function_tmp->GetParameter(0) << " +/- " << function_tmp->GetParError(0) << endl;
    }
    temp_profile->GetXaxis()->SetTitle(x_axis_name);
    temp_profile->GetYaxis()->SetTitle(y_axis_name);
    temp_profile->SetTitle(profile_title);
    temp_profile->Draw();
    function_tmp->Draw("same");
    if (fit)
    {
        TLegend *leg = new TLegend(0.1, 0.7, 0.48, 0.9);
        leg->AddEntry((TObject *)0, TString::Format("y-intercept = %.3e +/- %.3e", function_tmp->GetParameter(0), function_tmp->GetParError(0)).Data(), "");
        leg->Draw();
    }
    c1.Write();

    std::vector<float> temp_vector_TProfile;
    temp_vector_TProfile.push_back(function_tmp->GetParameter(0));
    temp_vector_TProfile.push_back(function_tmp->GetParError(0));

    return temp_vector_TProfile;
}

void plotting_TGraph_Helper(const char *title, const char *x_title, const char *y_title, int npts, std::vector<TString> x_labels, float *y1, float *y1_err, float *y2, float *y2_err, float *y3, float *y3_err, float *y4, float *y4_err, float *y5, float *y5_err, float *y6, float *y6_err)
{
    TGraphErrors *tmp_graph1 = new TGraphErrors(4);
    for (Int_t i = 0; i < 4; i++)
    {                                           // Loop over all entries
        tmp_graph1->SetPoint(i, i + 1., y1[i]); // Set The point itself
        tmp_graph1->SetPointError(i, 0, y1_err[i]);
        tmp_graph1->GetXaxis()->SetBinLabel(tmp_graph1->GetXaxis()->FindBin(i + 1.), x_labels.at(i).Data()); // Find out which bin on the x-axis the point corresponds to and set the bin label
        cout << "x_labels.at(i).Data() = " << x_labels.at(i).Data() << endl;
    }

    if (y2 != NULL)
    {
        TGraphErrors *tmp_graph2 = new TGraphErrors(4);
        for (Int_t i = 0; i < 4; i++)
        {                                           // Loop over all entries
            tmp_graph2->SetPoint(i, i + 1., y2[i]); // Set The point itself
            tmp_graph2->SetPointError(i, 0, y2_err[i]);
            tmp_graph2->GetXaxis()->SetBinLabel(tmp_graph2->GetXaxis()->FindBin(i + 1.), x_labels.at(i).Data()); // Find out which bin on the x-axis the point corresponds to and set the bin label
        }
    }
    if (y3 != NULL)
    {
        TGraphErrors *tmp_graph3 = new TGraphErrors(4);
        for (Int_t i = 0; i < 4; i++)
        {                                           // Loop over all entries
            tmp_graph3->SetPoint(i, i + 1., y3[i]); // Set The point itself
            tmp_graph3->SetPointError(i, 0, y3_err[i]);
            tmp_graph3->GetXaxis()->SetBinLabel(tmp_graph3->GetXaxis()->FindBin(i + 1.), x_labels.at(i).Data()); // Find out which bin on the x-axis the point corresponds to and set the bin label
        }
    }
    if (y4 != NULL)
    {
        TGraphErrors *tmp_graph4 = new TGraphErrors(4);
        for (Int_t i = 0; i < 4; i++)
        {                                           // Loop over all entries
            tmp_graph4->SetPoint(i, i + 1., y4[i]); // Set The point itself
            tmp_graph4->SetPointError(i, 0, y4_err[i]);
            tmp_graph4->GetXaxis()->SetBinLabel(tmp_graph4->GetXaxis()->FindBin(i + 1.), x_labels.at(i).Data()); // Find out which bin on the x-axis the point corresponds to and set the bin label
        }
    }
    if (y5 != NULL)
    {
        TGraphErrors *tmp_graph5 = new TGraphErrors(4);
        for (Int_t i = 0; i < 4; i++)
        {                                           // Loop over all entries
            tmp_graph5->SetPoint(i, i + 1., y5[i]); // Set The point itself
            tmp_graph5->SetPointError(i, 0, y5_err[i]);
            tmp_graph5->GetXaxis()->SetBinLabel(tmp_graph5->GetXaxis()->FindBin(i + 1.), x_labels.at(i).Data()); // Find out which bin on the x-axis the point corresponds to and set the bin label
        }
    }
    if (y6 != NULL)
    {
        TGraphErrors *tmp_graph6 = new TGraphErrors(4);
        for (Int_t i = 0; i < 4; i++)
        {                                           // Loop over all entries
            tmp_graph6->SetPoint(i, i + 1., y6[i]); // Set The point itself
            tmp_graph6->SetPointError(i, 0, y6_err[i]);
            tmp_graph6->GetXaxis()->SetBinLabel(tmp_graph6->GetXaxis()->FindBin(i + 1.), x_labels.at(i).Data()); // Find out which bin on the x-axis the point corresponds to and set the bin label
        }
    }

    outputfile->cd();
    TCanvas c1(title, title, 800, 800);
    tmp_graph1->SetMarkerStyle(kFullCircle);
    tmp_graph1->SetMarkerColor(kBlack);
    tmp_graph1->SetLineColor(kBlack);
    tmp_graph1->SetTitle(title);
    tmp_graph1->GetXaxis()->SetTitle(x_title);
    tmp_graph1->GetYaxis()->SetTitle(y_title);
    tmp_graph1->SetMarkerSize(3);
    tmp_graph1->Draw("AP");
    if (y2 != NULL)
    {
        // cout << "plotting second" << endl;
        tmp_graph2->SetMarkerStyle(kFullCircle);
        tmp_graph2->SetMarkerColor(kRed);
        tmp_graph2->SetLineColor(kRed);
        tmp_graph2->SetMarkerSize(3);
        tmp_graph2->Draw("SAMES P");
    }
    if (y3 != NULL)
    {
        tmp_graph3->SetMarkerStyle(kFullCircle);
        tmp_graph3->SetMarkerColor(kBlue);
        tmp_graph3->SetLineColor(kBlue);
        tmp_graph3->SetMarkerSize(3);
        tmp_graph3->Draw("SAMES P");
    }
    if (y4 != NULL)
    {
        tmp_graph4->SetMarkerStyle(kFullCircle);
        tmp_graph4->SetMarkerColor(kGreen);
        tmp_graph4->SetLineColor(kGreen);
        tmp_graph4->SetMarkerSize(3);
        tmp_graph4->Draw("SAMES P");
    }
    if (y5 != NULL)
    {
        tmp_graph5->SetMarkerStyle(kFullSquare);
        tmp_graph5->SetMarkerColor(kViolet);
        tmp_graph5->SetLineColor(kViolet);
        tmp_graph5->SetMarkerSize(3);
        tmp_graph5->Draw("SAMES P");
    }
    if (y6 != NULL)
    {
        tmp_graph6->SetMarkerStyle(kFullSquare);
        tmp_graph6->SetMarkerColor(44);
        tmp_graph6->SetLineColor(44);
        tmp_graph6->SetMarkerSize(3);
        tmp_graph6->Draw("SAMES P");
    }
    TLine l1(0.75, 0, 4.25, 0);
    l1.Draw();
    c1.Write();
}