#include "TFile.h"
#include "TCanvas.h"
#include "TProfile.h"

using namespace std;

TFile *input_file, *output_file;

void plot_plots_tgt(const char *name, const char *plot1, const char *plot2 = NULL, const char *plot3 = NULL);

void plot_cme_results(){
    input_file = new TFile("output.root");
    output_file = new TFile("drawn.root", "recreate");

    // plot_plots_tgt("v2_q2_Canvas", "all_v2_Q2_rho_v2_Q2", "all_v2_Q2_global_v2_Q2");
    // plot_plots_tgt("v2_q2_Canvas", "all_v2pair_Q2_rho_v2_Q2", "all_v2pair_Q2_global_v2_Q2");
    // plot_plots_tgt("v2_q2_Canvas", "all_v2_pairQ2_rho_v2_pairQ2", "all_v2_pairQ2_global_v2_pairQ2");
    // plot_plots_tgt("v2_q2_Canvas", "all_v2pair_pairQ2_rho_v2_pairQ2", "all_v2pair_pairQ2_global_v2_pairQ2");

    // plot_plots_tgt("v2_pairq2_Canvas", "all_v2_Q2_primordial_Gamma_os_Q2", "all_v2_Q2_both_prim_Gamma_os_Q2");
    // plot_plots_tgt("v2_pairq2_Canvas", "all_v2pair_Q2_primordial_Gamma_os_Q2", "all_v2pair_Q2_both_prim_Gamma_os_Q2");
    // plot_plots_tgt("v2_pairq2_Canvas", "all_v2_pairQ2_primordial_Gamma_os_pairQ2", "all_v2_pairQ2_both_prim_Gamma_os_pairQ2");
    // plot_plots_tgt("v2_pairq2_Canvas", "all_v2pair_pairQ2_primordial_Gamma_os_pairQ2", "all_v2pair_pairQ2_both_prim_Gamma_os_pairQ2");

    // plot_plots_tgt("v2_pairq2_Canvas", "all_v2_Q2_both_glob_Gamma_os_Q2", "all_v2_Q2_one_prim_Gamma_os_Q2");
    // plot_plots_tgt("v2_pairq2_Canvas", "all_v2pair_Q2_both_glob_Gamma_os_Q2", "all_v2pair_Q2_one_prim_Gamma_os_Q2");
    // plot_plots_tgt("v2_pairq2_Canvas", "all_v2_pairQ2_both_glob_Gamma_os_pairQ2", "all_v2_pairQ2_one_prim_Gamma_os_pairQ2");
    // plot_plots_tgt("v2_pairq2_Canvas", "all_v2pair_pairQ2_both_glob_Gamma_os_pairQ2", "all_v2pair_pairQ2_one_prim_Gamma_os_pairQ2");

    // plot_plots_tgt("v2_pairq2_Canvas", "all_v2_Q2_primordial_Gamma_os_Q2", "all_v2_Q2_both_glob_Gamma_os_Q2");
    plot_plots_tgt("v2_pairq2_Canvas", "all_v2pair_Q2_primordial_Gamma_os_Q2", "all_v2pair_Q2_both_glob_Gamma_os_Q2");
    // plot_plots_tgt("v2_pairq2_Canvas", "all_v2_pairQ2_primordial_Gamma_os_pairQ2", "all_v2_pairQ2_both_glob_Gamma_os_pairQ2");
    // plot_plots_tgt("v2_pairq2_Canvas", "all_v2pair_pairQ2_primordial_Gamma_os_pairQ2", "all_v2pair_pairQ2_both_glob_Gamma_os_pairQ2");

    // plot_plots_tgt("v2_pairq2_Canvas", "rho_v2_pt_profile_profile", "primordial_v2_pt_profile_profile", "global_v2_pt_profile_profile");
    

    input_file->Close();
    output_file->Close();
}

void plot_plots_tgt(const char *name, const char *plot1, const char *plot2, const char *plot3){

    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);
    gStyle->SetOptDate(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetPadColor(0);
    gStyle->SetPadBorderSize(0);
    gStyle->SetCanvasBorderSize(0);
    gStyle->SetPadGridX(false);
    gStyle->SetPadGridY(false);
    gStyle->SetLabelSize(.05, "X");
    gStyle->SetLabelSize(.05, "Y");
    gStyle->SetTitleSize(.06, "X");
    gStyle->SetTitleSize(.06, "Y");

    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(kWhite);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetTitleBorderSize(0);

    TGraphErrors *temp_profile1 = (TGraphErrors *) input_file->Get(plot1)->Clone();
    if(plot2 != NULL) {TGraphErrors *temp_profile2 = (TGraphErrors *) input_file->Get(plot2)->Clone();}
    if(plot3 != NULL) TGraphErrors *temp_profile3 = (TGraphErrors *) input_file->Get(plot3)->Clone();

    // TProfile *temp_profile1 = (TProfile *) input_file->Get(plot1)->Clone();
    // if(plot2 != NULL) {TProfile *temp_profile2 = (TProfile *) input_file->Get(plot2)->Clone();}
    // if(plot3 != NULL) TProfile *temp_profile3 = (TProfile *) input_file->Get(plot3)->Clone();

    cout << "after loading plots" << endl;
    float max_x = 0.12;
    TF1 *function_tmp = new TF1("function_tmp", "[0] + [1] * x + [2] * x * x", 0, max_x);
    TF1 *function_tmp2 = new TF1("function_tmp2", "[0] + [1] * x + [2] * x * x", 0, max_x);
    TF1 *function_tmp3 = new TF1("function_tmp3", "[0] + [1] * x + [2] * x * x", 0, max_x);
    TLegend *leg = new TLegend(0.1, 0.7, 0.48, 0.9);
    
    cout << "after functions" << endl;

    output_file->cd();
    TCanvas *c1 = new TCanvas(name, name, 40, 40, 740, 580);
    c1->cd();

    cout << "after TCanvas" << endl;

    // c1->SetTopMargin(0.06);
    c1->SetTopMargin(0.08);
    c1->SetRightMargin(0.03);
    c1->SetBottomMargin(0.17);
    c1->SetLeftMargin(0.11);
    c1->Draw();

    gPad->SetGridx(0);
    gPad->SetGridy(0);
    gPad->SetTickx();
    gPad->SetTicky();

    cout << "after gPad" << endl;

    temp_profile1->Fit(function_tmp, "QRN");
    // TLegendEntry* l1 = leg->AddEntry((TObject *)0, TString::Format("All Pions v_{2}(q^{2}) y-intercept = %.3e +/- %.3e %%", function_tmp->GetParameter(0) * 100.0, function_tmp->GetParError(0) * 100.0).Data(), "");
    // TLegendEntry* l1 = leg->AddEntry((TObject *)0, TString::Format("All Pions v_{2}(q^{2}_{pair}) y-intercept = %.3e +/- %.3e %%", function_tmp->GetParameter(0) * 100.0, function_tmp->GetParError(0) * 100.0).Data(), "");
    // l1->SetTextColor(kBlack);

    cout << "fitted" << endl;
    
    // temp_profile1->SetTitle("Pions v_{2}(q^{2})");
    // temp_profile1->SetTitle("Pions v_{2}(q^{2}_{pair})");
    
    // temp_profile1->GetXaxis()->SetTitle("Pion-Pair q^{2}_{pair}");
    // temp_profile1->GetXaxis()->SetTitle("Single Pion q^{2}");
    // temp_profile1->GetXaxis()->SetTitle("single v_{2}{single q^{2}_{2}}");
    temp_profile1->GetXaxis()->SetTitle("pair v_{2}{single q^{2}_{2}}");
    // temp_profile1->GetXaxis()->SetTitle("single v_{2}{pair q^{2}_{2}}");
    // temp_profile1->GetXaxis()->SetTitle("pair v_{2}{pair q^{2}_{2}}");
    // temp_profile1->GetXaxis()->SetTitle("pT (GeV/c)");
    // temp_profile1->GetXaxis()->SetLimits(0.0, 0.22);
    
    temp_profile1->GetXaxis()->SetLimits(0.0, max_x+0.001);
    temp_profile1->GetYaxis()->SetRangeUser(-0.0006, 0.0035);
    // temp_profile1->GetYaxis()->SetRangeUser(-0.09, 0.25);
    TGaxis::SetMaxDigits(3);
    temp_profile1->GetYaxis()->SetTitle("#Delta#gamma_{112}  ");
    // temp_profile1->GetYaxis()->SetTitle("v_{2}");
    
    temp_profile1->GetYaxis()->SetTitleOffset(0.67);
    temp_profile1->GetYaxis()->SetTitleSize(0.075);
    temp_profile1->GetXaxis()->SetTitleSize(0.08);
    temp_profile1->GetXaxis()->SetTitleOffset(0.90);

    
    temp_profile1->GetXaxis()->SetNdivisions(6);
    temp_profile1->GetYaxis()->SetNdivisions(605);
    double lsize = temp_profile1->GetXaxis()->GetLabelSize();
    // temp_profile1->GetYaxis()->SetLabelSize(lsize * 1.1);
    // temp_profile1->GetXaxis()->SetLabelSize(lsize * 1.2);
    temp_profile1->GetYaxis()->SetLabelSize(0.06);
    temp_profile1->GetXaxis()->SetLabelSize(0.06);

    temp_profile1->SetMarkerStyle(kOpenCircle);
    int marker_size = 2.5;
    temp_profile1->SetMarkerSize(marker_size);
    temp_profile1->SetMarkerColor(4);
    temp_profile1->SetLineColor(4);
    temp_profile1->SetFillColor(4);
    temp_profile1->SetLineStyle(1);
    temp_profile1->SetLineWidth(2);
    
    temp_profile1->Draw("ape1");
    function_tmp->SetLineColor(4);
    function_tmp->Draw("same");

    // TLatex *tex = new TLatex(0.01, 0.18, Form("y-int(#rho) = (%3.2f#pm%3.2f)#times10^{-3}", 1000 * function_tmp->GetParameter(0), 1000 * function_tmp->GetParError(0)));
    TLatex *tex = new TLatex(0.037, 0.0005, Form("y-int(all pairs) = (%3.2f#pm%3.2f)#times10^{-5}", 100000 * function_tmp->GetParameter(0), 100000 * function_tmp->GetParError(0)));
    // TLatex *tex = new TLatex(0.003, 0.000025, Form("y-int(#Delta#gamma_{112, 2D}) = (%3.2f#pm%3.2f)#times10^{-5}", 100000 * function_tmp->GetParameter(0), 100000 * function_tmp->GetParError(0)));
    // TLatex *tex = new TLatex(0.003, 0.000025, Form("y-int(#rho v_{2}(p_{T})) = (%3.2f#pm%3.2f)#times10^{-5}", 100000 * function_tmp->GetParameter(0), 100000 * function_tmp->GetParError(0)));
    // TLatex *tex = new TLatex(0.003, 0.000025, Form("#rho v_{2}(p_{T})", 100000 * function_tmp->GetParameter(0), 100000 * function_tmp->GetParError(0)));
    tex->SetTextSize(0.065);
    tex->SetTextColor(4);
    tex->Draw();

    TLatex *letter = new TLatex(0.112, -0.0002, "(d)");
    // TLatex *letter = new TLatex(0.135, 0.01, "(d)");
    letter->SetTextSize(0.08);
    letter->SetTextColor(kBlack);
    letter->Draw();

    TLatex *title_toy = new TLatex(0.045, 0.0026, "Toy Model");
    // TLatex *title_toy = new TLatex(0.05, 0.2, "Toy Model");
    title_toy->SetTextSize(0.095);
    title_toy->SetTextColor(kBlack);
    // title_toy->Draw();

    cout << "after plot1" << endl;

    if(plot2 != NULL){
        temp_profile2->Fit(function_tmp2, "QRN");
        temp_profile2->SetMarkerStyle(kOpenDiamond);
        temp_profile2->SetMarkerSize(marker_size);
        temp_profile2->SetMarkerColor(2);
        temp_profile2->SetLineColor(2);
        temp_profile2->SetFillColor(2);
        temp_profile2->SetLineStyle(1);
        temp_profile2->SetLineWidth(2);
        temp_profile2->Draw("pe1 same");
        // temp_profile2->SetMarkerColor(kRed);
        // temp_profile2->SetLineColor(kRed);
        // temp_profile2->Draw("same");
        // TLegendEntry* l2 = leg->AddEntry((TObject *)0, TString::Format("Primordial Pions v_{2}(q^{2}) y-intercept = %.3e +/- %.3e %%", function_tmp2->GetParameter(0) * 100.0, function_tmp2->GetParError(0) * 100.0).Data(), "");
        // TLegendEntry* l2 = leg->AddEntry((TObject *)0, TString::Format("Primordial Pions v_{2}(q^{2}_{pair}) y-intercept = %.3e +/- %.3e %%", function_tmp2->GetParameter(0) * 100.0, function_tmp2->GetParError(0) * 100.0).Data(), "");
        // l2->SetTextColor(kRed);
        function_tmp2->SetLineColor(2);
        function_tmp2->Draw("same");
        // TLatex *tex2 = new TLatex(0.01, -0.07, Form("y-int(decay pion) = (%3.2f#pm%3.2f)#times10^{-3}", 1000 * function_tmp2->GetParameter(0), 1000 * function_tmp2->GetParError(0)));
        // TLatex *tex2 = new TLatex(0.003, 0.000025, Form("y-int(#Delta#gamma_{112, 2P}) = (%3.2f#pm%3.2f)#times10^{-5}", 100000 * function_tmp2->GetParameter(0), 100000 * function_tmp2->GetParError(0)));
        TLatex *tex2 = new TLatex(0.026, 0.0008, Form("y-int(both decay) = (%3.2f#pm%3.2f)#times10^{-5}", 100000 * function_tmp2->GetParameter(0), 100000 * function_tmp2->GetParError(0)));
        // TLatex *tex2 = new TLatex(0.003, 0.000025, Form("y-int(#Delta#gamma_{112, 1P1D}) = (%3.2f#pm%3.2f)#times10^{-5}", 100000 * function_tmp2->GetParameter(0), 100000 * function_tmp2->GetParError(0)));
        // TLatex *tex2 = new TLatex(0.003, 0.000025, Form("primordial pions v_{2}(p_{T})", 100000 * function_tmp2->GetParameter(0), 100000 * function_tmp2->GetParError(0)));
        tex2->SetTextSize(0.065);
        tex2->SetTextColor(2);
        tex2->Draw();
    }

    if(plot3 != NULL){
        temp_profile3->Fit(function_tmp3, "QRN");
        temp_profile3->SetMarkerStyle(kOpenSquare);
        temp_profile3->SetMarkerSize(marker_size);
        temp_profile3->SetMarkerColor(6);
        temp_profile3->SetLineColor(6);
        temp_profile3->SetFillColor(6);
        temp_profile3->SetLineStyle(1);
        temp_profile3->SetLineWidth(2);
        temp_profile3->Draw("pe1 same");

        // TLegendEntry* l3 = leg->AddEntry((TObject *)0, TString::Format("Decay Pions v_{2}(q^{2}) y-intercept = %.3e +/- %.3e %%", function_tmp3->GetParameter(0) * 100.0, function_tmp3->GetParError(0) * 100.0).Data(), "");
        // TLegendEntry* l3 = leg->AddEntry((TObject *)0, TString::Format("Decay Pions v_{2}(q^{2}_{pair}) y-intercept = %.3e +/- %.3e %%", function_tmp3->GetParameter(0) * 100.0, function_tmp3->GetParError(0) * 100.0).Data(), "");
        // l3->SetTextColor(kBlue);
        TLatex *tex3 = new TLatex(0.003, 0.000025, Form("decay pions v_{2}(p_{T})", 100000 * function_tmp3->GetParameter(0), 100000 * function_tmp3->GetParError(0)));
        tex3->SetTextSize(0.055);
        tex3->SetTextColor(6);
        tex3->Draw();
        function_tmp3->SetLineColor(6);
        // function_tmp3->Draw("same");
    }
    leg->SetBorderSize(0);
    // leg->Draw();

    TLine *line1 = new TLine(0, 0, max_x, 0);
    // TLine *line1 = new TLine(0, 0, 10, 0);
    line1->SetLineColor(kBlack);
    line1->Draw("same");

    c1->Write();
}