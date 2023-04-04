#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include "TH1.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"



void z_vs_phi()
{

// Set up a range of theta values to plot
    const int Ntheta = 180;
    const double theta_min = 0;
    const double theta_max = 180;
    const double dtheta = (theta_max - theta_min) / (Ntheta - 1);
    double theta[Ntheta];
    for (int i = 0; i < Ntheta; i++) {
        theta[i] = theta_min + i*dtheta;
    }


// Create canvases and graph to display the graphs and histogram
//TCanvas *c = new TCanvas("c", "Functions Graph", 800, 600);
TGraph *gr1 = new TGraph();
TGraph *gr2 = new TGraph();
TF1 *f1 = new TF1();
TGraph2D* graph = new TGraph2D();


std::ifstream infile("data/data190_200.dat"); // Open input file

//reading data from input file
std::vector<double> xpos, ypos, zpos, time, ratio, phi, rho, max_bin_coord_1, max_bin_coord_2, zz;
Double_t xi, yi, zi, ti;

    while (infile >> xi >> yi >> zi >> ti) {
        xpos.push_back(xi);
        ypos.push_back(yi);
        zpos.push_back(zi);
	time.push_back(ti);

        
        }


    for (int i =0; i<zpos.size(); i++){
    
    graph->SetPoint(graph->GetN(), xpos[i], ypos[i], zpos[i]);

    double aa = atan2(ypos[i], xpos[i]);
    phi.push_back(aa);
    gr1->SetPoint(gr1->GetN(), phi[i], zpos[i]);
    gr1->SetMarkerSize(2);

 	}



// Create a multi-graph to hold the family of curves
    TMultiGraph* mg = new TMultiGraph("mg", "rho vs theta");

    // Loop over z and phi values and generate a curve for each one
    for (int i = 0; i < zpos.size(); i++) {
        TGraph* graph = new TGraph(Ntheta);
        for (int j = 0; j < Ntheta; j++) {
            double rho = phi[i]*(180/3.14)*cos(theta[j]*TMath::DegToRad()) + zpos[i]*sin(theta[j]*TMath::DegToRad());
            graph->SetPoint(j, theta[j], rho);
        }
        graph->SetLineColor(i+1);
        graph->SetLineWidth(2);
        mg->Add(graph);
    }


    
    // Create a 2D histogram to bin the curves into
    const int Nbins = 180;
    TH2D* hist = new TH2D("hist", "binning", Nbins, theta_min, theta_max, 1000, -500, 500);

    // Loop over each curve in the multi-graph and add it to the histogram
    for (int i = 0; i < mg->GetListOfGraphs()->GetEntries(); i++) {
        TGraph* graph = (TGraph*)mg->GetListOfGraphs()->At(i);
        for (int j = 0; j < graph->GetN(); j++) {
            double x, y;
            graph->GetPoint(j, x, y);
            hist->Fill(x, y);
        }
    }


    // Set up a canvas and draw the histogram
    TCanvas* canvas1 = new TCanvas("canvas1", "binning");
    hist->Draw("COLZ");



    // Create a graph for the original z vs phi points
    TGraph* data = new TGraph(zpos.size(), &phi[0], &zpos[0]);
    data->SetMarkerStyle(20);
    data->SetMarkerSize(0.5);
    data->SetLineColor(kBlack);
    data->SetLineWidth(2);

    // Set up a canvas and draw the multi-graph and data graph on it
    TCanvas* canvas = new TCanvas("canvas", "graph");
    mg->Draw("AL");
    

    // Set up axis labels and legend
    mg->GetXaxis()->SetTitle("#theta");
    mg->GetYaxis()->SetTitle("#rho");
    canvas->BuildLegend(0.8, 0.8, 1.0, 1.0);

    canvas->Update();


//=====================================

//to find the bin counts of each bins
  double maxCounts = -1;
  int xBin, yBin;
  std::vector<double> x_coord, y_coord, bin_counts;
  for (int i = 1; i <= hist->GetNbinsX(); i++) {
    for (int j = 1; j <= hist->GetNbinsY(); j++) {
      double counts = hist->GetBinContent(i, j);
        maxCounts = counts;
        xBin = i;
        yBin = j;

	x_coord.push_back(hist->GetXaxis()->GetBinCenter(xBin));
	y_coord.push_back(hist->GetYaxis()->GetBinCenter(yBin));
	bin_counts.push_back(maxCounts);
    }
  }



//To sort the bin counts using multimap(taking bin count as the key) 

    std::multimap<double, std::pair<double, double>> m;

    for (int i = 0; i < x_coord.size(); i++) {
        m.insert({bin_counts[i], {x_coord[i], y_coord[i]}});
    }

    int count = 0;
    for (auto it = m.rbegin(); it != m.rend(); it++) {
        if (count == 10) break;
        cout << it->second.first << "	" << it->second.second << "	" << it->first << endl;
	max_bin_coord_1.push_back(it->second.first);
	max_bin_coord_2.push_back(it->second.second);
        count++;
    }





    // Create a multi-graph to hold the curves(reversing the hough transformation)
    TMultiGraph* mg1 = new TMultiGraph("mg1", "z vs phi");

    // Loop over z and phi values and generate a curve for each one
    for (int i = 0; i < 10/*max_bin_coord_1.size()*/; i++) {

	Double_t aaa=1/sin(max_bin_coord_1[i]*TMath::DegToRad());
	Double_t bbb=1/tan(max_bin_coord_1[i]*TMath::DegToRad());
	Double_t rho1=max_bin_coord_2[i];

        TGraph* graph1 = new TGraph(phi.size());
        for (int j = 0; j < phi.size(); j++) {
            Double_t zz = aaa*rho1 - phi[j]*bbb;
            graph1->SetPoint(j, phi[j], zz);
        }
        graph1->SetLineColor(i+1);
        graph1->SetLineWidth(4);
        mg1->Add(graph1);
    }

    // Set up a canvas and draw the multi-graph and data graph on it
    TCanvas* canvas11 = new TCanvas("canvas11", "Curves11 and Data");

    //legend
    mg1->GetXaxis()->SetTitle("#phi");
    mg1->GetYaxis()->SetTitle("z");
    canvas11->BuildLegend(0.8, 0.8, 1.0, 1.0);

    canvas11->Update();
 
    mg1->Add(gr1, "AP*");
    gr1->SetMarkerSize(3);
    mg1->GetYaxis()->SetRangeUser(-240, 240);
    mg1->Draw("ALP");


//Drawing and plotting part




}
