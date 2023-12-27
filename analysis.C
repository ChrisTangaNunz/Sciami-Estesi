#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TGraphErrors.h"  


void plotHistogram(const char* filePath, TFile* outputFile, const char* outputFileName){
    // Leggi l'ultimo valore dal file per determinare totalDuration
    std::ifstream lastValueFile(filePath);
    double lastValue = 0.0;

    if (lastValueFile.is_open()) {
        lastValueFile >> std::ws; // Ignora gli spazi bianchi iniziali
        while (lastValueFile >> lastValue) {
            // Continua a leggere fino all'ultimo valore
        }
        lastValueFile.close();
    } else {
        std::cerr << "Errore nell'apertura del file: " << filePath << std::endl;
        return;
    }

    const double binWidth = 10.0;
    const double totalDuration = lastValue; // Utilizza l'ultimo valore come totalDuration
    const int numBins = static_cast<int>(totalDuration / binWidth);

    // Crea un istogramma 
    TH1F *histogram = new TH1F("histogram", "Istogramma dei tempi;Tempi [s];Eventi", numBins, 0, totalDuration);

    // Imposta lo stile dell'istogramma
    histogram->SetLineColor(kBlue);
    histogram->SetFillColor(kBlue);
    histogram->SetLineWidth(2);

    // Leggi i dati dal file e riempi l'istogramma
    std::ifstream inputFile(filePath);
    if (!inputFile.is_open()) {
        std::cerr << "Errore nell'apertura del file: " << filePath << std::endl;
        return;
    }

    double time;
    while (inputFile >> time) {
        histogram->Fill(time);
    }

    // Chiudi il file
    inputFile.close();

    // Crea un canvas e disegna l'istogramma
    TCanvas *canvas = new TCanvas("canvas", "Istogramma", 800, 600);
    histogram->Draw();

    canvas->SaveAs(Form("/home/chris/SciamiEstesi/23_11_2023/plots/%s_1D.png", outputFileName));
    outputFile->cd();
    histogram->Write();
    // Calcola il rate e crea un istogramma 2D
    TH2F *histogram2D = new TH2F("histogram2D", "Istogramma 2D;Tempi [s];Rate [Hz]", 100, 0, totalDuration, 100, 0, 30);
    TH1F *rateHistogram = new TH1F("rateHistogram", "Istogramma del Rate;Conteggi;Occorrenze", 100, 0, 200);

    // Crea vettori per i dati del grafico a errori
    std::vector<double> times;
    std::vector<double> rates;
    std::vector<double> errors;


    // Riempi l'istogramma 2D con il rate
    for (int i = 1; i <= histogram->GetNbinsX(); ++i) {
        double binContent = histogram->GetBinContent(i);
        double binWidth = histogram->GetBinWidth(i);
        double rate = binContent / binWidth;
        double error = std::sqrt(binContent) / binWidth;  



        histogram2D->Fill(histogram->GetBinCenter(i), rate);

        times.push_back(histogram->GetBinCenter(i));
        rates.push_back(rate);
        errors.push_back(error);
        // Riempi l'istogramma del rate con il binContent (non normalizzato)
        rateHistogram->Fill(binContent);
    }

    

    // Crea un secondo canvas e disegna l'istogramma 2D
    TCanvas *canvas2D = new TCanvas("canvas2D", "Istogramma 2D", 800, 600);
    histogram2D->Draw("colz");

    // Imposta il rate come etichetta dell'asse Z
    histogram2D->SetZTitle("Rate [Hz]");

    // Posiziona il box statistics in modo che non copra la scala dell'asse Z
    gPad->Update();
    TPaveStats *stats = (TPaveStats*)histogram2D->FindObject("stats");
    stats->SetX1NDC(0.78);
    stats->SetX2NDC(0.98);
    stats->SetY1NDC(0.78);
    stats->SetY2NDC(0.98);

    canvas2D->SaveAs(Form("/home/chris/SciamiEstesi/23_11_2023/plots/%s_2D.png", outputFileName));

    // Salva i risultati su un file ROOT
    outputFile->cd();
    histogram2D->Write();

    // Crea un secondo canvas per il grafico a errori
    TCanvas *canvasError = new TCanvas("canvasError", "Grafico Rate vs Tempi", 800, 600);

    // Crea il grafico a errori
    TGraphErrors *graphErrors = new TGraphErrors(times.size(), &times[0], &rates[0], 0, &errors[0]);

    // Imposta gli stili del grafico a errori
    graphErrors->SetMarkerStyle(20);
    graphErrors->SetMarkerSize(1.5);
    graphErrors->SetTitle("Rate vs Tempi;Tempi [s];Rate [Hz]");

    // Disegna il grafico a errori
    graphErrors->Draw("AP");

    canvasError->SaveAs(Form("/home/chris/SciamiEstesi/23_11_2023/plots/%s_Rate.png", outputFileName));

    // Fai il fit della distribuzione poissoniana
    TF1 *poissonFit = new TF1("poissonFit", "[0]*TMath::Poisson(x, [1])", 0, 200);
    poissonFit->SetParameters(10, 10);  // Parametri iniziali per il fit

    // Esegui il fit dell'istogramma del rate
    rateHistogram->Fit("poissonFit", "L");

    // Estrai i parametri del fit
    double fitMean = poissonFit->GetParameter(1);
    double fitMeanError = poissonFit->GetParError(1);
    // Stampa i risultati del fit
    std::cout << "Fit Mean: " << fitMean << " +/- " << fitMeanError << std::endl;
    // Salva i risultati del fit su un file ROOT
    outputFile->cd();
    poissonFit->Write("poissonFit");

    // Crea un terzo canvas per l'istogramma del rate
    TCanvas *canvasRate = new TCanvas("canvasRate", "Istogramma del Rate", 800, 600);

    // Disegna l'istogramma del rate
    rateHistogram->Draw();

    canvasRate->SaveAs(Form("/home/chris/SciamiEstesi/23_11_2023/plots/%s_RateHistogram.png", outputFileName));

    // Salva i risultati su un file ROOT
    outputFile->cd();
    rateHistogram->Write();


}

void calculateTimeDifferences(const char* filePath, std::vector<double>& timeDifferences) {
    std::ifstream inputFile(filePath);
    if (!inputFile.is_open()) {
        std::cerr << "Errore nell'apertura del file: " << filePath << std::endl;
        return;
    }

    double time;
    std::vector<double> eventTimes;

    while (inputFile >> time) {
        eventTimes.push_back(time);
    }

    // Calcola le differenze di tempo e riempi il vettore timeDifferences
    for (size_t i = 1; i < eventTimes.size(); ++i) {
        double deltaTime = eventTimes[i] - eventTimes[i - 1];
        timeDifferences.push_back(deltaTime);
    }

    inputFile.close();
}

void plotTimeDifferences(const char* filePath, TFile* outputFile, const char* outputFileName) {
    std::vector<double> timeDifferences;
    calculateTimeDifferences(filePath, timeDifferences);

    // Crea un istogramma delle differenze di tempo
    TH1F *timeDifferenceHistogram = new TH1F("timeDifferenceHistogram", "Differenze di tempo;Differenze di tempo [s];Occorrenze", 100, 0, 1.5);

    // Riempie l'istogramma con le differenze di tempo
    for (const auto& deltaTime : timeDifferences) {
        timeDifferenceHistogram->Fill(deltaTime);
    }

    // Aggiorna le opzioni di visualizzazione della box di statistica
    timeDifferenceHistogram->SetStats(kTRUE);  // Abilita la visualizzazione della box di statistica

    // Aggiungi il chi quadro alla box di statistica
    gStyle->SetOptFit(1111);  // Visualizza chi quadro nella box di statistica


    // Fai il fit dell'istogramma delle differenze di tempo con una funzione esponenziale
    TF1 *exponentialFit = new TF1("exponentialFit", "[0]*exp(-[1]*x)", 0, 1.5);
    exponentialFit->SetParameters(100, 1);  // Parametri iniziali per il fit

    timeDifferenceHistogram->Fit("exponentialFit", "L");

    // Estrai i parametri del fit
    double fitAmplitude = exponentialFit->GetParameter(0);
    double fitDecay = exponentialFit->GetParameter(1);

    // Stampa i risultati del fit
    std::cout << "Fit Amplitude: " << fitAmplitude << std::endl;
    std::cout << "Fit Decay: " << fitDecay << std::endl;
    



    // Salva i risultati del fit su un file ROOT
    outputFile->cd();
    exponentialFit->Write("exponentialFit");
    // Crea un canvas e disegna l'istogramma delle differenze di tempo
    TCanvas *canvasTimeDifferences = new TCanvas("canvasTimeDifferences", "Istogramma delle differenze di tempo", 800, 600);
    timeDifferenceHistogram->Draw();

    canvasTimeDifferences->SaveAs(Form("/home/chris/SciamiEstesi/23_11_2023/plots/%s_TimeDifferences.png", outputFileName));
    // Crea un secondo canvas con scala logaritmica per l'asse delle y
    TCanvas *canvasLogScale = new TCanvas("canvasLogScale", "Istogramma con scala logaritmica", 800, 600);
    gPad->SetLogy();
    timeDifferenceHistogram->Draw();

    // Fai il fit con una funzione esponenziale
    TF1 *expFitLogScale = new TF1("expFitLogScale", "[0]*exp(-[1]*x)", 0, 3);
    expFitLogScale->SetParameters(fitAmplitude, fitDecay);  // Usa i parametri del fit precedente come iniziali

    timeDifferenceHistogram->Fit("expFitLogScale", "L");

    // Estrai i parametri del nuovo fit
    double fitAmplitudeLogScale = expFitLogScale->GetParameter(0);
    double fitDecayLogScale = expFitLogScale->GetParameter(1);

    // Stampa i risultati del nuovo fit
    std::cout << "Fit Amplitude (log scale): " << fitAmplitudeLogScale << std::endl;
    std::cout << "Fit Decay (log scale): " << fitDecayLogScale << std::endl;




    // Salva i risultati del nuovo fit su un file ROOT
    outputFile->cd();
    expFitLogScale->Write("expFitLogScale");

    canvasLogScale->SaveAs(Form("/home/chris/SciamiEstesi/23_11_2023/plots/%s_TimeDifferences_LogScale.png", outputFileName));

    // Salva i risultati su un file ROOT
    outputFile->cd();
    timeDifferenceHistogram->Write();

}

void plotTimeSeriesEfficiencies(const char* filePath, TFile* outputFile, const char* outputFileName) {

}


int analysis() {
    TFile *outputFile = new TFile("output_combined.root", "RECREATE");

    const char* filePath1 = "/home/chris/SciamiEstesi/23_11_2023/data/clean_data/tempi_corretti_canale_1_231123.txt";
    const char* filePath2 = "/home/chris/SciamiEstesi/23_11_2023/data/clean_data/tempi_corretti_canale_2_231123.txt";
    const char* filePath3 = "/home/chris/SciamiEstesi/23_11_2023/data/clean_data/tempi_corretti_canale_3_231123.txt";

    plotHistogram(filePath1, outputFile, "output_canale_1");
    plotHistogram(filePath2, outputFile, "output_canale_2");
    plotHistogram(filePath3, outputFile, "output_canale_3");
    // Aggiungi la creazione di istogrammi delle differenze di tempo per ciascun file
    plotTimeDifferences(filePath1, outputFile, "output_canale_1");
    plotTimeDifferences(filePath2, outputFile, "output_canale_2");
    plotTimeDifferences(filePath3, outputFile, "output_canale_3");
    outputFile->Close();

    return 0;
}
