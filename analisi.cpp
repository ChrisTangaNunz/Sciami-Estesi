#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

std::vector<double> loadFiles(const char* filePath){
    // Vettore per memorizzare i tempi
    std::vector<double> times; 

    // Apre il file in modalità di lettura
    std::ifstream file(filePath);

    // Verifica se il file è aperto con successo
    if (!file.is_open()) {
        std::cerr << "Errore nell'apertura del file: " << filePath << std::endl;
        return times;
    }

    double time;
    while(file >> time){
        times.push_back(time);
    }

    file.close();

    return times;
}



void createErrorGraph(const TH1F* histogram, const char* outputFileName, const char* channel) {
    // Crea vettori per i dati del grafico a errori
    std::vector<double> timesVec;
    std::vector<double> rates;
    std::vector<double> errors;

    for (int i = 1; i <= histogram->GetNbinsX(); ++i) {
        double binContent = histogram->GetBinContent(i);
        double binWidth = histogram->GetBinWidth(i);
        double rate = binContent / binWidth;
        double error = std::sqrt(binContent) / binWidth;

        timesVec.push_back(histogram->GetBinCenter(i));
        rates.push_back(rate);
        errors.push_back(error);
    }

    // Crea il grafico a errori
    TGraphErrors graphErrors(timesVec.size(), &timesVec[0], &rates[0], 0, &errors[0]);

    // Imposta gli stili del grafico a errori
    graphErrors.SetMarkerStyle(20);
    graphErrors.SetMarkerSize(1.5);
    graphErrors.SetTitle("Rate vs Tempi;Tempi [s];Rate [Hz]");

    // Crea un canvas per il grafico a errori
    TCanvas canvasError("canvasError", "Grafico Rate vs Tempi", 800, 600);
    graphErrors.Draw("AP");

    // Salva il canvas
    canvasError.SaveAs(Form("/home/chris/SciamiEstesi/23_11_2023/plots/channel_%s/%s_Rate.png", channel,outputFileName));
}


void plot_histogram2D(const TH1F* histogram, const char* outputFileName, const char* channel){

    std::vector<double> times2D;
    std::vector<double> rates2D;
    std::vector<double> errors2D;
    TH2F *histogram2D = new TH2F("histogram2D", "Istogramma 2D;Tempi [s];Rate [Hz]", 100, 0, histogram->GetXaxis()->GetBinUpEdge(histogram->GetNbinsX()), 100, 0, 30);


    // Riempi l'istogramma 2D con il rate
    for (int i = 1; i <= histogram->GetNbinsX(); ++i) {
        double binContent = histogram->GetBinContent(i);
        double binWidth = histogram->GetBinWidth(i);
        double rate = binContent / binWidth;

        histogram2D->Fill(histogram->GetBinCenter(i), rate);

        times2D.push_back(histogram->GetBinCenter(i));
        rates2D.push_back(rate);
  

    }
    TCanvas *canvas2D = new TCanvas("canvas2D", "Istogramma 2D", 800, 600);
    histogram2D->Draw("colz");



    // Posiziona il box statistics in modo che non copra la scala dell'asse Z
    gPad->Update();
    TPaveStats *stats = (TPaveStats*)histogram2D->FindObject("stats");
    stats->SetX1NDC(0.68);
    stats->SetX2NDC(0.88);
    stats->SetY1NDC(0.78);
    stats->SetY2NDC(0.98);

    canvas2D->SaveAs(Form("/home/chris/SciamiEstesi/23_11_2023/plots/channel_%s/%s_2D.png",channel, outputFileName));


}

void plotHistogram_Rate(const TH1F* histogram,  const char* outputFileName, const char* channel) {

    TH1F *rateHistogram = new TH1F("rateHistogram", "Istogramma del Rate;Conteggi;Occorrenze", 100, 0, 200);
    for (int i = 1; i <= histogram->GetNbinsX(); ++i) {
        double binContent = histogram->GetBinContent(i);
        rateHistogram->Fill(binContent);
    }
    // Fai il fit della distribuzione poissoniana
    TF1 *poissonFit = new TF1("poissonFit", "[0]*TMath::Poisson(x, [1])", 0, 200);
    poissonFit->SetParameters(10, 10);  // Parametri iniziali per il fit

    // Esegui il fit dell'istogramma del rate
    rateHistogram->Fit("poissonFit", "L");
    // Crea un terzo canvas per l'istogramma del rate
    TCanvas *canvasRate = new TCanvas("canvasRate", "Istogramma del Rate", 800, 600);
    // Disegna l'istogramma del rate
    rateHistogram->Draw();

    canvasRate->SaveAs(Form("/home/chris/SciamiEstesi/23_11_2023/plots/channel_%s/%s_RateHistogram.png", channel, outputFileName));


    // Estrai i parametri del fit
    double fitMean = poissonFit->GetParameter(1);
    double fitMeanError = poissonFit->GetParError(1);
    // Stampa i risultati del fit
    std::cout << "Fit Mean: " << fitMean << " +/- " << fitMeanError << std::endl;

}


void plot_histogram(const char* filePath, const char* outputFileName, const char* channel) {
    std::vector<double> time = loadFiles(filePath);
    if (time.empty()) {
        std::cerr << "Nessun dato da caricare. Uscita." << std::endl;
        return;
    }

    double lastValue = time.back();
    const double binWidth = 10;
    const int numBins = static_cast<int>(std::ceil(lastValue / binWidth));

    TH1F histogram("histogram", "Histogram", numBins, 0., lastValue);
    // Imposta lo stile dell'istogramma
    histogram.SetLineColor(kBlue);
    histogram.SetFillColor(kBlue);
    histogram.SetLineWidth(2);

    for (size_t i = 0; i < time.size(); ++i) {
        histogram.Fill(time[i]);
    }

    // Creazione di un canvas per visualizzare l'istogramma
    TCanvas canvasH1("canvas", "Canvas", 800, 600);
    histogram.Draw();

    // Mostra il canvas
    canvasH1.Draw();
    canvasH1.SaveAs(Form("/home/chris/SciamiEstesi/23_11_2023/plots/channel_%s/%s_1D.png", channel, outputFileName));

    createErrorGraph(&histogram, outputFileName, channel);
    plot_histogram2D(&histogram, outputFileName, channel);
    plotHistogram_Rate(&histogram, outputFileName, channel);
}


void calculateTimeDifferences(const char* filePath, std::vector<double>& timeDifferences) {
    std::vector<double> time = loadFiles(filePath);
    if (time.empty()) {
        std::cerr << "Nessun dato da caricare. Uscita." << std::endl;
        return;
    }


    // Calcola le differenze di tempo e riempi il vettore timeDifferences
    for (size_t i = 1; i < time.size(); ++i) {
        double deltaTime = time[i] - time[i - 1];
        timeDifferences.push_back(deltaTime);
    }

}



void plotTimeDifferences(const char* filePath,  const char* outputFileName, const char* channel) {
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
    
    // Crea un canvas e disegna l'istogramma delle differenze di tempo
    TCanvas *canvasTimeDifferences = new TCanvas("canvasTimeDifferences", "Istogramma delle differenze di tempo", 800, 600);
    timeDifferenceHistogram->Draw();

    canvasTimeDifferences->SaveAs(Form("/home/chris/SciamiEstesi/23_11_2023/plots/channel_%s/%s_TimeDifferences.png", channel, outputFileName));
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


    canvasLogScale->SaveAs(Form("/home/chris/SciamiEstesi/23_11_2023/plots/channel_%s/%s_TimeDifferences_LogScale.png", channel, outputFileName));

}

int analisi(){
    //CANALE1
    const char* filePath1 = "/home/chris/SciamiEstesi/23_11_2023/data/clean_data/tempi_corretti_canale_1_231123.txt"; 
    plot_histogram(filePath1, "output_canale1", "1");
    plotTimeDifferences(filePath1, "output_canale_1", "1");
    //CANALE2
    const char* filePath2 = "/home/chris/SciamiEstesi/23_11_2023/data/clean_data/tempi_corretti_canale_2_231123.txt"; 
    plot_histogram(filePath2, "output_canale2", "2");
    plotTimeDifferences(filePath2, "output_canale_2", "2");
    //CANALE3
    const char* filePath3 = "/home/chris/SciamiEstesi/23_11_2023/data/clean_data/tempi_corretti_canale_3_231123.txt"; 
    plot_histogram(filePath3, "output_canale3", "3");
    plotTimeDifferences(filePath3, "output_canale_3", "3");


    return 0;
}
