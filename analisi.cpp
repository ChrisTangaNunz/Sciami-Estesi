#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include <TPaveStats.h>
#include <algorithm> 


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



void createErrorGraph(const TH1F* histogram, const char* channel, double efficiency) {
    // Crea vettori per i dati del grafico a errori
    std::vector<double> timesVec;
    std::vector<double> bins_Content;
    std::vector<double> rates;
    std::vector<double> errors;
    std::vector<double> timeErrors;


    for (int i = 1; i <= histogram->GetNbinsX(); ++i) {
        double binContent = histogram->GetBinContent(i);
        double binWidth = histogram->GetBinWidth(i);
        double rate = binContent / binWidth;
        double error = std::sqrt(binContent) / binWidth;
        double err_x = binWidth / 2;

        timesVec.push_back(histogram->GetBinCenter(i)/3600);
        bins_Content.push_back(binContent);
        rates.push_back(rate);
        errors.push_back(error);
        timeErrors.push_back(err_x/3600);

    }
    
    // Moltiplica i valori sull'asse y per 1000 durante la creazione del grafico
    /*

    for (size_t i = 0; i < rates.size(); ++i) {
        rates[i] *= 1000;
        errors[i] *= 1000;
    }

    */  
    //CORREZIONE CON EFFICIENZA
        for (size_t i = 0; i < rates.size(); ++i) {
        rates[i] *= efficiency;
        errors[i] *= efficiency;
    }
    // Crea il grafico a errori
    TGraphErrors graphErrors(timesVec.size(), &timesVec[0], &rates[0], &timeErrors[0], &errors[0]);
    // Esegui il fit con una retta
    //TF1* fitFunction = new TF1("fitFunction", "pol1", timesVec[0], timesVec.back());  // Usa una funzione lineare (pol1)
    //graphErrors.Fit(fitFunction, "R");

    // Imposta gli stili del grafico a errori
    graphErrors.SetMarkerStyle(20);
    graphErrors.SetMarkerSize(0.01);
    graphErrors.SetTitle("Rate vs Tempi; Tempi [h] in bins di 3600 s; Rate [Hz]"); //Aggiusta mHz se moltiplichi i dato per un fattore 1000
    // Calcola la media dei valori sull'asse y del grafico
    double meanY = 0.0;
    double cont_bin = 0.0;
    double sumSquaredErrors = 0.0;

    for (size_t i = 0; i < rates.size(); ++i) {
        meanY += rates[i];
        cont_bin += bins_Content[i];
        double squaredError = errors[i] * errors[i];
        sumSquaredErrors += squaredError;
    }
    meanY /= rates.size();
    //double err_meanY=std::sqrt(cont_bin)/3600;
    double err_meanY = std::sqrt(sumSquaredErrors)/ errors.size();


    // Crea un canvas per il grafico a errori
    TCanvas canvasError("canvasError", "Grafico Rate vs Tempi", 800, 600);
    canvasError.cd();  // Sposta il contesto all'interno del canvas
    //fitFunction->SetParameters(0.0, meanY); // coefficiente angolare e intercetta

    graphErrors.Draw("AP");
    // Aggiungi una linea tratteggiata corrispondente al rate medio
    TLine* meanLine = new TLine(graphErrors.GetXaxis()->GetXmin(), meanY, graphErrors.GetXaxis()->GetXmax(), meanY);
    meanLine->SetLineStyle(2);  // Imposta la linea come tratteggiata
    meanLine->SetLineColor(kRed);  // Imposta il colore della linea
    meanLine->Draw("same");  // Disegna la linea sullo stesso grafico


    // Aggiungi il box statistico
    TPaveStats* stats = new TPaveStats(0.7, 0.8, 0.9, 0.9, "NDC");
    stats->SetName("stats");
    stats->SetBorderSize(1);
    stats->SetFillColor(0);
    stats->SetTextAlign(12);
    stats->SetTextSize(0.025);
    stats->AddText(Form(("Rate medio [Hz]:")));
    stats->AddText(Form("%.3f $\\pm$ %.3f$", meanY, err_meanY));    //ricorda di moltiplicare *1000 quando fai le coincidenze tra telescopi per portare il rate in mHz 
    //stats->AddText(Form("Intercept: %.3f", fitFunction->GetParameter(1)));    
    //stats->AddText(Form("Slope: %.3f", fitFunction->GetParameter(0)));    

    canvasError.Modified();  // Necessario per rendere effettive le modifiche al canvas
    canvasError.Update();    // Aggiorna il canvas
    stats->Draw("same");     // Disegna il box statistico sopra il grafico

    // Salva il canvas
    canvasError.SaveAs(Form("/home/chris/SciamiEstesi/23_11_2023/plots/channel_%s/Rate_corrected_with_efficiency.png", channel));
}

void plot_histogram2D(const TH1F* histogram, const char* channel){

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
    // Posiziona il box statistics in modo che non copra la scala dell'asse Z per l'istogramma 2D
    gPad->Update();
    TPaveStats *stats2D = (TPaveStats*)histogram2D->FindObject("stats");
    stats2D->SetX1NDC(0.68);
    stats2D->SetX2NDC(0.88);
    stats2D->SetY1NDC(0.78);
    stats2D->SetY2NDC(0.98);

    canvas2D->SaveAs(Form("/home/chris/SciamiEstesi/23_11_2023/plots/channel_%s/histogram2D.png", channel));



}

void plotHistogram_Rate(const TH1F* histogram, const char* channel) {
  
    TH1F *rateHistogram = new TH1F("rateHistogram", "Istogramma del Rate;Conteggi in #Deltat di 10 s;Occorrenze", 100, 0, 200);
    double sum = 0.0;
    double sumSquared = 0.0;

    for (int i = 1; i <= histogram->GetNbinsX(); ++i) {
        double binContent = histogram->GetBinContent(i);
        //binContent *= efficiency; Ricordatii di aggiungere , double efficiency tra gli argomenenti della funzione
        sum += binContent;
        sumSquared += binContent * binContent;
        rateHistogram->Fill(binContent);
    }
    // Calcola la media dei dati 
    double mean = sum / rateHistogram->GetEntries();
    //Calcola l'errore sulla media dei dati
    double errorOnMean = mean / sqrt(rateHistogram->GetEntries());
    // Calcola la deviazione standard dei dati
    double stdev = sqrt((sumSquared / rateHistogram->GetEntries() - mean * mean));
    //Calcola l'errore della deviazione standard dei dati
    double err_stdev = stdev / sqrt(2 * rateHistogram->GetEntries());


    //Imposta correttamente il range dell'asse x
    rateHistogram->GetXaxis()->SetRangeUser(100, 205);

    // Fai il fit della distribuzione poissoniana
    TF1 *poissonFit = new TF1("poissonFit", "[0]*TMath::Poisson(x, [1])", 0, 200);
    poissonFit->SetParameters(61310, 45);  // Parametri iniziali per il fit
    // Rinomina i parametri del fit
    poissonFit->SetParName(0, "N");
    poissonFit->SetParName(1, "#mu");
    // Esegui il fit dell'istogramma del rate
    rateHistogram->Fit("poissonFit", "L");

    // Crea un terzo canvas per l'istogramma del rate
    TCanvas *canvasRate = new TCanvas("canvasRate", "Istogramma del Rate", 800, 600);
    // Disegna l'istogramma del rate
    rateHistogram->Draw();
    // Abilita la visualizzazione delle informazioni del fit nel box statistico
    gPad->Update(); // Aggiorna il pad prima di ottenere il riferimento alla funzione di fit
    TPaveStats *stats = (TPaveStats*)rateHistogram->FindObject("stats");
    stats->SetOptFit(1111);
    // Imposta la precisione desiderata per la stampa dei parametri del fit
    stats->SetOptStat(1111); // 1111 consente di visualizzare la funzione, i parametri e le statistiche con le opzioni desiderate
    gStyle->SetOptStat(1111); // Assicurati che la formattazione sia applicata anche a livello di stile globale

    // Imposta la precisione dei parametri del fit nel formato "N.F" (N cifre totali, F cifre decimali)
    stats->SetFitFormat("5.4f"); // Modifica 5 e 4 secondo le tue esigenze di precisione

    canvasRate->SaveAs(Form("/home/chris/SciamiEstesi/23_11_2023/plots/channel_%s/RateHistogram.png", channel));


    // Estrai i parametri del fit
    double fitMean = poissonFit->GetParameter(1);
    double fitMeanError = poissonFit->GetParError(1);
    double fitVariance = sqrt(fitMean);
    double erroritVariance = fitMeanError/(2*(sqrt(fitMean)));

    std::cout<< "Histogram RATE CHANNEL: "<< channel << std::endl;
    // Stampa i risultati del fit.
    std::cout << "Data Mean: " << mean << " +/- " << errorOnMean << std::endl;
    std::cout << "Fit Mean #mu: " << fitMean << " +/- " << fitMeanError << std::endl;
    std::cout << "Fit Variance: " << fitVariance << " +/- " << erroritVariance << std::endl;
    std::cout << "Standard Deviation: " << stdev << " +/- " << err_stdev << std::endl;

}


void plot_histogram(const char* filePath, const char* channel, double efficiency) {
    std::vector<double> time = loadFiles(filePath);
    if (time.empty()) {
        std::cerr << "Nessun dato da caricare. Uscita." << std::endl;
        return;
    }
/*
    // Convert times from seconds to hours
    for (size_t i = 0; i < time.size(); ++i) {
        time[i] /= 3600.0;
    }
*/
    double lastValue = time.back();
    const double binWidth = 10.0;
    const int numBins = static_cast<int>(std::ceil(lastValue / binWidth));

    TH1F histogram("histogram", "Istogramma dei tempi in cui un evento viene registrato", numBins, 0., lastValue);
    // Imposta lo stile dell'istogramma
    //histogram.SetLineColor(kBlue);
    //histogram.SetFillColor(kBlue);
    histogram.SetLineWidth(2);
    histogram.SetFillStyle(4050);


    for (size_t i = 0; i < time.size(); ++i) {
        histogram.Fill(time[i]);
    }
    histogram.GetYaxis()->SetTitle("Eventi");
    histogram.GetXaxis()->SetTitle("Tempi [h] in bins di 1 h");

    // Creazione di un canvas per visualizzare l'istogramma
    TCanvas canvasH1("canvas", "Canvas", 800, 600);
    histogram.Draw();

    // Mostra il canvas
    canvasH1.Draw();
    canvasH1.SaveAs(Form("/home/chris/SciamiEstesi/23_11_2023/plots/channel_%s/histogram1D.png", channel));

    createErrorGraph(&histogram, channel, efficiency);
    plot_histogram2D(&histogram, channel);
    plotHistogram_Rate(&histogram, channel); //Ricordati di aggiungere efficency tra gli argomenti della funzione nel caso
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



void plotTimeDifferences(const char* filePath, const char* channel) {
    std::vector<double> timeDifferences;
    calculateTimeDifferences(filePath, timeDifferences);

    // Crea un istogramma delle differenze di tempo
    TH1F *timeDifferenceHistogram = new TH1F("timeDifferenceHistogram", "Differenze di tempo;Differenze di tempo [s];Occorrenze", 100, 0, 1.2);

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
    exponentialFit->SetParName(0, "N");
    exponentialFit->SetParName(1, "#lambda");
    timeDifferenceHistogram->Fit("exponentialFit", "L");

    // Estrai i parametri del fit
    double fitAmplitude = exponentialFit->GetParameter(0);
    double fitDecay = exponentialFit->GetParameter(1);
    double fitDecayError = exponentialFit->GetParError(1);


    // Stampa i risultati del fit
    std::cout << "Fit Amplitude: " << fitAmplitude << std::endl;
    std::cout << "Fit Decay: " << fitDecay << std::endl;
    
    // Crea un canvas e disegna l'istogramma delle differenze di tempo
    TCanvas *canvasTimeDifferences = new TCanvas("canvasTimeDifferences", "Istogramma delle differenze di tempo", 800, 600);
    timeDifferenceHistogram->Draw();

    // Abilita la visualizzazione delle informazioni del fit nel box statistico
    gPad->Update(); // Aggiorna il pad prima di ottenere il riferimento alla funzione di fit
    TPaveStats *statsNonLog = (TPaveStats*)timeDifferenceHistogram->FindObject("stats");
    statsNonLog->SetOptFit(1111);
    // Imposta la precisione desiderata per la stampa dei parametri del fit
    statsNonLog->SetOptStat(1111); // 1111 consente di visualizzare la funzione, i parametri e le statistiche con le opzioni desiderate
    gStyle->SetOptStat(1111); // Assicurati che la formattazione sia applicata anche a livello di stile globale

    // Imposta la precisione dei parametri del fit nel formato "N.F" (N cifre totali, F cifre decimali)
    statsNonLog->SetFitFormat("5.4f"); // Modifica 5 e 4 secondo le tue esigenze di precisione

    canvasTimeDifferences->SaveAs(Form("/home/chris/SciamiEstesi/23_11_2023/plots/channel_%s/TimeDifferences.png", channel));
    // Crea un secondo canvas con scala logaritmica per l'asse delle y
    TCanvas *canvasLogScale = new TCanvas("canvasLogScale", "Istogramma con scala logaritmica", 800, 600);
    gPad->SetLogy();
    timeDifferenceHistogram->Draw();
    // Abilita la visualizzazione delle informazioni del fit nel box statistico
    gPad->Update(); // Aggiorna il pad prima di ottenere il riferimento alla funzione di fit
    TPaveStats *stats = (TPaveStats*)timeDifferenceHistogram->FindObject("stats");
    stats->SetOptFit(1111);
    // Imposta la precisione desiderata per la stampa dei parametri del fit
    stats->SetOptStat(1111); // 1111 consente di visualizzare la funzione, i parametri e le statistiche con le opzioni desiderate
    gStyle->SetOptStat(1111); // Assicurati che la formattazione sia applicata anche a livello di stile globale

    // Imposta la precisione dei parametri del fit nel formato "N.F" (N cifre totali, F cifre decimali)
    stats->SetFitFormat("5.4f"); // Modifica 5 e 4 secondo le tue esigenze di precisione

    // Fai il fit con una funzione esponenziale
    TF1 *expFitLogScale = new TF1("expFitLogScale", "[0]*exp(-[1]*x)", 0, 3);
    expFitLogScale->SetParameters(fitAmplitude, fitDecay);  // Usa i parametri del fit precedente come iniziali
    expFitLogScale->SetParName(0, "N");
    expFitLogScale->SetParName(1, "#lambda");
    timeDifferenceHistogram->Fit("expFitLogScale", "L");

    // Estrai i parametri del nuovo fit
    double fitAmplitudeLogScale = expFitLogScale->GetParameter(0);
    double fitDecayLogScale = expFitLogScale->GetParameter(1);
    double fitDecayLogScaleError = expFitLogScale->GetParError(1);
    double tauFit = 1 / fitDecayLogScale;
    double error_tauFit  = fitDecayLogScaleError / std::pow(fitDecayLogScale, 2);

    // Stampa i risultati del nuovo fit
    std::cout << "Risultati Distribuzione temporale del Delta t tra due eventi successi nel CANALE: " << channel << std::endl;
    std::cout << "Fit Decay (log scale) #lambda: " << fitDecayLogScale << " +/- " << fitDecayLogScaleError << std::endl;
    std::cout << "Tau fit: " << tauFit << " +/- " << error_tauFit << std::endl;



    canvasLogScale->SaveAs(Form("/home/chris/SciamiEstesi/23_11_2023/plots/channel_%s/TimeDifferences_LogScale.png", channel));

}



void calculateAndPlotDiffCoinc(const char* fileName, const char* channel) {
    // Apri il file
    std::ifstream file(fileName);

    if (!file.is_open()) {
        std::cerr << "Errore nell'apertura del file " << fileName << std::endl;
        return;
    }

    // Leggi i dati da entrambe le colonne
    std::vector<double> col1;
    std::vector<double> col2;
    double value1, value2;

    while (file >> value1 >> value2) {
        col1.push_back(value1);
        col2.push_back(value2);
    }

    file.close();

    // Calcola la differenza tra le colonne
    std::vector<double> differences;
    for (size_t i = 0; i < col1.size(); ++i) {
    differences.push_back(std::abs(col1[i] - col2[i])*(std::pow(10, 9)));
    }
    /*
    std::cout << "Differenze:" << std::endl;
    for (size_t i = 0; i < differences.size(); ++i) {
        std::cout << differences[i] << std::endl;
    }
    */
    // Crea un istogramma delle differenze
    TH1F* histogramCoinc = new TH1F("histogramCoinc", "Istogramma coincidenze", 100, 0, 100);
    
    histogramCoinc->SetTitle("Istogramma coincidenze due telescopi; Differenza di Tempi tra eventi successivi [ns]; Occorrenze");
    // Riempie l'istogramma con le differenze
    for (size_t i = 0; i < differences.size(); ++i) {
        histogramCoinc->Fill(differences[i]);
    }

    // Crea un canvas
    TCanvas* canvas = new TCanvas("canvas", "Canvas", 800, 600);

    // Disegna l'istogramma
    histogramCoinc->Draw("HIST");

    // Salva il canvas
    canvas->SaveAs(Form("/home/chris/SciamiEstesi/23_11_2023/plots/channel_%s/DiffCoincDoppie.png", channel));

    delete histogramCoinc;
    delete canvas;

}


void plot_RatevsParameter(const char* filePath, const char* parameterName, const char* channel) {
    // Disabilita la casella delle statistiche di default
    //gStyle->SetOptFit(0);
    gStyle->SetOptFit(0111); 
    // Regola le dimensioni della casella delle statistiche
    gStyle->SetStatW(0.1);
    gStyle->SetStatH(0.1);
   
    // Carica i dati dal file
    std::ifstream infile(filePath);

    if (!infile.is_open()) {
        std::cerr << "Errore nell'apertura del file: " << filePath << std::endl;
        return;
    }

    std::vector<double> parameter_values, rates, y_errors;

    double parameter_value, rate;
    while (infile >> parameter_value >> rate) {
        parameter_values.push_back(parameter_value);
        rates.push_back(rate);
        y_errors.push_back(std::sqrt(rate*300)/300);
    }

    infile.close();
    // Creazione di un TCanvas per il plot
    TCanvas *canvas = new TCanvas("Canvas", Form("Rates vs %s", parameterName), 800, 600);
    canvas->SetGrid();
    

    // Creazione di un TGraphErrors per lo scatter plot con rate NON corretti
    TGraphErrors *graph = new TGraphErrors(parameter_values.size(), &parameter_values[0], &rates[0], 0, &y_errors[0]);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.25);
    graph->SetMarkerColor(kBlue);
    graph->SetTitle(Form("Rates vs %s", parameterName));
    graph->GetYaxis()->SetTitle("Rate [Hz]");
    // Imposta il titolo sull'asse x in base al tipo di parametro
    if (strcmp(parameterName, "temeperatura") == 0) {
    graph->GetXaxis()->SetTitle(Form("Temperature [%cC]", 0xB0));
    } else if (strcmp(parameterName, "pressione") == 0) {
        graph->GetXaxis()->SetTitle("Pressure [mbar]");
    } else if (strcmp(parameterName, "umidità") == 0) {
        graph->GetXaxis()->SetTitle("Humidity [%]");
    } else {
        // Se il parametro non corrisponde a nessuno dei precedenti, imposta un titolo generico
        graph->GetXaxis()->SetTitle(Form("%s", parameterName));
    }
    graph->GetYaxis()->SetRangeUser(14.4, 16.3);


    // Disegna lo scatter plot
    graph->Draw("AP");

    // Creazione di una funzione lineare per il fit
    TF1 *linearFit = new TF1("linearFit", "[0] + [1]*x", *std::min_element(parameter_values.begin(), parameter_values.end()), *std::max_element(parameter_values.begin(), parameter_values.end()));
    // Imposta i nomi dei parametri
    linearFit->SetParName(0, "Intercept");
    linearFit->SetParName(1, "Slope");
    // Esegui il fit
    graph->Fit("linearFit", "R");
    // Estrai i parametri del fit
    double intercept = linearFit->GetParameter(0);
    double slope = linearFit->GetParameter(1);


    // Salva il canvas come immagine
    canvas->SaveAs(Form("/home/chris/SciamiEstesi/23_11_2023/plots/parametri_atmosferici/%s/channel_%s/channel%s.png", parameterName, channel, channel));
    // Rilascia la memoria
    delete linearFit;
    delete graph;
    delete canvas;

}



void plot_CorrectedRatevsParameter(const char* filePath, const char* parameterName, const char* channel, double efficiency) {
    // Disabilita la casella delle statistiche di default
    //gStyle->SetOptFit(0);
    gStyle->SetOptFit(0111); 
    // Regola le dimensioni della casella delle statistiche
    gStyle->SetStatW(0.1);
    gStyle->SetStatH(0.1);
   
    // Carica i dati dal file
    std::ifstream infile(filePath);

    if (!infile.is_open()) {
        std::cerr << "Errore nell'apertura del file: " << filePath << std::endl;
        return;
    }

    std::vector<double> parameter_values, rates, y_errors;

    double parameter_value, rate;
    while (infile >> parameter_value >> rate) {
        parameter_values.push_back(parameter_value);
        rates.push_back(rate);
        y_errors.push_back(std::sqrt(rate*300)/300);
    }

    infile.close();
    //CORREZIONE CON EFFICIENZA
    std::vector<double> rates_corrected, errors_corrected;

 
        for (size_t i = 0; i < rates.size(); ++i) {
        double corrected_rate = rates[i] * efficiency;
        double corrected_error = y_errors[i] * efficiency;
        rates_corrected.push_back(corrected_rate);
        errors_corrected.push_back(corrected_error);
    }

    // Creazione di un TCanvas per il plot
    TCanvas *canvas = new TCanvas("Canvas", Form("Rates vs %s", parameterName), 800, 600);
    canvas->SetGrid();
    

    // Creazione di un TGraphErrors per lo scatter plot con rate NON corretti
    TGraphErrors *graph = new TGraphErrors(parameter_values.size(), &parameter_values[0], &rates_corrected[0], 0, &errors_corrected[0]);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.25);
    graph->SetMarkerColor(kBlue);
    graph->SetTitle(Form("Rates corretti vs %s", parameterName));
    graph->GetYaxis()->SetTitle("Rate [Hz]");
    // Imposta il titolo sull'asse x in base al tipo di parametro
    if (strcmp(parameterName, "temeperatura") == 0) {
    graph->GetXaxis()->SetTitle(Form("Temperature [%cC]", 0xB0));
    } else if (strcmp(parameterName, "pressione") == 0) {
        graph->GetXaxis()->SetTitle("Pressure [mbar]");
    } else if (strcmp(parameterName, "umidità") == 0) {
        graph->GetXaxis()->SetTitle("Humidity [%]");
    } else {
        // Se il parametro non corrisponde a nessuno dei precedenti, imposta un titolo generico
        graph->GetXaxis()->SetTitle(Form("%s", parameterName));
    }
    graph->GetYaxis()->SetRangeUser(4.1, 5);


    // Disegna lo scatter plot
    graph->Draw("AP");

    // Creazione di una funzione lineare per il fit
    TF1 *linearFit = new TF1("linearFit", "[0] + [1]*x", *std::min_element(parameter_values.begin(), parameter_values.end()), *std::max_element(parameter_values.begin(), parameter_values.end()));
    // Imposta i nomi dei parametri
    linearFit->SetParName(0, "Intercept");
    linearFit->SetParName(1, "Slope");
    // Esegui il fit
    graph->Fit("linearFit", "R");
    // Estrai i parametri del fit
    double intercept = linearFit->GetParameter(0);
    double slope = linearFit->GetParameter(1);
    double error_slope = linearFit->GetParError(1);

    // Calcola la distanza sigma dallo zero
    double sigma_distance = slope / error_slope; 
    std::cout << "PARAMETRO: " << parameterName << std::endl;
    std::cout << "Distanza sigma dallo zero: " << std::fixed << std::setprecision(3) << sigma_distance << std::endl;

    // Salva il canvas come immagine
    canvas->SaveAs(Form("/home/chris/SciamiEstesi/23_11_2023/plots/parametri_atmosferici/%s/channel_%s/channel_CorrectedRates%s.png", parameterName, channel, channel));
    // Rilascia la memoria
    delete linearFit;
    delete graph;
    delete canvas;


}



struct TelescopeInfo {
    const char* filePath;
    const char* channel;
    double efficiency;
};

struct AtmosphericParameterInfo {
    const char* filePath;
    const char* parameterName;
    const char* channel;
    double efficiency;
};

void plot_AtmosphericParameter(const AtmosphericParameterInfo& info) {
    plot_CorrectedRatevsParameter(info.filePath, info.parameterName, info.channel, info.efficiency);
    //plot_RatevsParameter(info.filePath, info.parameterName, info.channel);

}


int analisi(){
    TelescopeInfo telescopes[] = {
        {"/home/chris/SciamiEstesi/23_11_2023/data/clean_data/tempi_corretti_canale_1_231123.txt", "1", 0.35},
        {"/home/chris/SciamiEstesi/23_11_2023/data/clean_data/tempi_corretti_canale_2_231123.txt", "2", 0.48},
        {"/home/chris/SciamiEstesi/23_11_2023/data/clean_data/tempi_corretti_canale_3_231123.txt", "3", 0.75}
    };

    for (const auto& telescope : telescopes) {
        plot_histogram(telescope.filePath, telescope.channel, telescope.efficiency);
        plotTimeDifferences(telescope.filePath, telescope.channel);
    }

    const char* coincidences[] = {
        "/home/chris/SciamiEstesi/23_11_2023/data/concidence_data/coincidenze_triple_231123.txt",
        "/home/chris/SciamiEstesi/23_11_2023/data/concidence_data/coincidenze_doppie23_231123.txt",
        "/home/chris/SciamiEstesi/23_11_2023/data/concidence_data/coincidenze_doppie13_231123.txt"
    };

    const char* channels[] = {"123", "23", "13"};
    double efficiencies[] = {0.2, 0.3, 0.4};

    for (size_t i = 0; i < 3; ++i) {
        plot_histogram(coincidences[i], channels[i], efficiencies[i]);
        calculateAndPlotDiffCoinc(coincidences[i], channels[i]);
    }



    AtmosphericParameterInfo temperatureInfos[] = {
        {"/home/chris/SciamiEstesi/23_11_2023/data/parametri_atmosferici/temperatura/canale_1_231123/Rates_and_temperatura_canale_1_231123.txt", "temperatura", "1", 0.35},
        {"/home/chris/SciamiEstesi/23_11_2023/data/parametri_atmosferici/temperatura/canale_2_231123/Rates_and_temperatura_canale_2_231123.txt", "temperatura", "2", 0.48},
        {"/home/chris/SciamiEstesi/23_11_2023/data/parametri_atmosferici/temperatura/canale_3_231123/Rates_and_temperatura_canale_3_231123.txt", "temperatura", "3", 0.75}
    };

    AtmosphericParameterInfo pressureInfos[] = {
        {"/home/chris/SciamiEstesi/23_11_2023/data/parametri_atmosferici/pressione/canale_1_231123/Rates_and_pressione_canale_1_231123.txt", "pressione", "1", 0.35},
        {"/home/chris/SciamiEstesi/23_11_2023/data/parametri_atmosferici/pressione/canale_2_231123/Rates_and_pressione_canale_2_231123.txt", "pressione", "2", 0.48},
        {"/home/chris/SciamiEstesi/23_11_2023/data/parametri_atmosferici/pressione/canale_3_231123/Rates_and_pressione_canale_3_231123.txt", "pressione", "3", 0.75}
    };

    AtmosphericParameterInfo humidityInfos[] = {
        {"/home/chris/SciamiEstesi/23_11_2023/data/parametri_atmosferici/umidità/canale_1_231123/Rates_and_umidità_canale_1_231123.txt", "umidità", "1", 0.35},
        {"/home/chris/SciamiEstesi/23_11_2023/data/parametri_atmosferici/umidità/canale_2_231123/Rates_and_umidità_canale_2_231123.txt", "umidità", "2", 0.48},
        {"/home/chris/SciamiEstesi/23_11_2023/data/parametri_atmosferici/umidità/canale_3_231123/Rates_and_umidità_canale_3_231123.txt", "umidità", "3", 0.75}
    };

    for (const auto& info : temperatureInfos) {
        plot_AtmosphericParameter(info);
    }

    for (const auto& info : pressureInfos) {
        plot_AtmosphericParameter(info);
    }

    for (const auto& info : humidityInfos) {
        plot_AtmosphericParameter(info);
    }

    std::cout << "TERMINAZIONE SCRIPT" << std::endl;

    return 0;
}
