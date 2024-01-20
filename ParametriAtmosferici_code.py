import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime, timedelta
from tqdm import tqdm
import pdb 
import matplotlib.dates as mdates


def process_atmospheric_parameter(parametro, channel_date, colonna):
    
    # Carica i dati dei tempi
    file_path_tempi_TOT = '/home/chris/SciamiEstesiCOPIA/23_11_2023/data/clean_data/tempi_corretti__TOTALE.txt'
    df_tempi_TOT = pd.read_csv(file_path_tempi_TOT, header=None, names=['tempi'])
    ultimo_valore = df_tempi_TOT.iloc[-1, 0]


    # Data e ora di fine acquisizione
    data_ora_fine_acquisizione = datetime(2023, 11, 27, 4, 52, 38)

    # Durata totale dell'acquisizione (in secondi)
    durata_acquisizione = ultimo_valore

    # Calcola il timestamp di fine acquisizione
    timestamp_fine_acquisizione = int(data_ora_fine_acquisizione.timestamp())

    # Calcola il timestamp di inizio acquisizione
    timestamp_inizio_acquisizione = timestamp_fine_acquisizione - durata_acquisizione

    # Converti i timestamp in oggetti datetime
    data_ora_inizio_acquisizione = datetime.utcfromtimestamp(timestamp_inizio_acquisizione)
    data_ora_inizio_acquisizione = data_ora_inizio_acquisizione.replace(microsecond=0)

    # Stampare i risultati
    print("Inizio acquisizione:", data_ora_inizio_acquisizione)
    print("Fine acquisizione:", data_ora_fine_acquisizione)


    # Crea una colonna con i timestamp utilizzando l'ora esatta di inizio acquisizione
    file_path_tempi_channel = f'/home/chris/SciamiEstesiCOPIA/23_11_2023/data/clean_data/tempi_corretti_{channel_date}.txt'
    df_tempi_channel = pd.read_csv(file_path_tempi_channel, header=None, names=['tempi'])
    df_tempi_channel['timestamp'] = data_ora_inizio_acquisizione + pd.to_timedelta(df_tempi_channel['tempi'], unit='s')
    #pdb.set_trace() per debugging

    # Filtro per limitare i dati di df_tempi prima dell2 22:55 del 26 novembre 2023 perché i dati della temperatura arrivano fino alle 22:55 del 26.11.2023
    df_tempi_channel = df_tempi_channel[df_tempi_channel['timestamp'] < datetime(2023, 11, 26, 22, 55)]

    # Carica i dati dal file
    data_parametri = np.loadtxt("/home/chris/SciamiEstesiCOPIA/23_11_2023/data/parametri_atmosferici/2023-11-27_downld08.txt", usecols=(0, 1, colonna), dtype=str)

    # Creazione di un DataFrame 'df_parametri' utilizzando i dati dalla variabile 'data_parametri'
    # Le colonne del DataFrame saranno etichettate come 'data', 'ora', e 'temperatura'
    df_parametri = pd.DataFrame(data_parametri, columns=['data', 'ora', parametro])


    # Combina le colonne 'data' e 'ora' in una colonna datetime 'timestamp'
    df_parametri['timestamp'] = pd.to_datetime(df_parametri['data'] + ' ' + df_parametri['ora'], format='%d/%m/%y %H:%M')

    '''
    # Salva il DataFrame in un file CSV, alla fine serve solo come debug
    df_parametri.to_csv('/home/chris/Tempi.txt', index=False)


    # Elimina colonne 'data' e 'ora' se non sono più necessarie
    df_parametri = df_parametri.drop(['data', 'ora'], axis=1)
    '''
    # Ordina i DataFrame per i timestamp
    df_tempi_channel= df_tempi_channel.sort_values('timestamp')
    df_parametri = df_parametri.sort_values('timestamp')

    # Esegui unione basata sul timestamp più vicino (precedente)
    df_combined = pd.merge_asof(df_tempi_channel, df_parametri, on='timestamp', direction='backward')

    '''
    # Salva il risultato in un nuovo file CSV
    output_file_path = '/home/chris/Tempi_e_Parametri_Completi.txt'
    with tqdm(total=len(df_combined), desc="Progresso") as pbar:
        df_combined.to_csv(output_file_path, index=False)
        pbar.update(len(df_combined))
    '''

    # Dati
    df_combined['timestamp'] = pd.to_datetime(df_combined['timestamp'])

    #tempi = df_combined['tempi']
    #temperature = df_combined['temperatura']
    # Imposta l'inizio del primo bin come l'ora di inizio acquisizione
    start_time = df_combined['timestamp'].min()
    max_time = df_combined['timestamp'].max()

    # Calcola l'intervallo di 5 minuti
    interval = pd.Timedelta(hours=1)



    # Aggiungi un bin aggiuntivo per coprire il periodo dall'ora massima fino alla fine dell'acquisizione
    bins = pd.date_range(start=start_time, end=max(max_time, data_ora_fine_acquisizione), freq=interval)


    # Calcola i valori centrali di ogni bin
    bin_centers = bins[:-1] + (bins[1] - bins[0]) 

    # Trova il valore più vicino in df_combined['timestamp'] per ogni valore centrale di bin
    nearest_index = np.searchsorted(df_combined['timestamp'].values, bin_centers) - 1
    # Assicurati che gli indici siano all'interno del range di df_combined
    nearest_index = np.clip(nearest_index, 0, len(df_combined) - 1)
    # Seleziona i timestamp più vicini, le ore associate e le temperature associate
    nearest_times = df_combined['timestamp'].values[nearest_index]
    nearest_ora_values = df_combined['ora'].values[nearest_index]
    parametro_values = df_combined[parametro].values[nearest_index]


    parametro_ass = []

    # Stampa i timestamp più vicini e le relative temperature
    for bin_center, ora, par in zip(bin_centers, nearest_ora_values, parametro_values):
        print(f"Timestamp più vicino: {bin_center}, Ora associata: {ora}, {parametro} : {par}")
        # Aggiungi la temperatura associata al vettore
        parametro_ass.append(par)
    '''
    # Crea l'istogramma
    plt.hist(df_combined['timestamp'], bins=bins, edgecolor='black')

    '''
    # Calcola i rate
    bin_counts, bin_edges = np.histogram(df_combined['timestamp'], bins=bins)

    # Calcola la larghezza di ciascun bin
    bin_widths = np.diff(bin_edges)

    # Converti bin_counts in tipo float prima di eseguire la divisione
    bin_counts_float = bin_counts.astype(float)

    # Converti bin_widths in secondi prima di eseguire la divisione
    bin_widths_seconds = bin_widths.astype('timedelta64[s]').astype(float)

    # Calcola i rate come rapporto tra il contenuto del bin e la larghezza del bin
    rates = bin_counts_float / bin_widths_seconds
    parametro_ass = pd.to_numeric(parametro_ass, errors='coerce')
    rate_medio = np.mean(rates)
    print(rate_medio)

    # Creare un DataFrame con rates e temperature_ass
    df_export = pd.DataFrame({parametro: parametro_ass, 'Rate': rates})

    # Salva il DataFrame in un file di testo
    output_txt_path = f'/home/chris/SciamiEstesiCOPIA/23_11_2023/data/parametri_atmosferici/{parametro}/{channel_date}/Rates_and_{parametro}_{channel_date}.txt'
    df_export.to_csv(output_txt_path, sep='\t', index=False, header=False)

parametro_da_processare = 'umidità'
channel_date = 'canale_3_231123'
colonna = 5 #colonna 2 = temperatura; colonna 5 = umidità; colonna 16 = pressione
process_atmospheric_parameter(parametro_da_processare, channel_date, colonna)

