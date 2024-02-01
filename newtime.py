import matplotlib.pyplot as plt
import numpy as np
import time


def correzione_dati(tempi):
    tempi_corretti = [tempi[0]]

    for i in range(1, len(tempi)):
        if tempi[i] > tempi[i - 1]:
            tempi_corretti.append(tempi[i])
        else:
            diff = tempi[i - 1] - tempi[i]
            for j in range(i, len(tempi)):
                tempi[j] += diff
            tempi_corretti.append(tempi[i])

    tempi[:] = tempi_corretti
    
    
    
def correzione_dati2(tempi):
    tempi_corretti = [tempi[0]]
    corr = 0
    delta_corr= 0 
    for i in range(1, len(tempi)):
        if tempi[i] < tempi[i - 1]:
            diff = tempi[i - 1] - tempi[i]
            corr=corr+diff+delta_corr
              
        tempi_corretti.append(tempi[i]+corr)

    tempi[:] = tempi_corretti
    

# Caricamento dati
runnumber=231123
channel, original_time = np.loadtxt("/home/chris/SciamiEstesi/23_11_2023/data/231123test2.dat", unpack=True)



time = np.copy(original_time)


correzione_dati2(time)


# Separazione dei tempi in base al canale
tempi_canale_1 = time[channel == 1]
tempi_canale_2 = time[channel == 2]
tempi_canale_3 = time[channel == 3]
'''
np.savetxt('//home/chris/SciamiEstesi/23_11_2023/data/clean_data/tempi_corretti__pyTOTALE.txt', time)

# Salvataggio dei tempi corretti per ogni canale
np.savetxt(f'/home/chris/SciamiEstesi/23_11_2023/data/clean_data/tempi_corretti_canale_1_{runnumber}.txt', tempi_canale_1, fmt='%.8f')
np.savetxt(f'/home/chris/SciamiEstesi/23_11_2023/data/clean_data/tempi_corretti_canale_2_{runnumber}.txt', tempi_canale_2, fmt='%.8f')
np.savetxt(f'/home/chris/SciamiEstesi/23_11_2023/data/clean_data/tempi_corretti_canale_3_{runnumber}.txt', tempi_canale_3, fmt='%.8f')
'''
print('Dimesione array dei tempi non corretti:', len(original_time))
# Plot dell'original_time
plt.plot(original_time, label='Original Time')
plt.xlabel('Indice')
plt.ylabel('Tempi Non Corretti')
plt.xlim(0, 0.2e6)  # Imposta i limiti sull'asse x

plt.title('Original Time vs Index')
plt.savefig(f'/home/chris/SciamiEstesi/23_11_2023/plots/grafico_tempi_non_corretti.png')

plt.scatter(range(len(tempi_canale_1)), tempi_canale_1, label='Canale 1')
plt.scatter(range(len(tempi_canale_2)), tempi_canale_2, label='Canale 2')
plt.scatter(range(len(tempi_canale_3)), tempi_canale_3, label='Canale 3')

plt.xlabel('Indice')
plt.ylabel('Tempi Corretti')
plt.legend()
plt.title('Tempi Corretti vs Indice')
plt.xlim(0, len(tempi_canale_1))  # Imposta i limiti sull'asse x


plt.savefig(f'/home/chris/SciamiEstesi/23_11_2023/plots/grafico_tempi_Corretti.png')

plt.show()