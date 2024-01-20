import pandas as pd
import time
import numpy as np 


# Registra il timestamp iniziale
start_time = time.time()



def coincidenza_doppia(tel1, tel2):
    coincidenze = []

    idx1, idx2 = 0, 0
    len1, len2= len(tel1), len(tel2)

    while idx1 < len1 and idx2 < len2:
        time1, time2= tel1[idx1], tel2[idx2]
        min_time = min(time1, time2)
        max_time = max(time1, time2)

        delta_1 = abs(max_time - min_time)

        if delta_1 <= 200e-9 and delta_1 >= 150e-9:
            coincidenze.append((time1, time2))

        # Incrementa l'indice del telescopio con il tempo più basso
        if min_time == time1:
            idx1 += 1
        elif min_time == time2:
            idx2 += 1


    return coincidenze


def coincidenze(telescopio_master, telescopio1, telescopio2):
    j=0
    k=0
    len1=len(telescopio1)
    len2=len(telescopio2)
    t1_max = 100e-9
    t1_min = 50e-9
    t2_max = 200e-9
    t2_min = 150e-9

    coincidenze=[]
    
    for i in range(len(telescopio_master)):
        while ( telescopio1[j] < telescopio_master[i]+ t1_min  and j<len1): 
            j+=1
        if telescopio1[j] < telescopio_master[i]+ t1_max:
            while ( telescopio2[k] < telescopio_master[i]+ t2_min  and k<len2): 
                k+=1
            if telescopio2[k] < telescopio_master[i]+ t2_max:
                coincidenze.append((telescopio_master[i], telescopio1[j], telescopio2[k]))
                j+=1  
                k+=1
    return coincidenze
           
           


     

def coincidenze_doppie( telescopio_master, telescopio2):
    j=0
    len1=len(telescopio2)
    t1_max = 100e-9
    t1_min = 50e-9
    t2_max = 200e-9
    t2_min = 150e-9
    t3_max = 150e-9
    t3_min = 80e-9 


    coincidenze=[]
    
    for i in range(len(telescopio_master)):
        while (telescopio2[j] < telescopio_master[i]+ t3_min  and j<len1): 
            j+=1            
        if telescopio2[j] < telescopio_master[i]+ t3_max:
            coincidenze.append((telescopio_master[i], telescopio2[j]))
            j+=1
               
    return coincidenze
             
                
 
                
                                
            

                  
    





tempi_canale_1 = np.loadtxt("/home/chris/SciamiEstesiCOPIA/23_11_2023/data/clean_data/tempi_corretti_canale_1_231123.txt")
tempi_canale_2 = np.loadtxt("/home/chris/SciamiEstesiCOPIA/23_11_2023/data/clean_data/tempi_corretti_canale_2_231123.txt")
tempi_canale_3 = np.loadtxt("/home/chris/SciamiEstesiCOPIA/23_11_2023/data/clean_data/tempi_corretti_canale_3_231123.txt")

#v = coincidenze(tempi_canale_3, tempi_canale_1, tempi_canale_2)
doppie = coincidenze_doppie(tempi_canale_2, tempi_canale_1)
triple = coincidenze(tempi_canale_3, tempi_canale_1, tempi_canale_2)

#np.savetxt('/home/chris/SciamiEstesiCOPIA/SciamiCoincidenza/coincidenze_triple.txt', v, fmt='%.8f')
runnumber=231123
np.savetxt(f'/home/chris/SciamiEstesiCOPIA/23_11_2023/data/concidence_data/coincidenze_doppie12_{runnumber}.txt', doppie, fmt='%.8f')
triple_array = np.array(triple)  # Converti la lista di tuple in un array numpy
np.savetxt(f'/home/chris/SciamiEstesiCOPIA/23_11_2023/data/concidence_data/coincidenze_triple_{runnumber}.txt', triple_array[:, 0], fmt='%.8f')


# Registra il timestamp finale
end_time = time.time()

# Calcola la differenza e stampa il tempo di esecuzione
execution_time = end_time - start_time
print(f"Tempo di esecuzione: {execution_time} secondi")

print(len(tempi_canale_3))
print(len(tempi_canale_2))
print(len(tempi_canale_1))