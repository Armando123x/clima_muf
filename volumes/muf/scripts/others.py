import os
import re
import numpy
import traceback
import subprocess

from copy import deepcopy
from multiprocessing import Pool
from dateutil.relativedelta import relativedelta 

PATH_IRI = "/iri2020/"

def PREPARE_TXT(lat,lon,date,hour):

    filename = os.path.join(PATH_IRI,'input.txt')

    buffer = numpy.empty(10,dtype=object)
  
    buffer[0] = '0 {:2.2f} {:2.2f}'.format(lat,lon)
    buffer[1] = '{:4d} {:02d}{:02d} 1 {}'.format(date.year,
                                                 date.month,
                                                 date.day,
                                                 hour)
    buffer[2] = '0'
    buffer[3] = '1'
    buffer[4] = '60 1000 10'
    buffer[5] = '0'
    buffer[6] = '1000'
    buffer[7] = '0'
    buffer[8] = '0'
    buffer[9] = '0'

    try:
        with open(filename,'w') as file:

            file.writelines(line +'\n' for line in buffer)
    except Exception as e:

        error = traceback.format_exc()

        print(f"Se produjo un error al generar el archivo input.txt: {error}.")

    return 

def RUN_IRI_PROC():

    ejecutable = str(os.path.join(PATH_IRI,"iri"))
    
 
    # thread = Pool(1,maxtasksperchild=1000)

    current_path = os.getcwd()

    os.chdir(PATH_IRI)


    resultado = subprocess.run([ejecutable],shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    salida = resultado.stdout
    error = resultado.stderr

    os.chdir(current_path)
 
 


    return 

def GET_VALUES_IRI():

    filename = os.path.join(PATH_IRI,'output.txt')

    exists = os.path.isfile(filename)

    ident = 36

    keys = ['km', 'Ne/cm-3','Ne/NmF2', 'Tn/K','Ti/K','Te/K','O+','N+','H+','He+','O2+','NO+', 'Clust','TEC','t/%']

    

    payload = list()


    if exists:

        with open(filename,'r') as file:

            buffer = numpy.array(file.read().splitlines(),dtype=object)

            buffer = buffer[ident:]


            patron = r'[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?|[-+]?\d+'

            for linea in buffer:
              
                numeros = re.findall(patron, linea)
                numeros = [float(numero) for numero in numeros]

                if len(numeros)>2:
                    
                    struct = dict.fromkeys(keys,0)

                    for key,value in zip(keys,numeros):

                        struct[key] = value
                    
                    payload.append(deepcopy(struct))
        
    else:
        print("Archivo 'output.txt' de IRI2020 no existe")
        
    payload = numpy.array(payload,dtype=object)


    return payload 

               
               
def vectorize(dict_,key):
    #return a array with all values that belong key from dict 
    return (numpy.array([data[key] for data in dict_ ]))


def CHECK_VALUES(payload):

    buffer = list()

    for key in ['Clust','TEC','t/%']:

        temporal = vectorize(payload,key)
        
        buffer.append(numpy.nanmean(temporal))

    
    buffer = numpy.nanmean(buffer)

    if int(buffer) ==-1:
        return False
    
    return True

def RUN_IRI(lat,lon,date,hour):

    while 1:

        PREPARE_TXT(lat,lon,date,hour)

        RUN_IRI_PROC()

        payload = GET_VALUES_IRI()

        # Revisar si el promedio no es -1 

        flag = CHECK_VALUES(payload)

        ##############################

        ne = vectorize(payload,'Ne/cm-3')
    
        heights = vectorize(payload,'km')

        if flag:

            return ne,heights,date
        
        date = date + relativedelta(days=-1)

def JULDAY(datetime):
    
    a = (14 - datetime.month) // 12
    y = datetime.year + 4800 - a
    m = datetime.month + 12 * a - 3
    jdn = datetime.day + (153 * m + 2) // 5 + 365 * y + y // 4 - y // 100 + y // 400 - 32045
    
    
    return jdn

def READ_SAO_FILE(IU,linea):
        '''
        Metodo que realiza la lectura del archivo .SAO

        Parameters
        ----------
        IU : TYPE
            DESCRIPTION.
        linea : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''   
        IERR = 1
    
        #DATA FILE INDEX
        IDFI = [0]*80
        for i in range(40):
            IDFI[i] = int(linea[i*3:i*3+3])
        
        linea = IU.readline()
        for i in range(40):
            IDFI[i+40] = int(linea[i*3:i*3+3])
        
        if(not(IDFI[0] > 0 and IDFI[0] < 17)):
            IERR = 1
            return IERR 
    
        IERR = 0
    
        #GEOPHYSICAL CONSTANTS
        GCONST = IDFI[0]*[0]
        if(IDFI[0] > 0):
            linea = IU.readline()
            for i in range(IDFI[0]):
                GCONST[i] = float(linea[i*7:(i+1)*7])
    
        #SYSTEM DESCRIPTION AND OPERATOR'S MESSAGE
        SYSDES = ''
        if(IDFI[1] > 0):
            linea = IU.readline()
            SYSDES = linea[0:120]
    
        #TIME STAMP AND SOUNDER SETTINGS
        IPREF = IDFI[2]*['']
        if(IDFI[2] > 0):
            linea = IU.readline()
            for i in range(IDFI[2]):
                IPREF[i] = linea[i:i+1]
                
        #SCALED IONOSPHERIC CHARACTERISTICS
        SCALED = IDFI[3]*[0]
        if(IDFI[3] > 0):
            for i in range(int(IDFI[3]/15)):
                linea = IU.readline()
                for j in range(15):
                    SCALED[j+15*i] = float(linea[j*8:(j+1)*8])
            if(IDFI[3]%15 != 0):
                linea = IU.readline()
                for i in range(IDFI[3]%15):
                    SCALED[int(IDFI[3]/15)*15 + i] = float(linea[i*8:(i+1)*8])
                    
        #ANALYSIS FLAGS
        IAF = IDFI[4]*[0]
        if(IDFI[4] > 0):
            linea = IU.readline()
            for i in range(IDFI[4]):
                IAF[i] = linea[i*2:(i+1)*2]
    
        #DOPPLER TRANSLATION TABLE
        DTT = IDFI[5]*[0]
        if(IDFI[5]>0):
            linea = IU.readline()
            for i in range(IDFI[5]):
                DTT[i] = float(linea[i*7:(i+1)*7])
    
        #O-TRACE POINTS - F2 LAYER
        #VIRTUAL HEIGHTS
        OTF = []
        if(IDFI[6] > 0):
            OTF = IDFI[6]*[0]
            for i in range(int(IDFI[6]/15)):
                linea = IU.readline()
                for j in range(15):
                    OTF[j+15*i] = float(linea[j*8:(j+1)*8])
            if(IDFI[6]%15 != 0):
                linea = IU.readline()
                for i in range(IDFI[6]%15):
                    OTF[int(IDFI[6]/15)*15 + i] = float(linea[i*8:(i+1)*8])
    
        #TRUE HEIGHTS
        OTHF = []
        if(IDFI[7] > 0):
            OTHF = IDFI[7]*[0]
            for i in range(int(IDFI[7]/15)):
                linea = IU.readline()
                for j in range(15):
                    OTHF[j+15*i] = float(linea[j*8:(j+1)*8])
            if(IDFI[7]%15 != 0):
                linea = IU.readline()
                for i in range(IDFI[7]%15):
                    OTHF[int(IDFI[7]/15)*15 + i] = float(linea[i*8:(i+1)*8])
    
        #AMPLITUDES
        IOAF = []
        if(IDFI[8] > 0):
            IOAF = IDFI[8]*[0]
            for i in range(int(IDFI[8]/40)):
                linea = IU.readline()
                for j in range(40):
                    IOAF[j + 40*i] = int(linea[j*3:(j+1)*3])
            if(IDFI[8]%40 != 0):
                linea = IU.readline()
                for i in range(IDFI[8]%40):
                    IOAF[int(IDFI[8]/40)*40+i] = int(linea[i*3:(i+1)*3])
    
        #DOPPLER NUMBERS
        IODF = []
        if(IDFI[9] > 0):
            IODF = IDFI[9]*[0]
            linea = IU.readline()
            for i in range(IDFI[9]):
                IODF[i] = linea[i:i+1]
    
        #FREQUENCIES
        FTOF = []
        if(IDFI[10] > 0):
            FTOF = IDFI[10]*[0]
            for i in range(int(IDFI[10]/15)):
                linea = IU.readline()
                for j in range(15):
                    FTOF[j+15*i] = float(linea[j*8:(j+1)*8])
            if(IDFI[10]%15 != 0):
                linea = IU.readline()
                for i in range(IDFI[10]%15):
                    FTOF[int(IDFI[10]/15)*15 + i] = float(linea[i*8:(i+1)*8])
    
        #O-TRACE POINTS - F1 LAYER
        #VIRTUAL HEIGHTS
        OTF1 = []
        if(IDFI[11] > 0):
            OTF1 = IDFI[11]*[0]
            for i in range(int(IDFI[11]/15)):
                linea = IU.readline()
                for j in range(15):
                    OTF1[j+15*i] = float(linea[j*8:(j+1)*8])
            if(IDFI[11]%15 != 0):
                linea = IU.readline()
                for i in range(IDFI[11]%15):
                    OTF1[int(IDFI[11]/15)*15 + i] = float(linea[i*8:(i+1)*8])
    
        #TRUE HEIGHTS
        OTHF1 = []
        if(IDFI[12] > 0):
            OTHF1 = IDFI[12]*[0]
            for i in range(int(IDFI[12]/15)):
                linea = IU.readline()
                for j in range(15):
                    OTHF1[j+15*i] = float(linea[j*8:(j+1)*8])
            if(IDFI[12]%15 != 0):
                linea = IU.readline()
                for i in range(IDFI[12]%15):
                    OTHF1[int(IDFI[12]/15)*15 + i] = float(linea[i*8:(i+1)*8])
    
        #AMPLITUDES
        IOAF1 = []
        if(IDFI[13] > 0):
            IOAF1 = IDFI[13]*[0]
            for i in range(int(IDFI[13]/40)):
                linea = IU.readline()
                for j in range(40):
                    IOAF1[j + 40*i] = int(linea[j*3:(j+1)*3])
            if(IDFI[13]%40 != 0):
                linea = IU.readline()
                for i in range(IDFI[13]%40):
                    IOAF1[int(IDFI[13]/40)*40+i] = int(linea[i*3:(i+1)*3])
    
        #DOPPLER NUMBERS
        IODF1 = []
        if(IDFI[14] > 0):
            IODF1 = IDFI[14]*[0]
            linea = IU.readline()
            for i in range(IDFI[14]):
                IODF1[i] = linea[i:i+1]
    
        #FREQUENCIES
        FTOF1 = []
        if(IDFI[15] > 0):
            FTOF1 = IDFI[15]*[0]
            for i in range(int(IDFI[15]/15)):
                linea = IU.readline()
                for j in range(15):
                    FTOF1[j+15*i] = float(linea[j*8:(j+1)*8])
            if(IDFI[15]%15 != 0):
                linea = IU.readline()
                for i in range(IDFI[15]%15):
                    FTOF1[int(IDFI[15]/15)*15 + i] = float(linea[i*8:(i+1)*8])
    
        #O-TRACE POINTS - E LAYER
        #VIRTUAL HEIGHTS
        OTE = []
        if(IDFI[16] > 0):
            OTE = IDFI[16]*[0]
            for i in range(int(IDFI[16]/15)):
                linea = IU.readline()
                for j in range(15):
                    OTE[j+15*i] = float(linea[j*8:(j+1)*8])
            if(IDFI[16]%15 != 0):
                linea = IU.readline()
                for i in range(IDFI[16]%15):
                    OTE[int(IDFI[16]/15)*15 + i] = float(linea[i*8:(i+1)*8])
    
        #TRUE HEIGHTS
        OTHE = []
        if(IDFI[17] > 0):
            OTHE = IDFI[17]*[0]
            for i in range(int(IDFI[17]/15)):
                linea = IU.readline()
                for j in range(15):
                    OTHE[j+15*i] = float(linea[j*8:(j+1)*8])
            if(IDFI[17]%15 != 0):
                linea = IU.readline()
                for i in range(IDFI[17]%15):
                    OTHE[int(IDFI[17]/15)*15 + i] = float(linea[i*8:(i+1)*8])
    
        #AMPLITUDES
        IOAE = []
        if(IDFI[18] > 0):
            IOAE = IDFI[18]*[0]
            for i in range(int(IDFI[18]/40)):
                linea = IU.readline()
                for j in range(40):
                    IOAE[j + 40*i] = int(linea[j*3:(j+1)*3])
            if(IDFI[18]%40 != 0):
                linea = IU.readline()
                for i in range(IDFI[18]%40):
                    IOAE[int(IDFI[18]/40)*40+i] = int(linea[i*3:(i+1)*3])
    
        #DOPPLER NUMBERS
        IODE = []
        if(IDFI[19] > 0):
            IODE = IDFI[19]*[0]
            linea = IU.readline()
            for i in range(IDFI[19]):
                IODE[i] = linea[i:i+1]
    
        #FREQUENCIES
        FTOE = []
        if(IDFI[20] > 0):
            FTOE = IDFI[20]*[0]
            for i in range(int(IDFI[20]/15)):
                linea = IU.readline()
                for j in range(15):
                    FTOE[j+15*i] = float(linea[j*8:(j+1)*8])
            if(IDFI[20]%15 != 0):
                linea = IU.readline()
                for i in range(IDFI[20]%15):
                    FTOE[int(IDFI[20]/15)*15 + i] = float(linea[i*8:(i+1)*8])
    
        #X-TRACE POINTS - F2 LAYER
        #VIRTUAL HEIGHTS
        XTF = []
        if(IDFI[21] > 0):
            XTF = IDFI[21]*[0]
            for i in range(int(IDFI[21]/15)):
                linea = IU.readline()
                for j in range(15):
                    XTF[j+15*i] = float(linea[j*8:(j+1)*8])
            if(IDFI[21]%15 != 0):
                linea = IU.readline()
                for i in range(IDFI[21]%15):
                    XTF[int(IDFI[21]/15)*15 + i] = float(linea[i*8:(i+1)*8])
    
        #AMPLITUDES
        IXAF = []
        if(IDFI[22] > 0):
            IXAF = IDFI[22]*[0]
            for i in range(int(IDFI[22]/40)):
                linea = IU.readline()
                for j in range(40):
                    IXAF[j + 40*i] = int(linea[j*3:(j+1)*3])
            if(IDFI[22]%40 != 0):
                linea = IU.readline()
                for i in range(IDFI[22]%40):
                    IXAF[int(IDFI[22]/40)*40+i] = int(linea[i*3:(i+1)*3])
    
        #DOPPLER NUMBERS
        IXDF = []
        if(IDFI[23] > 0):
            IXDF = IDFI[23]*[0]
            linea = IU.readline()
            for i in range(IDFI[23]):
                IXDF[i] = int(linea[i:i+1])
    
        #FREQUENCIES
        FTXF = []
        if(IDFI[24] > 0):
            FTXF = IDFI[24]*[0]
            for i in range(int(IDFI[24]/15)):
                linea = IU.readline()
                for j in range(15):
                    FTXF[j+15*i] = float(linea[j*8:(j+1)*8])
            if(IDFI[24]%15 != 0):
                linea = IU.readline()
                for i in range(IDFI[24]%15):
                    FTXF[int(IDFI[24]/15)*15 + i] = float(linea[i*8:(i+1)*8])
            
        #X-TRACE POINTS - F1 LAYER
        #VIRTUAL HEIGHTS
        XTF1 = []
        if(IDFI[25] > 0):
            XTF1 = IDFI[25]*[0]
            for i in range(int(IDFI[25]/15)):
                linea = IU.readline()
                for j in range(15):
                    XTF1[j+15*i] = float(linea[j*8:(j+1)*8])
            if(IDFI[25]%15 != 0):
                linea = IU.readline()
                for i in range(IDFI[25]%15):
                    XTF1[int(IDFI[25]/15)*15 + i] = float(linea[i*8:(i+1)*8])
    
        #AMPLITUDES
        IXAF1 = []
        if(IDFI[26] > 0):
            IXAF1 = IDFI[26]*[0]
            for i in range(int(IDFI[26]/40)):
                linea = IU.readline()
                for j in range(40):
                    IXAF1[j + 40*i] = int(linea[j*3:(j+1)*3])
            if(IDFI[26]%40 != 0):
                linea = IU.readline()
                for i in range(IDFI[26]%40):
                    IXAF1[int(IDFI[26]/40)*40+i] = int(linea[i*3:(i+1)*3])
    
        #DOPPLER NUMBERS
        IXDF1 = []
        if(IDFI[27] > 0):
            IXDF1 = IDFI[27]*[0]
            linea = IU.readline()
            for i in range(IDFI[27]):
                IXDF1[i] = int(linea[i:i+1])
    
        #FREQUENCIES
        FTXF1 = []
        if(IDFI[28] > 0):
            FTXF1 = IDFI[28]*[0]
            for i in range(int(IDFI[28]/15)):
                linea = IU.readline()
                for j in range(15):
                    FTXF1[j+15*i] = float(linea[j*8:(j+1)*8])
            if(IDFI[28]%15 != 0):
                linea = IU.readline()
                for i in range(IDFI[28]%15):
                    FTXF1[int(IDFI[28]/15)*15 + i] = float(linea[i*8:(i+1)*8])
    
        #X-TRACE POINTS - E LAYER
        #VIRTUAL HEIGHTS
        XTE = []
        if(IDFI[29] > 0):
            XTE = IDFI[29]*[0]
            for i in range(int(IDFI[29]/15)):
                linea = IU.readline()
                for j in range(15):
                    XTE[j+15*i] = float(linea[j*8:(j+1)*8])
            if(IDFI[29]%15 != 0):
                linea = IU.readline()
                for i in range(IDFI[29]%15):
                    XTE[int(IDFI[29]/15)*15 + i] = float(linea[i*8:(i+1)*8])
    
        #AMPLITUDES
        IXAE = []
        if(IDFI[30] > 0):
            IXAE = IDFI[30]*[0]
            for i in range(int(IDFI[30]/40)):
                linea = IU.readline()
                for j in range(40):
                    IXAE[j + 40*i] = int(linea[j*3:(j+1)*3])
            if(IDFI[30]%40 != 0):
                linea = IU.readline()
                for i in range(IDFI[30]%40):
                    IXAE[int(IDFI[30]/40)*40+i] = int(linea[i*3:(i+1)*3])
    
        #DOPPLER NUMBERS
        IXDE = []
        if(IDFI[31] > 0):
            IXDE = IDFI[31]*[0]
            linea = IU.readline()
            for i in range(IDFI[31]):
                IXDE[i] = int(linea[i:i+1])
    
        #FREQUENCIES
        FTXE = []
        if(IDFI[32] > 0):
            FTXE = IDFI[32]*[0]
            for i in range(int(IDFI[32]/15)):
                linea = IU.readline()
                for j in range(15):
                    FTXE[j+15*i] = float(linea[j*8:(j+1)*8])
            if(IDFI[32]%15 != 0):
                linea = IU.readline()
                for i in range(IDFI[32]%15):
                    FTXE[int(IDFI[32]/15)*15 + i] = float(linea[i*8:(i+1)*8])
    
        #MEDIAN AMPLITUDES OF F ECHOES
        MEDF = []
        if(IDFI[33] > 0):
            MEDF = IDFI[33]*[0]
            for i in range(int(IDFI[33]/40)):
                linea = IU.readline()
                for j in range(40):
                    MEDF[j + 40*i] = int(linea[j*3:(j+1)*3])
            if(IDFI[33]%40 != 0):
                linea = IU.readline()
                for i in range(IDFI[33]%40):
                    MEDF[int(IDFI[33]/40)*40+i] = int(linea[i*3:(i+1)*3])
    
        #MEDIAN AMPLITUDES OF E ECHOES
        MEDE = []
        if(IDFI[34] > 0):
            MEDE = IDFI[34]*[0]
            for i in range(int(IDFI[34]/40)):
                linea = IU.readline()
                for j in range(40):
                    MEDE[j + 40*i] = int(linea[j*3:(j+1)*3])
            if(IDFI[34]%40 != 0):
                linea = IU.readline()
                for i in range(IDFI[34]%40):
                    MEDE[int(IDFI[34]/40)*40+i] = int(linea[i*3:(i+1)*3])
    
        #MEDIAN AMPLITUDES OF ES ECHOES
        MEDES = []
        if(IDFI[35] > 0):
            MEDES = IDFI[35]*[0]
            for i in range(int(IDFI[35]/40)):
                linea = IU.readline()
                for j in range(40):
                    MEDES[j + 40*i] = int(linea[j*3:(j+1)*3])
            if(IDFI[35]%40 != 0):
                linea = IU.readline()
                for i in range(IDFI[35]%40):
                    MEDES[int(IDFI[35]/40)*40+i] = int(linea[i*3:(i+1)*3])
    
        #TRUE HEIGHTS COEFFICIENTS F2 LAYER UMLCAR METHOD
        THF2 = []
        if(IDFI[36] > 0):
            THF2 = IDFI[36]*[0]
            for i in range(int(IDFI[36]/10)):
                linea = IU.readline()
                for j in range(10):
                    THF2[j + 10*i] = float(linea[j*11:(j+1)*11])
            if(IDFI[36]%10 != 0):
                linea = IU.readline()
                for i in range(IDFI[36]%10):
                    THF2[int(IDFI[36]/10)*10+i] = float(linea[i*11:(i+1)*11])
    
        #TRUE HEIGHTS COEFFICIENTS F1 LAYER UMLCAR METHOD            
        THF1 = []
        if(IDFI[37] > 0):
            THF1 = IDFI[37]*[0]
            for i in range(int(IDFI[37]/10)):
                linea = IU.readline()
                for j in range(10):
                    THF1[j + 10*i] = float(linea[j*11:(j+1)*11])
            if(IDFI[37]%10 != 0):
                linea = IU.readline()
                for i in range(IDFI[37]%10):
                    THF1[int(IDFI[37]/10)*10+i] = float(linea[i*11:(i+1)*11])
    
        #TRUE HEIGHTS COEFFICIENTS E LAYER UMLCAR METHOD
        THE = []
        if(IDFI[38] > 0):
            THE = IDFI[38]*[0]
            for i in range(int(IDFI[38]/10)):
                linea = IU.readline()
                for j in range(10):
                    THE[j + 10*i] = float(linea[j*11:(j+1)*11])
            if(IDFI[38]%10 != 0):
                linea = IU.readline()
                for i in range(IDFI[38]%10):
                    THE[int(IDFI[38]/10)*10+i] = float(linea[i*11:(i+1)*11])
    
        #QUAZI-PARABOLIC SEGMENTS FITTED TO THE PROFILE
        QPCOEF = []
        if(IDFI[39] > 0):
            QPCOEF = IDFI[39]*[0]
            for i in range(int(IDFI[39]/6)):
                linea = IU.readline()
                for j in range(6):
                    QPCOEF[j + 6*i] = float(linea[j*20:(j+1)*20])
            if(IDFI[39]%6 != 0):
                linea = IU.readline()
                for i in range(IDFI[39]%6):
                    QPCOEF[int(IDFI[39]/6)*6+i] = float(linea[i*20:(i+1)*20])
    
        #EDIT FLAGS - CHARACTERISTICS
        IEDF = []
        if(IDFI[40] > 0):
            IEDF = IDFI[40]*[0]
            linea = IU.readline()
            for i in range(IDFI[40]):
                IEDF[i] = int(linea[i:i+1])
    
        #VALLEY DESCRIPTION - W,D UMLCAR MODEL
        THVAL = []
        if(IDFI[41] > 0):
            THVAL = IDFI[41]*[0]
            for i in range(int(IDFI[41]/10)):
                linea = IU.readline()
                for j in range(10):
                    THVAL[j + 10*i] = float(linea[j*11:(j+1)*11])
            if(IDFI[41]%10 != 0):
                linea = IU.readline()
                for i in range(IDFI[41]%10):
                    THVAL[int(IDFI[41]/10)*10+i] = float(linea[i*11:(i+1)*11])
    
        #O-TRACE POINTS - ES LAYER
        #VIRTUAL HEIGHTS
        OTSE = []
        if(IDFI[42] > 0):
            OTSE = IDFI[42]*[0]
            for i in range(int(IDFI[42]/15)):
                linea = IU.readline()
                for j in range(15):
                    OTSE[j+15*i] = float(linea[j*8:(j+1)*8])
            if(IDFI[42]%15 != 0):
                linea = IU.readline()
                for i in range(IDFI[42]%15):
                    OTSE[int(IDFI[42]/15)*15 + i] = float(linea[i*8:(i+1)*8])
    
        #AMPLITUDES
        IOASE = []
        if(IDFI[43] > 0):
            IOASE = IDFI[43]*[0]
            for i in range(int(IDFI[43]/40)):
                linea = IU.readline()
                for j in range(40):
                    IOASE[j + 40*i] = int(linea[j*3:(j+1)*3])
            if(IDFI[43]%40 != 0):
                linea = IU.readline()
                for i in range(IDFI[43]%40):
                    IOASE[int(IDFI[43]/40)*40+i] = int(linea[i*3:(i+1)*3])
    
        #DOPPLER NUMBERS
        IODSE = []
        if(IDFI[44] > 0):
            IODSE = IDFI[44]*[0]
            linea = IU.readline()
            for i in range(IDFI[44]):
                IODSE[i] = int(linea[i:i+1])
    
        #FREQUENCIES
        FTOSE = []
        if(IDFI[45] > 0):
            FTOSE = IDFI[45]*[0]
            for i in range(int(IDFI[45]/15)):
                linea = IU.readline()
                for j in range(15):
                    FTOSE[j+15*i] = float(linea[j*8:(j+1)*8])
            if(IDFI[45]%15 != 0):
                linea = IU.readline()
                for i in range(IDFI[45]%15):
                    FTOSE[int(IDFI[45]/15)*15 + i] = float(linea[i*8:(i+1)*8])
    
        #O-TRACE POINTS - E AURORAL LAYER
        #VIRTUAL HEIGHTS
        OTSE = []
        IOTSE = []
        if(IDFI[46] > 0):
            OTSE = IDFI[46]*[0]
            IOTSE = IDFI[46]*[0]
            for i in range(int(IDFI[46]/15)):
                linea = IU.readline()
                for j in range(15):
                    OTSE[j+15*i] = float(linea[j*8:(j+1)*8])
                    IOTSE[j+15*i] = int(OTSE[j+15*i])
            if(IDFI[46]%15 != 0):
                linea = IU.readline()
                for i in range(IDFI[46]%15):
                    OTSE[(IDFI[46]/15)*15 + i] = float(linea[i*8:(i+1)*8])
                    IOTSE[int(IDFI[46]/15)*15 + i] = int(OTSE[(IDFI[46]/15)*15 + i])
    
        #AMPLITUDES
        IOASE = []
        if(IDFI[47] > 0):
            IOASE = IDFI[47]*[0]
            for i in range(int(IDFI[47]/40)):
                linea = IU.readline()
                for j in range(40):
                    IOASE[j + 40*i] = int(linea[j*3:(j+1)*3])
            if(IDFI[47]%40 != 0):
                linea = IU.readline()
                for i in range(IDFI[47]%40):
                    IOASE[int(IDFI[47]/40)*40+i] = int(linea[i*3:(i+1)*3])
    
        #DOPPLER NUMBERS
        IODSE = []
        if(IDFI[48] > 0):
            IODSE = IDFI[48]*[0]
            linea = IU.readline()
            for i in range(IDFI[48]):
                IODSE[i] = int(linea[i:i+1])
    
        #FREQUENCIES
        FTOSE = []
        if(IDFI[49] > 0):
            FTOSE = IDFI[49]*[0]
            for i in range(int(IDFI[49]/15)):
                linea = IU.readline()
                for j in range(15):
                    FTOSE[j+15*i] = float(linea[j*8:(j+1)*8])
            if(IDFI[49]%15 != 0):
                linea = IU.readline()
                for i in range(IDFI[49]%15):
                    FTOSE[int(IDFI[49]/15)*15 + i] = float(linea[i*8:(i+1)*8])
    
        #TRUE HEIGHT PROFILE
        #TRUE HEIGHTS
        HTAB = []
        if(IDFI[50] > 0):
            HTAB = IDFI[50]*[0]
            for i in range(int(IDFI[50]/15)):
                linea = IU.readline()
                for j in range(15):
                    HTAB[j+15*i] = float(linea[j*8:(j+1)*8])
            if(IDFI[50]%15 != 0):
                linea = IU.readline()
                for i in range(IDFI[50]%15):
                    HTAB[int(IDFI[50]/15)*15 + i] = float(linea[i*8:(i+1)*8])
    
        #PLASMA FREQUENCIES
        FTAB = []
        if(IDFI[51] > 0):
            FTAB = IDFI[51]*[0]
            for i in range(int(IDFI[51]/15)):
                linea = IU.readline()
                for j in range(15):
                    FTAB[j+15*i] = float(linea[j*8:(j+1)*8])
            if(IDFI[51]%15 != 0):
                linea = IU.readline()
                for i in range(IDFI[51]%15):
                    FTAB[int(IDFI[51]/15)*15 + i] = float(linea[i*8:(i+1)*8])
    
        #ELECTRON DENSITIES[e/cm^3]
        NTAB = [] 
        if(IDFI[52] > 0):
            NTAB = IDFI[52]*[0]
            for i in range(int(IDFI[52]/15)):
                linea = IU.readline()
                for j in range(15):
                    NTAB[j+15*i] = float(linea[j*8:(j+1)*8])
            if(IDFI[52]%15 != 0):
                linea = IU.readline()
                for i in range(IDFI[52]%15):
                    NTAB[int(IDFI[52]/15)*15 + i] = float(linea[i*8:(i+1)*8])
    
        #URSI QUALIFYING AND DESCRIPTIVE LETTERS
        #QUALIFYING LETTERS
        QL = []
        if(IDFI[53] > 0):
            QL = IDFI[53]*['']
            linea = IU.readline()
            for i in range(IDFI[53]):
                QL[i] = linea[i:i+1]
    
        #DESCRIPTIVE LETTERS
        DL = []
        if(IDFI[54] > 0):
            DL = IDFI[54]*['']
            linea = IU.readline()
            for i in range(IDFI[54]):
                DL[i] = linea[i:i+1]
    
        #EDIT FLAGS - TRACES AND PROFILE
        IEDFTP = []
        if(IDFI[55] > 0):
            IEDFTP = IDFI[55]*[0]
            linea = IU.readline()
            for i in range(IDFI[55]):
                IEDFTP[i] = int(linea[i:i+1])
    
        IREAD = 1
    
        #IOTF, IOTHF, IOAFB, IOTF1, IOTHF1, IOAF1B, IOAF1B, IOTE, IOTHE, IXTF,
        #IXTF1, IXTE, OTFB, OTF1B, HTABB, FTABB 
        return IDFI, GCONST, SYSDES, IPREF, SCALED, IAF, DTT, OTF, IODF, FTOF, IODF1, \
               FTOF1, IOAE, IODE, FTOE, IXAF, IXDF, FTXF, IXAF1, IXDF1, FTXF1, \
               IXAE, IXDE, FTXE, MEDF, MEDE, MEDES, THF2, THF1, THE, QPCOEF, THVAL, \
               IEDF, IOTSE, IOASE, IODSE, FTOSE, OTHF1, OTE, OTHE, XTF, XTF1, XTE, \
               OTSE, HTAB, FTAB, QL, DL, IEDFTP
