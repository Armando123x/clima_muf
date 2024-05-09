#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 13 11:11:07 2022

@author: arx
"""
import os 
import gc 
import sys 
import pytz
import math 
import scipy 
import numpy
import subprocess
import numpy as np  
import datetime as dt
import matplotlib.pyplot as plt

from math import *
from io import BytesIO
from pathlib import Path
from multiprocessing import Pool,cpu_count
from .others import *
from scipy.optimize import curve_fit
from matplotlib.offsetbox import AnchoredText
from dateutil.relativedelta import relativedelta
from scipy.interpolate import InterpolatedUnivariateSpline



num_cpus = cpu_count()


class plot_IRI ():
    '''
    plot_IRI creates an object to run the model IRI thought Pyglow package
    and Basemap 
    
    Paramaters
    ----------
    latitude, longitude= Float values of a point desired.
    version_IRI= Four possible values -> 2007, 2012, 2016,2020
    
    IRI Notes
    ---------
    -> IRI actually run with the new version 2020 written in Fortran. 


    '''
    
    
    def __config_dict(self):
        self.values = dict()

        self.values['Chachapoyas']   =   numpy.array([-6.23169,-77.86903])
        self.values['Huaraz']        =   numpy.array([-9.52779, -77.52778])
        self.values['Abancay']        =   numpy.array([-13.63389,-72.88139])
        self.values['Arequipa']        =   numpy.array([-16.39889,-71.535])
        self.values['Ayacucho']        =   numpy.array([-13.15878,-74.22321])
        self.values['Cajamarca']        =   numpy.array([-7.1637,-78.50027])
        self.values['Cuzco']        =   numpy.array([-13.52264,-71.96734])
        self.values['Huancavelica']        =   numpy.array([-12.78261, -74.97266])
        self.values['Huanuco']        =   numpy.array([-9.93062, -76.24223])
        self.values['Ica']        =   numpy.array([-14.06777, -75.72861])
        self.values['Huancayo']        =   numpy.array([ -12.06513, -75.20486])
        self.values['Trujillo']        =   numpy.array([ -8.11599, -79.02998])
        self.values['Chiclayo']        =   numpy.array([ -6.77137, -79.84088])
        self.values['Iquitos']        =   numpy.array([ -3.74912, -73.25383])
        self.values['Puerto_M']        =   numpy.array([ -12.59331, -69.18913])
        self.values['Moquegua']        =   numpy.array([-17.19832, -70.93567])
        self.values['Cerro_de_Pasco']        =   numpy.array([-10.66748, -76.25668])
        self.values['Piura']            =   numpy.array([ -5.19449,  -80.63282])
        self.values['Puno']             =   numpy.array([-15.8422,  -70.0199])
        self.values['Moyobamba']        =   numpy.array([-6.03416,  -76.97168])
        self.values['Tacna']            =   numpy.array([-18.01465,  -70.25362])
        self.values['Tumbes']           =   numpy.array([-3.56694,  -80.45153])
        self.values['Pucallpa']         =   numpy.array([-8.37915,  -74.55387]) 
        self.values['Jicamarca']        =   numpy.array([-11.7393521,   -76.704258])
        self.values['San_Martin']       =   numpy.array([ -6.51389,   -76.7408])
    
    def __init__(self,version_IRI):

        
        self.zona_horaria = pytz.timezone('Etc/GMT+5')


        self.lat1=None
        self.lon1=None
        self.version_IRI=version_IRI
        
        self.time_utc = dt.datetime.now() + relativedelta(hours=-5)
        self.time_now = self.time_utc + relativedelta(hours=5)

        self.__config_dict()
        
        self.DIR_FILE=os.path.dirname(__file__)
        
        self.get_date()
        
        
        self.get_mdi()
        self.start_values()
        self.set_saofile()
        
        
        

    
    def start_values(self,):
        '''
        Inicia variables de cálculo para las correciones que se comparten
        entre metodos. 

        Returns
        -------
        None.

        '''
        
        self.ra = 0
        self.ho = 0
        self.ym = 0
        self.dsk = 0
        self.xo = 0
        self.set_center_point()
        
    def set_center_point(self,lat=-11.98, lon=-76.87):
        '''
        Función que permite establecer otro punto de origen
        de los datos proporcionados en el archivos .SAO. 
        Por defecto toma las coordenadas de Jicamarca. 

        Parameters
        ----------
        lat : TYPE, latitud
            DESCRIPTION. The default is -11.98.
        lon : TYPE, longitud
            DESCRIPTION. The default is -76.87.

        Returns
        -------
        None.

        '''
        self.latj=lat
        self.lonj=lon
    
    def run(self,):
        '''
        Metodo que obtiene los valores de la densidad de electrones
        e invoca la función para la obtencion de valores MUF.
        

        Returns
        -------
        None.

        '''

        keys_name = self.values.keys()

        for key in keys_name:

            self.lat1,self.lon1 = self.values[key][0],self.values[key][1]
        
            _,_,_,fof2o=self.get_valueIRI(self.lat1, self.lon1)
            _,_,_,fof2j=self.get_valueIRI(self.latj, self.lonj)
            
            dfiridps = self.get_correction_MUF(fof2j)
            mmufper=self.get_MUF(lat1=self.lat1,lon1=self.lon1,dfiridps=dfiridps)


            
            self.plot_MUF(mmufper=mmufper,
                        lat1=self.lat1,
                        lon1=self.lon1,
                        fof2o=fof2o,
                        dfiridps=dfiridps,
                        key=key)


            self.value_=fof2o+dfiridps
        
    def to_seconds(self,date):
        '''
        Función que calcula la cantidad de segundos transcurridos a la fecha 
        proporcionada.

        Input
        -----
        date = Datetime object

        Returns
        -------
        sec : int
            Cantidad de segundos transcurridos UNIX (1970)

        '''
        
        fecha_referencia = dt.datetime(1970,1,1).replace(tzinfo=None)

        diff = date.replace(tzinfo=None) - fecha_referencia

        segundos = int(diff.total_seconds())

        return segundos
    

    def get_value(self):
        return self.value_        
            
    def plot_MUF(self,mmufper,lat1,lon1,fof2o,dfiridps,key):
        '''
        
        Método que plotea los valores MUF en el mapa.
        
         Parameters
         ----------
         mmufper : TYPE
            DESCRIPTION. Matriz (22,17)dim con valores MUF
        lat1 : TYPE
            DESCRIPTION. Latitud de punto de transmision.
        lon1 : TYPE
            DESCRIPTION. Longitud de punto de transmision
        fof2o : TYPE
            DESCRIPTION. Frecuencia en el punto de transmision
        dfiridps : TYPE
            DESCRIPTION. Valor de correcion

         Returns
         -------
         None.
 
        '''
 
            
           
        
        #abrimos mapa 
        
        pmap = "/images/mapa.png"
      
    
        img = plt.imread(pmap)
        
        fig = plt.figure(figsize=(9,11),dpi=120)
        ax=fig.add_subplot(1,1,1)
        
        ax.imshow(img, extent=[-84, -68, -20, 1])
        
        
        values = dict()
        
        #latitud,longitud
        values['Chachapoyas']   =   numpy.array([-6.23169,-77.86903])
        values['Huaraz']        =   numpy.array([-9.52779, -77.52778])
        values['Abancay']        =   numpy.array([-13.63389,-72.88139])
        values['Arequipa']        =   numpy.array([-16.39889,-71.535])
        values['Ayacucho']        =   numpy.array([-13.15878,-74.22321])
        values['Cajamarca']        =   numpy.array([-7.1637,-78.50027])
        values['Cuzco']        =   numpy.array([-13.52264,-71.96734])
        values['Huancavelica']        =   numpy.array([-12.78261, -74.97266])
        values['Huanuco']        =   numpy.array([-9.93062, -76.24223])
        values['Ica']        =   numpy.array([-14.06777, -75.72861])
        values['Huancayo']        =   numpy.array([ -12.06513, -75.20486])
        values['Trujillo']        =   numpy.array([ -8.11599, -79.02998])
        values['Chiclayo']        =   numpy.array([ -6.77137, -79.84088])
        values['Iquitos']        =   numpy.array([ -3.74912, -73.25383])
        values['Puerto M.']        =   numpy.array([ -12.59331, -69.18913])
        values['Moquegua']        =   numpy.array([-17.19832, -70.93567])
        values['Cerro de Pasco']        =   numpy.array([-10.66748, -76.25668])
        values['Piura']            =   numpy.array([ -5.19449,  -80.63282])
        values['Puno']             =   numpy.array([-15.8422,  -70.0199])
        values['Moyobamba']        =   numpy.array([-6.03416,  -76.97168])
        values['Tacna']            =   numpy.array([-18.01465,  -70.25362])
        values['Tumbes']           =   numpy.array([-3.56694,  -80.45153])
        values['Pucallpa']         =   numpy.array([-8.37915,  -74.55387]) 

        values['Jicamarca']        =   numpy.array([-11.7393521,   -76.704258])
        values['San Martin']       =   numpy.array([ -6.51389,   -76.7408])


        ciudades = values.keys()
        latitudes= numpy.array([x[0] for x in values.values() ])
        longitudes= numpy.array([x[1] for x in values.values() ])

        latitudes2=latitudes+0.3 
        longitudes2=longitudes-0.5
        
 
            
        clats = -20.0 + np.linspace(0,21,22)
        clons = -84.0 + np.linspace(0,16,17)    
        levmuf = np.linspace(0,59,60)/2.0
        levmuf1 = np.linspace(0,29,30)


        CS = plt.contour(clons, clats, mmufper,13,linewidths=0.6, colors = 'black', alpha = .77, levels = levmuf1)

        plt.clabel(CS, fontsize = 7, inline = 1)
        plt.contour(clons, clats, mmufper,5, colors ='black', alpha = .45, linewidths=0.6, 
                    linestyles = 'dashdot',levels=levmuf)
         
        
        
        plt.text(-80,-18, 'OCEANO PACIFICO')
        # lon1,lat1=m(lon1,lat1)
        plt.text(lon1-1, lat1-1, '%5.2f' % (fof2o+dfiridps), color = 'black', alpha = 1,
                fontsize=10,fontweight=970)
        
        plt.scatter(lon1,lat1,
                        color='#102C54',s=35,alpha=0.85,marker="x"
                        )
        
        chain = "MUF - IRI Ver.2020 {}/{:02d}/{:02d} LT {:02d}:{:02d} (UTC-5)".format( self.time_now.year,
                                                                                      int(self.time_now.month),
                                                                                      int(self.time_now.day),
                                                                                      int(self.time_now.hour),
                                                                                      int(0))
        

        
        plt.title(chain)
         
        
        
        self.buf=None
        
        name_fig = f"{key}.png"
        path_img = os.path.join("/images",name_fig)

        fig.savefig(path_img,format="png")

        print(f"Imagen ha sido generada: {path_img}")
     

        fig.clf()
        plt.close('all')


          
          
    def get_buf(self,):
        '''
        Retorna el buffer de bits que contiene el mapa con valores MUF

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        return self.buf
    
    def get_saofile(self,):
        
        '''
        Retorna el nombre del archivo SAO
        '''
        
        return self.saofile
        
    def set_saofile(self,year=2009,mdi=170):
        
        '''
        Permite establecer un archivo SAO diferente
        '''
        
        self.saofile = self.DIR_FILE +'/IRI_files/JI91J_%4d%03d.SAO' %((2005), abs(mdi))
        
    
    def get_MUF(self,lat1,lon1,dfiridps):
        '''
        Obtiene valores MUF en todo el territorio.

        Parameters
        ----------
        lat1 : Latitud
        lon1 : Longitud
        dfiridps : Valor de correción

        Returns
        -------
        mmufper : Matriz MUF.

        '''
       
        mmufper = np.zeros((22,17))
        latmap = [-20,1]
        lonmap = [-84,-68]
        lonigrid = lon1 - int(lon1 - lonmap[0])
        latigrid = lat1 - int(lat1 - latmap[0])
        clat = []
        clon = []
        
        
        for i in range(22):
            clat.append(latigrid + float(i))
        for i in range(17):
            clon.append(lonigrid + float(i))

        for ilon in range(17):
            for jlat in range(22):

                mmufper[jlat,ilon] = self.process_iri_muf(ilon,jlat,dfiridps,lonigrid,latigrid,lat1,lon1)
 
 

        return mmufper
 
                
            

    
    def process_iri_muf(self, ilon,jlat,dfiridps,lonigrid,latigrid,lat1,lon1):
                

        lon2 = lonigrid + ilon
        lat2 = latigrid + jlat
        loni = (lon1+lon2)/2.0
        lati = (lat1+lat2)/2.0

        fof2,hl,_,_ = self.get_valueIRI(lati, loni)
        
        
        ver1 = fof2.shape
        ver2 = hl.shape
        
        
        nmf2p = numpy.nanmax(fof2)
        fof2p = numpy.sqrt((nmf2p)*1e6/1.24*1e-10)
        hmf2p = []
        indhm = []

        #Mejorar 
        for i in range((fof2.shape[0])):
            if(fof2[i] == nmf2p):
                hmf2p.append(hl[i])
                indhm.append(i)


        NhmCs = [ 135.0]
        self.NhmCs1 = max(fof2)
        self.NhmCs2 = max(hmf2p)
        ind1 = indhm[0]-9
        ind2 = indhm[0]+1



        nelfit, chisq = curve_fit(self.parbfunc, 
                                    np.array(hl[ind1:ind2]), 
                                    np.array(fof2[ind1:ind2]),
                                    p0=NhmCs)
        NhmCs = nelfit

        hpbex = []
        for i in range(40):
            
            hpbex.append(i*10.0+80.0)
            
        hpbex = np.array(hpbex)

        nelpb = self.parbfunc(hpbex,NhmCs[0])
    
            
        hop = np.interp(0.0,nelpb[0:21],hpbex[0:21])
        
        self.ym = max(hmf2p) - hop
        ym=self.ym
        
        self.ra = 6370.0
        ra=self.ra
        
        self.ho = hop
        ho=self.ho
        
        self.dsk = self.map_2points(lat1,lon1,lat2,lon2)
        dsk=self.dsk
        
        self.xo = 1.8
        xo=self.xo
        
        self.philim1 = math.acos((xo - math.sqrt(xo*xo-4*xo*xo*ym/(ra+ho)*\
                                            (1-xo*xo*ym/(ra+ho))))/2/\
                            (xo*xo*ym/(ra+ho)))
        
    
            
        self.philim2 = math.atan(math.sqrt(ra/2.0/ho))
        philim2=self.philim2
        
        xlim = (-1*math.cos(philim2)+math.sqrt(math.cos(philim2)**2+4*\
                                                (ym*(1-math.cos(philim2)**2)\
                                                /(ra+ho))))/2/(ym*(1-math.cos(philim2)**2)/(ra+ho))
        
    
        
        mufdsk= self.Mufx2(dsk,xlim-0.06)*(fof2p+dfiridps)


        return mufdsk
                

    def get_correction_MUF(self,fof2j=0):
        '''
        Genera la correcion del MUF mediante el archivo SAO

        Parameters
        ----------
        fof2j : TYPE, Frecuencia

        Returns
        -------
        dfiridps : Valor de correccion

        '''
        
        dps_time, dps_scaled, dps_otf, dpsTimev, gconst = self.get_sao(self.saofile)
        nde = len(dps_time)
        
        nde = 0
        dfiridps = 0
        foF2_dps = 0
        tdps=0
        itdps = 0
        
        if(nde > 1):

            for j in range(nde) :
                if (j == 0) :
                    if (self.utnow < (dps_time[j]+dps_time[j+1])/2.) :
                        dfiridps = dps_scaled[0][0] - fof2j
                        foF2_dps = dps_scaled[0][0]
                        tdps = dps_time[j]
                        itdps = j
                if (j>0) and (j< (nde-1)) :
                    x1 = (dps_time[j-1]+dps_time[j])/2.
                    x2 = (dps_time[j]+dps_time[j+1])/2.
                    if (self.utnow >= x1) and (self.utnow < x2) :
                        dfiridps = dps_scaled[j][0] - fof2j
                        foF2_dps = dps_scaled[j][0] 
                        tdps = dps_time[j]
                        itdps = j
                if (j== (nde-1)) :
                    if (self.utnow >= (dps_time[j]+dps_time[j-1])/2.) :
                        dfiridps = dps_scaled[j][0] - fof2j
                        foF2_dps = dps_scaled[j][0]
                        tdps = dps_time[j]
                        itdps = j
 
        if(nde == 0):
            dfiridps = 0
          
        if nde == 1 : 
            if abs(self.utnow- dps_time[0])< 3600 :
                dfiridps = dps_scaled[0][0] - fof2j

            else :
                dfiridps = 0
  
        return dfiridps

        
        
        
    def get_sao(self,fileopen):
        '''
        Método para la obtención de valores contenidos en el .SAO file

        Parameters
        ----------
        fileopen : Dirección del archivo .SAO
        Returns
        -------
        timed : TYPE
            DESCRIPTION.
        scaledb : TYPE
            DESCRIPTION.
        OTF : TYPE
            DESCRIPTION.
        timev : TYPE
            DESCRIPTION.
        GCONST : TYPE
            DESCRIPTION.

        '''
        
        heightd = 60*[0.0]
        for i in range(60):
            heightd[i] = 15.0*i + 100.0
        
        timed = []
        timev = []
        scaledb = []
        PLASFRE = []
        OTFB = []
        OTF = []
        GCONST = []
        
        self.set_saofile()
        fileopen=self.saofile

        try :

            file1 = open(fileopen, 'r')

        except IOError :

            print ('No hay SAO en la dirección: {}'.format(fileopen))
            return timed, scaledb, OTF, timev, GCONST
        else: 

            while True:
                linea = file1.readline()
                if linea == '':
                    break
                (IDFI, GCONST, SYSDES, IPREF, SCALED, IAF, DTT, OTF, IODF, FTOF, IODF1, \
                FTOF1, IOAE, IODE, FTOE, IXAF, IXDF, FTXF, IXAF1, IXDF1, FTXF1, \
                IXAE, IXDE, FTXE, MEDF, MEDE, MEDES, THF2, THF1, THE, QPCOEF, THVAL, \
                IEDF, IOTSE, IOASE, IODSE, FTOSE, OTHF1, OTE, OTHE, XTF, XTF1, XTE, \
                OTSE, HTAB, FTAB, QL, DL, IEDFTP) = READ_SAO_FILE(file1, linea)
                #IOTF, IOTHF, IOAFB, IOTF1, IOTHF1, IOAF1B, IOAF1B, IOTE, IOTHE, IXTF,
                #IXTF1, IXTE, OTFB, OTF1B, HTABB, FTABB 

                #Time of Measurements
                yyyy = int(IPREF[2]+IPREF[3]+IPREF[4]+IPREF[5])
                doy = int(IPREF[6]+IPREF[7]+IPREF[8])
                mmm = int(IPREF[9]+IPREF[10])
                dom = int(IPREF[11]+IPREF[12])
                hh = int(IPREF[13]+IPREF[14])
                mm = int(IPREF[15]+IPREF[16])
                ss = int(IPREF[17]+IPREF[18]+IPREF[19])

                date = dt.datetime(yyyy,doy,hh,mm,ss)
                timeb = self.to_seconds(date)
                timed.append(timeb)

                #Scaled Ionogram Parameters
                if(IDFI[3] > 0):
                    
                    for i in range(len(SCALED)):
                        if(SCALED[i] == 9999.0): 
                            SCALED[i] = float('nan')
                    if(len(SCALED) < 50):
                        for i in range(50 - len(SCALED)):
                            SCALED.append(float('nan'))
                    scaledb.append(SCALED)
                #N(h) Tabulation
                h = []
                if(IDFI[50] > 0):
                    if(max(HTAB) < max(heightd)):
                        for i in range(len(heightd)):
                            if(heightd[i] <= max(HTAB)): h.append(heightd[i])
                    else:
                        h = heightd

                    plasfreb = InterpolatedUnivariateSpline(HTAB,FTAB)(h)
                    temp_plasfre = []
                    for i in range(len(plasfreb)):
                        temp_plasfre.append(plasfreb[i])
                    if(len(plasfreb) < 60):
                        for i in range(60 - len(temp_plasfre)): 
                            temp_plasfre.append(float('nan'))
                    PLASFRE.append(temp_plasfre)
                #O-trace F2 Points
                #Virtual Height
                if(IDFI[6] > 0):
                    for i in range(len(OTF)):
                        if(OTF[i] == 9999.0): OTF[i] = float('nan')
                    if(len(OTF) < 120): 
                        for i in range(120 - len(OTF)):
                            OTF.append(float('nan'))
                    timev.append(timeb)
                    OTFB.append(OTFB)    
            file1.close()
            
            return timed, scaledb, OTF, timev, GCONST
    
    def set_date(self,year=None,month=None,day=None,hour=None,mm=None):
        '''
        Método que permite configurar fecha y hora distinta a la del sistema

        Parameters
        ----------
        year : Año
        month : Mes
        day : Dia

        Returns
        -------
        None.

        '''

        if hour is None:
            hour = 0
        if mm is None:
            mm = 0 


    
        self.time_now = dt.datetime(year,month,day,hour,mm,tzinfo=self.zona_horaria)
        self.time_utc = self.time_now + relativedelta(hours=+5)
 
        self.generate_dt()
        self.get_mdi()
        self.get_hourmuf()
        self.utnow=self.to_seconds(self.time_now)
    
    def get_date(self):
        '''
        Metodo que obtiene la hora actual del sistema por defecto.

        Returns
        -------
        None.

        '''
 
        self.time_now = dt.datetime.now(dt.timezone.utc)
        self.time_utc = dt.datetime.now(pytz.utc)
        
        self.utnow = self.to_seconds(self.time_now)
        
        self.time_now =dt.datetime.now()
        
 
            
 
 
        
        self.get_hourmuf()
        self.generate_dt()
        
        
    def get_mdi(self):
        '''
        Obtiene el valor mdi

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        
        self.mdi=self.get_doy(self.time_utc)
        
        return self.mdi
    

   
    
    def generate_dt(self,):
        '''
        Metodo que genera el tipo datetime.

        Returns
        -------
        None.

        '''
        pass
        #seteamos el time_utc al iri date 
        # self.iri.set_date(self.time_utc)
    
    def get_doy(self,date):
        '''
        Obtiene el valor doy o numero de días de la fecha proporcionada.

        Parameters
        ----------
        year : year
        month : month
        day : day

        Returns
        -------
        doy : TYPE
            DESCRIPTION.

        '''
        date_init = dt.datetime(date.year,1,1) 


        doy = JULDAY(date)-JULDAY(date_init) + 1
    
        return doy
        
    def get_hourmuf(self):
        '''
        Hora MUF

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        self.hourmuf=self.time_utc.hour+int((self.time_utc.minute+int(self.time_utc.second/60))/60.0)
        

        return self.hourmuf
        
    def get_nh (self,lat,lon):
        
        pass
 
    def get_valueIRI (self,lat,lon):
        '''
        Obtiene diferentes parametros IRI

        Parameters
        ----------
        lat : TYPE
            DESCRIPTION.
        lon : TYPE
            DESCRIPTION.

        Returns
        -------
        fof2 : TYPE
            DESCRIPTION.

        '''
        
        ne,h1=self.get_values_IRI( lat=lat, 
                                  lon=lon)
        
        ne=np.array(ne)
        h1=np.array(h1)
        
 

        nmf2o=np.nanmax(ne)
        fof2 = np.sqrt(float(nmf2o)*1.0e6/1.24*1.0e-10)
        
        return ne,h1,nmf2o,fof2
    
    def get_Nef2(self, lat, lon):
        
        '''
        Obtiene vector de valores Ne
        '''
        ne,h1=self.get_values_IRI(version=self.version_IRI, 
                             lat=lat, 
                             lon=lon,
                             dn=self.time_utc,
                             year=self.year,
                             mdi=self.mdi,
                             hourmuf=self.hourmuf)
        
        return ne,h1
    
    def get_values_IRI(self, lat,lon,heights=None ):
        
        
            ne,heights,date = RUN_IRI(lat,lon,self.time_utc,self.hourmuf)

            self.time_utc = date
            
            return ne,heights
    
 
     
        
        
    def bisrt(self,func,x1,x2,xacc=1.0e-8):

        

        if (func=='DerDrphixy'):
            yi=self.DerDrphixy(x1)
            yf=self.DerDrphixy(x2)
            
        elif (func=='distfx'):
            yi=self.distfx(x1)
            yf=self.distfx(x2)          
     

        if(yi*yf < 0):
            if(yi < 0):
                rtb = x1
                dx = x2 - x1
            else:
                rtb = x2
                dx = x1 - x2
            i = 0
            ymb = yi
            
  
            while (i <= 50 and abs(dx) > xacc and ymb != 0.0):
                dx = dx/2
                xm = rtb + dx
       
                if (func=='DerDrphixy'):
                    ymb=self.DerDrphixy(xm)
                elif (func=='distfx'):
                    ymb=self.distfx(xm)

                if ymb<= 0: 
                    rtb = xm
            return rtb
        else:
            return float('nan')
             

    def parbfunc(self,X,c):
        '''
        Funcion de ajuste
        '''
        a = self.NhmCs1
        b = self.NhmCs2
        t = (X - b)/c
        F = a*(1.0 - t**2)
        return F

    def map_2points(self,lat1,lon1,lat2,lon2):
        '''
        Calcula la distancia de dos puntos

        Parameters
        ----------
        lat1 : TYPE
            DESCRIPTION.
        lon1 : TYPE
            DESCRIPTION.
        lat2 : TYPE
            DESCRIPTION.
        lon2 : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        R = 6378.2064        #Clarke 1866 radius of Earth (km)
        dLat = math.radians(lat2-lat1)
        dLon = math.radians(lon2-lon1)
        a = math.sin(dLat/2.0) * math.sin(dLat/2.0) + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * \
            math.sin(dLon/2.0) * math.sin(dLon/2)
        c = 2.0* math.atan2(math.sqrt(a), math.sqrt(1.0-a))
        return R*c

    def DerDrphixy(self,phii):
       x = self.xo
       ym=self.ym
       ra=self.ra
       ho=self.ho
       
       #print x
       numr = (1 - x*x*ym/(ra + ho)*math.sin(phii)**2 + x*math.cos(phii))
       denor = (1 - x*x*ym/(ra + ho)*math.sin(phii)**2 - x*math.cos(phii))
       numr1 = (1 - ym/(ra+ho)*x*x*(1+math.cos(phii)**2))
       denor1 = ((1-x*x*ym/(ra+ho)*math.sin(phii)**2)**2 - (x*math.cos(phii))**2)
       argLn = numr/denor
       argrc = 1.0/(math.tan(phii)**2) - 2.0*float(ho)/float(ra)
       #print numr, denor, numr1, denor1, argLn, argrc
       if(argrc > 0) and (argLn >= 0):
           F = ra/(ra+ho)*x*ym*(math.log(numr/denor)*math.cos(phii)-2*x* \
                                math.sin(phii)**2*(numr1/denor1)) + 2*ra/ \
                                (math.sin(phii)**2)*(1/math.tan(phii)/math.sqrt(\
                                    1.0/math.tan(phii)**2-2.0*float(ho)/float(ra))-1)
       else:
           F = float('nan')
       #print F
       return F
   
    
    def distfx(self,x):
       self.xo = x
       xo=self.xo
       x = self.xo
       ym=self.ym
       ra=self.ra
       ho=self.ho
       #Se repiten calculos
       self.philim1 = math.acos((xo-math.sqrt(xo*xo-4*xo*xo*ym/(ra+ho)*(1.0-xo*xo*ym/(ra+ho))))\
                           /2/(xo*xo*ym/(ra+ho)))
       self.philim2 = math.atan(math.sqrt(ra/2/ho))
       yminf = self.bisrt(func='DerDrphixy',x1=self.philim1+0.001,x2=self.philim2-0.001)
       distfx_result = self.DistF2xy(yminf)
       return distfx_result    
   
    
    def Mufx2(self,x, y):
     
       dsk = x
       x1 = 1.001
       x2 = y
       if dsk < 80: 
           mufbs = 1.0
           return mufbs
       mufbs = self.bisrt('distfx',x1,x2)

       return mufbs
   
    def DistF2xy(self,phii):
        
       dsk=self.dsk
       x = self.xo 
       ra=self.ra
       ho=self.ho
       ym=self.ym
       
       
       numr = (1 - x*x*ym/(ra + ho)*math.sin(phii)**2 + x*math.cos(phii))
       denor = (1 - x*x*ym/(ra + ho)*math.sin(phii)**2 - x*math.cos(phii))
       argLn = numr/denor
       argrc = 1.0/(math.tan(phii)**2)-2.0*float(ho)/float( ra)
       
       if(argrc > 0) and (argLn >= 0):
           F = ra/(ra + ho)*math.sin(phii)*x*ym*math.log(numr/denor) + \
               2.0*float(ra)*(1.0/math.tan(phii))-2.0*float(ra)*   \
               math.sqrt(1.0/(math.tan(phii)**2)-2.0*float(ho)/float(ra))-dsk
       
       else:
        
            F= float('nan')      #to check use math.isnan(F)
       
       return float(F)
   
    
 


class IRI_old():
    '''
    Clase que permite calcular valores IRI para la version 
    2007 en Fortran. 
    Es necesario descargar los archivos IRI del sitio web 
    y modificar los directorios respectivos.
    
    --------
    Nota
    --------
    No se recomienda usar esta version por ser antigua, en su lugar 
    use la versión 2015 mediante Pyglow o de la codificación original
    en Fortran. 
    '''
    
    def __init__(self, lat,lon,year,mdi,hourmuf):
        self.lat=lat
        self.lon=lon
        self.year=year
        self.mdi=mdi
        self.hourmuf=hourmuf
        self.fileinp='salidaIri.txt'
        self.processing()()
        
    def inpfiliri2(self,coord = 0, lat = -11.98, lon = 283.13, year = 2004,\
                   md = 1109, ut = 1, hour = 11.875, hx = 0, htec_max = 1000.0,\
                   var = 1, vbeg = 80.0, vend = 175.0,stp = 5.0, choice = 1):


        file10 = open('./iri/iri/iniri.txt', 'w')
        file10.write('     '+str(coord)+'     '+str(lat)+'     '+str(lon)+'\n')
        file10.write('     '+str(year)+'     '+str(md)+'     '+str(ut)+'     '+\
                     str(hour)+'\n')
        file10.write('     '+str(hx)+'\n')
        file10.write('     '+str(htec_max)+'\n')
        file10.write('     '+str(var)+'\n')
        file10.write('     '+str(vbeg)+'     '+str(vend)+'     '+str(stp)+'\n')
        file10.write('     '+str(choice)+'\n')
        
        

        if(choice != 0):
            #Compute Ne, T, Ni? (enter: t,t,t  if you want all)
            Neop = 't'
            Top = 't'
            Niop = 't'
            file10.write(Neop+' '+ Top+' '+ Niop+'\n')
            if(Neop == 't'):
                #LAY version: t=standard ver., f=LAY version. {t}
                file10.write('t\n')
                #Ne Topside: t=IRI-2001, f=new options {t}
                file10.write('f\n')
                #Ne Topside: t=IRI01_corrt, f=NeQuick {t}
                file10.write('f\n')
                #Ne Topside: t=F10.7<188, f=unlimited {t}
                file10.write('t\n')
                #foF2 model: t=CCIR, f=URSI-88 {standard:f}
                file10.write('t\n')
                #foF2: t=with storm model, f=without {t}
                file10.write('t\n')
                #F2 peak density or foF2: t=model,f=user input {t}
                file10.write('t\n')
                #F2 peak height or M3000F2: t=model, f=user input {t}
                file10.write('t\n')
                #Bottomside thickness B0: t=Table-option, f=Gulyaeva {t}
                file10.write('t\n')
                #F1 peak density or foF1: t=model,f=user input {t}
                file10.write('t\n')
                Layvers = 0
                #F1 peak height: t=model, f=user input {t}
                if(Layvers):
                    file10.write('t\n')
                #F1: t=with probability model, f=without   {t}
                file10.write('t\n')
                #F1: t=standard probability, f=with L condition {t}
                file10.write('t\n')
                # E peak density or foE: t=model, f=user input {t}
                file10.write('t\n')
                #E peak height: t=model, f=user input {t}'
                file10.write('t\n')
                #D: t=old model, f=new options {t}'
                file10.write('t\n')
            if Top == 't':
                #Te(Ne) model: t=not used, f=correlation is used. {t}
                file10.write('t\n')
                #Te: t=AEROS/ISIS model, f=InterKosmos model {f}
                file10.write('t\n')
            if Niop == 't':
                #Ion comp. model: t=DS78/DY85, f=DS95/TTS05 {f}'
                file10.write('f\n')
                #Ni: t=ion composition in %, f=ion densities in cm-3 {t}
                file10.write('t\n')
            #Equat. Vert. Ion Drift: t=computed, f=not computed {t}
            file10.write('t\n')
            #Spread-F probability: t=computed, f=not computed {t}
            file10.write('t\n')
            #Sunspot index: t=from file, f=user input.  {t}'
            file10.write('t\n')
            #Ionospheric index: t=from file, f=user input. {t}'
            file10.write('t\n')
            #F10.7D Index: t=from file, f=user input {t}'
            file10.write('t\n')
            #UT/LT computation: t=no date change, f=ut_lt subroutine {t}
            file10.write('t\n')
            #Message output unit: t=(UNIT=6), f=(UNIT=12). {t}'
            file10.write('t\n')
        if(hx < 0):
            #Three additional output parameters (number:1-48)'
            #NmF2','hmF2','NmF1','hmF1','NmE','hmE','NmD','hmD',
            #'h05','B0','NVmin','hVtop','Tpeak','hTpek','T300','T400','T600',
            #'T1400','T3000','T120','Ti450','hTeTi','sza','sndec','dip',
            #'dipla','modip','dela','Srise','Sset','seasn','nseas','Rz12',
            #'cov','B1','M3000','TEC','TECtp','IG12','F1prb','F107d',
            #'C1','daynr','vdrft','foF2r','F1noL','F1+L','sp_F'
            # or 0,0,0 for default (F1 probability[40], equ. vert. ion drift[44], foF2_st/_quiet[45])'
            file10.write('     '+str(0)+'     '+str(0)+'     '+str(0)+'\n')
        #Enter 0 to exit or 1 to generate another profile?
        file10.write(str(0))
        file10.close()
        
 
        os.spawnl(os.P_WAIT,"/home/arx/MUF_IRI_MODULE/iri/iri/iri", 'iri')
        return_code = subprocess.call(["/home/arx/MUF_IRI_MODULE/iri/iri/iri", "iri"])
    
    def leeprofiri07(self,):
        file10 = open(self.fileinp, 'r')
        self.hl = []
        mtec = []
        nh = 0
        self.nmf2 = []
        hmf2 = []

        linea = ' '
        for i in range(1,30):
            linea = file10.readline()

        for linea in file10:
            hax = linea.split()[0]
            nmf2ax = linea.split()[1]
            self.hl.append(float(hax))
            nh = nh+1
            #mtec.append(tecax)
            self.nmf2.append(float(nmf2ax))
            #hmf2.append(hmf2ax)
        file10.close()

    def set_fileinp(self,string):
        self.fileinp=string
    
    
    def get_values(self):
        self.processing()
        return self.nmf2, self.hl
    
    def processing(self):
        self.inpfiliri2(0, self.lat, self.lon, self.year, self.mdi, 1, 
                        self.hourmuf, 0, 1000.0, 1, 80, 
                            1000.0, 10.0, 1)
        self.leeprofiri07()


def get_options_MUF(option):
    '''
    Funcion que retorna valores de latitud y longitud de acuerdo a la opción elegida
    en la pagina WEB.

    Parameters
    ----------
    option : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.
    TYPE
        DESCRIPTION.

    '''
    
    latitude=[-13.64, -16.39,13.16,-7.15,-12.06,-10.68,-13.51,-6.23,-6.77,
              -12.79,-12.06,-9.93,-9.53,-14.07,-3.73,-11.95,-12.04,-17.19,
              -6.03,-5.20,-8.37,-12.59,-15.84,-18.01,-8.11,-3.57]
    
    longitude=[-72.88,-71.54,-74.22,-78.51,-77.14,-76.26,-71.98,-77.85,
               -79.85,-74.97,-75.21,-76.23,-77.53,-75.72,-73.24,-76.87,
               -77.03,-70.93,-76.97,-80.63,-74.53,-69.18,-70.02,-70.25,
               -79.03,-80.46]
    
    return latitude[option-1],longitude[option-1]        


# import time

# # Registra el tiempo de inicio
# tiempo_inicio = time.time()

# l,l_=get_options_MUF(7)
# obj=plot_IRI(latitude=l,
#           longitude=l_,
#           version_IRI=2016)


# year=2021
# month=12
# day=2
# hour=12
# mm=0

# obj.set_date(year=year,
#              month=month,
#              day=day,
#              hour=hour,
#              mm=mm)
# obj.set_saofile(year=2012,mdi=338)
# obj.run()
# del obj

# tiempo_fin = time.time()
# tiempo_ejecucion = tiempo_fin - tiempo_inicio
# print("tiempo de ejec",tiempo_ejecucion)
