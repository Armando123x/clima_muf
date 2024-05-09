import urllib.request
import time 
import os 
import shutil

 
def download(url,ptmp):
    while 1:
        try:
            urllib.request.urlretrieve(url, ptmp)
    
        except Exception as e:
            print(f"Error al descargar el archivo: {e}")
            if attempts==4:
                break 

            attempts +=1 
            time.sleep(10)
        else:
            shutil.copyfile(ptmp,pfile)
            try:
                os.remove(ptmp)
            except: 
                pass
            print("Exito al actualizar archivo ")
            break 

        finally:

            try:
                os.remove(ptmp)
            except: 
                pass


url = "https://chain-new.chain-project.net/echaim_downloads/apf107.dat"
path_iri = '/iri2020'
nombre_archivo = "apf107.dat"   
nombre_archivo_tmp = "tmp.dat"

pfile = os.path.join(path_iri,nombre_archivo)
ptmp = os.path.join(path_iri,nombre_archivo_tmp)

attempts = 0 


download(url,ptmp)




url = "https://irimodel.org/indices/ig_rz.dat"
path_iri = '/iri2020'
nombre_archivo = "ig_rz.dat"   
nombre_archivo_tmp = "tmp.dat"

pfile = os.path.join(path_iri,nombre_archivo)
ptmp = os.path.join(path_iri,nombre_archivo_tmp)

attempts = 0 


download(url,ptmp)