import Metashape
from osgeo import osr
from osgeo import ogr
path = 'C:/Users/jkarl/Downloads/Matador_Emond'
doc = Metashape.Document()
doc.open(path+"/photogrammetry/projectPS.psx")
chunk = doc.chunk