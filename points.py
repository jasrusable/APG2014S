import math

def decrad (dms):   
    deg = float(dms[1])
    min = float(dms[2])  
    sec = float(dms[3])  
    dd = deg + (min / 60) + (sec / 3600)
    if dms[0] == '-':
       dd = -dd
    return math.radians(dd)

class Point (object):

	lattitude = []

	longitude = []

	height = False

	@property
	def phi(self):
		return decrad(self.lattitude)

	@property
	def lambda_(self):
		return decrad(self.longitude)