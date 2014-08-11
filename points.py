import math

def decrad (dms):   
    deg = float(dms[0])
    min = float(dms[1])  
    sec = float(dms[2])
    if deg < 0:
       deg = abs(deg)
       return math.radians(-(deg + (min / 60) + (sec / 3600)))
    else:
    	return math.radians(deg + (min / 60) + (sec / 3600))

class Point (object):
	name = ''

	lattitude = []

	longitude = []

	height = None

	phi = None

	lambda_ = None

	def __init__ (self, name = None, phi = None, lambda_ = None, lattitude = None, longitude = None, height = None):
		if lattitude:
			self.phi = decrad(lattitude)
		if longitude:
			self.lambda_ = decrad(longitude)
		if phi or phi == 0:
			self.phi = (phi)
		if lambda_ or lambda_ == 0:
			self.lambda_ = (lambda_)


	
