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

	height = None

	phi = None

	lambda_ = None

	def __init__ (self, phi = None, lambda_ = None, lattitude = None, longitude = None, height = None):
		if lattitude:
			self.phi = decrad(lattitude)
		if longitude:
			self.lambda_ = decrad(longitude)
		if phi or phi == 0:
			self.phi = (phi)
		if lambda_ or lambda_ == 0:
			self.lambda_ = (lambda_)


	
