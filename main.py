import math
from points import Point

a = float(6378137)
b = float(6356752.31)
f = 1/298.257223563

ctn = Point()
ctn.lattitude = ['-',33,57,05.163]
ctn.longitude = ['',18,28,06.7862]
ctn.height = 83.614

lgbn = Point()
lgbn.lattitude = ['-', 32,58,20.9551]
lgbn.longitude = ['',18,9,27.9488]
lgbn.height = 64.178

# point = [[lat],[long],height]

def decrad (dms):   
    deg = float(dms[1])
    min = float(dms[2])  
    sec = float(dms[3])  
    dd = deg + (min / 60) + (sec / 3600)
    if dms[0] == '-':
       dd = -dd
    return math.radians(dd)

def get_flattening():
    return (a - b) / a

def get_first_eccentricity():
    return math.sqrt((math.pow(a,2) - math.pow(b,2))/math.pow(a,2))

def get_second_eccentricity():
    return math.sqrt((math.pow(a,2) - math.pow(b,2))/math.pow(b,2))

def get_m(phi):
    e = get_first_eccentricity()
    num =  a * (1 - math.pow(e,2))
    denom = math.pow((1 - math.pow(e,2) * math.pow(math.sin(phi),2)), 3 / 2)
    return num / denom

def get_n(phi):
    e = get_first_eccentricity()
    num = a
    denom = math.pow((1 - math.pow(e,2) * math.pow(math.sin(phi), 2)), 1 / 2)
    return num / denom

def get_mean(a, b):
    return (a + b) / 2

def get_delta(a, b):
    return b - a

def get_ssina(mean_phi, delta_phi, delta_lambda):
    m = get_m(mean_phi)
    n = get_n(mean_phi)
    n2 = math.pow(get_second_eccentricity(), 2) * math.pow(math.cos(mean_phi),2)
    t2 = math.pow(math.tan(mean_phi), 2)
    btemp1 = n * math.cos(mean_phi) * delta_lambda
    btemp2 = (n * math.cos(mean_phi) / 24) * delta_lambda
    bbracket1 = (1 + n2 - 9*n2*t2)
    botherpart = (math.pow(delta_phi,2) / (1 + n2))
    blastpart = t2 * math.pow(delta_lambda, 2) * math.pow(math.cos(mean_phi), 2)
    return btemp1 + btemp2 * ( bbracket1 * botherpart - blastpart)

def get_scosa(mean_phi, delta_phi, delta_lambda):
    m = get_m(mean_phi)
    n = get_n(mean_phi)
    n2 = math.pow(get_second_eccentricity(), 2) * math.pow(math.cos(mean_phi),2)
    t2 = math.pow(math.tan(mean_phi), 2)
    atemp1 = m * delta_phi
    atemp2 = (m/24) * delta_phi
    abracket1 = (3 * n2) - (3 * n2 * t2) + 3*(n2*n2)
    aotherpart = math.pow(delta_phi,2)/1+n2
    abracket2 = 2 + (3 * t2) + (2 * n2)
    aotherotherpart = math.pow(delta_lambda,2)*math.pow(mean_phi,2)
    return atemp1 + atemp2 * ( abracket1 * aotherpart - abracket2 * aotherotherpart ) 

def get_s(p1, p2):
    delta_phi =  p2.phi - p1.phi
    delta_lambda = p2.lambda_ - p1.lambda_
    mean_phi = get_mean(p1.phi, p2.phi)
    scosa = get_scosa(mean_phi, delta_phi, delta_lambda)
    ssina = get_ssina(mean_phi, delta_phi, delta_lambda)
    return math.sqrt(math.pow(ssina, 2) + math.pow(scosa, 2))

def get_mean_a(p1, p2):
    delta_phi =  p2.phi - p1.phi
    delta_lambda = p2.lambda_ - p1.lambda_
    mean_phi = get_mean(p1.phi, p2.phi)
    scosa = get_scosa(mean_phi, delta_phi, delta_lambda)
    ssina = get_ssina(mean_phi, delta_phi, delta_lambda)
    return math.atan(ssina / scosa)

def get_delta_a(p1, p2):
	mean_phi = get_mean(p1.phi, p2.phi)
	n = get_n(mean_phi)
	a = get_mean_a(p1, p2)
	s = get_s(p1, p2)
	t = math.tan(mean_phi)
	t2 = math.pow(math.tan(mean_phi), 2)
	n2 = math.pow(get_second_eccentricity(), 2) * math.pow(math.cos(mean_phi),2)
	first = (s/n) * t * math.sin(a)
	second = (math.pow(s,3)/(24 * math.pow(n, 3))) * t
	third = (2 + t2 + (2 * n2)) * math.pow(a, 3)
	fourth = (2 + (7 * n2) + (9 * n2 * t2) + (5 * n2 * n2)) * math.sin(a) * math.pow(math.cos(a), 2)
	return first + second * (third + fourth)

def get_a1(p1, p2):
	delta_a = get_delta_a(p1, p2)
	return get_mean_a(p1, p2) - (0.5 * delta_a)

def get_a2(p1, p2):
	delta_a = get_delta_a(p1, p2)
	return get_mean_a(p1, p2) + (0.5 * delta_a)

print math.degrees(get_a1(ctn, lgbn)) + 360
print math.degrees(get_a2(ctn, lgbn)) + 360 - 180
print get_s(ctn, lgbn)