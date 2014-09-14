import math
from points import Point

a = float(6378137)
b = float(6356752.31)
f = 1/298.257223563

points = {}

def read_file(path):
    with open(path) as f:
        content = f.readlines()
        for line in content:
            part = line.rstrip('\n').split(' ')
            points[part[0]] = Point(lattitude=part[1].split(','), longitude=part[2].split(','), height=part[3])

read_file('coords.txt')

def write_file():
    file = open("new_coords.txt", "w")
    for pointname in points:
        file.write(str(pointname) + ' ' + str(points[pointname].phi) + ' ' + str(points[pointname].lambda_) + ' ' + str(points[pointname].height))
        file.write("\n")
    file.close()

write_file()

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

def to_positive_radian(radian):
    if radian < 0:
        return radian + 2 * math.pi
    else:
        return radian

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

def get_delta_a(p1, p2, s = None):
    if not s:
        s = get_s(p1, p2)
    mean_phi = get_mean(p1.phi, p2.phi)
    n = get_n(mean_phi)
    a = get_mean_a(p1, p2)
    t = math.tan(mean_phi)
    t2 = math.pow(math.tan(mean_phi), 2)
    n2 = math.pow(get_second_eccentricity(), 2) * math.pow(math.cos(mean_phi), 2)
    first = (s/n) * t * math.sin(a)
    second = (math.pow(s,3)/(24 * math.pow(n, 3))) * t
    third = (2 + t2 + (2 * n2)) * math.pow(a, 3)
    fourth = (2 + (7 * n2) + (9 * n2 * t2) + (5 * n2 * n2)) * math.sin(a) * math.pow(math.cos(a), 2)
    return first + second * (third + fourth)

def get_delta_phi(p1, p2, s = None):
    if not s:
        s = get_s(p1, p2)
    a = get_mean_a(p1, p2)
    mean_phi = get_mean(p1.phi, p2.phi)
    m = get_m(mean_phi)
    t2 = math.pow(math.tan(mean_phi), 2)
    n2 = math.pow(get_second_eccentricity(), 2) * math.pow(math.cos(mean_phi),2)
    first = (s / m) * math.cos(a)
    second = (math.pow(s, 3) / (24 * math.pow(m, 3))) * (1 / math.pow((1 + n2), 2))
    third = (2 + (3 * t2) + (2 * n2)) * math.pow(math.sin(a), 2) * math.cos(a)
    fourth = 3 * ((n2 * -1) + (3 * n2 * t2) - (n2 * n2)) * math.pow(math.cos(a), 3)
    return first + second * (third + fourth)

def get_delta_lambda(p1, p2):
    mean_phi = get_mean(p1.phi, p2.phi)
    s = get_s(p1, p2)
    n = get_n(mean_phi)
    a = get_mean_a(p1, p2)
    t2 = math.pow(math.tan(mean_phi), 2)
    n2 = math.pow(get_second_eccentricity(), 2) * math.pow(math.cos(mean_phi), 2)
    first = (s / (n * math.cos(mean_phi))) * math.sin(a)
    second = (math.pow(s, 3) / (24 * math.pow(n, 3) * math.cos(mean_phi)))
    third = t2 * math.pow(math.sin(a), 3)
    fourth = (-1 - n2 + (9 * n2 * t2)) * math.sin(a) * math.pow(math.cos(a), 2)
    return first + second * (third + fourth)

def get_a1(p1, p2):
	delta_a = get_delta_a(p1, p2)
	return get_mean_a(p1, p2) - (0.5 * delta_a)

def get_a2(p1, p2):
	delta_a = get_delta_a(p1, p2)
	return get_mean_a(p1, p2) + (0.5 * delta_a)





def solve_inverse(p1, p2):
    s = get_s(p1, p2)
    a1 = to_positive_radian(get_a1(p1, p2))
    a2 = to_positive_radian(get_a2(p1, p2))
    # Sort out 180 stuff
    if p2.phi > p1.phi:
        a2 = to_positive_radian(a2 - math.pi)

    return math.degrees(a1), math.degrees(a2), s

def solve_direct(p1, s, a1):
    p2 = Point(phi = -0.6, lambda_ = 0)
    a2 = 0
    dec = 4
    while round((p2.phi - p1.phi), dec) != round(get_delta_phi(p1, p2, s), dec):
        p2.phi = p2.phi + math.pow(10, -dec)
        
    print p2.phi
    dec = 4
    while round((a2 - a1), dec) != round(get_delta_a(p1, p2, s), dec):
        a2 = a2 + math.pow(10, -dec)
    print a2


#print solve_direct(points['ctn'], 10, 2)
#print solve_inverse(points['ctn'], points['lgbn'])






####### Assignment 2 #######

def get_ellipsoidal_coords2(point):
    get_cartesian_coords(point)
    n = get_n(point.phi)
    p = (n + point.height) * math.cos(point.phi)
    # p = math.sqrt(math.pow(x, 2) + math.pow(y, 2))
    e2 = get_second_eccentricity()
    h = point.height
    z = point.z
    lambda_ = math.pow(math.tan(point.y / point.x), -1)
    #phi = math.pow(math.tan((z/p)*(1-(e2 * n)/(N + h))), -1)

    theta = math.pow(math.tan((z/p)*(a/b)), -1)
    # Direct solution
    phi = math.pow((z + e2 * b * math.pow(math.sin(theta), 3))/ (p - e2 * a * math.pow(math.cos(theta, 3))))

    return phi, lambda_, h

def get_cartesian_coords(point):
    cos = math.cos
    sin = math.sin
    tan = math.tan
    n = get_n(point.phi)
    e2 = get_second_eccentricity()
    h = point.height
    x = (n + h) * cos(point.phi) * cos(point.lambda_)
    y = (n + h) * cos(point.phi) * sin(point.lambda_)
    z = (n * (1 - e2) + h) * sin(point.phi)
    return x, y, z

def get_ellipsoidal_coords(point):
    x,y,z = get_cartesian_coords(point)
    cos = math.cos
    sin = math.sin
    tan = math.tan
    h = point.height
    n = get_n(point.phi)
    e2 = get_second_eccentricity()
    p = (n + h) * cos(point.phi)

    theta = tan((z/p)*(a/b))**-1

    lambda_ = tan(y / x)**-1
    phi = tan((z + (e2 * b * (sin(theta)**3)))/(p - e2 * a * (cos(theta))**3))**-1
    return point.phi, point.lambda_, point.height

def get_local_geodetic(point):
    x,y,z = get_cartesian_coords(point)
    cos = math.cos
    sin = math.sin
    tan = math.tan
    import numpy
    matrix = numpy.matrix
    ry = point.phi - math.pi / 2
    rz = point.lambda_ - math.pi
    rx = 0
    Py = matrix([[1, 0, 0],
                [0, -1, 0],
                [0, 0, 1],
                ])
    Rx = matrix([[1, 0, 0],
                [0, cos(rx), sin(rx)],
                [0, -sin(rx), cos(rx)],
                ])
    Ry = matrix([[cos(ry), 0, -sin(ry)],
                [0, 1, 0],
                [sin(ry), 0, cos(ry)],
                ])
    Rz = matrix([[cos(rz), sin(rz), 0],
                [-sin(rz), cos(rz), 0],
                [0, 0, 1],
                ])
    Rgs = matrix([[x], [y], [z]
                ])
    rlgs = Py * Ry * Rz * Rgs
    return rlgs

#print get_cartesian_coords(points['CTWN'])
#print get_ellipsoidal_coords(points['CTWN'])
print get_local_geodetic(points['CTWN'])