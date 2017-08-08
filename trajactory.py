from poliastro.twobody.propagation import propagate
from poliastro.plotting import plot, OrbitPlotter
from poliastro.twobody import Orbit
from poliastro.bodies import Sun, Earth
from poliastro.twobody import angles as ang
from poliastro.iod.izzo import lambert
from astropy import units as u
from astropy import time
from numpy.linalg import norm
import matplotlib.pyplot as plt
import numpy as np
import poliastro
import astropy

#create time and time step
start_date = time.Time('2018-01-01 00:00', scale = 'utc')
step = time.TimeDelta(1.0, format='jd')
end_date = time.Time('2018-12-31 00:00', scale = 'utc')

#create orbit with start_date
#for asteroid
a = (1.11+0.90)/2*u.AU
ecc = 0.104*u.one
inc = 7.77*u.deg
raan = 66.51*u.deg
argp = 307.23*u.deg
nu = np.rad2deg(ang.M_to_nu(np.deg2rad(297.53*u.deg), ecc))
asteroid = Orbit.from_classical(Sun,a,ecc,inc,raan,argp,nu, start_date)
asteroid = propagate(asteroid, 150*3600*24 * u.second)

#for earth
a = 1*u.AU
ecc = 0.0167086*u.one
nu = np.rad2deg(ang.M_to_nu(np.deg2rad(358.617*u.deg), ecc))
inc = 7.155*u.deg
raan = -11.26064*u.deg
argp = 114.207*u.deg
earth = Orbit.from_classical(Sun,a,ecc,inc,raan,argp,nu, start_date)
op = OrbitPlotter()
op.plot(asteroid, label = 'asteroid')
op.plot(earth, label = 'earth')
plt.show()

day = 3600*24 * u.second
hour = 3600*u.second
minute = 60*u.second
tol = 1

pool = {}


def init():
    global flight, old_dv, new_dv, temp, itera, flag
    flight = 15*3600*24 * u.second
    old_dv  = 100000
    new_dv = 100
    temp = old_dv
    flag = True
    itera = 0

while start_date != end_date:
    init()
    while np.absolute(norm(old_dv) - norm(new_dv)) > tol:
        old_dv = temp
        new_ast = propagate(asteroid, flight)
        result = [x for x in lambert(Sun.k, earth.rv()[0], new_ast.rv()[0], flight)]
        dv = result[0][1] - result[0][0]
        new_dv = np.absolute(norm(result[0][0] - result[0][1]))
        flight += hour
        temp = new_dv
        if norm(result[0][1]) > 2000 or norm(result[0][0]) > 2000:
            flag = False
            break
        itera += 1
    if flag:
        pool[new_dv] = [flight, result[0][0], result[0][1]]
    init()
    while np.absolute(norm(old_dv) - norm(new_dv)) > tol:
        old_dv = temp
        new_ast = propagate(asteroid, flight)
        result = [x for x in lambert(Sun.k, earth.rv()[0], new_ast.rv()[0], flight), start_date]
        dv = result[0][1] - result[0][0]
        new_dv = np.absolute(norm(result[0][0] - result[0][1]))
        flight -= hour
        temp = new_dv
        if norm(result[0][1]) > 2000 or norm(result[0][0]) > 2000:
            flag = False
            break
        itera += 1
    if flag:
        pool[new_dv] = [flight, result[0][0], result[0][1],start_date]
    start_date += step
    asteroid = propagate(asteroid, day)
    earth = propagate(earth, day)

with open('outte', 'w+') as out:
    out.truncate()
    for key, value in pool.items():
        out.write('{}\n{}\n'.format(key,value))
