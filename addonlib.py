#coding: utf-8
import numpy as np
from scipy.signal import butter, lfilter


DEFAULT_STATS = {
    'network': "",
    'location': "LOC",
    "calib": 1.0,
}

# SQL queries
# id, date_e, time_e, lat, lon, energy, timep, times, filename
CREATE_TABLE_SQL = """\
CREATE TABLE data (
id AUTOINCREMENT,
datee varchar(15),
timee varchar(15),
lat Single,
lon Single,
energy Single,
timep varchar(15),
times varchar(15),
filename varchar(100),
distance Single,
snr Single,
azab Single,
azba Single,
profil Integer
)\
"""


INSERT_DATA_SQL = """\
INSERT INTO data (datee, timee, lat, lon, energy, timep, times, filename, snr, distance, azab, azba) \
VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)\
"""


SELECT_DATA_SQL = """\
SELECT id, datee, timee, energy, timep, times, filename, distance, snr, azab, azba
FROM data
WHERE (energy BETWEEN ? AND ?)
AND (distance <= ?)\
"""


UPDATE_EVENT_BY_ID_SQL = """\
UPDATE data
SET SNR=?
WHERE id=?\
"""

# make saving freqs/windows sql queries
WINDOWS = ("P", "S", "C")
CHANNELS = ("N", "E", "Z")
ATTRS = [wind+"_"+ch for ch in CHANNELS for wind in WINDOWS]

CREATE_TABLE_CALC_SQL = """\
CREATE TABLE calc (
calc_id Integer,
Freq Single,
P_N Double,
S_N Double,
C_N Double,
P_E Double,
S_E Double,
C_E Double,
P_Z Double,
S_Z Double,
C_Z Double,
LN_P Double,
LN_S Double)
"""

#"""""".join([wind+"_"+ch+" Double" for ch in CHANNELS for wind in WINDOWS]) +\

INSERT_CALC_SQL = """\
INSERT INTO calc (calc_id, Freq, {}, {}, {}, {}, {}, {}, {}, {}, {}, LN_P, LN_S)
VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
""".format(*ATTRS)


UPDATE_PROFILES_SQL = """\
UPDATE data
SET profil=? WHERE azAB between ? and ?\
"""

# final query: view result
VIEW_RESULT_SQL = """\
SELECT
    data.id, data.datee, data.timee, data.energy, data.lat, data.lon, data.distance, data.snr, data.azab, \
    calc.LN_P, calc.LN_S, calc.Freq, data.profil
FROM data INNER JOIN calc ON data.id = calc.calc_id
WHERE
    ((data.energy) > ?)
    AND ((data.distance) Between ? And ?)
    AND ((data.snr)>?)
    AND ((calc.Freq)=?)
    AND ((data.profil) BETWEEN 1 and 10)\
"""


CREATE_Q_TABLE = """\
CREATE TABLE Q (
id AUTOINCREMENT,
freq Single,
profil Integer,
N Integer,
mp Double,
ms Double,
Qp Double,
Qs Double,
dist_min Single,
dist_max Single
)\
"""


#===


def calcVincentyInverse(lat1, lon1, lat2, lon2):
    """ Computes the distance between two geographic points on the WGS84
    ellipsoid and the forward and backward azimuths between these points.
    """
    # Check inputs
    if lat1 > 90 or lat1 < -90:
        msg = "Latitude of Point 1 out of bounds! (-90 <= lat1 <=90)"
        raise ValueError(msg)
    while lon1 > 180:
        lon1 -= 360
    while lon1 < -180:
        lon1 += 360
    if lat2 > 90 or lat2 < -90:
        msg = "Latitude of Point 2 out of bounds! (-90 <= lat2 <=90)"
        raise ValueError(msg)
    while lon2 > 180:
        lon2 -= 360
    while lon2 < -180:
        lon2 += 360
    # Data on the WGS84 reference ellipsoid:
    a = 6378137.0          # semimajor axis in m
    f = 1 / 298.257223563  # flattening
    b = a * (1 - f)        # semiminor axis
    #
    if (abs(lat1 - lat2) < 1e-8) and (abs(lon1 - lon2) < 1e-8):
        return 0.0, 0.0, 0.0
    # convert latitudes and longitudes to radians:
    lat1 = lat1 * 2.0 * np.pi / 360.
    lon1 = lon1 * 2.0 * np.pi / 360.
    lat2 = lat2 * 2.0 * np.pi / 360.
    lon2 = lon2 * 2.0 * np.pi / 360.
    # 
    TanU1 = (1 - f) * np.tan(lat1)
    TanU2 = (1 - f) * np.tan(lat2)
    #
    U1 = np.arctan(TanU1)
    U2 = np.arctan(TanU2)
    #
    dlon = lon2 - lon1
    last_dlon = -4000000.0  # an impossible value
    omega = dlon
    #
    # Iterate until no significant change in dlon or iterlimit has been reached
    iterlimit = 100
    try:
        while (last_dlon < -3000000.0 or dlon != 0 and
               abs((last_dlon - dlon) / dlon) > 1.0e-9):
            sqr_sin_sigma = pow(np.cos(U2) * np.sin(dlon), 2) + \
                pow((np.cos(U1) * np.sin(U2) - np.sin(U1) *
                     np.cos(U2) * np.cos(dlon)), 2)
            Sin_sigma = np.sqrt(sqr_sin_sigma)
            Cos_sigma = np.sin(U1) * np.sin(U2) + np.cos(U1) * \
                np.cos(U2) * np.cos(dlon)
            sigma = np.arctan2(Sin_sigma, Cos_sigma)
            Sin_alpha = np.cos(U1) * np.cos(U2) * np.sin(dlon) / \
                np.sin(sigma)
            alpha = np.arcsin(Sin_alpha)
            Cos2sigma_m = np.cos(sigma) - \
                (2 * np.sin(U1) * np.sin(U2) / pow(np.cos(alpha), 2))
            C = (f / 16) * pow(np.cos(alpha), 2) * \
                (4 + f * (4 - 3 * pow(np.cos(alpha), 2)))
            last_dlon = dlon
            dlon = omega + (1 - C) * f * np.sin(alpha) * \
                (sigma + C * np.sin(sigma) *
                    (Cos2sigma_m + C * np.cos(sigma) *
                        (-1 + 2 * pow(Cos2sigma_m, 2))))

            u2 = pow(np.cos(alpha), 2) * (a * a - b * b) / (b * b)
            A = 1 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
            B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
            delta_sigma = B * Sin_sigma * \
                (Cos2sigma_m + (B / 4) *
                    (Cos_sigma * (-1 + 2 * pow(Cos2sigma_m, 2)) - (B / 6) *
                        Cos2sigma_m * (-3 + 4 * sqr_sin_sigma) *
                        (-3 + 4 * pow(Cos2sigma_m, 2))))

            dist = b * A * (sigma - delta_sigma)
            alpha12 = np.arctan2(
                (np.cos(U2) * np.sin(dlon)),
                (np.cos(U1) * np.sin(U2) - np.sin(U1) * np.cos(U2) *
                 np.cos(dlon)))
            alpha21 = np.arctan2(
                (np.cos(U1) * np.sin(dlon)),
                (-1 * np.sin(U1) * np.cos(U2) + np.cos(U1) * np.sin(U2) *
                 np.cos(dlon)))
            iterlimit -= 1
            if iterlimit < 0:
                # iteration limit reached
                raise StopIteration
    except ValueError:
        raise StopIteration

    if alpha12 < 0.0:
        alpha12 = alpha12 + (2.0 * np.pi)
    if alpha12 > (2.0 * np.pi):
        alpha12 = alpha12 - (2.0 * np.pi)

    alpha21 = alpha21 + np.pi

    if alpha21 < 0.0:
        alpha21 = alpha21 + (2.0 * np.pi)
    if alpha21 > (2.0 * np.pi):
        alpha21 = alpha21 - (2.0 * np.pi)

    # convert to degrees:
    alpha12 = alpha12 * 360 / (2.0 * np.pi)
    alpha21 = alpha21 * 360 / (2.0 * np.pi)

    return dist, alpha12, alpha21


#from obspy.core.util.geodetics import gps2DistAzimuth
def gps2DistAzimuth(lat1, lon1, lat2, lon2):
    """ Computes the distance between two geographic points on the WGS84
    example, distance between Moscow and UUD:
    lat1, lon1, lat2, lon2
    55.75222, 37.61556, 51.87000, 107.66000
    """
    try:
        values = calcVincentyInverse(lat1, lon1, lat2, lon2)
        if np.alltrue(np.isnan(values)):
            raise StopIteration
        return values
    except StopIteration:
        return (20004314.5, 0.0, 0.0)
    except ValueError as e:
        raise e


def locations2degrees(lat1, long1, lat2, long2):
    """ Function to calculate the great circle distance between two points """
    # Convert to radians
    lat1 = np.radians(lat1)
    lat2 = np.radians(lat2)
    long1 = np.radians(long1)
    long2 = np.radians(long2)
    long_diff = long2 - long1
    gd = np.degrees(np.arctan2(
                    np.sqrt(
                        np.power(np.cos(lat2) * np.sin(long_diff), 2) +
                        np.power(np.cos(lat1) * np.sin(lat2) - np.sin(lat1) *
                            np.cos(lat2) * np.cos(long_diff), 2)),
                    np.sin(lat1) * np.sin(lat2) + np.cos(lat1) * np.cos(lat2) *
                    np.cos(long_diff)))
    return gd


def rotate_NE_RT(n, e, ba):
    """
    Rotates horizontal components N and E in Radial Transversal Component.
    The angle is given as the back-azimuth, that is defined as the angle
    measured between the vector pointing from the station to the source and
    the vector pointing from the station to the North.
    """
    if len(n) != len(e):
        raise TypeError("North and East component have different length.")
    if ba < 0 or ba > 360:
        raise ValueError("Back Azimuth should be between 0 and 360 degrees.")
    r = e * np.sin((ba + 180) * 2 * np.pi / 360) + n * np.cos((ba + 180) * 2 * np.pi / 360)
    t = e * np.cos((ba + 180) * 2 * np.pi / 360) - n * np.sin((ba + 180) * 2 * np.pi / 360)
    return r, t
#========================

def calc_spectrum(y, dt):
    """ calc spectrum """
    _len = y.size
    half = int( np.fix(_len / 2.) )
    freqs = np.fft.fftfreq(_len, d=dt)
    freqs = freqs[:half]
    pz = np.fft.fft(y)
    pz = np.fft.fftshift(pz)
    ampsp = np.abs(pz)
    ampsp = ampsp / _len
    if _len % 2 == 0: ampsp = ampsp[half:]
    else: ampsp = ampsp[half:-1]
    ampsp = 2 * ampsp
    return freqs, ampsp


def butter_bandpass_filter(data, lowcut, highcut, fs, order=4):
    """ bandpass Butterworth filter """
    nyquist = 0.5 * fs
    Wn = np.array([lowcut, highcut]) / nyquist
    b, a = butter(order, Wn, btype='bandpass')
    return lfilter(b, a, data)


def calc_seconds_from_time(_timeStr):
    """ вычисляет время в секундах с начала дня из текстового времени """
    hour, minute, second = map(float, _timeStr.split(":"))
    return hour * 3600 + minute * 60 + second
