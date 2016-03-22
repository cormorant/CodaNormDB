#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
#
#  CodaNormdb.py
#  
#  Copyright 2014 petr <petr@SEISMOGRAMMA>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
"""
Программа расчета среднеквадратических амплитуд прямых волн и коды сейсмограмм.
Входной формат данных - GSE2 (опционально, добавить MSEED, исходный Байкал-5).

Расчет коды осуществляется по 2-м алгоритмам:
    - алгоритм нормализации.
    - по амплитудным спектрам коды, волн P, S.
"""
APP_NAME = "CodaNormdb"
__version__="0.0.1"
COMPANY_NAME = 'GIN SB RAS'

import os
import sys
import datetime
import string
from obspy.gse2.core import isGSE2, readGSE2
from obspy import UTCDateTime, Trace, Stream
import numpy as np
import ConfigParser
import pyodbc
import xlwt

# grapihcs
import matplotlib.pyplot as plt
from itertools import cycle
COLORS = cycle(("g", "c", "r"))


from BaikalFile import BaikalFile


from addonlib import (DEFAULT_STATS, gps2DistAzimuth, rotate_NE_RT,
    calc_spectrum, butter_bandpass_filter, calc_seconds_from_time,
    CREATE_TABLE_SQL, INSERT_DATA_SQL, SELECT_DATA_SQL, UPDATE_EVENT_BY_ID_SQL,
    CREATE_TABLE_CALC_SQL, INSERT_CALC_SQL, UPDATE_PROFILES_SQL, VIEW_RESULT_SQL)


def module_path():
    if hasattr(sys, "frozen"):
        return os.path.dirname(unicode(sys.executable, sys.getfilesystemencoding( )))
    return os.path.dirname(unicode(__file__, sys.getfilesystemencoding( )))


# get current dir, may vary if run from EXE file
CurrDir = module_path()
#print CurrDir
CONFIG_FILENAME = os.path.join(CurrDir, "codadb.conf")

# path to database folder
DB_FOLDER = os.path.join(CurrDir, "db")

# to seisan data
PATH = os.path.join(CurrDir, "seisan")

# path to data
#PATH = "seisan/hrm_all"
SEISOBR_PATH = "D:/Work/seis/seisobr"

# порядок каналов в файлах Байкал
CHANNELS = ("N", "E", "Z")


# profiles azimutes from station HRM
PROFILES = {
    1: (245, 300), # 40
    2: (320, 360), # 40
    3: (20, 50),   # 25
    # east
    4: (80, 117),  # 35
    # south
    5: (117, 235),
}


LIMITS = {
}


def calculate(Freq, f1, f2, stream, secondsE, secondsP, secondsS, azBA):
    """ расчет параметров коды """
    # load settings
    PLOT = Settings["plot"]
    SD = Settings["sd"]
    KOEF = Settings["koef"]
    rotate = Settings["rotate"] # rotate or not
    # rotate here traces 1 and 2 (N, E)
    if rotate:
        r, t = rotate_NE_RT(stream[0].data, stream[1].data, azBA)
        stream[0].data, stream[1].data = r, t
    # get results for this filename at this freqs
    if PLOT:
        fig, axes = plt.subplots(figsize=(20, 12), nrows=stream.count(), sharex=True)
        fig.suptitle("File {}. Freq {} ({}-{})".format(stream.__str__().splitlines()[1],
            Freq, f1, f2), fontsize=14)
    else:
        axes = np.arange(stream.count())
    # будем считать
    result = []
    for ax, trace in zip(axes, stream.traces):
        # delete mean value
        trace.data -= trace.data.mean()
        # sampling rate, 100 Hz normally
        sr = trace.stats.sampling_rate
        y = butter_bandpass_filter(trace.data, f1, f2, sr)
        if PLOT:
            color = COLORS.next()
            ax.plot(trace.times(), trace.data, color, label=trace.stats.channel, alpha=0.5, lw=0.5)
            # plot filtered
            ax.plot(trace.times(), y, "k", lw=0.75)
            #ax.set_ylim(y.min(), y.max())
            ax.legend()
        # за 0 берем время начала файла
        secondsStart = calc_seconds_from_time(trace.stats.starttime.time.strftime("%H:%M:%S.%f"))
        # уменьшим времена в секундах относительно начала файла
        timeE = secondsE - secondsStart
        timeP = secondsP - secondsStart
        timeS = secondsS - secondsStart
        # calc time of start Coda
        timeCoda = int(np.ceil(KOEF *  (timeS - timeE)) + timeE)
        if PLOT:
            # mark time of Event, P and S time by vertical lines, start -- end
            ax.axvline(x=timeE, linestyle="--", color="y") # Event
            ax.axvline(x=timeP, linestyle="--", color="r") # P
            ax.axvline(x=timeP+SD, linestyle="--", color="r") # P
            ax.axvline(x=timeS, linestyle="--", color="r") # S
            ax.axvline(x=timeS+SD, linestyle="--", color="r") # S
            # mark coda
            ax.axvline(x=timeCoda, linestyle="-", color="r") # coda
            ax.axvline(x=timeCoda+SD, linestyle="-", color="r") # coda
        #=== все вычисления по данной компоненте
        #= window for P
        Pindex1 = int( np.ceil(timeP * sr) )
        Pindex2 = int( Pindex1 + SD * sr )
        Pwindow = y[Pindex1:Pindex2] # P
        #= window for S
        Sindex1 = int( np.ceil(timeS*sr) )
        Sindex2 = int( Sindex1 + SD * sr )
        Swindow = y[Sindex1:Sindex2] # S
        #!!!TODO: проверка что окна не пересекаются (P, S window)
        if Pindex2 > Sindex1:
            # пересекаются
            print("Error: windows for P and S windows intersects! Skipping...")
            return
        # C - кода
        #coda_index1 = int(np.ceil(KOEF * (timeS - timeE)) * sr + np.ceil(timeE * sr))
        coda_index1 = int(timeCoda * sr)
        coda_index2 = int(coda_index1 + SD * sr)
        # вырезаный по индексам массив с кодой
        coda_window = y[coda_index1:coda_index2]
        # check if all size equal
        if not Pwindow.size == Swindow.size == coda_window.size:
            print("Error: Size of windows P and S mismatch!")
            return
        if not coda_window.size:
            print("Error: Not enough time for coda to cut!")
            return
        # check also if size of windows P, S and coda differs
        elif coda_window.size < Pwindow.size:
            print("Error: Size of coda-window is smaller then for P (S).")
            return
        else:
            pass
        #=== spectrums or just normalize
        # считать средне-квадр. значения по спектру
        if "SPECTR" in Settings["algorithm"].upper():
            _freqs, valuesP = calc_spectrum(Pwindow, trace.stats.delta) # P
            _freqs, valuesS = calc_spectrum(Swindow, trace.stats.delta) # S
            _freqs, valuesC = calc_spectrum(coda_window, trace.stats.delta) # C
        # или по исходным данным
        else:
            valuesP, valuesS, valuesC = Pwindow, Swindow, coda_window
        # calc result on windows P, S, coda
        for values in (valuesP, valuesS, valuesC):
            # Среднее квадратическое
            result.append( np.sqrt( np.sum(np.power(values, 2).astype(np.float64))
                / values.size ) )
    # plot details
    if PLOT:
        # plot text with calculated values (P, S, coda)
        ax.text(x=timeP+2, y=y.max(), s="%.4f"%result[-3]) # P
        ax.text(x=timeS+2, y=y.max(), s="%.4f"%result[-2]) # S
        ax.text(x=timeCoda+2, y=y.max(), s="%.4f"%result[-1]) # coda
        #fig.tight_layout(h_pad=0.01)
        ax.set_xlim(timeE - 10, timeCoda + SD + 20)
        plt.show()
        outfilename = "png/{}_{}.png".format(stream[0].stats.starttime, Freq).replace(":", "-")
        #plt.savefig(outfilename)
        plt.close()
    # the end
    return result


def read_config_file(config_filename):
    config = ConfigParser.SafeConfigParser()
    config.read(config_filename)
    # read main options
    section = "main"
    if not config.has_section(section):
        print("No section '{}' in config file {}! Exiting.".format(section,
            config_filename))
        sys.exit(0)
    # read all config in dictionary (default type is str)
    Settings = dict( (k, v) for k, v in config.items(section) )
    # read INT value
    Settings['sd'] = config.getint(section, "sd")
    # read float values
    for item in ("station_lat", "station_lon", 'koef', "vp", "vs"):
        Settings[item] = config.getfloat(section, item)
    # freqs
    Settings["freqs"] = map(float, config.get(section, 'freqs').split())
    # plot, rotate? (bool)
    for item in ("plot", "rotate", "calc"):
        Settings[item] = config.getboolean(section, item)
    #TODO: assert ing CALCQ not in (P, S)
    return Settings


def calc_signal_noise_ratio(stream, secE, secP, secS):
    """ calc singal-to-noise ratio """
    # get_window for noise (before avent 5 seconds)
    trace = stream[-1]
    SD = 5
    sr = trace.stats.sampling_rate
    # уменьшим времена в секундах относительно начала файла
    secondsStart = calc_seconds_from_time(trace.stats.starttime.time.strftime("%H:%M:%S.%f"))
    timeE = secE - secondsStart
    timeP = secP - secondsStart
    timeS = secS - secondsStart
    #?!! What to do if start of file (secondsStart) as after time of Event (timeE)
    if secondsStart > secE:
        #return -1
        # тогда взять другой кусок
        timeP = timeS
    #=== noise
    index1 = int( (timeP - SD - 2) * sr ) # get only 2 seconds for example
    # до события
    index2 = int(timeP * sr)
    noise_window = trace.data[index1:index2]
    if noise_window.size < 1:
        #just get 2 seconds for noise window from start of file
        noise_window = trace.data[:int(2 * sr)]
    # AVG SQRT calc
    noise_A_sqrt = np.sqrt( np.sum(np.power(noise_window.astype(np.float64), 2))
        / noise_window.size )
    #=== S window
    #= window for S
    Sindex1 = int( timeS *sr )
    Sindex2 = Sindex1 + int( (SD - 2) * sr )
    Swindow = trace.data[Sindex1:Sindex2]
    S_A_sqrt = np.sqrt( np.sum(np.power(Swindow.astype(np.float64), 2))
        / Swindow.size )
    return S_A_sqrt / noise_A_sqrt


def create_stream_from_baikal_file(bf, use_coefs=False):
    """ получить stream (XX-5) """
    # header
    header = bf.MainHeader
    # datetime
    date = datetime.datetime(header["year"], header["month"], header["day"])
    delta = datetime.timedelta(seconds=header["to"])
    dt = date + delta
    # make utc datetime
    utcdatetime = UTCDateTime(dt, precision=3)
    # названия каналов в выходном файле (N, E, Z). Грубые каналы отбрасываются
    data_traces = bf.traces.astype(np.float64)[:3]
    # умножить на коэф-ты
    if use_coefs:
        for _i, ch_header in enumerate(bf.ChannelHeaders):
            koef = ch_header['koef_chan']
            data_traces[_i] = data_traces[_i] * koef
    traces = []
    for channel, data in zip(CHANNELS, data_traces):
        # подготовить заголовок
        stats = DEFAULT_STATS.copy()
        stats.update({
            "station": header['station'].upper()[:3],
            'channel': channel,
            'sampling_rate': round( 1./header["dt"] ),
            "delta": header["dt"],
            "npts": data.size,
            'starttime': utcdatetime,
        })
        # создать трассу
        trace = Trace(data=data, header=stats)
        # объединять все в одну трассу
        traces.append(trace)
    # create Stream
    stream = Stream(traces)
    return stream


def setup_db_conn(mdbfile, verbose=False):
    conn_str = "Driver={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=%s;" % mdbfile
    try:
        conn = pyodbc.connect(conn_str)
        cursor = conn.cursor()
    except pyodbc.Error as msg:
        print("Error acessing to mdb file: %s" % msg)
        return None, None
    else:
        if verbose: print("Connected succesfully to database %s" % mdbfile)
        return conn, cursor


def load_data_into_db(data, drop_table=False):
    """ загрузка данных из текстового файла в БД """
    if drop_table:
        # drop table data
        try:
            cursor.execute("DROP TABLE data")# or just delete all from table?
        except pyodbc.ProgrammingError as msg:# table doesn't exists
            print(msg[1].decode("cp1251"))
    # создать таблицу с исходными данными
    try:
        cursor.execute(CREATE_TABLE_SQL)
    except pyodbc.Error as e:
        print(e[1].decode("cp1251"))
    else:
        conn.commit()
    # перед загрузкой данных, выполним расчет расстояний и азимутов
    lat1, lon1 = Settings["station_lat"], Settings["station_lon"]
    print("Loading data..."),
    # заполнить таблицу data
    for line in data:
        params = list(line)
        params[0] = str(params[0])
        # calc distances
        lat, lon = params[2:4]
        dist, azAB, azBA = gps2DistAzimuth(lat1, lon1, lat, lon)
        dist /= 1000.
        # paramsplus SNR (0) + etc
        params += [-1, dist, azAB, azBA]
        # "DATE_E TIME_E LAT LON K TIME_P TIME_S FILENAME"
        try:
            cursor.execute(INSERT_DATA_SQL, list(params))
        except pyodbc.Error as msg:
            print(msg)
        else:
            conn.commit()
    print("succesful.")


def main(**kwargs):
    """ выполняя запрос согласно условиям (К, азимуты и т.д.), получаем список для обработки """
    # выполнить запрос согласно заданным условиям
    params = [ kwargs[k] for k in "K_min K_max max_distance".split() ]
    try:
        cursor.execute(SELECT_DATA_SQL, params)
    except pyodbc.Error as e:
        print(e)
        return
    else:
        result = cursor.fetchall()
    #===
    Freqs, Limits = Settings["freqs"][::2], Settings["freqs"][1::2]
    assert len(Freqs) == len(Limits), "Number of Freqs must be = number of Limits. Check config file!"
    # получили списки для обработки
    # timee, energy, timep, times, filename, distance, snr, azab
    for _id, datee, timee, energy, timep, times, filename, distance, snr, azab, azba in result:
        # открывать исходный файл
        if filename.startswith("20"):
            filename = os.path.join(PATH, filename)
            if not isGSE2(filename):
                print("File {} is not in GSE2 format!")
                raw_input()
                continue
            #    trace.data = trace.data.astype(np.float64) * calib
            stream = readGSE2(filename)
            #TODO: проверка порядка каналов: должно быть N, E, Z
            # проверка количества каналов, д.б 3
            if not len(stream) == 3:
                print("Warning: {} channel(s) in file {}. Must be 3!".format(len(stream), filename))
            # calib - is coefficient
            # умножать не будем, у нас всё равно нормируются значения (P/Z)
            for trace in stream:
                if trace.stats.calib == 0.: trace.stats.calib = 1.
        else:
            # файл может быть в формате Байкал
            filename = os.path.join(SEISOBR_PATH, filename.strip("/"))
            if not os.path.exists(filename):
                print("File %s not found...:(" % filename)
                raw_input()
                continue
            bf = BaikalFile(filename)
            if not bf.valid:
                print("Invalid or corrupt Baikal-5 file %s!" % filename)
            stream = create_stream_from_baikal_file(bf)
        #== convert time from str to time in seconds (time of Event in seconds since start of day)
        secondsE, secondsP, secondsS = map(calc_seconds_from_time, (timee, timep, times))
        # проверки
        if secondsP < secondsE:
            print("Time P ({}) must be after time of Event ({}) for {}".format(secondsP, secondsE, _id))
            continue
        #= считать соотношение сигнал/шум
        SNR = calc_signal_noise_ratio(stream, secondsE, secondsP, secondsS)
        if SNR == 0:
            print _id, datee, timee, energy, timep, times, filename, SNR
            raw_input()
        # nice output
        s = " ".join(map(str, [_id, datee, timee, energy, distance, SNR]))
        sys.stdout.write("\r" + s)
        sys.stdout.flush()
        # save calculated SNR value into db with id
        cursor.execute(UPDATE_EVENT_BY_ID_SQL, (SNR, _id))
        conn.commit()
        # check SNR
        SNR_min = kwargs.get("SNR_min", -1)
        if SNR <= SNR_min: continue
        # для всех частот считать...
        # порядок записи частот: частота порог
        for Freq, plus_min in zip(Freqs, Limits):
            # подсчитать границы от центральной частоты, low and high freq corners
            f1, f2 = Freq - plus_min, Freq + plus_min
            #=== calculations
            result = calculate(Freq, f1, f2, stream,
                secondsE, secondsP, secondsS, azba)
            # if no values return, skip other freqs...
            if result is None:
                if Settings["plot"]: plt.close()
                # не надо считать по другим частотам, если по 1-й неудачно
                break
            #===
            # next step...
            P_NS, S_NS, C_NS, P_EW, S_EW, C_EW, P_Z, S_Z, C_Z = result
            # считать логарифм
            # LN( Ap / (Az * DIST_KM))
            P = np.log([P_Z / (C_Z * distance)])[0]
            # S value
            S = np.log([S_NS / (C_NS * distance)])[0]
            # save it into db
            insert_params = [_id, Freq] + result + [P, S]
            try:
                cursor.execute(INSERT_CALC_SQL, list(insert_params))
            except pyodbc.Error as e:
                print e[1].decode("cp1251")
            conn.commit()
            #=== end calculations
    print


def write_output_excel_file(**kwargs):
    """ writing output data """
    Freqs = Settings['freqs'][::2]
    # собирать значения и сохранять в файл Эксель
    WorkBook = xlwt.Workbook(encoding="utf8")
    # add sheets with Freq as name
    sheets = [WorkBook.add_sheet(str(freq)) for freq in Freqs]
    #
    headers = "id datee timee energy lat lon distance snr azab LN_P LN_S Freq profil".split()
    # заголовки на каждом листе
    for sheet in sheets:
        for col, header in enumerate(headers):
            sheet.write(0, col, header)
    # записать остальные заголовки
    for freq, sheet in zip(Freqs, sheets):
        the_rest = [np.pi, freq, .01, Settings['vp'], xlwt.Formula("N1*O1/(P1*Q1)"),
            .01, Settings['vs'], xlwt.Formula("N1*O1/(S1*T1)")]
        # дописать служебные данные (Pi freq...)
        for col, item in enumerate(the_rest):
            sheet.write(0, col+len(headers), item)
    #===
    # for each freq
    freq_num = 1
    for freq, sheet in zip(Freqs, sheets):
        # выполнить запрос и записать результат
        # params: energy > ? distance Between ? AND ? AND snr>? AND Freq=?
        params = [kwargs[k] for k in "K_min min_distance max_distance SNR_min".split()]
        params += [freq]
        try:
            cursor.execute(VIEW_RESULT_SQL, params)
        except pyodbc.Error as e:
            print(e)
            return
        else:#TODO: fetch only some items; chain or smth
            result = cursor.fetchall()
        #=== iter over result, and writing
        row = 1
        # for every profil, do
        profil_result = dict([(i, [[], [], []]) for i in range(1, max(PROFILES.keys())+1)])
        for items in result:
            # parse items
            _id, _, _, _, _lat, _lon, distance, _, _, LN_P, LN_S, _freq, profil = items
            if profil in range(1, 7):
                profil_result[profil] [0] += [distance]
                profil_result[profil] [1] += [LN_P]
                profil_result[profil] [2] += [LN_S]
            # write em all in a row
            for col, item in enumerate(items):
                sheet.write(row, col, item)
            row += 1
        #===
        for profil in profil_result.keys():
            x, y1, y2 = profil_result[profil]
            x = np.array(x)
            # there are limits for every profil
            #limit_min, limit_max = LIMITS[profil]
            # find where to take values
            #indexes = np.where( (x >= limit_min) & (x <= limit_max))
            #x = x[indexes]
            # для каждого профиля строить график
            fig, axes = plt.subplots(nrows=2)
            for ax, y, name, V in zip(axes, [y1, y2], ["P", "S"], [Settings["vp"], Settings["vs"]]):
                # calc regression (polyfit)
                y = np.array(y)
                #y = y[indexes]
                m, b = np.polyfit(x, y, 1)
                Y = m * x + b
                ax.set_title("Profil {} ({}). Freq = {}".format(profil, name, freq))
                ax.plot(x, y, 'ro', label="x-y")
                ax.plot(x, Y, "-k")
                # add equation label
                ax.legend([r"$y=%.5f x + %.5f$" % (m, b)], framealpha=0.5)
                # Далее рассчитать сами значения добротности Q
                Q = -1 * np.pi * freq / (V * m)
                print freq, profil, x.size,
                print "%.5f" % m,
                print("Q{0}={1:.3f}".format(name, Q))
            #plt.show()
            outimagename = "png/profil_{}__{}_{}.png".format(profil, freq_num, freq)
            plt.savefig(outimagename)
            plt.close()
        freq_num += 1
    # save results
    # filename must consist StationName, algorithm etc
    outfilename = "{station}_{algorithm}__{sd}s.xls".format(**Settings)
    try:
        WorkBook.save(outfilename)
    except IOError as e:
        print(e)


if __name__ == '__main__':
    #===
    # Для каждой записи получить:
    # времена события (задается во входном каталоге)
    # - время вступления волны P для данной станции
    # - время вступления волны S для данной станции
    # - имя файла для обработки, откуда считывать данные
    #=== считывание файла настроек
    try:
        Settings = read_config_file(CONFIG_FILENAME)
    except BaseException as e:
        print(e)
        sys.exit(1)
    # get name if input file
    InputFile = Settings["input_file"]
    if not os.path.exists(InputFile):
        print("Input file not found!")
        sys.exit(1)
    #=== loading data file
    # загрузить входной файл (guess types in file, hoping for the best)
    data = np.genfromtxt(InputFile, dtype=None, autostrip=True, names=True,
        loose=True, invalid_raise=False)
    # имена столбцов (в верхнем регистре)
    names = map(string.upper, data.dtype.names)
    # проверить что все нужные столбцы в файле есть
    for col in ("TIME_E", "TIME_P", "TIME_S", "FILENAME", "LAT", "LON"):
        if not col in names:
            print("Coulmn {} not found in input file. Fix it!".format(col))
            sys.exit(1)
    #=== database and other preparations 
    # connect to database
    database_filename = os.path.join(DB_FOLDER, Settings["database"])
    conn, cursor = setup_db_conn(database_filename)
    # загрузить данные в БД
    if Settings["calc"]:
        #try:
        load_data_into_db(data)
        #except pyodbc.Error as e:
        #    print(e)
        # also lets create (or clean) tables for results
        for query in ("DROP TABLE calc", CREATE_TABLE_CALC_SQL):
            try:
                cursor.execute(query)
            except pyodbc.ProgrammingError:# table exists (or not)
                pass
            conn.commit()
    #===
    # ограничения:
    Restrictions = {
        # К
        'K_min': 7,
        'K_max': 18,
        # Azimuth (AB)
        #'azAB_min': 250,
        #'azAB_max': 280,
        # maximal distance
        'min_distance': 30,
        'max_distance': 350,
        # Signal-to-NoiseRatio
        "SNR_min": 2,
    }
    # main func according given K, az, etc
    try:
        if Settings["calc"]:
            # search with restrictions
            main(**Restrictions)
            # update profile values
            for profil, values in PROFILES.items():
                az_min, az_max = values
                # execute update query
                update_params = (profil, az_min, az_max)
                cursor.execute(UPDATE_PROFILES_SQL, update_params)
                conn.commit()
            # also set 0 where None value
            cursor.execute("UPDATE data SET profil=0 WHERE NOT profil BETWEEN 1 AND 10;")
            conn.commit()
        # also write output Excel file
        write_output_excel_file(**Restrictions)
    finally:
        conn.close()
