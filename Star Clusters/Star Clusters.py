from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import sep
import astropy.units as u
from pylab import *
from dust_extinction.averages import GCC09_MWAvg
from scipy import optimize
import pandas as pd
import gzip
import random
from mpl_toolkits.mplot3d import Axes3D


def clean(array):
    for i in range(array.shape[0]):
        array[i, 256] = (array[i, 256 + 1] + array[i, 256 - 1]) / 2
        array[i, 783] = (array[i, 783 + 2] + array[i, 783 - 2]) / 2
        array[i, 784] = (array[i, 784 + 2] + array[i, 784 - 2]) / 2
        array[i, 1002] = (array[i, 1002 + 2] + array[i, 1002 - 2]) / 2


def getData(runs, science=False):
    """"
    Makes a list of all the inputed files and removed dead pixels.
    """
    file = [];
    filedata = []  # make the lists
    # get all the files
    for i in runs:
        file.append(fits.open('fits_files/d{}.fits'.format(i)))

    # get all the files' data
    for i in range(len(file) - 1):
        filedata.append(file[i][0].data)

    filedata = np.array(filedata)  # make it a np array
    filedata = np.delete(filedata, np.s_[1024:], axis=2)  # remove covered pixels

    if not science:
        masterfile = np.median(filedata, axis=0)
        clean(masterfile)
        return masterfile
    else:
        return filedata


def getScience(scienceruns, master_flat, master_bias, runtime):
    sciencedata = getData(scienceruns, science=True)

    clean_flat = master_flat - master_bias
    normalized_flat = clean_flat / (np.sum(clean_flat) / clean_flat.size)
    clean(normalized_flat)
    clean_science = sciencedata - master_bias
    persec_science = clean_science / runtime
    master_science = np.median(persec_science, axis=0)
    calibrated_science = master_science / normalized_flat

    return calibrated_science


def plots(s, f, indecies=[]):
    name = 'calibrated_{}_{}'.format(s, f)
    run = runs['calibrated_{}_{}'.format(s, f)]
    m, std = np.mean(run), np.std(run)
    plt.imshow(run, interpolation='nearest', vmin=m - std, vmax=m + std, origin='lower')
    plt.colorbar()
    plt.xlabel('X', fontsize=15)
    plt.ylabel('Y', fontsize=15)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.title('Calibrated {} {}'.format(s, f), fontsize=20)
    # plt.savefig('{}.pdf'.format(name), dpi=300, bbox_inches='tight')
    # plt.show()

    # Showing what got identified
    data_sub = runs['bkg_sub_{}_{}'.format(s, f)]  # data - background or bkg_sub files
    objects = runs['objects_{}_{}'.format(s, f)]  # objects_ files

    # plot background-subtracted image
    fig, ax = plt.subplots()
    m, std = np.mean(data_sub), np.std(data_sub)
    im = ax.imshow(data_sub, interpolation='nearest', vmin=m - std, vmax=m + std, origin='lower')
    plt.colorbar(im)
    plt.title(name, fontsize=20)

    # plot an ellipse for each object
    index = runs['index_{}_{}'.format(s, f)]
    nindex = [x for x in range(len(objects)) if x not in index]
    if indecies == []:
        if s == 'Calib': r = 35
        if s == 'M67': r = 22
        if s == 'M3': r = 8
        for i in range(len(objects)):
            e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                        width=6 * objects['a'][i],
                        height=6 * objects['b'][i],
                        # width = r,
                        # height = r,
                        angle=objects['theta'][i] * 180. / np.pi)
            e.set_facecolor('none')
            e.set_edgecolor('red')
            ax.add_artist(e)
        # for i in nindex:
        #     e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
        #                 width=6 * objects['a'][i],
        #                 height=6 * objects['b'][i],
        #                 angle=objects['theta'][i] * 180. / np.pi)
        #     e.set_facecolor('none')
        #     e.set_edgecolor('purple')
        #     ax.add_artist(e)
    else:
        j = 0
        for i in indecies:
            e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                        width=6 * objects['a'][i],
                        height=6 * objects['b'][i],
                        angle=objects['theta'][i] * 180. / np.pi)
            e.set_facecolor('none')
            if j % 5 == 0: e.set_edgecolor('red')
            if j % 5 == 1: e.set_edgecolor('blue')
            if j % 5 == 2: e.set_edgecolor('orange')
            if j % 5 == 3: e.set_edgecolor('green')
            if j % 5 == 4: e.set_edgecolor('purple')
            j += 1
            ax.add_artist(e)
    plt.title('Objects {} {}'.format(s, f), fontsize=20)
    plt.xlabel('X', fontsize=15)
    plt.ylabel('Y', fontsize=15)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    # plt.savefig('objects_{}_{}.pdf'.format(s, f), dpi=300, bbox_inches='tight')
    # plt.show()


def makedic():
    # setting all the run numbers
    global runs
    runs = {
        # Calibration Runs
        'bias_runs': np.arange(100, 110),
        'flat_runs_B': np.arange(110, 115),
        'flat_runs_V': np.arange(115, 120),
        'flat_runs_R': np.arange(120, 125),

        # Calib Runs
        'Calib_runs_B': np.arange(126, 131),
        'Calib_runs_V': np.arange(131, 136),
        'Calib_runs_R': np.arange(136, 141),

        # M3 Science Runs
        'M3_runs_B': np.arange(197, 200),
        'M3_runs_V': np.arange(200, 203),
        'M3_runs_R': np.arange(203, 208),

        # M67 Science Runs
        'M67_runs_B': np.arange(141, 146),
        'M67_runs_V': np.arange(146, 151),
        'M67_runs_R': np.arange(151, 156),

        # Calib Apparent Magnitudes and Error
        'appMag_Calib_B': 13.056,
        'appMag_Calib_V': 13.327,
        'appMag_Calib_R': 13.456,

        'appMagErr_Calib_B': 0.002,
        'appMagErr_Calib_V': 0.002,
        'appMagErr_Calib_R': 0.002,

        # Extinction magnitudes
        'atm_ext_B': 0.309,  # magnitude / 1 airmass
        'atm_ext_V': 0.185,
        'atm_ext_R': 0.142,

        'dust_ext_Calib': 0.148,  # transmission per V
        'dust_ext_M3': 0.083,
        'dust_ext_M67': 0.165,

        # Extinction Factors
        'dust_ext_fact_B': 1.296,  # factor to get to the correct wavelengths
        'dust_ext_fact_V': 1.042,
        'dust_ext_fact_R': 0.8617,

        'atm_ext_fact_Calib': 1.224571943283,  # airmasses per source; factor needed to get to correct magnitude
        'atm_ext_fact_M3': 1.038031101227,
        'atm_ext_fact_M67': 1.288133382797,

        # Distances and Error
        'mas_M3': 0.110,  # miliarcsecond
        'mas_M67': 1.1325,

        'masErr_M3': 0.010,
        'masErr_M67': 0.0011,
    }
    runs['master_bias'] = getData(runs['bias_runs'])  # only need this once


def getIndex(s, f):
    objects = runs['objects_{}_{}'.format(s, f)]
    name = 'index_{}_{}'.format(s, f)
    flux = 'flux_{}_{}'.format(s, f)
    runs[name] = []

    if s == 'Calib':  # this just gets the calibration star
        a, b = objects['x'], objects['y']
        index = np.where((a > 450) & (a < 460) & (b > 520) & (b < 540))[0][0]
        runs[name].append(index)
    else:
        for i in range(len(objects)):
            if runs[flux][i] >= 0:  # goes through only good stars to make sure flux is not negative
                xdiff = np.abs(objects['xmax'][i] - objects['xmin'][i])  # spread in x
                ydiff = np.abs(objects['ymax'][i] - objects['ymin'][i])  # spread in y
                a, b = objects['a'][i], objects['b'][i]
                circ = np.abs(a - b)  # test if circular as expected
                secx = objects['x2'][i]  # second detection in x
                secy = objects['y2'][i]  # second detection in y
                npix = objects['npix'][i]  # number of pixels
                # these conditions work well to pick out good stars only
                if s == 'M3' and (xdiff < 44 and ydiff < 44 and circ < 2 and (secx < 22 or secy < 22)):
                    runs[name].append(i)
                if s == 'M67' and (npix > 100 and circ < 2):
                    runs[name].append(i)


def exp(s, f):
    runs['master_flat_{}'.format(f)] = getData(runs['flat_runs_{}'.format(f)])
    if s == 'Calib':
        if f == 'B' or f == 'V':
            rt = 80  # Calib B-80 V-80 R-120 M3 B-250 V-250 R-150  M67 B-120 V-60 R-60
        else:
            rt = 120
    if s == 'M3':
        if f == 'B' or f == 'V':
            rt = 250
        else:
            rt = 150
    if s == 'M67':
        if f == 'V' or f == 'R':
            rt = 60
        else:
            rt = 120

    name = 'calibrated_{}_{}'.format(s, f)  # the source-filter combo we're looking at
    rl = runs['{}_runs_{}'.format(s, f)]  # the runlist for that source
    mflat = runs['master_flat_{}'.format(f)]  # the master flat for the right filter
    mbias = runs['master_bias']

    # get and clean science data
    runs[name] = getScience(rl, mflat, mbias, rt)
    clean(runs[name])

    # get the background
    bkg = 'bkg_{}_{}'.format(s, f)
    runs[bkg] = sep.Background(runs[name])

    # subtract the background from the image
    bkg_sub = 'bkg_sub_{}_{}'.format(s, f)
    runs[bkg_sub] = runs[name] - runs[bkg]

    # detect objects
    objects = 'objects_{}_{}'.format(s, f)
    runs[objects] = sep.extract(runs[bkg_sub], 2, err=runs[bkg].globalrms)
    # get flux and flux err
    flux = 'flux_{}_{}'.format(s, f)
    fluxerr = 'fluxerr_{}_{}'.format(s, f)
    if s == 'Calib': radius = 35
    if s == 'M67': radius = 22
    if s == 'M3': radius = 8
    runs[flux], runs[fluxerr], flag = sep.sum_circle(runs[bkg_sub], runs[objects]['x'], runs[objects]['y'],
                                                     radius, err=runs[bkg].globalrms, gain=1.0)

    getIndex(s, f)  # picks out the good stars
    index = runs['index_{}_{}'.format(s, f)]
    # rm, j = [], 0
    # for i in index:
    #     # remove the objects if the flux is less than 0
    #     if runs[flux][i] < 0: # goes through only good stars to make sure flux is not negative
    #         rm.append(i)
    #         j += 1
    # print(s, f, j)
    # index = [x for x in index if x not in rm] # removes index for all fluxes below 0
    # print('objects lost', len(runs[objects]) - len(index), 'of', len(runs[objects]), '. Remaining', len(index))
    runs[flux] = runs[flux][index]
    runs[fluxerr] = runs[fluxerr][index]
    ##### After this point, index only needs to be used if objects is directly used ######
    # plots(s,f)


def filterWL(wavelengths, transmissions):
    plt.plot(wavelengths, transmissions)
    y = plt.gca().get_ylim()
    plt.axhline(y[1] / 2, color='red')
    plt.show()


def quaderr(d1, d2):
    return np.sqrt(d1 ** 2 + d2 ** 2)


def getMag(s, f):
    flux = runs['flux_{}_{}'.format(s, f)]
    fluxerr = runs['fluxerr_{}_{}'.format(s, f)]
    aefac = runs['atm_ext_fact_{}'.format(s)];
    ae = runs['atm_ext_{}'.format(f)]
    defac = runs['dust_ext_fact_{}'.format(f)];
    de = runs['dust_ext_{}'.format(s)]
    appMag = 'appMag_{}_{}'.format(s, f)
    appMagErr = 'appMagErr_{}_{}'.format(s, f)

    # print(s,f, 'ATM', ae * aefac)
    # print(s,f, 'Dust', de * defac)

    if s == 'Calib':  # get the zero point magnitude to use for the others
        zpMag = 'zpMag_{}'.format(f)
        zpMagErr = 'zpMagErr_{}'.format(f)
        runs[zpMag] = runs[appMag] + 2.5 * np.log10(flux) + ae * aefac + de * defac
        runs[zpMagErr] = quaderr(runs[appMagErr], 2.5 * fluxerr / (flux * np.log(10)))
        # print(f, 'zpmag', runs[zpMag])
        # print(f, 'zpmagerr', runs[zpMagErr])
        # print(s,f, 'flux', flux)
        # print(s,f, 'fluxerr', fluxerr)
    else:  # get the apparent and absolute magnitudes for other s,f combos
        zpMag = runs['zpMag_{}'.format(f)]
        zpMagErr = runs['zpMagErr_{}'.format(f)]

        runs[appMag] = np.zeros(len(flux))
        runs[appMagErr] = np.zeros(len(flux))
        # get the apparent magnitude for each star and puts it in to an array with key appMag_s_f
        runs[appMag] = - 2.5 * np.log10(flux) - ae * aefac - de * defac + zpMag
        runs[appMagErr] = quaderr(zpMagErr, 2.5 * fluxerr / (flux * np.log(10)))

        # get the absolute magnitude and put it into an array with key absMag_s_f
        absMag = 'absMag_{}_{}'.format(s, f)
        absMagErr = 'absMagErr_{}_{}'.format(s, f)
        runs[absMag] = np.zeros(len(flux))
        runs[absMagErr] = np.zeros(len(fluxerr))

        mas = runs['mas_{}'.format(s)]
        masErr = runs['masErr_{}'.format(s)]
        d = 1 / (mas * 0.001)
        print(s, d)
        runs[absMag] = runs[appMag] - 5 * np.log10(d) + 5
        derr = masErr / (mas ** 2 * 0.001)
        # derr = 1 / (masErr * 0.001)
        runs[absMagErr] = quaderr(runs[appMagErr], 5 * derr / (d * np.log(10)))


def HR(s, f1, f2, isofit=False, file=''):
    # div = np.min(runs['objects_{}_{}'.format(s, f1)]['npix'])/2  # gets the psf of the smallest stars
    div = 5
    # print('div',div)
    # good indices
    index1 = runs['index_{}_{}'.format(s, f1)]
    index2 = runs['index_{}_{}'.format(s, f2)]
    # filter objects
    objects1 = runs['objects_{}_{}'.format(s, f1)][index1]
    objects2 = runs['objects_{}_{}'.format(s, f2)][index2]
    # filter magnitudes
    absMag1 = runs['absMag_{}_{}'.format(s, f1)]
    absMag2 = runs['absMag_{}_{}'.format(s, f2)]
    # filter magnitude errs
    absMagErr1 = runs['absMagErr_{}_{}'.format(s, f1)]
    absMagErr2 = runs['absMagErr_{}_{}'.format(s, f2)]

    absMagDiff = []  # difference of the starts will go here. This is the x-axis.
    absMagUsed1 = []  # the absolute magnitude in B filter. This is the y-axis.

    absMagErrUsed1 = []  # magnitude error in the 1st filter. This is y-axis error
    absMagErrUsed2 = []  # the quaderr of these is the x-axis error

    indices1, indices2 = [], []  # used for plotting different colored circles

    # gets the coords of the stars' centers
    x1, y1 = objects1['x'], objects1['y']
    x2, y2 = objects2['x'], objects2['y']
    for i in range(len(objects1)):
        for j in range(len(objects2)):
            # compares stars between pictures and claims they are the same if they are within div pixels of each other
            if x1[i] - div < x2[j] < x1[i] + div and y1[i] - div < y2[j] < y1[i] + div:
                diff = absMag1[i] - absMag2[j]  # takes the diff in magnitudes of the stars

                absMagDiff.append(diff)  # adds it to the x-list
                absMagUsed1.append(absMag1[i])  # adds the mag to the y-list
                # gets the corresponding errs
                absMagErrUsed1.append(absMagErr1[i])
                absMagErrUsed2.append(absMagErr2[j])
                # gets the indices
                indices1.append(i)
                indices2.append(j)

    # print(indices1)
    indices1 = [index1[i] for i in indices1]
    indices2 = [index2[i] for i in indices2]
    # makes them np arrays
    absMagDiff = np.array(absMagDiff)  # x-axis
    absMagUsed1 = np.array(absMagUsed1)  # y-axis
    absMagErrUsed1 = np.array(absMagErrUsed1)
    absMagErrUsed2 = np.array(absMagErrUsed2)
    absMagDiffErr = quaderr(absMagErrUsed1, absMagErrUsed2)
    # plots(s,f1,indices1)
    # plots(s,f2,indices2)

    fig, ax = plt.subplots()

    if s == 'M3':
        # print(np.where(absMagDiff<-1))
        a = [939, 1016, 643]
        absMagDiff = np.delete(absMagDiff, a)
        absMagDiffErr = np.delete(absMagDiffErr, a)
        absMagUsed1 = np.delete(absMagUsed1, a)
        absMagErrUsed1 = np.delete(absMagErrUsed1, a)

    plt.scatter(absMagDiff, absMagUsed1, color='black', zorder=10)
    plt.errorbar(absMagDiff, absMagUsed1, xerr=absMagDiffErr, yerr=absMagErrUsed1, ls='none', zorder=0)
    plt.title("Hertzsprung-Russell Diagram for {}".format(s), fontsize=20)
    if isofit:
        xaxis, yaxis, age = iso(f1, f2, file);
        for i in range(len(xaxis)):
            if s == 'M3' and ((10 < age[i] < 10.23 and file == '001') or (
                    10.1 >= age[i] >= 10 and file == '0004')):  # these are the best for M3 #10.05
                plt.plot(xaxis[i], yaxis[i], zorder=20, label='10^{} yr'.format(age[i]))
            if s == 'M67' and ((9.5 < age[i] < 9.66 and file == '019m') or (
                    9.4 < age[i] < 9.7 and file == '030')):  # these are the best for M67 9.75
                plt.plot(xaxis[i], yaxis[i], zorder=20, label='10^{} yr'.format(age[i]))
        plt.title("HR diagram for {} with Isocrone {}".format(s, file))
        plt.legend()
    plt.gca().invert_yaxis()
    plt.xlabel("{} - {}".format(f1, f2), fontsize=15)
    plt.ylabel("Absolute Magnitude in {}".format(f1), fontsize=15)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    if isofit:
        plt.savefig("Hertzsprung-Russell Diagram for {} with Isocrone {}.pdf".format(s, file), dpi=300,
                    bbox_inches='tight')
    else:
        plt.savefig("HR for {} werr.pdf".format(s), dpi=300, bbox_inches='tight')
    plt.show()


def iso(f1, f2, M):  # B,V,R = 8,9,10
    file = "C:/Users/antho/OneDrive/Desktop/135/Experiment 1/Isocrones/isoc_z{}.dat.gz".format(M)
    # file = "C:/Users/antho/OneDrive/Desktop/135/Experiment 1/Isocrones/isochrones.tar.gz"
    data = np.loadtxt(file)

    if f1 == 'B': i = 8
    if f1 == 'V': i = 9
    if f1 == 'R': i = 10

    if f2 == 'B': k = 8
    if f2 == 'V': k = 9
    if f2 == 'R': k = 10

    u, index = np.unique(data[:, 0], return_index=True)
    # u, index = np.unique(data, return_index=True)
    # print(u)
    # print(index)
    mag, magDiff, age = [], [], []
    for j in range(len(index)):
        c = index[j]
        if j != len(index) - 1:
            n = index[j + 1] - 1
        else:
            n = len(index)
        age.append(data[c, 0])
        # print(age[j])
        magDiff.append(data[c:n, i] - data[c:n, k])
        mag.append(data[c:n, i])
    return magDiff, mag, age  # xaxis, yaxis


def extra(num_walks, num_steps, step_length=0.5, radius=69.634e9):
    final_x = np.zeros(num_walks)  # array of final displacement
    for iwalk in range(num_walks):
        current_x = 0  # starts at the origin
        for istep in range(num_steps):
            current_x += random.choice([-step_length, step_length])  # randomly goes forward or backward
        final_x[iwalk] = current_x  # set the final displacement

    plt.hist(final_x, bins=30)
    plt.xlabel('Displacement (cm)', fontsize=15)
    plt.ylabel('Number of Walks Ending Here', fontsize=15)
    plt.title('Displacement of Random Walk', fontsize=20)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.savefig("1d random walk.pdf", dpi=300, bbox_inches='tight')
    plt.show()

    ave = np.average(final_x)  # take the average
    rms = np.sqrt(np.sum(final_x ** 2) / num_walks)

    # print("With {} trials, {} steps per trial, and a step of {} cm, the average displacement is {} cm and the rms displacements is.".format(num_walks, num_steps, step_length, ave))
    print('Average Displacement = {} cm'.format(ave))
    print('RMS Displacement = {} cm'.format(rms))
    print()

    td_ave = ave / np.sqrt(3)
    td_rms = rms / np.sqrt(3)
    c = 3e8
    tot_steps = radius / np.abs(td_ave) * num_steps
    print(tot_steps / 1e12, "trillion steps")
    tot_time = tot_steps * 0.5 / c
    print("total time to escape is {} hours".format(tot_time / 60 / 60))

    # Escaping the Sun
    # for iwalk in range(1):
    #     current_x = 0
    #     current_y = 0
    #     current_z = 0
    #     xs,ys,zs = [], [], []
    #     i = 0
    #     limit = 1000000
    #     while i < limit:
    #         xs.append(current_x)
    #         ys.append(current_y)
    #         zs.append(current_z)
    #         theta = random.uniform(0, np.pi)
    #         phi = random.uniform(0, 2. * np.pi)
    #         x_step = step_length * np.sin(theta) * np.cos(phi)
    #         y_step = step_length * np.sin(theta) * np.sin(phi)
    #         z_step = step_length * np.cos(theta)
    #
    #         current_x += x_step
    #         current_y += y_step
    #         current_z += z_step
    #         i += 1
    #
    #     # Create the 3D plot
    #     fig = plt.figure()
    #     ax = fig.add_subplot(111, projection='3d')
    #     ax.plot(xs, ys, zs)
    #     ax.scatter(0,0,0,color='red',s=50)
    #     ax.scatter(xs[-1], ys[-1], zs[-1], color='orange',s=50)
    #
    #     # Set labels and title
    #     ax.set_xlabel('X (cm)',fontsize=15)
    #     ax.set_ylabel('Y (cm)',fontsize=15)
    #     ax.set_zlabel('Z (cm)',fontsize=15)
    #     ax.set_title('Random Walk of Photons in the Sun',fontsize=20)
    #     plt.savefig("3d sun walk.pdf", dpi=300, bbox_inches='tight')
    #     plt.show()
    #     print(xs[-1], ys[-1], zs[-1])

    # # Escaping the Sun
    # for iwalk in range(1):
    #     current_x = 0
    #     current_y = 0
    #     current_z = 0
    #     r = np.sqrt(current_z**2 + current_y**2 + current_x**2)
    #     i = 0
    #     while r < radius:
    #         theta = random.uniform(0, np.pi)
    #         phi = random.uniform(0, 2. * np.pi)
    #         x_step = step_length * np.sin(theta) * np.cos(phi)
    #         y_step = step_length * np.sin(theta) * np.sin(phi)
    #         z_step = step_length * np.cos(theta)
    #
    #         current_x += x_step
    #         current_y += y_step
    #         current_z += z_step
    #         r = np.sqrt(current_z ** 2 + current_y ** 2 + current_x ** 2)
    #         i += 1
    #
    #         if i % 1e6 == 0:
    #             print(i, r)
    #
    #     print('After {} steps, the photon escaped the sun!'.format(i))


def dustExtinction():
    # generate the curves and plot them
    x = np.arange(0.3, 10.0, 0.1) / u.micron
    ext_model = GCC09_MWAvg()
    plt.plot(1. / x, ext_model(x))
    plt.plot(0.550, 1, 'o', label='A($/lambda$) = A(V)')
    plt.axvline(.432626834, color='blue')
    plt.axvline(.5324672, color='green')
    plt.axvline(.624324, color='red')
    title('Normalized Extinction')
    xlabel('$/lambda$ [$/mu m$]')
    ylabel('$A(/lambda)/A(V)$')
    plt.xscale('log')
    plt.xlim(0.09, 4.0)
    legend(loc='best')
    # tight_layout()
    plt.show()


if __name__ == "__main__":
    ########## Done With this Stuff #############
    # get wavelengths of filter
    # B = np.arange(3000,5300,100); B_trans = np.array([.01,.01,.01,.01,.01,.01,.02,.09,.25,.4,.51,.58,.63,.66,.66,.63,.58,.47,.33,.19,.1,.04,.01])
    # V = np.arange(4700,6400,100); V_trans = np.array([.01,.02,.38,.79,.88,.88,.85,.79,.72,.61,.5,.37,.24,.15,.09,.05,.02])
    # R = np.arange(5400,7600,100); R_trans = np.array([.01,.02,.15,.47,.71,.73,.72,.71,.68,.64,.59,.54,.49,.44,.39,.34,.29,.24,.21,.15,.12,.11])
    #
    # filterWL(B, B_trans) # 4326.26834
    # filterWL(V, V_trans) # 5324.672
    # filterWL(R, R_trans) # 6243.24
    # dustExtinction()

    ########### Current Needs ##############
    sources = ['Calib', 'M3', 'M67']
    filters = ['B', 'V', 'R']
    makedic()

    for s in sources:
        for f in filters:
            # print(s,f)
            exp(s, f)
            getMag(s, f)

    # HR(sources[1],filters[0],filters[2]) # M3 2268, 2007, 2142 -> B-R
    # HR(sources[2], filters[1], filters[2]) # M67 191, 209, 192 -> B-V

    # Zs = ['0','0001','0004','001','004','008','019','030','019nov',
    #       '008s','019s','040s','070s','008a','019a','040a','070a','004m','008m','019m']
    # Zs = ['001','0004'] # M3
    # Zs = ['019m','030'] # M67
    # for z in Zs:
    # iso(filters[0], filters[2],z)
    # print(z)
    # HR(sources[1], filters[0], filters[2],True,z) # '0004' for 10.0 or 10.05 Gyr
    # HR(sources[2], filters[1], filters[2],True,z)

    ########## Extra Credit ##############
    # extra(1000, 10000)