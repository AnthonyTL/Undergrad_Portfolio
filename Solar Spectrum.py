from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.units as u
from astropy import units as u
import numpy as np
import matplotlib as mpl
from glob import glob
from scipy.optimize import curve_fit
from pathlib import Path

# Set the default text font size
plt.rc('font', size=16)
# Set the axes title font size
plt.rc('axes', titlesize=20)
# Set the axes labels font size
plt.rc('axes', labelsize=20)
# Set the font size for x tick labels
plt.rc('xtick', labelsize=15)
# Set the font size for y tick labels
plt.rc('ytick', labelsize=15)
# Set the legend font size
plt.rc('legend', fontsize=18)
# Set the font size of the figure title
plt.rc('figure', titlesize=30)


def set_size(width=700, fraction=1, subplots=(1, 1)):
    """Set figure dimensions to avoid scaling in LaTeX.

    Parameters
    ----------
    width: float or string
            Document width in points, or string of predined document type
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy
    subplots: array-like, optional
            The number of rows and columns of subplots.
    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """

    # Width of figure (in pts)
    fig_width_pt = width * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    # https://disq.us/p/2940ij3
    golden_ratio = (5**.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])
    if subplots == (1,1):
        return(12,12)
    else:
        return (20,30)
    # return (fig_width_in, fig_height_in)

def makeFiles():
    global LoRes, HiRes, wl, array, Micrometer, tit, lab, tik, foundele

    tit = 30  # title sizes
    lab = 20  # axis label sizes
    tik = 13  # tick sizes

    Micrometer = [4.0, 4.4, 4.8, 5.2, 5.6, 6.0, 6.4, 6.8]
    
    sunSearch = """
    393.3682\t2.0253\tCa II
    440.4761\t0.0898\tFe I
    394.4016\t0.0488\tAl I	
    441.5135\t0.0417\tFe I
    396.1535\t0.0621\tAl I
    452.8627\t0.0275\tFe I
    396.8492\t1.5467\tCa II
    455.4036\t0.0159\tBa II
    404.5825\t0.1174\tFe I
    470.3003\t0.0326\tMg I
    406.3605\t0.0787\tFe I
    486.1342\t0.3680\tH
    407.1749\t0.0723\tFe I
    489.1502\t0.0312\tFe I
    407.7724\t0.0428\tSr II		
    492.0514\t0.0471\tFe I
    410.1748\t0.3133\tH		    
    495.7613\t0.0696\tFe I
    413.2067\t0.0404\tFe I		
    516.7327\t0.0935\tMg I
    414.3878\t0.0466\tFe I		
    517.2698\t0.1259\tMg I
    416.7277\t0.0200\tMg I		
    518.3619\t0.1584\tMg I
    420.2040\t0.0326\tFe I		
    525.0216\t0.0062\tFe I
    422.6740\t0.1476\tCa I		
    526.9550\t0.0478\tFe I
    423.5949\t0.0385\tFe I		
    532.8051\t0.0375\tFe I
    425.0130\t0.0342\tFe I		
    552.8418\t0.0293\tMg I
    425.0797\t0.0400\tFe I		
    588.9973\t0.0752\tNa I (D2)
    425.4346\t0.0393\tCr I		
    589.5940\t0.0564\tNa I (D1)
    426.0486\t0.0595\tFe I		
    610.2727\t0.0135\tCa I
    427.1774\t0.0756\tFe I		
    612.2226\t0.0222\tCa I
    432.5775\t0.0793\tFe I		
    616.2180\t0.0222\tCa O
    434.0475\t0.2855\tH		    
    630.2499\t0.0083\tFe I
    438.3557\t0.1008\tFe I		
    656.2808\t0.1020\tH 
    """
    array = [line.strip().split('\t')[:3] for line in sunSearch.strip().split('\n')] # gets the first column
    array = sorted(array, key=lambda x: float(x[0])) # orders highest to lowest

    LoRes = {"Sun" : "Fits/Team_SCN_Sun_5.6mm_Low.00000101.FIT",
             "Hg"  : "Fits/Team_SCN_HG_5.6mm_Low.00000100.FIT"}

    HiRes = {"Sun" : {4.0 : 'Fits/Team_SCN_Sun_4.0mm_High.00000103.FIT',
                      4.4 : 'Fits/Team_SCN_Sun_4.4mm_High.00000105.FIT',
                      4.8 : 'Fits/Team_SCN_Sun_4.8mm_High.00000107.FIT',
                      5.2 : 'Fits/Team_SCN_Sun_5.2mm_High.00000109.FIT',
                      5.6 : 'Fits/Team_SCN_Sun_5.6mm_High.00000113.FIT',
                      6.0 : 'Fits/Team_SCN_Sun_6.0mm_High.00000116.FIT',
                      6.4 : 'Fits/Team_SCN_Sun_6.4mm_High.00000120.FIT',
                      6.8 : 'Fits/Team_SCN_Sun_6.8mm_High.00000122.FIT',
                      },
             "Hg"  : {4.0 : 'Fits/Team_SCN_HG_4.0mm_High.00000102.NOAUTODARK.FIT',
                      4.4 : 'Fits/Team_SCN_HG_4.4mm_High.00000104.FIT',
                      4.8 : 'Fits/Team_SCN_HG_4.8mm_High.00000106.FIT',
                      5.2 : 'Fits/Team_SCN_HG_5.2mm_High.00000108.FIT',
                      5.6 : 'Fits/Team_SCN_HG_5.6mm_High.00000112.FIT',
                      6.0 : 'Fits/Team_SCN_HG_6.0mm_High.00000115.FIT',
                      6.4 : 'Fits/Team_SCN_HG_6.4mm_High.00000117.FIT',
                      6.8 : 'Fits/Team_SCN_HG_6.8mm_High.00000121.FIT',}}

    wl = {
        # wavelengths identified with Hg lamp
        # pixel : wavelength
        # 360 - 440
        1 : {261 : 404.656,
            290 : 407.783,
            552 : 435.833},
        # 400 - 480
        2 : {81  : 404.656,
            111 : 407.783,
            373 : 435.833,},
        # 440-520
        3 : {12 : 435.833,},
        # 480 - 560
        4 : {675 : 546.074,},
        # 520 - 640
        5 : {307 : 546.074,
            598 : 576.960,
            618 : 579.066,},
        # 560 - 640
        6: {229 : 576.960,
            249 : 579.066,},
        # 600-680
        7 : {108 : 603.213,
             119 : 604.322,
             440 : 638.471,
             470 : 641.631,},
             #717 : 671.643},
        # 640 - 720
        8 : {615 : 696.543,
            712 : 706.722,}
    }

    foundele = {
        # found elements in the sun from array
        # j : [ks]
        # plot number : identified spectral line ordinal
        1 : [1],
        2 : [1],
        3 : [1],
        4 : [1],
        5 : [1],
        6 : [1],
        7 : [1],
        8 : [1],

    }


def spectrum(img_data, row_min, row_max):
    # Flip the spectrum so that wavelength increases going from left to right
    return np.flip(np.average(img_data[row_min:row_max, :], axis=0))


def getSpectrum(source, img=False, tenrows=False, ss=False, bb=False):
    # a and b calculations for low res
    p1, l1 = 306, 546.074
    p2, l2 = 383, 579.066
    a = (l1 - l2) / (p1 - p2)
    b = l1 - p1 * (l1 - l2) / (p1 - p2)
    if source == 'Hg': x = 1
    else: x = 2
    sigpix, sigl = getErr(a, x)
    print('Low Res: a = {}  b = {}'.format(a,b))
    print('Pix err = {}  Lambda err = {}'.format(sigpix, sigl))

    hdu = fits.open(LoRes[source])

    image_hdu = hdu[0]
    spec_data = image_hdu.data
    if img:  # shows and saves the full image
        plt.figure(figsize=set_size(fraction=1))
        print('Image Shape:', spec_data.shape)
        plt.imshow(spec_data)
        plt.title('Hg Low Res Image')
        plt.xlabel('Pixels')
        plt.ylabel('Pixels')
        plt.plot()
        plt.savefig('Low Res {} Full.pdf'.format(source))
        plt.show()


    img_cut = spec_data[250:260]
    if tenrows:  # shows and saves the 10 pixel slice of the image
        plt.figure(figsize=set_size(fraction=1))
        plt.imshow(img_cut)
        plt.axis('off')
        plt.savefig('Low Res {} Strip.pdf'.format(source))
        plt.show()


    if source == 'Hg' or ss or bb:
        fig, ax = plt.subplots(1, 1, figsize=set_size())
        spectrum = np.average(img_cut, axis=0)  # gets the image data
        # applies correct ticks and labels
        xlabels, xticks = getLabels(a, b, len(spectrum))
        ax.set_xticks(xticks)
        ax.set_xticklabels(xlabels)
        plt.xlabel('Wavelength (nm)')
        # plt.xlabel('Pixel')
        plt.ylabel('Intensity (counts)')

    if source == 'Hg': # shows and save the low reolution Hg spectrum
        spec = np.flip(np.average(spec_data[250:260, :], axis=0))
        for i in range(1, 9):
            for key in wl[i]:
                lam = wl[i][key]
                pix = wltopix(lam, a, b)
                ax.vlines(pix, np.min(spec), np.max(spec), color='green', zorder=20, label=lam)
                ax.fill_betweenx([np.min(spec), np.max(spec)], pix - sigpix, pix + sigpix, color='gray', alpha=0.3, zorder=0)
        ax.axvline(x=306, color='red', zorder=20)  # 2 Identified lines
        ax.axvline(x=383, color='red', zorder=20)
        ax.plot(spec, zorder=100)
        plt.title('Mercury Low Res Spectrum')
        plt.savefig('Low Res Hg.pdf')
        plt.show()


    if source == 'Sun':
        if ss: # plots the spectrum with lines and errors
            for i in range(len(array)):
                lam = float(array[i][0])
                pix = wltopix(lam, a, b)
                ax.vlines(pix, np.min(spectrum), np.max(spectrum), color='red', zorder=20)
                # ax.vlines(x=(lam - b) / a, color='green', zorder=20)
                ax.fill_betweenx([np.min(lam), np.max(lam)], pix - sigpix, pix + sigpix, color='gray', alpha=0.3,
                                 zorder=0)

            plt.ylim([5000,18500])
            plt.plot(spectrum)
            plt.title('Solar Low Res Spectrum')
            plt.savefig('Low res Sun.pdf')
            plt.show()

        
        if bb:  # plots black body
            xdata = a * np.arange(0,len(spectrum)) + b
            ydata = spectrum
            plt.plot(xdata, ydata)
            popt, pcov = curve_fit(blackbody, xdata, ydata, p0=[1e10, 5800])
            plt.plot(xdata, blackbody(xdata, *popt), 'r-',
                     label='Fit: A = {:.2e} +- {:.2e}   T = {:.2e} +- {:.2e}'.format(popt[0], pcov[0,0], popt[1], pcov[1,1]))
            plt.legend()
            plt.title("Sun's Black Body Spectrum")
            plt.xlabel("Wavelength (nm)")
            plt.ylabel("intensity (counts)")
            plt.savefig('Blackbody.pdf')
            plt.show()

            print('blackbody errors', pcov[0, 0], pcov[1, 1])


def plotSpectra(numbers, test7=False, exact=False, check=False, model=False):
    Hg_spectra = []
    Sun_spectra = []
    for microm in Micrometer:
        Hg_hdu = fits.open(HiRes['Hg'][microm])
        Sun_hdu = fits.open(HiRes['Sun'][microm])
        Hg_spec_data = Hg_hdu[0].data
        Sun_spec_data = Sun_hdu[0].data
        Hg_spectra.append(spectrum(Hg_spec_data, 250, 260))
        Sun_spectra.append(spectrum(Sun_spec_data, 250, 260))
        # at this point we have the 1D strip for each

    ########################  START HG  ########################
    fig, axs = plt.subplots(len(numbers), 1, figsize=set_size(subplots=(len(numbers), 1)))
    Hg_spectra = [Hg_spectra[i - 1] for i in numbers]
    if len(numbers) == 1: axs = [axs]
    for d, ax, j in zip(Hg_spectra, axs, numbers):
        print(j)
        getAB(j)
        a = wl[j]['a']
        b = wl[j]['b']

        xlabels, xticks = getLabels(a, b, len(d))
        ax.set_xticks(xticks)
        ax.set_xticklabels(xlabels)  # plots ticks

        ax.plot(d, zorder=100)  # plots data
        ax.set_yscale('log')  # sets to log scale in y

        if test7 and j == 7:  # used for calibration using argon lines
            AgWl = np.array([6032.13, 6043.22, 6384.71, 6416.31, 6604.28, 6677.28, 6716.43, 6752.83])
            AgWl = (AgWl - 6000) / 1.07
            ax.vlines(AgWl, np.min(d), np.max(d), color='red')
        else:
            for pix in wl[j]: # plots the lines at the peaks
                if pix not in ('a', 'b'): # key is the pixel number
                    sigpix, sigl = getErr(a,1)
                    st = list(wl[j].keys())[0]
                    nd = list(wl[j].keys())[-3] # these are the lines used for calibration
                    lam = wl[j][pix]
                    if pix in (st, nd):
                        ax.vlines(pix, np.min(d), np.max(d), color='red', zorder=20)
                        ax.fill_betweenx([np.min(d), np.max(d)], pix - sigpix, pix + sigpix, color='gray', alpha=0.3)
                    else:
                        if exact: # plot the lines exactly where I found them
                            c = 'red'
                        if check: # plots the middle lines according to a and b
                            pix = wltopix(lam, a, b)
                            c = 'blue'
                        ax.vlines(pix, np.min(d), np.max(d), color=c, zorder=10)
                        ax.fill_betweenx([np.min(d), np.max(d)], pix - sigpix, pix + sigpix, color='gray', alpha=0.3)

                    print('Pixel: {} +- {}'.format(pix, sigpix))
                    print('Wavle: {} +- {}'.format(lam, sigl))

    plt.suptitle("Mercury High Res Spectra")
    plt.subplots_adjust(hspace=0.5)
    plt.xlabel('Wavelength (nm)')
    fig.text(0.07, 0.5, 'Intensity (counts)', va='center', rotation='vertical')

    # builds name for the pdf
    name = str(numbers)
    if test7: name = name + ' test7'
    if exact: name = name + ' exact'
    if check: name = name + ' check'
    plt.savefig('Hg High Res {}.pdf'.format(name))

    plt.show()
    #########################  END HG  #########################

    ###################################################################################################################

    ########################  START SUN  ########################
    fig, axs = plt.subplots(len(numbers), 1, figsize=set_size(subplots=(len(numbers), 1)))
    Sun_spectra = [Sun_spectra[i-1] for i in numbers]
    if len(numbers) == 1: axs = [axs]
    for d, ax, j in zip(Sun_spectra, axs, numbers):
        a = wl[j]['a']
        b = wl[j]['b']
        sigpix, sigl = getErr(a,2)
        print(j,a ,b)
        print('errors', sigpix,sigl)

        xlabels, xticks = getLabels(a, b, len(d))
        ax.set_xticks(xticks)
        ax.set_xticklabels(xlabels)  # sets ticks
        ax.set_yscale('log')
        ax.plot(d, zorder=100)  # Plot the lines

        # wavelength ranges
        minn = b
        maxx = a * len(d) + b

        k = 1
        for i in range(len(array)):  # plots the lines if they are in the right range
            testwl = float(array[i][0])
            c = 'red'
            status = 'X'
            
            if minn < testwl < maxx:  # if the element is in the right part of the spectrum
                pix = (testwl - b) / a  # get the pixel it should be at

                round_int = lambda x: int(np.round(x))
                checkpix = round_int(pix)
                low = round_int(pix - sigpix)
                high = round_int(pix + sigpix)

                # print(low, high)
                p=low
                for u in range(low,high+1):
                    current = d[p]
                    if d[i] < current:
                        p = u
                obs = pixtowl(low, a, b)
                drangeleft = [d[i] for i in range(low, checkpix+1)]
                drangeright = [d[i] for i in range(checkpix, high+1)]
                # ratio = d[checkpix] / max(drange)
                ratioleft = d[checkpix] / max(drangeleft)
                ratioright = d[checkpix] / max(drangeright)
                # print(ratioright, ratioleft)
                if ratioleft <= 0.9 and ratioright <= 0.9 and low > 0:  # if it's the lowest within 5 pixels on each side
                    c = 'green'
                    status = '0'

                ax.vlines(pix, np.min(d), np.max(d), color=c, zorder=10)
                ax.fill_betweenx([np.min(d), np.max(d)], pix - sigpix, pix + sigpix, color='gray', alpha=0.3)

                # print('{:.2f}\t{:.2f}\t{}\t{}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}'.format(j, k, array[i][-1], status, obs, testwl, sigl, abs(testwl-obs), abs((obs - testwl) / testwl * 100)))
                # print('Pixel: {} +- {}'.format(pix, sigpix))
                # print('Wavle: {} +- {}'.format(testwl, sigl))
                k += 1
        print('---------------------------')

    # title labels
    plt.suptitle("Sun High Res Spectra")
    # plt.subplots_adjust(hspace=0.5)
    plt.xlabel('Wavelength (nm)')
    fig.text(0.07, 0.5, 'Intensity (counts)', va='center', rotation='vertical')

    # builds name for the pdf
    name = str(numbers)
    plt.savefig('Sun High Res {}.pdf'.format(name))

    plt.show()

    # Comparison Part
    if model:
        plots = [1, 6]
        for j in plots:
            a = wl[j]['a']
            b = wl[j]['b']
            d = Sun_spectra[j-1]  # Sun_spectra counts from 0-7, j counts from 1-8

            minn = b
            maxx = a * len(d) + b

            directory = Path('Compare')
            # Loop over files in the directory
            for file in directory.iterdir():
                if file.is_file():
                    plt.figure(figsize=set_size(fraction=1))
                    xdata, ydata = compare(file.name, minn, maxx) # wavelength, flux

                    xlabels, xticks = getLabels(a, b, len(d))
                    plt.xticks(xticks, xlabels)

                    plt.title('{:.0f} - {:.0f} nm  {}'.format(minn, maxx, file.name[:-5]))
                    plt.xlabel('Wavelength (nm)')
                    plt.ylabel('Normalized Intensity')

                    if j == 1: shift = 1
                    else: shift = 0.5
                    y_data = ydata / np.median(ydata) + shift  # normalize *2 plot and raise by a little
                    x_data = wltopix(xdata, a, b)  # convert wavelength to pixel
                    d = d / np.median(d)  # normalize sun


                    plt.plot(np.arange(0, len(d)), d, label='Sun', zorder=100)  # plot sun
                    plt.plot(x_data, y_data, label='Model', zorder=100)  # plot compare

                    plt.savefig('{} {}.pdf'.format(j, file.name[:2]))

                    plt.show()
    ########################  END SUN  ########################

def getAB(j):
    keys = list(wl[j].keys())  # the pixel number
    if len(keys) != 1:
        pixels = [keys[i] for i in [0, -1]]  # gets the two furthest points apart
        lam = [wl[j][p] for p in pixels]  # gets the associated wavelength
        tempa = (lam[0] - lam[1]) / (pixels[0] - pixels[1])
        tempb = lam[0] - pixels[0] * (lam[0] - lam[1]) / (pixels[0] - pixels[1])
        wl[j]['a'], wl[j]['b'] = tempa, tempb  # adds to dictionary
    else:
        pixel = int(keys[0]) # gets the 1 pixel
        lam = wl[j][pixel] # gets the associated wavelength
        tempa = wl[j-1]['a']
        tempb = lam - tempa * pixel
        wl[j]['a'], wl[j]['b'] = tempa, tempb  # adds to dictionary


def getLabels(a, b, length, steps=85):
    # converts pixels to wavelegths for x-axis labels
    xlabels = a * np.arange(length, step=steps) + b  # gets wavelengths to label x axis with
    xlabels = np.append(xlabels, a * length + b)
    xlabels = np.round(xlabels, decimals=0) # rounds to 3 decimals
    xticks = np.arange(length, step=steps)  # places to put labels
    xticks = np.append(xticks, length)
    return xlabels, xticks


def blackbody(lam, A, T):
    h = 6.62607015e-34 # J s
    kb = 1.380649e-23 # m^2 kg s^-2 K^-1
    c = 3e8 # m/s
    f = c / (lam / 1e9)
    return 2 * np.pi * h * f**3 / c**2 * A / (np.exp(h*f/(kb*T)) - 1)


def compare(file, minn, maxx):
    # gets the wavelengths and fluxes between the min and max nm ranges from the given file
    hdu = fits.open('Compare/{}'.format(file))
    data = hdu[1].data
    xdata, ydata = [], []
    for i in range(len(data)):
        val = 10**data[i][0] / 10  # convert from log scale and from A to nm
        if minn <= val <= maxx:
            xdata.append(val)
            ydata.append(data[i][1])
    xdata = np.array(xdata)
    ydata = np.array(ydata)
    return xdata, ydata


def pixtowl(p ,a, b):
    return a * p + b


def wltopix(l, a, b):
    return (l - b) / a


def getErr(a,x):
    delpix = 0.5 * a  # nm from pixel uncertainty
    dell = 0.5        # nm from micrometer uncertainty
    lerr = np.sqrt(x * delpix**2 + dell**2 )  # error in wavelength
    pixerr = wltopix(dell, a, 0)  # error in pixel
    return pixerr, lerr


if __name__ == "__main__":
    makeFiles() # loads the files into a dictionary
    
    source = 'Sun' # Hg or sun, just Hg will show its spectrum
    img = False  # shows the collected image
    tenrows = False # shows the reduced strip
    ss = False # shows the sun spectrum
    bb = False  # shows the black body, bb or ss must be true for something to happen
    
    getSpectrum(source, img, tenrows, ss, bb)


    # numbers = [1,2,3,4,5,6,7,8]  # the spectra to plot; cannot send only 1 at a time
    # numbers = [7]
    numbers = [1,2,3,4]
    # numbers = [5,6,7,8]

    # numbers = [1,2]
    # numbers = [3, 4]
    # numbers = [5, 6]
    # numbers = [7, 8]

    #  these are for Hg only; only 1 should be true at a time
    test7 = False  # plots lines on the 7th spectrum for testing; only use True if number = [7]
    exact = False   # plots all lines exactly where they should be
    check = True  # plots all middle lines using a and b; one of these must be true
    ###############

    # these are for Sun only; can run one or both
    model = False     # shows all the model comparison plots
    ###############

    plotSpectra(numbers, test7, exact, check, model)
