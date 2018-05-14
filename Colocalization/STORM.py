from __future__ import print_function
import numpy as _np
from scipy import optimize as _optimize
from scipy.signal import fftconvolve as _fftconvolve
import random as _random
from scipy.spatial import distance as _distance
import copy as _copy
import struct as _struct
import ImageProcessing as _ip
import numba
import threading as _threading
from sys import stdout as _stdout


mlist_dtype = [('c', 'i4'),
               ('x', 'f4'),
               ('y', 'f4'),
               ('xc', 'f4'),
               ('yc', 'f4'),
               ('h', 'f4'),
               ('A', 'f4'),
               ('w', 'f4'),
               ('phi', 'f4'),
               ('ax', 'f4'),
               ('b', 'f4'),
               ('I', 'f4'),
               ('f', 'i4'),
               ('len', 'i4'),
               ('link', 'i4'),
               ('valid', 'i4'),
               ('z', 'f4'),
               ('zc', 'f4')]

mlist_column_number = {'c': 0, 'xc': 3, 'yc': 4}


def point_set_distance_distribution(mlist_set, channel_set, mlist_points,
                                    channel_points, movie_size, pixel_size,
                                    bin_size, r_max, nonspec=0, verbose=False):
    movie_size = float(movie_size)
    pixel_size = float(pixel_size)
    bin_size = float(bin_size)
    r_max = float(r_max)
    nonspec = float(nonspec)

    zoom = pixel_size/bin_size

    if verbose:
        print('Generating pixel-based super-resolution images...')
    img_set = generate_density_map(mlist_set, channel_set, movie_size, zoom)
    img_points = generate_density_map(mlist_points, channel_points, movie_size,
                                      zoom)
    if nonspec:
        img_0 = generate_density_map(mlist_points, 0, movie_size, zoom)
        img_points -= nonspec*img_0
    img_set = _ip.pad(img_set)
    img_points = _ip.pad(img_points)

    rho = bin_size*_ip.radial_coordinate(img_set.shape)
    n_bins = int(_np.ceil(r_max/bin_size))
    r = bin_size*_np.arange(n_bins+1)
    
    cum_ps_dist_hist = _np.zeros(n_bins+1)
    cum_area = _np.zeros(n_bins+1)
    
    for i in xrange(1, n_bins+1):
        if verbose:
            _stdout.write('\rCalculating point-distribution for distance {}/{}.'.format(i, n_bins+1))
        disk = rho <= r[i]
        img_area = _ip.conv(img_set, disk, padding=False) > 0.5
        cum_ps_dist_hist[i] = _np.sum(img_points*img_area)
        cum_area[i] = _np.sum(img_area)
    
    if verbose:
        _stdout.write('\n')
   
    XY = (movie_size*pixel_size)**2
    N_points = _np.sum(img_points)
    area = bin_size*bin_size*_np.diff(cum_area)
    
    psdd = XY*_np.diff(cum_ps_dist_hist)/(area*N_points)
        
    return r[:-1], psdd, area


def pair_distance_distribution(mlist_A, channel_A, mlist_B, channel_B,
                               movie_size, pixel_size, bin_size, r_max,
                               nonspec_A=0, nonspec_B=0, same=False, n_threads=1,
                               verbose=False):
    movie_size = float(movie_size)
    pixel_size = float(pixel_size)
    bin_size = float(bin_size)
    r_max = float(r_max)
    mols_A = mlist_A[mlist_A['c']==channel_A]
    mols_B = mlist_B[mlist_B['c']==channel_B]
    NA = len(mols_A)
    NB = len(mols_B)
    bin_size_pxl = bin_size/pixel_size
    r_max_pxl = r_max/pixel_size
    if verbose:
        print("Correlating {}x{} molecules.".format(NA,NB))
        print("Calculating distance histogram...")
    dists_A_B = distance_histogram(mols_A, mols_B, bin_size_pxl, r_max_pxl, same,
                                   n_threads, verbose)
    if nonspec_A:
        if verbose:
            print("Calculating image A to non-specific distance histogram...")
        mols_A0 = mlist_A[mlist_A['c']==0]
        dists_A0_B = distance_histogram(mols_A0, mols_B, bin_size_pxl, r_max_pxl,
                                        False, n_threads, verbose)
        dists_A_B -= nonspec_A*dists_A0_B
        NA -= nonspec_A*len(mols_A0)
    if nonspec_B:
        if verbose:
            print("Calculating image B to non-specific distance histogram...")
        mols_B0 = mlist_B[mlist_B['c']==0]
        dists_A_B0 = distance_histogram(mols_A, mols_B0, bin_size_pxl, r_max_pxl,
                                        False, n_threads, verbose)
        dists_A_B -= nonspec_B*dists_A_B0
        NB -= nonspec_B*len(mols_B0)
    if nonspec_A and nonspec_B:
        if verbose:
            print("Calculating non-specific auto-distance histogram...")
        dists_A0_B0 = distance_histogram(mols_A0, mols_B0, bin_size_pxl, r_max_pxl,
                                         True, n_threads, verbose)
        dists_A_B += nonspec_A*nonspec_B*dists_A0_B0
    XY = (movie_size*pixel_size)**2
    n_bins = _np.ceil(r_max/bin_size)
    r = bin_size*_np.arange(n_bins)+bin_size
    A = _np.pi*bin_size*(2*r+bin_size)
    pdd = XY*dists_A_B/(NA*NB*A)
    return r-bin_size, pdd


def distance_histogram(mols_A, mols_B, bin_size, r_max, same=False, n_threads=1,
                       verbose=False):
    n_bins = int(_np.ceil(r_max/bin_size))
    if n_threads>1:
        NA = len(mols_A)
        mols_per_thread = NA/n_threads
        start = mols_per_thread*_np.arange(n_threads)
        end = start + mols_per_thread
        end[-1] = NA
        results = [_np.zeros(n_bins) for _ in range(n_threads)]
        numba_func = numba.jit(_numba_distance_histogram, nopython=True,
                               nogil=True)
        args = [(mols_A['xc'][s:e],
                 mols_A['yc'][s:e],
                 mols_B['xc'],
                 mols_B['yc'],
                 bin_size,
                 n_bins,
                 same,
                 False,
                 result) for s,e,result in zip(start,end,results)]
        threads = [_threading.Thread(target=numba_func, args=_) for _ in args]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()

        return sum(results)
    else:
        numba_func = numba.jit(_numba_distance_histogram, nopython=True)
        return numba_func(mols_A['xc'], mols_A['yc'],
                          mols_B['xc'], mols_B['yc'],
                          bin_size, n_bins, same, verbose,
                          _np.zeros(n_bins))


def _numba_distance_histogram(x1, y1, x2, y2, bin_size, n_bins, same, verbose, dist_hist):
    i_percent = int(len(x1)/100)
    if i_percent==0:
        i_percent = 1
    for i in range(len(x1)):
        if verbose and i%i_percent==0:
            print(i/i_percent)
        for j in range(len(x2)):
            if (not same) or (i!=j):
                d = _np.sqrt((x1[i]-x2[j])**2 + (y1[i]-y2[j])**2)
                bin = int(d/bin_size)
                if bin<n_bins:
                    dist_hist[bin] += 1
    return dist_hist


def pair_distance_distribution_fft(mlistA, channelA, mlistB, channelB, movie_size, pixel_size, bin_size, nonspecA=0, nonspecB=0, verbose=False):
    movie_size = float(movie_size)
    pixel_size = float(pixel_size)
    bin_size = float(bin_size)
    zoom = pixel_size/bin_size
    mapA = generate_density_map(mlistA, channelA, movie_size, zoom)
    mapB = generate_density_map(mlistB, channelB, movie_size, zoom)
    if nonspecA:
        if verbose:
            print('Removing non-specific blinking from image A.')
        mapA0 = generate_density_map(mlistA, 0, movie_size, zoom)
        mapA -= nonspecA*mapA0
    if nonspecB:
        if verbose:
            print('Removing non-specific blinking from image B.')
        mapB0 = generate_density_map(mlistB, 0, movie_size, zoom)
        mapB -= nonspecB*mapB0
    if verbose:
        print('Calculating image correlation.')
    xcorr = _ip.xcorr(mapA, mapB)
    XY = (movie_size*pixel_size)**2
    A_sr_pxl = bin_size**2
    NA = _np.sum(mapA)
    NB = _np.sum(mapB)
    pdd = XY*_ip.radial_average(xcorr)/(NA*NB*A_sr_pxl)
    r = bin_size*_np.arange(len(pdd))
    return r, pdd


def generate_density_map(mlist, channel, movie_size, zoom):
	molecules = mlist['c']==channel
	xc_sr_pixel = _np.floor(zoom*(mlist['xc'][molecules]-0.5))
	yc_sr_pixel = _np.floor(zoom*(mlist['yc'][molecules]-0.5))
	n_pixels = _np.ceil(movie_size*zoom)
	valid = (xc_sr_pixel>=0) & (xc_sr_pixel<n_pixels) & (yc_sr_pixel>=0) & (yc_sr_pixel<n_pixels)
	coords_sr_pixel = zip(yc_sr_pixel[valid], xc_sr_pixel[valid])
	density_map = _np.zeros((n_pixels, n_pixels))
	for coord in coords_sr_pixel:
		density_map[coord] += 1
	return density_map

def default_mlist(n_molecules):
	# Integer ones
	oi = _np.ones(n_molecules, dtype='i4')
	# Float zeros
	zf = _np.zeros(n_molecules, dtype='f4')
	# Float ones
	of = _np.ones(n_molecules, dtype='f4')
	mlist = _np.array(zip(oi,zf,zf,zf,zf,100*of,zf,2*of,zf,of,100*of,1000*of,oi,oi,-oi,oi,zf,zf), dtype=mlist_dtype)
	return mlist

def load_mlist(filename, columns=None):
	usecols = None
	if columns:
		usecols = [mlist_column_number[_] for _ in columns]
	return _np.genfromtxt(filename, dtype=mlist_dtype, skip_header=1, usecols=usecols)

def write_mlist(filename, mlist):
	with open(filename, 'w') as output:
		output.write('Cas{}\tX\tY\tXc\tYc\tHeight\tArea\tWidth\tPhi\tAx\tBG\tI\tFrame\tLength\tLink\tValid\tZ\tZc\n'.format(len(mlist)))
		for m in mlist:
			output.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(m['c'],
																										   m['x'],
																										   m['y'],
																										   m['xc'],
																										   m['yc'],
																										   m['h'],
																										   m['A'],
																										   m['w'],
																										   m['phi'],
																										   m['ax'],
																										   m['b'],
																										   m['I'],
																										   m['f'],
																										   m['len'],
																										   m['link'],
																										   m['valid'],
																										   m['z'],
																										   m['zc']))





def create_dummy_spe(filename, frames, shape):

	with open(filename, 'wb') as handle:
		handle.write(_struct.pack(4100*'x'))

		handle.seek(108)
		handle.write(_struct.pack('i', 3))

		handle.seek(42)
		handle.write(_struct.pack('I', shape[0]))

		handle.seek(656)
		handle.write(_struct.pack('I', shape[1]))

		handle.seek(1446)
		handle.write(_struct.pack('q', frames))

		handle.seek(4100)
		for frame in xrange(frames):
			handle.write(_np.zeros(shape, dtype=_np.uint16).tostring())



def localization_precision(s, N, b, a, EM=True, FWHM=True):
    '''
    Reference: Mortensen et al., Nat. Meth., 2010.

    Parameters
    ----------
    s: standard deviation of PSF
    N: number of photons
    b: number of background photons per pixel
    a: pixel size
    EM: specify whether camera has EM gain.
    FWHM: return localization precision in FWHM, otherwise in standard deviation.
    '''
    s2 = _np.float32(s)**2
    N = _np.float32(N)
    b2 = _np.float32(b)   # We don't need power of 2, see Ref.
    a2 = _np.float32(a)**2
    sa2 = s2 + a2/12
    var = sa2 * ( 16/9. + 8*_np.pi*sa2*b2 / ( N*a2 ) ) / N
    if EM:
        var *= 2
    lp = _np.sqrt(var)
    if FWHM:
        lp *= 2.355
    return lp

def simulateLines(n_points_per_line, length, distance, center, angle):

    x = _np.linspace(-length/2., length/2., n_points_per_line) + center[0]
    y1 = center[1] + distance/2.
    y2 = center[1] - distance/2.

    positions = []

    for i in xrange(n_points_per_line):
        positions.append((x[i], y1, 0))
        positions.append((x[i], y2, 0))

    return positions



def simulateCircle(n_points, diameter=1, center=(0,0), rotation=0, label_probability=1):

    radius = diameter/2.
    angles = _np.linspace(0, 2*_np.pi, n_points+1) + rotation

    detected = _np.random.random_sample(n_points)
    indices = _np.where(detected < label_probability)

    x = radius * _np.cos(angles[indices]) + center[0]
    y = radius * _np.sin(angles[indices]) + center[1]

    return x, y


def simulateSTORMImage(positions, avg_n_switching_cycles,
        localization_precision, min_n_switching_cycles=0):

    x0, y0 = positions
    n_positions = len(x0)
    n_switching_cycles = _np.round(_np.random.exponential(
        avg_n_switching_cycles, n_positions))
    n_switching_cycles[n_switching_cycles<min_n_switching_cycles] = \
            min_n_switching_cycles
    x = _np.array([])
    y = _np.array([])
    for i in range(n_positions):
        x = _np.append(x, _np.random.normal(x0[i], localization_precision,
            n_switching_cycles[i]))
        y = _np.append(y, _np.random.normal(y0[i], localization_precision,
            n_switching_cycles[i]))
    return (x,y)


def radial_distribution_u(listA, listB, bin_size, r_max):

    coordsA = _np.vstack((listA.x, listA.y)).T
    coordsB = _np.vstack((listB.x, listB.y)).T
    distances = _distance.cdist(coordsA, coordsB)
    bins = _np.arange(0, r_max, bin_size)
    hist, bins = _np.histogram(distances, bins)
    return bins[:-1], hist


def radial_distribution(listA, listB, bin_size, r_max):

    bins_lower, distribution = radial_distribution_u(listA, listB, bin_size, r_max)
    return bins_lower, distribution/len(listA)


def pair_correlation(listA, listB, bin_size, r_max, X=1, Y=1):

    bins_lower, distribution = radial_distribution(listA, listB, bin_size, r_max)
    area = _np.pi*(2*bins_lower*bin_size + bin_size**2)
    return bins_lower, distribution*X*Y/(area*len(listA)*len(listB))



def STORM_correlation(listA, listB):
    '''
    This is a brute force looping solution.
    Here, listA and listB are (x,y) coordinate tuples, where
    x and y are arrays.
    '''
    xcorr_x = _np.array([])
    xcorr_y = _np.array([])
    for molA in zip(*listA):
        xcorr_x = _np.append(xcorr_x, molA[0] - listB[0])
        xcorr_y = _np.append(xcorr_y, molA[1] - listB[1])
    return (xcorr_x, xcorr_y)


def get_rand_n_photons(size, avg_n_photons, photons_dist_width,
        n_min_photons=0,
        n_max_photons=None):

    n_photons = _np.random.normal(avg_n_photons,
            photons_dist_width*avg_n_photons, size)

    while _np.any(n_min_photons > n_photons > n_max_photons):
        _ = n_min_photons > n_photons > n_max_photons
        n_photons[_] = _np.random.normal(avg_n_photons,
                photons_dist_width*avg_n_photons, _.sum())

    return n_photons



def evaluatePSF(size, position, width):

    nx, ny = size
    x0, y0 = position

    x, y = _np.mgrid[:nx, :ny]

    s = width/2.

    # Gaussian
    #return _np.exp( -((x-x0)**2 + (y-y0)**2) / (2*s**2) ) / (2*_np.pi*s**2)

    # Bo's experimental expression
    return ( _np.exp(-2 * ((x-x0)**2 + (y-y0)**2) / 1.72**2) \
        + 0.0208 * _np.exp(-2 * (_np.sqrt((x-x0)**2 + (y-y0)**2) - 2.45)**2 \
        / 1.10**2) ) / 5.088



def simulateRawMovie(size, mlist, avg_background=70,
        std_read_noise=0.5, em_gain=1, preamp_gain=1, ad_unit=1, baseline=100,
        pixel_size=150):

    photons = avg_background * _np.ones(size)

    for mol in mlist:
        photons[mol.frame] += mol.i * evaluatePSF(size[1:], (mol.x, mol.y),
                mol.w/float(pixel_size))

    read_noise = _np.random.normal(scale=std_read_noise, size=size)

    electrons = em_gain * preamp_gain * _np.random.poisson(photons) + read_noise

    counts = electrons/ad_unit + baseline

    return counts


def gaussian2d(x, y, x0, y0, w, h):
    return h*_np.exp(-( ( (x-x0)**2/(2*w**2) ) + ( (y-y0)**2/(2*w**2) ) ))


def multiple_gaussian2d(x, y, x0, y0, w, h):
    data = _np.zeros(x.shape)
    for i in xrange(len(x0)):
        data += gaussian2d(x, y, x0[i], y0[i], w[i], h[i])
    return data


def fit_gaussians_at_positions(data, positions, w0, h0, b0):

    n = len(positions)
    x, y = _np.meshgrid(range(data.shape[0]), range(data.shape[1]))
    p0 = n*[w0, h0] + [b0]
    print(p0)

    def errorfunction(p):

        x0 = [_[0] for _ in positions]
        y0 = [_[1] for _ in positions]
        w = p[0:-1:2]
        h = p[1::2]
        b = p[-1]
        return _np.ravel(multiple_gaussian2d(x, y, x0, y0, w, h, b) - data)

    return _optimize.leastsq(errorfunction, p0)
