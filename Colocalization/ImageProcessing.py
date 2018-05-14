import numpy as _np
from numpy import fft as _fft


def radial_coordinate(size):
    center_x = size[0]/2.0-0.5
    center_y = size[1]/2.0-0.5
    y, x = _np.indices(size)
    return _np.sqrt((x-center_x)**2 + (y-center_y)**2)


def radial_average(image):
    center_x = image.shape[0]/2-0.5
    center_y = image.shape[1]/2-0.5
    y, x = _np.indices(image.shape)
    r = _np.sqrt((x-center_x)**2 + (y-center_y)**2)
    r = r.astype(_np.int)
    return _np.bincount(r.ravel(), image.ravel())/_np.bincount(r.ravel())


def radial_sum(image):
    center_x = image.shape[0]/2-0.5
    center_y = image.shape[1]/2-0.5
    y, x = _np.indices(image.shape)
    r = _np.sqrt((x-center_x)**2 + (y-center_y)**2)
    r = r.astype(_np.int)
    return _np.bincount(r.ravel(), image.ravel())


def pad(image):
    return _np.pad(image, ((image.shape[0], image.shape[0]),
                           (image.shape[1], image.shape[1])), mode='constant')


def xcorr(image_A, image_B, padding=True):
    if padding:
        image_A = pad(image_A)
        image_B = pad(image_B)
    return _fft.fftshift(_np.real(_fft.ifft2(_fft.fft2(image_A) *
                                             _np.conj(_fft.fft2(image_B)))))

                                             
def conv(image_A, image_B, padding=True):
    if padding:
        image_A = pad(image_A)
        image_B = pad(image_B)
    return _fft.fftshift(_np.real(_fft.ifft2(_fft.fft2(image_A) *
                                             _fft.fft2(image_B))))