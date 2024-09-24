#!/usr/bin/env python3

import numpy as np

from scipy import integrate
from scipy.signal import butter, lfilter
from scipy.fft import fft
from scipy import signal
from scipy.interpolate import griddata


class Record:
    def __init__(self, array, time_vector, f_min, f_max, dt, sim_time, n_han_win=13):
        # IF THE DATA COMES FROM SW4, INPUT IT AS A NUMPY ARRAY. OTHERWISE, PLEASE PROVIDE THE LINK

        self.__sim_time = sim_time
        self.__dt = dt
        self.__time_vector = np.arange(int(sim_time / dt)) * dt
        self.__f_min = f_min
        self.__f_max = f_max
        self.__n_han_win = n_han_win

        self.__time_vector_raw = time_vector
        self.__vel = self.synthetic_record_preparation(array)

        self.__acc = np.diff(self.__vel, prepend=0) / dt
        self.__fas, self.__freq = self.fourier_transform()

    def synthetic_record_preparation(self, array):
        new_signal = griddata(self.__time_vector_raw, array, self.__time_vector, method='linear')
        new_signal[np.isnan(new_signal)] = 0
        return self.butter_bandpass_filter(new_signal)

    def butter_bandpass_filter(self, vector, order=3):
        nyq = 0.5 / self.__dt
        f_min_norm = self.__f_min / nyq
        f_max_norm = self.__f_max / nyq
        #print(f_min_norm, f_max_norm, self.__dt, nyq)
        b, a = butter(order, (f_min_norm, f_max_norm), btype='band')
        return lfilter(b, a, vector)

    def fourier_transform(self):
        han_win = np.hanning(self.__n_han_win) / np.sum(np.hanning(self.__n_han_win))
        taper_rec = signal.windows.tukey(self.vel.size, alpha=0.05)

        fft_x = np.abs(fft(self.__vel * taper_rec, norm='ortho')[:self.__vel.size // 2]) * 2
        freq = np.arange(fft_x.size) * (1 / self.__sim_time)

        return np.convolve(fft_x, han_win, mode='same'), freq

    def sdof_response(self, osc_freq, damping=0.05):
        def transfer_function(freq_vector, f0, damping):
            omega_0 = 2 * np.pi * f0
            omega_vector = freq_vector * (2 * np.pi)

            return -1 * omega_0 ** 2 / (omega_0 ** 2 - omega_vector ** 2 - 2j * omega_0 * omega_vector * damping)

        freq = np.fft.fftfreq(self.__time_vector.size, d=self.__dt)
        tf_func = transfer_function(freq, f0=osc_freq, damping=damping)

        sdof_resp = np.fft.ifft(np.fft.fft(self.__acc * signal.windows.tukey(self.__acc.size, alpha=0.02)) * tf_func)

        return sdof_resp

    @property
    def vel(self):
        return self.__vel

    @property
    def acc(self):
        return self.__acc

    @property
    def time_vector(self):
        return self.__time_vector

    @property
    def dt(self):
        return self.__dt

    @property
    def fas(self):
        return self.__fas

    @property
    def frequency_vector(self):
        return self.__freq

    @property
    def mean_fas(self):
        f_max_index = np.argmin(np.abs(self.__f_max - self.__freq))
        f_min_index = np.argmin(np.abs(self.__f_min - self.__freq))
        return self.__fas[f_min_index:f_max_index].mean()

    @property
    def pga(self):
        return np.max(np.abs(self.__acc))

    @property
    def pgv(self):
        return np.max(np.abs(self.__vel))

    @property
    def pgd(self):
        return np.max(integrate.cumtrapz(self.__vel, dx=self.__dt, initial=0))

    def duration(self, osc_freq, xf, xi, damping):
        #print("In duration.")
        sfod_waveforms = self.sdof_response(osc_freq, damping)
        integral = integrate.cumtrapz(sfod_waveforms ** 2, self.__time_vector)
        norm_ai = integral / integral[-1]

        ti = np.argmin(np.abs(norm_ai - xi / 100))
        tf = np.argmin(np.abs(norm_ai - xf / 100))

        #ret_val = self.__time_vector[tf] - self.__time_vector[ti]
        #print(type(ret_val))
        return self.__time_vector[tf] - self.__time_vector[ti]

    def arias_intensity(self):
        return integrate.trapz(self.__acc ** 2, dx=self.__dt) * np.pi / (2 * 9.8)

    def ai_as_function_of_time(self):
        return integrate.cumtrapz(self.__acc ** 2, dx=self.__dt, initial=0) * np.pi / (2 * 9.8)

    def energy_integral(self):
        return integrate.trapz(self.__vel ** 2, dx=self.__dt)

    def ei_as_function_of_time(self):
        return integrate.cumtrapz(self.__vel ** 2, dx=self.__dt, initial=0)
