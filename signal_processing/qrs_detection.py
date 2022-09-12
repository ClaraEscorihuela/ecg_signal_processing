"""
Algorithm for detecting qrs complexes, and q - r - s peaks in ecg data

Based on:
PAN, Jiapu; TOMPKINS, Willis J.
"A real-time QRS detection algorithm."
IEEE transactions on biomedical engineering, 1985, no 3, p. 230-236.

and,
D. S. Raju, M. S. Manikandan, and R. Barathram.
“An automated methodfor detecting systolic peaks from arterial blood pressure signals”
In Proceedings of the 2014 IEEE Students Technology Symposium, 2014, pp. 41–46
"""
import os
import sys
curr_dir = os.getcwd()+'\signal_processing'
sys.path.append(curr_dir)

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from itertools import groupby
from scipy import signal
import numpy as np


"""
@brief Denoise ecg signal with a buterworth bandpass filter, subtract the first and last tw seconds 
       to avoid noise at the beginning/ end of the signal and remove baseline
@param signal_ecg: signal to denoise
@param fs: sample frequency
@param w1: low cut off frequency
@param w2: high cut off frequency
@param tw: second to eliminate from the beginning/end of the signal
@return ecg: denoised signal
"""
def denoise_ecg(signal_ecg,fs = 1000, w1 = 6, w2=20, tw = 5):
    tw = tw
    # En el paper ponen entre 5 i 15, pero destroza maximo la senyal. Buscar otras frecuencias
    sos = signal.butter(N=10, Wn=[w1,w2], btype='bandpass', fs=fs, output='sos')
    ecg = np.zeros(signal_ecg.size)
    ecg[int(fs*tw/1000):len(signal_ecg)-int(fs*tw/1000)] = signal_ecg[int(fs*tw/1000):-int(fs*tw/1000)]
    ecg[int(fs*tw/1000):len(signal_ecg)-int(fs*tw/1000)] = signal.sosfiltfilt(sos, ecg[int(fs*tw/1000):len(signal_ecg)-int(fs*tw/1000)])
    return ecg

"""
@brief Filter to remove baseline and low frequencies
@param signal: signal to denoise
@param fs: sample frequency
@param fc: low cut off frequency
@return ecg: denoise signal
"""
def remove_basal(signal, fs, fc):
    n = signal.size
    m = int(fs/(2*fc))
    if m % 2 == 0:
        m += 1
    r = int((m-1)/2)
    y = np.zeros(signal.size)
    for i in range((r+1), (n-r)):
        sum = 0
        for s in range((i-r),(i+r)):
            sum += signal[s]

        # Estimated Baseline
        y[i]=sum/m
    ecg = signal - y

    return ecg

"""
@brief Gaussian filter and gaussian derivative filter: this filter aims at smoothing the signal or
       calculating the derivative respectively
@param dat: signal to filter
@param M: kernel's size
@param sigma2: variance (spread) of the kernel
@param der: Boolean to determine derivative (True)/ "normal"(False) gaussian filter
@param plot: Boolean to draw the filter's kernel
@return filtered: convolution of the designed filter with the signal
"""
def guassian_derivative_filter(data,M = 5,sigma2 = 0.4, der = True, plot = False):
    m = np.array(range(M))+1
    g = np.exp(-1*(np.power(m-(M/2),2))/(2*sigma2)) #Guassian Kernel

    if der:
        h = g[1:]-g[:-1]
    else:
        h = g  #Guassian derivative Kernel.

    if plot:
        plt.plot(h)
        plt.show()

    filtered = np.convolve(data, h, mode='same')
    return filtered

"""
@brief Power filter: This filter aims at elevating the derivative filter to the power of 2
@param data: signal to filter
@param N: number to power the signal
@return filtered: absolute value of the signal powered N
"""
def square_filter(data, N=2):
    return abs((data)**N)

"""
@brief Integral filter: This filter aims at calculating the integral of the power signal
@param data: signal to filter
@param N: kernel size
@return filtered: signal filter, position x(t) == sum(x(t),...,x(t-N))
"""
def integration_filter(data, N=5):
    energy = np.zeros(data.shape)
    for i in range(N, len(data)):
        energy[i] = np.sum(data[i-N:i])

    threshold_energy = energy > 0.005
    return energy, threshold_energy

"""
@brief Slope filter: This filter aims at detecting the qrs complex, by selecting the 
                     positive slope of the derivative's integral
@param data: signal to filter
@param min: minimum time to select the slope of the signal must be higher than th_grad 
@param fs: sample frequency
@param th_grad: value that the gradient of the signal must be higher
@return periods: list of list, where each list contains the index of the onset and the offset of the qrs complex
"""
def slope_duration(data, min=0.75, fs=1000, th_grad = 0.1):
    #print('th_grad =' , th_grad)
    pos_grad = np.gradient(data, 1 / fs) > th_grad
    i = 0
    min_lenght = min*fs
    periods = []
    for k, g in groupby(pos_grad):
        g = list(g)
        # Make sure the qrs complex time is higher than a phisiologic factor
        if (k ==True and len(g)>min_lenght):
            periods.append([i, i+len(g)])
        i += len(g)

    if len(periods) <= 1:
        th_grad -= 0.02
        #print('check this patient,because somthing is happening, th_grad =', th_grad)

        if th_grad > 0:
            periods = slope_duration(data, min=0.75, fs=1000, th_grad=th_grad)

        else:
            periods = True

    return periods

"""
@brief find the r peaks in the qrs complex
@param signal: signal to find the r peaks
@param periods: qrs periods to look for the r peak
@param min_refractor: time after a qrs complex, where another qrs complex 
                      is not physiological reasonable to raise
@param fs: sample frequency
@return r_p_ind, r_p, p: indexes for the r peaks, values of the r peaks, qrs periods
"""
def find_r(signal, periods, min_refractor = 0.01, fs=1000):
    t_refractory = min_refractor*fs
    r_peaks_ind = []
    r_peaks = []
    r_p_ind = []
    p = []

    # Make sure the qrs interval is in the physiologic range time
    for i, period in enumerate(periods):
        r_peaks_ind.append(period[0] + np.argmax(signal[period[0]:period[1]]))
        r_peaks.append(np.max(signal[period[0]:period[1]]))

    # Make sure rpeaks are separated with a reasonable phisiological time 200ms
    a = -1
    #for i in range(0,len(r_peaks_ind)-2):
    while a < (len(r_peaks_ind)-2):
        a += 1
        if (abs(periods[a][1] - periods[a+1][0])) >= t_refractory:
            r_p_ind.append(r_peaks_ind[a])
            p.append(periods[a])
        else:
            ind_max = np.argmax([r_peaks[a],r_peaks[a+1]])
            r_p_ind.append(r_peaks_ind[a + ind_max])
            p.append(periods[a + ind_max])
            a += 1


    r_p_ind, unique_ind = np.unique(np.array(r_p_ind), return_index=True)
    r_p = signal[r_p_ind]
    p = np.array(p)[unique_ind]

    # print(f'number of qrs detected before cleaning = {len(r_peaks_ind)}')
    # print(f'number of qrs detected aftere cleaning = {len(r_p_ind)}')

    return r_p_ind, r_p, p


"""
@brief Find q and s peaks 
@param signal: signal to find the q and s peaks
@param r_peaks_ind: index of r peaks
@param N: number of samples to look before/after the r peak to fins the s/q peak
@return q_peaks_ind, q_peaks, s_peaks_ind, s_peaks: indexes and values for the q and s peaks
"""
def find_qs(signal, r_peaks_ind, N =75):
    q_peaks_ind = []
    q_peaks = []
    s_peaks_ind = []
    s_peaks = []

    for r_ind in r_peaks_ind:
        q_peaks_ind.append(r_ind - (len(signal[r_ind-N:r_ind]) - (np.argmin(signal[r_ind-N:r_ind]))))
        q_peaks.append(np.min(signal[r_ind-N:r_ind]))

        s_peaks_ind.append(r_ind + np.argmin(signal[r_ind:r_ind+N]))
        s_peaks.append(np.min(signal[r_ind:r_ind+N]))

    return q_peaks_ind, q_peaks, s_peaks_ind, s_peaks

"""
@brief Correction for the qrs complex, make sure the complex contains the q and s peak, 
                                       and the onset and offset is below the avg of the signal 
@param periods: list of list with the qrs complexes
@param q_peaks: index for the q peaks
@param s_peaks: index for the s peaks
@return periods: list of list where each list contains the onset of the q, and the offset
                 of the s peak for each qrs complex
"""
def final_corrections(signal, periods, q_peaks, s_peaks, N = 30):

    # Make sure the qrs interval arrive to the minimum/maximum local point
    avg = np.mean(signal)
    print('avg = ', avg)
    for i, period in enumerate(periods):
        if q_peaks[i] < period[0]:
            period[0] = q_peaks[i]
        if s_peaks[i] > period[1]:
            period[1] = s_peaks[i]

        # Fix q peak onset
        if i!=0:
            values = signal[period[0]-N: period[0]] < avg
            below_zero = [[k,len(list(g))] for k, g in groupby(values)]
            below_zero = below_zero[-1]
            if below_zero[0] == True:
                period[0] = period[0]-below_zero[1]
            elif below_zero[0] == False:
                period[0] = period[0]+below_zero[1]

        # Fix s peak offset
        if i!=(len(periods)-2):
            values = signal[period[1]: period[1]+N] > avg
            above_zero = [[k,len(list(g))] for k, g in groupby(values)]
            above_zero = above_zero[0]
            if above_zero[0] == True:
                period[1] = period[1]-above_zero[1]
            elif above_zero[0] == False:
                period[1] = period[1]+above_zero[1]

    return periods


"""
@brief Process pipeline to detect qrs complexes: denoise of the ecg, derivative of the denoise signal,
       integration of the square signal, smooth of the integral and detection of its increasing slope
@param signal_ecg_1: signal to find the qrs complexes
@param fs: sample frequency
@param wind_int: kernel size for the integral filter 
@param m_gaus: kernel size for the smoothing filter
@param sigma_gaus: variance for the smoothing filter
@param th_grad: threshold for accepting the gradient of the integral 
@param t_qrs_min: minimum time for a qrs complex
@return periods: list of list, where each list contains the onset and the offset of the q/s peaks
                 for each qrs complex
"""
def processing(signal_ecg_1, fs = 1000, wind_int = 150, m_gaus = 20, sigma_gaus = 6,
               th_grad = 0.1, t_qrs_min = 0.6):

    # Passaband filter and remove baseline
    ecg_filter = denoise_ecg(signal_ecg_1,fs)

    # Remove baseline manually (chose btw one of those)
    # ecg_no_baseline = remove_basal(signal_ecg_2, fs, 5)

    # Gaussian derivative filter
    ecg_derivative = guassian_derivative_filter(ecg_filter)

    # Square threshold filter
    ecg_square = square_filter(ecg_derivative, 2)

    # Moving windows integration (s=150) # To take a windws of 150 ms == paper
    ecg_int, th_energy = integration_filter(ecg_square, wind_int)

    # QRS complex detection
    ecg_int_smooth = guassian_derivative_filter(ecg_int, M=m_gaus, sigma2=sigma_gaus, der=False)
    # pos_grad = np.gradient(ecg_int_smooth, 1 / fs) > th_grad
    periods = slope_duration(ecg_int_smooth, min=t_qrs_min, fs=fs, th_grad = th_grad)

    return periods

"""
@brief Process pipeline to detect q,r,s peaks in each qrs complex
@param signal: signal to find the q,r,s peaks
@param periods: list of list, where with the onset and the offset of the q/s peaks
                 for each qrs complex
@param t_ref: physiological time where a qrs complex cannot occur after another qrs complex 
@param fs: sample frequency
@param w1/w2: sut of low high frequency for denoising the signal

@return periods, r_peaks_ind, q_peaks_ind, s_peaks_ind: list of list, where each list contains the onset and the offset of the q/s peaks
                 for each qrs complex, list with indexes of r peaks, list with indexes q peaks, list with indexes of s peaks
"""
def peak_detection(signal, periods, t_ref = 0.2, fs=1000, w1 =2, w2=50):
    signal_filter = denoise_ecg(signal, fs=fs, w1=w1, w2=w2)

    indxs_pre = []
    # Q,R,S peak ind detection
    r_peaks_ind, _, periods = find_r(signal_filter, periods, min_refractor=t_ref, fs=fs)
    q_peaks_ind, _, s_peaks_ind, _ = find_qs(signal_filter, r_peaks_ind)

    periods_after = final_corrections(signal_filter,periods, q_peaks_ind, s_peaks_ind)

    return periods_after, r_peaks_ind, q_peaks_ind, s_peaks_ind

"""
@brief Draw al the signal for the process to detect the qrs complexs
"""
def plot_process(fs, ecg_filter, ecg_derivative, ecg_square, ecg_int_smooth,indxs):
    time = np.arange(0, 1/fs* (ecg_filter.size), 1/fs, dtype=float)
    plt.plot(time, ecg_filter, color = 'blue')
    plt.plot(time, ecg_derivative, color = 'red')
    plt.plot(time, ecg_square, color='green')
    plt.plot(time, ecg_int_smooth, color='black')
    #plt.plot(time, np.gradient(ecg_int_smooth, 1 / fs))
    plt.scatter(time[indxs], ecg_int_smooth[indxs], marker = '*', color='yellow')
    plt.show()

"""
@brief Plot of the ecg signal, the qrs complexs, and the q,r,s peaks detected
"""
def final_plot(signal_filter, fs, periods, r_peaks_ind, q_peaks_ind, s_peaks_ind):

    fig, ax = plt.subplots()
    time = np.arange(0, 1 / fs * (signal_filter.size), 1 / fs, dtype=float)
    ax.plot(time, signal_filter, color='black', linewidth=0.25, label='ECG filter')
    for period in periods:
        ax.plot(time[period[0]:period[1]], signal_filter[period[0]:period[1]], linewidth = 0.4, color='red',alpha = 0.5)
    ax.scatter(time[q_peaks_ind], signal_filter[q_peaks_ind], color='black', marker='*', s=15, label='Q peak')
    ax.scatter(time[r_peaks_ind], signal_filter[r_peaks_ind], color='blue', marker='*', s=15, label='R peak')
    ax.scatter(time[s_peaks_ind], signal_filter[s_peaks_ind], color='green', marker='*', s=15, label='S peak')

    plt.legend(loc = 'lower left')
    #plt.title(f'ECG patient {path_name[-3:]}', loc='left')
    #matplotlib.rcParams['axes.spines.right'] = False
    #matplotlib.rcParams['axes.spines.top'] = False
    plt.ylabel('mV')
    plt.xlabel('Time [sec]')

    plt.show()
