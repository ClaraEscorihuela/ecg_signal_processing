import os
import sys
curr_dir = os.getcwd()+'\signal_processing'
sys.path.append(curr_dir)

import analyze_data
import qrs_detection

import numpy as np

class ecg_waves():
    def __init__(self):
        self.waves = []

    def add_wave(self, wave):
        '''
        Function to save all waves from the same patient
        :param subject: class subject
        :return: list of subjects
        '''
        self.waves.append(wave)

class ecg_wave():
    def __init__(self, fs, signal_1, signal_2):
        self.signal_1 = signal_1  # derivation i, used to detect the qrs complex
        self.signal_2 = signal_2  # derivation selected by the user to detect the peaks and plot the ecg
        self.fs = fs
        self.period = 1/fs
        self.duration = self.period * signal_1.size

    def qrs_complex(self, fs=1000, wind_int=150, m_gaus=20, sigma_gaus=6, th_grad=0.1, t_qrs_min=0.06):
        '''
        Detect qrs complexes (refer to processing in signal_processing file)
        '''
        periods = qrs_detection.processing(self.signal_1, fs=fs, wind_int=wind_int, m_gaus=m_gaus, sigma_gaus=sigma_gaus, th_grad=th_grad, t_qrs_min=t_qrs_min)
        self.periods = periods


    def peaks(self, t_ref=0.2):
        '''
        Detect q-r-s peaks (refer to peaks in signal_processing file)
        '''
        periods, r_peaks_ind, q_peaks_ind, s_peaks_ind = qrs_detection.peak_detection(self.signal_2, self.periods, t_ref=t_ref,
                                                                                     fs=self.fs)
        self.periods = periods
        self.r_peaks_ind = r_peaks_ind
        self.q_peaks_ind = q_peaks_ind
        self.s_peaks_ind = s_peaks_ind
        self.q_onset_ind = [period[0] for period in periods]
        self.s_offset_ind = [period[1] for period in periods]

    def ecg_fix(self):
        '''
        Fix ecf for patient where the algorithm does not properly work
        '''
        self.r_peaks_ind = [0]


    def calculate_heart_rate(self):
        """
        Calculate heart rate, periodicity and variability per ecg
        Heart rate: Take the number of qrs and divide by the length of the signal (the last two beats are never detected,
        therefore, we need to sum two to the number of beats)
        Periodic beat: Calculate the mean for the time two consecutive r peaks happen
        Variaility beat: Calculate the std for the time two consecutive r peaks happen
        """
        if type(self.periods) != bool:
            self.heart_rate = 60 * (len(self.r_peaks_ind) + 2) / self.duration  # beats/60sec = beats/min
            self.periodic_beat = np.mean(np.diff(self.r_peaks_ind) * self.period)  # rate for having a beat
            self.variability_beat = np.std(np.diff(self.r_peaks_ind) * self.period)  # rate variability for beat

        else:
            # There has been a problem
            self.heart_rate = 0
            self.periodic_beat = 0
            self.variability_beat = 0

    def main_properties(self,fs=1000, wind_int=150, m_gaus=20, sigma_gaus=6, th_grad=0.1, t_qrs_min=0.06,  t_ref=0.2, plot = True):
        """
        Pipeline for all the functions in the wave object
        """

        self.qrs_complex( fs=fs, wind_int=wind_int, m_gaus=m_gaus, sigma_gaus=sigma_gaus, th_grad=th_grad, t_qrs_min=t_qrs_min)
        if type(self.periods) != bool:
            self.peaks(t_ref=t_ref)

            if plot:
                self.final_plot()

        self.calculate_heart_rate()


        print(f'\n Heart rate = {self.heart_rate} (beats/min) \n'
              f'Heart rate periodicity = {self.periodic_beat} (mean seconds between beats) \n'
              f'Heart rate variability = {self.variability_beat} (std seconds between beats)')

    def final_plot(self):
        '''
        Plot of the ecg with the qrs complexes, and the q r s peaks (refer to plot in signal_processing file)
        '''
        signal_filter = qrs_detection.denoise_ecg(self.signal_2, fs=1000, w1=2, w2=50)
        qrs_detection.final_plot(signal_filter, self.fs, self.periods, self.r_peaks_ind, self.q_peaks_ind, self.s_peaks_ind)










