### ecg_signal_processing

This repository contains code to detect the QRS complexes and the heart rate in ECG signals. The algorithm implemented is based on the Pan-Tompkins algorithm [1] to detect QRS complexes and the Raju et al. [2] algorithm to identify systolic peaks from the Arterial Blood Pressure signal. 

The methodology goes as follows: The signal is denoised by a Butterworth passband filter. The next step is differentiation, followed by squaring and integration with a moving windows filter. The resulting signal contains information about the slope and the duration of the QRS complex, therefore, as fig 5 in [1] shows, finding its raising intervals allows us to identify the QRS complexes. Because of physiological reasons, the minimum time for a valid QRS complex is 0.06s, moreover, there should be a 0.2s refractory period before the next one can occur. Finally, the Q, R, and S peaks are detected in each complex by selecting the local maximum and the two local minimums. The methods section presents a detailed explanation of the designed algorithm [1]. 

To evaluate its performance the number of QRS complexes detected in each signal has been compared to the number of QRS complexes detected while applying the neurokit2 library. The neurokit2 library allows to process ECG signals and detects the T, QRS, and P segments as well as their onset and offset points: (https://neuropsychology.github.io/NeuroKit/.) The overall performance of the algorithm shows a mean error of -3.96 +- 20.04, and an absolute error of 5.93 +- 19.24 respectively while comparing the detected peaks by eq. 1: 

mean error = (R peaks detected by designed algorithm - R peaks detected by neurokit2 library)/number of ECG (eq.1)
absolute mean error = |(R peaks detected by designed algorithm - R peaks detected by neurokit2 library)|/number of ECG signals (eq.1)

### Implementation
The following link contains a folder in google drive with ECG data from physionet, two Colab notebooks, a library called signal_processing with three .py files (analyze_data, ecg_detection, and analyze_data), and a csv file with the results (heart rate, heart rate periodicity, heart rate variability, number of peaks detected and number of peaks detected by neurokit2). To reproduce and try its implementation create an "acceso directo" of this folder in your folder "mi unidad" in your google drive. Finally, open the two Colab notebooks from Colab:

https://drive.google.com/drive/folders/1j8oZGBm5wskjtK8EsSeW6IM5YlW743xB?usp=sharing

Notebook statistical_analysis_database: It contains code to do a statistical analysis of the demographic data from the database
Notebook detect_qrs: It contains an example to detect the QRS complexes, the heart rate and its periodicity and variability, and a comparison of the results between the implemented algorithm and the neurokit2 library. 

### Methodology
The implemented algorithm follows two main steps: The QRS complexes detection and the Q, R, and S peaks identification

#### **1. QRS Complex detection:** 

The algorithm uses the derivation 'i' signal to detect the position of the QRS complexes. It contains 5 steps:

  **1.1. Denoised signal:** This substracts the first *tw* seconds of the signal to avoid noisy patterns from the sensor (ex. movement while recording the ECG). Then, it applies a Butterworth passband filter to eliminate frequencies lower than 5 and higher than 20 Hz. This emphasizes the QRS energy and reduces the influence of baseline and T-wave interference [1]. 
       
       Function: denoise_ecg in the qrs_detection.py file 
  
  **1.2. Derivative:** After filtering, the signal is differentiated to provide QRS complex slope information [2]. The differentiator is a gaussian derivative filter with a kernel size = 5 and a sigma value equal to 0.5. These parameters have been found experimentally. 
       
       Function: guassian_derivative_filter in the qrs_detection.py file
 
 **1.3.  Square Function:** After differentiation, the signal is square point by point to make all points positive and to emphasize the output of the derivative.
       
       Function: square_filter in the qrs_detection.py file
 
 **1.4.  Moving window integration:** The integration aims at finding the location of the QRS complex. As figure 5 in [1] shows the raising period of the integration corresponds to the duration of the QRS complex. The window's size of the integration is extremely important, it needs to include all the samples of the QRS complex. As normal QRS complexes range between 0.06 and 0.1 seconds [1] the kernel size has been set to 0.15s to assure anormal/bigger QRS complexes are included in the entire window. Therefore, the kernel size has been set to 150 samples due to the sampling frequency in this example is equal to 1000Hz.  
 
       Function: integration_filter in the qrs_detection.py file
       
  **1.5. Slope detection:** The detection of the increasing slope is crucial to identify the QRS complexes. To identify the slope of the signal its gradient must be higher than a threshold equal to 0.075 (determined empirically) during a period equal to or longer than 0.06 seconds. ECG signals are very sensitive to this threshold, therefore a recursive approach has been implemented, to reduce this value in case the algorithm is not capable to detect any QRS complex. 
       
       Function: slope_duration in the qrs_detection.py file
       
#### **2. QRS Complex adjustment, Q-R-S peak detection:** 

The second part of the algorithm aims at adjusting the QRS complexes detected and finding the Q-R-S peaks in each one of them. The user can select the derivation to apply this algorithm. The algorithm will use the QRS complexes found below to detect the peaks and adjust the QRS complex to the selected signal. This can be done due to all derivations are well aligned.

  **2.1. R peak:** The algorithm assigns the r peak to the maximum value in each of the QRS complexes. It also makes sure that no QRS complexes happen in a lower period than 0.2s (N=200 samples)
      
      Function: find_r in the qrs_detection.py file
      
  **2.2. Q and S peaks:** Given the position of the R peak, the algorithm looks for a local minimum to the left and to the right of the R index to find the Q and the S peak respectively. 
      
      Function: find_qs in the qrs_detection.py file
  
  **2.3. QRS compelexs adjustment:** To improve the onset and the offset detection of the QRS complex, the algorithm makes sure that all the values between the last point before the Q point with a higher value than the average filtered signal, and the Q point, as well as the values between the S point and the last point with a lower average value of the signal after it belongs to this complex. 
      
      Function: final corrections
      
      
 Finally, the heart rate, and its periodicity and variability has been calculated as follows: 
     Heart rate = (number of r peaks + 1) / (number of samples * 1/sample frequency)
     Periodicity = average(time(r_peak[i+1]-r_peak[i]))
     Variability = std(time(r_peak[i+1]-r_peak[i]))
     
 * The heart rate contains a +1 factor, due to the algorithm cannot detect the last peak of the signal. 
 

### Evaluation 
To evaluate the performance of the algorithm the number of r peaks detected by the neurokit2 library has been set as a reference point. The csv file present in the folder shows the results for the heart rate, the periodicity, the variability, and the number of peaks detected for the designed algorithm and the neurokit2 library. The algorithm does not work correctly (cannot detect any peak) for 17 out of the 531 signals (ecg0, ecg1 and ecg2 of patient number 10, ecg3 patient 29, ecg1 patient 38, ecg3 patient 97, ecg0 patient 148, ecg3 patient 180, ecg1 patient 183, ecg 1 patient 189, ecg0 patient 198, ecg0 patient 215 and ecg0 patient 220, ecg0, ecg1, ecg2 patient 278 and ecg 0 patient 285). The main reason is that the duration of the slope in the integral signal is shorter than 0.06 seconds, consequently, no QRS periods can be detected. Including these 17 signals, the algorithm presents a mean and an absolute mean error equal to -8.53 +- 31.73, and 10.41 +- 31.10 peaks of difference.

However, while excluding these signals, the mean and absolute mean error are reduced to -3.96 +- 20.04 and 5.93 +- 19.24 respectively. Besides, the designed algorithm cannot detect the last r peak of the signal, therefore adding this value will reduce by one more point the different number of peaks. This shows promising results and a promising research line to continue improving this algorithm.

### Improvements
This section proposes possible improvements and curiosities to this challenge:

1. Artifact detection: Patient number 5 shows values over and below the physiological range for an ECG signal. It would have been interesting to apply an artifact detection filter to detect values out of the normal range and correct them by signal processing techniques such as interpolation

2. Improving the QRS detection of the 17 signals where the algorithm does not work correctly

3. Analysis of the relation of heart rate, periodicity, and variability regarding the reason for admission and demographic data like gender or smoker patients

4. Apply different techniques to calculate heart rate, such as calculating the heart rate in a window of 1 second over the entire signal and averaging the results





[1] PAN, Jiapu; TOMPKINS, Willis J. "A real-time QRS detection algorithm." IEEE transactions on biomedical engineering, 1985, no 3, p. 230-236.
[2] D. S. Raju, M. S. Manikandan, and R. Barathram. “An automated method for detecting systolic peaks from arterial blood pressure signals” In Proceedings of the 2014 IEEE Students Technology Symposium, 2014, pp. 41–46
