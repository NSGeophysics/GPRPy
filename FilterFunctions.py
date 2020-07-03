from scipy.signal import butter,filtfilt
from scipy import signal,fftpack
from scipy.signal import freqz
import numpy as np
import matplotlib.pyplot as plt

def nextpow2(i):
    n = 1
    while n < i: 
        n *= 2
    return n 

def butter_bandpass(lowcut, highcut, fs, order):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='bandpass', analog=False)
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    #w, h = freqz(b, a, worN=2000)
    #plt.plot((fs * 0.5 / np.pi) * w, abs(h), label="order = %d" % order)
    y = filtfilt(b,a,data)
    return y

def butter_lowpass(highcut, fs, order):
    nyq = 0.5 * fs
    high = highcut / nyq
    b, a = butter(order, high, btype='lowpass', analog=False)
    return b, a

def butter_lowpass_filter(data, highcut, fs, order):
    b, a = butter_lowpass(highcut, fs, order=order)
    y = filtfilt(b,a,data)
    return y

def butter_highpass(lowcut, fs, order):
    nyq = 0.5 * fs
    low = lowcut / nyq
    b, a = butter(order, low, btype='highpass', analog=False)
    return b, a

def butter_highpass_filter(data, lowcut, fs, order):
    b, a = butter_highpass(lowcut, fs, order=order)
    y = filtfilt(b,a,data)
    return y

def filtre_hautbas_calcul(var_filtre, data, cut, fs, ordre):
    if var_filtre == 1:
        datain_filtered = butter_lowpass_filter(data, cut, fs, ordre)
    elif var_filtre == 2:
        datain_filtered = butter_highpass_filter(data, cut, fs, ordre)
    num = len(datain_filtered)
    N = nextpow2(num)
    A_filt = fftpack.fft(datain_filtered, N)
    A_shift_filt = fftpack.fftshift(A_filt)
    Xpos_filt = A_shift_filt[N//2:]
    return datain_filtered, Xpos_filt, N

def filtre_passband_calcul(twtt,dataset, lowcut, highcut, fs, ordre, numero_trace, trace_pour_afficher):
    data = dataset[numero_trace]
    datain_filtered = butter_bandpass_filter(data, lowcut, highcut, fs, ordre)
    num=len(datain_filtered)
    N = nextpow2(num)
    A_filt = fftpack.fft(datain_filtered, N)
    A_shift_filt = fftpack.fftshift(A_filt)
    Xpos_filt = A_shift_filt[N//2:]
    
    if numero_trace == trace_pour_afficher :
        
        b, a = butter_bandpass(lowcut, highcut, fs, ordre)
        w, h = freqz(b, a, worN=2000)
        
        plt.subplot(231)
        plt.plot(twtt, data.tolist()[0])
        plt.xlabel("Temps (ns)")
        plt.ylabel("Amplitude")
        plt.title('Signal')
        
        plt.subplot(232)
        plt.plot((fs * 0.5 / np.pi) * w, abs(h), label="order = %d" % ordre)
        plt.xlabel("Fréquence (GHz)")
        plt.ylabel("Gain (dB)")
        plt.title('Filtre de Butterworth')
        
        plt.subplot(233)
        plt.plot(twtt, datain_filtered.tolist()[0])
        plt.xlabel("Temps (ns)")
        plt.ylabel("Amplitude")
        plt.title('Signal filtré')
    
        plt.waitforbuttonpress(0) # this will wait for indefinite time
        plt.close()

    
    return datain_filtered, Xpos_filt, N
