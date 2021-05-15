# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 15:35:16 2020

@author: skukm
"""

import serial 
from pyfirmata import Arduino
import numpy as np
import matplotlib.pyplot as plt


ard=serial.Serial('COM4',9600)

data=list() #input data
for i in range(2850):
    a=ard.readline().decode('utf-8')
    a=float(a)
    print(a)
    data.append(a)

data=np.array(data)
 
from biosppy.signals import ecg

#Rpeaks, filtering
def peaks_to_series(signal,peaks):
    pks=np.zeros(len(signal))
    pks[peaks]=1
    return pks

def filter_signal(signal, sampling_rate):
    order=int(0.3*sampling_rate)
    filtered, _,_,=ecg.st.filter_signal(signal=signal,
                                        ftype='FIR',
                                        band='bandpass',
                                        order=order,
                                        frequency=[3,45],
                                        sampling_rate=sampling_rate)
    rpeaks,=ecg.hamilton_segmenter(signal=filtered, sampling_rate=sampling_rate)
    rpeaks_corrected,=ecg.correct_rpeaks(signal=filtered,
                                         rpeaks=rpeaks,
                                         sampling_rate=sampling_rate,
                                         tol=0.05)
    rpeaks=peaks_to_series(signal, rpeaks)
    rpeaks_corrected=peaks_to_series(signal, rpeaks_corrected)
    return filtered, rpeaks_corrected

#Q,R,S 특징점 검출
def QRS(ecg,rpeaks):
    R=np.array(0)
    for i in range(len(rpeaks)):
        if rpeaks[i]==1:
            R=np.append(R,i)
    R=np.array(R)
    Q=np.zeros(len(R))
    S=np.zeros(len(R))
    #Q
    for i in range(1,len(R)): 
        q=[]
        for j in range(round((R[i-1]+R[i])/2), R[i]):    
            if ecg[j] < ecg[j+1]:
                q.append(j)
        Q[i]=(q[np.argmin(ecg[q])])
     #S   
    for i in range(1,len(R)): 
        s=[]
        if i==len(R)-1: 
            for j in range(R[i], len(ecg)-1): 
                if ecg[j] < ecg[j+1]:
                    s.append(j)
            S[i]=(s[np.argmin(ecg[s])])
            break

        for j in range(R[i], round((R[i]+R[i+1])/2)): 
            if ecg[j] < ecg[j+1]:
                s.append(j)
        S[i]=(s[np.argmin(ecg[s])])
    return Q,R,S

Fs=190
T=1/Fs
t=np.arange(0,len(data)/Fs,T)
ECG,rpeaks=filter_signal(data, Fs) #Rpeaks
Q,R,S=QRS(ECG,rpeaks)

#qrx x축
x_qrs=np.array([Q,R,S],dtype=int)

#qrs y축
def y_qrs(x_qrs):
    q,r,s=np.array(()),np.array(()),np.array(())
    q=[np.append(q,ECG[i]) for i in x_qrs[0]]
    r=[np.append(r,ECG[i]) for i in x_qrs[1]]
    s=[np.append(s,ECG[i]) for i in x_qrs[2]]
    y_qrs=np.array((q,r,s),dtype=float).reshape(3,len(R))
    return y_qrs

y_qrs=y_qrs(x_qrs)
x_qrs=np.transpose(x_qrs).reshape(len(R)*3)
y_qrs=np.transpose(y_qrs).reshape(len(R)*3)

# plt.plot(t,ECG)
# plt.plot(t,y_qrs[0,][1:16])
# plt.plot(t,y_qrs[1,][1:16])
# plt.plot(t,y_qrs[2,][1:16])

#특정 벡터 생성
#vector_F[:,i] x_qrs,y_qrs i=3,j=2

vector_F=np.transpose(np.vstack((x_qrs,y_qrs))).reshape(len(R),3,2)
n=np.array([list(rpeaks[:i*Fs]).count(1) for i in range(int(len(ECG)*T+1))])

def aver_vector(Weight=np.ones(n[-1]), **args):
    vector=args["vector"]
    AQRSF=np.zeros(3*2).reshape(3,2)
    for p in range(1,n[-1]):
        AQRSF=AQRSF+(vector[p]*Weight[p])
    AQRSF=1/n[-1]*AQRSF
    return AQRSF

def Weight(**args):
    AQRSF=args["AQRSF"]
    vector=args["vector"]

    abs_ARQSF_1=np.linalg.norm(AQRSF[:,0])
    abs_ARQSF_2=np.linalg.norm(AQRSF[:,1])
    X_simil=np.zeros(n[-1])
    Y_simil=np.zeros(n[-1])
    
    for p in range(1,n[-1]):
        abs_vector_F_1=np.linalg.norm(vector[p][:,0])
        abs_vector_F_2=np.linalg.norm(vector[p][:,1])   
        X_simil[p]=np.dot(AQRSF[:,0],vector[p][:,0])/(abs_ARQSF_1*abs_vector_F_1)
        Y_simil[p]=np.dot(AQRSF[:,1],vector[p][:,1])/(abs_ARQSF_2*abs_vector_F_2)
        
    error=np.array((X_simil,Y_simil))
    error=1-error
    W=(X_simil+Y_simil)/2    
    return W, error

def Feature_vector(vector_f):
    AQRSF=aver_vector(vector=vector_f)
    W,error=Weight(AQRSF=AQRSF,vector=vector_f)
    AQRSF_W=aver_vector(Weight=W, vector=vector_f)
    return AQRSF_W

# AQRSF=aver_vector(vector=vector_F)
# W,error=Weight(AQRSF=AQRSF, vector=vector_F)
# AQRSF_W=aver_vector(Weight=W, vector=vector_F)

User_vector=Feature_vector(vector_F)
User=User_vector.reshape(6,order='F')
f= np.loadtxt(r'C:\Users\skukm\Downloads\data.txt')
if f.size==0:    
    f = open("C:/Users/skukm/Downloads/data.txt","w")
    for i in range(len(User)):
        f.write(str(User[i]))
        f.write("\n")
    f.close()
f= np.loadtxt(r'C:\Users\skukm\Downloads\data.txt')


def one_hot_encode(sequence, n_features):
    encoding=list()
    for value in sequence:
        vector=[0 for _ in range(n_features)]
        vector[int(value)] =1
        encoding.append(vector)
    return np.array(encoding)

n_features =2000
encoding0=one_hot_encode(User, n_features)
encoding1=one_hot_encode(f, n_features)

fff=np.array([encoding0,encoding1])
timestep = fff.shape[1]
n_features = fff.shape[2]

from sklearn.preprocessing import LabelEncoder
symptom=np.array((["diff"],["same"]))
label_encoder=LabelEncoder()
integer_encoded = label_encoder.fit_transform(symptom)
y=one_hot_encode(integer_encoded,n_features)

from tensorflow.keras.layers import LSTM, Dense
from tensorflow.keras.models import Sequential
model = Sequential()
model.add(LSTM(128, input_shape=(timestep, n_features)))
model.add(Dense(n_features, activation='softmax'))
model.compile(loss='categorical_crossentropy',optimizer= 'adam',metrics=['acc'])
model.summary()

X =fff
for i in range(1000):
    model.fit(X,y,epochs=1, verbose=0)

encoded =one_hot_encode(User, n_features)
X=encoded.reshape(1,len(User),n_features)
yhat = model.predict(X)
from numpy import argmax
inverted_yhat=label_encoder.inverse_transform([argmax(yhat)])
print(inverted_yhat)


ard.close()
