# -*- coding: utf-8 -*-
"""
Created on 2018-03-24

@author: Harry S. Reichard
"""
# import csv
import math
import matplotlib.pyplot as plt
import numpy as np
import cmath  # Complex math.

# from matplotlib.ticker import FuncFormatter
# import pylab
# from scipy import signal
# from copy import copy
# from scipy.interpolate import spline


def AllPassZTransformCoefficients(w0, Q, Ts, mode):
    # Returns a0, a1, a2, where (note negative exponents of z)
    # F(z) = a0 + a1 * z^-1 + a2 * z^-2
    w0 *= Ts  # Normalize by sampling interval Ts.
    k1 = math.cos(w0)
    k2 = math.sin(w0)
    k3 = math.cos(2.0 * w0)
    # k4 = math.sin(2.0 * w0)  # Not used.
    k5 = -1.0 / math.tan(w0)
    k6 = k2 - k1 * k5
    k7 = -k5
    k8 = -4.0 * Q / w0
    k9 = k8 - 2.0
    k10 = 2.0 * k1 * (k8 - 1.0)
    k11 = 2.0 * k1 * (k8 + 1.0)
    k12 = 2.0 * k3 * k8
    k13 = k8 + 2.0
    # k14 = -1.0
    k15 = -k6 / k7
    k16 = k8 + k9 * k15 * k15 + k10 * k15
    k17 = -2.0 * k9 * k15 - k10 + k11 + k12 * k15
    k18 = k9 - k12 + k13
    a0 = 1.0  # Note:  a0 = 1.
    a1 = k17 * k17 - 4.0 * k16 * k18
    if a1 >= 0.0:
        a1 = math.sqrt(a1)
    else:
        a1 = 0.0
    a1 = (-k17 - a1) / (2.0 * k16)  # Select negative root.
    a2 = -1.0 + k15 * a1
    """
    print('k1 =', k1)
    print('k2 =', k2)
    print('k3 =', k3)
    print('k4 =', k4)
    print('k5 =', k5)
    print('k6 =', k6)
    print('k7 =', k7)
    print('k8 =', k8)
    print('k9 =', k9)
    print('k10', k10)
    print('k11', k11)
    print('k12', k12)
    print('k13', k13)
    print('k14', k14)
    print('k15', k15)
    print('k16', k16)
    print('k17', k17)
    print('k18', k18)
    print('a0', a0)
    print('a1', a1)
    print('a2', a2)
    print(' ')
    """
    return a0, a1, a2


def AmplitudeVsFrequency(ctbq, dtbq,  w0, Q, Ts, mode, frequencyrange, df):

    f = np.arange(df, frequencyrange, df)
    y1 = []
    y2 = []
    [y1.append(abs(ctbq.FrequencyResponse(item))) for item in f]
    [y2.append(abs(dtbq.FrequencyResponse(item))) for item in f]

    fig, ax = plt.subplots()

    # X-axis
    ax.set_xlabel('Frequency (Hz)')
    ax.xaxis.grid(which='major')
    ax.tick_params(axis='x', which='major', bottom='on')
    # ax.set_xlim()

    # Y-axis
    ax.set_ylabel('Amplitude')
    ax.yaxis.grid(which='both')
    ax.tick_params(axis='y', which='minor', left='on')
    # ax.set_ylim(0.0, +1.05)

    # Curve(s)
    ax.plot(f, y1, 'r:', label='Continuous')
    ax.plot(f, y2, 'r', label='Discrete')
    # ax.semilogy(x, y, 'r', label='')
    ax.legend(loc='upper right', bbox_to_anchor=(1, 0.84), framealpha=1)

    # Caption
    line0 = mode
    line1 = ''
    if line1 == '':
        fig.subplots_adjust(top=0.92)  # Single line title.
    else:
        line0 += '\n' + line1  # Two-line title.
        fig.subplots_adjust(top=0.86)
    fig.suptitle(line0, fontsize=10, fontweight='bold')

    # Text box
    text = 'fo = ' + '{:3.1f}'.format(w0 / (2.0 * math.pi)) + ' Hz'
    text += '\nQ = ' + '{:2.1f}'.format(Q)
    # ax.text(0.97, 0.95, text, fontsize=10, horizontalalignment='right',
    #         verticalalignment='top', transform=ax.transAxes,
    #         bbox=dict(facecolor='white', alpha=1.0))
    # plt.legend()
    plt.show()
    print('Response at fo =', abs(dtbq.FrequencyResponse(w0 / (2.0 * math.pi))))


class ContinuousTimeBiquad:

    def __init__(self, w0, Q, Ts, mode):
        # Response vs frequency only.  Ts is disregarded.
        self.w0 = w0

        self.d0 = +1.0
        self.d1 = +1.0 / (Q * w0)
        self.d2 = +1.0 / (w0 * w0)

        if mode == 'Allpass':
            self.n0 = +self.d0
            self.n1 = -self.d1
            self.n2 = self.d2
        elif mode == 'Bandpass':
            self.n0 = 0.0
            self.n1 = self.d1
            self.n2 = 0.0
        elif mode == 'Highpass':
            self.n0 = 0.0
            self.n1 = 0.0
            self.n2 = self.d2 / Q
        elif mode == 'Lowpass':  # ***Single pole***
            self.n0 = 1.0
            self.n1 = 0.0
            self.n2 = 0.0
            self.d0 = 1.0
            self.d1 = 1.0 / w0
            self.d2 = 0.0
        elif mode == 'Notch':
            self.n0 = self.d0
            self.n1 = 0.0
            self.n2 = self.d2

    def FrequencyResponse(self, f):
        s = 0.0 + 2.0j * math.pi * f
        y = self.n0 + self.n1 * s + self.n2 * s * s
        y /= self.d0 + self.d1 * s + self.d2 * s * s
        return y

    def Next(self, t):
        return 1.0 - math.exp(-self.w0 * t)
# end class ContinuousTimeBiquad


def DiscreteTimeBiquad_1(w0, Q, Ts, mode):
    # Returns
    # n0, n1, n2 numerator coefficients of z^-1.
    # d0, d1, d1 denominator coefficients of z^-1.

    d0, d1, d2 = AllPassZTransformCoefficients(w0, Q, Ts, mode)
    # d0, d1, d2 = LaplaceToZ(1.0, 1.0 / Q, 1.0, w0 * Ts)

    z = cmath.exp(w0 * Ts * 1.0j)
    # print('F(jw0) =', np.angle(d0 + d1 / z + d2 / (z * z), deg=True))
    print('F(jw0) =', np.angle(d0 / (z * z) + d1 / z + d2, deg=True))

    if mode == 'Allpass':
        n0 = +d2
        n1 = d1
        n2 = d0
    elif mode == 'Bandpass':
        n0 = (d0 - d2) / 2.0
        n1 = 0.0
        n2 = (d2 - d0) / 2.0
    elif mode == 'Lowpass':  # ***single pole***
        n0 = 1 - math.exp(-w0 * Ts)
        n1 = 0.0
        n2 = 0.0
        d0 = 1.0
        d1 = -math.exp(-w0 * Ts)
        d2 = 0
    elif mode == 'Lowpasx':
        n0 = (-Q * d1 + 1.0 + d1 + d2) / (2.0 * Q)
        n1 = 0.0
        n2 = n0 + d1
    elif mode == 'Notch':
        n0 = (d0 + d2) / 2.0
        n1 = d1
        n2 = (d2 + d0) / 2.0
    elif mode == 'Highpass':
        n0 = 2.0 * (d0 - Ts * Ts) / Q
        n1 = 2.0 * d1 / Q
        n2 = 2.0 * d2 / Q

    return n0, n1, n2, d0, d1, d2


class DiscreteTimeBiquad_2:

    def __init__(self, n0, n1, n2, d0, d1, d2, Ts):
        # n0, n1, n2 numerator coefficients of z^-1.
        # d0, d1, d1 denominator coefficients of z^-1.
        self.n0 = n0
        self.n1 = n1
        self.n2 = n2
        self.d0 = d0
        self.d1 = d1
        self.d2 = d2
        self.Ts = Ts

        self.y1 = 0.0  # Signal at z^-1 node.
        self.y2 = 0.0  # Signal at z^-2 node.

    def Next(self, x):
        y0 = 1.0 / self.d0 * (x - self.d1 * self.y1 - self.d2 * self.y2)
        y = self.n0 * y0 + self.n1 * self.y1 + self.n2 * self.y2
        self.y2 = self.y1
        self.y1 = y0
        return y

    def FrequencyResponse(self, f):
        # ζ (lower case zeta) = z^-1.
        ζ = cmath.exp(-2.0j * math.pi * f * self.Ts)
        y = self.n0 + self.n1 * ζ + self.n2 * ζ * ζ
        y /= self.d0 + self.d1 * ζ + self.d2 * ζ * ζ
        return y
# End class DiscreteTimeBiquad_2


def KroneckerImpulse(start, end, Ts):
    armed = True
    t = []
    y = []
    for i in range(0, math.ceil((end - start) / Ts)):
        t.append(start + i * Ts)
        if armed:
            if t[i] >= 0:
                y.append(1.0)
                armed = False
        else:
                y.append(0.0)
    return np.array(t), np.array(y)


def LaplaceToZ(A0, A1, A2, Ts, Mode):
    # Convert a quadratic in s [A0 + A1 * s + A2 * s^2] to a corresponding 
    #quadratic in z^-1 [a0 + a1 * z^-1 + a2 * z^-2] in accordance with mode
    # (difference equation, step invariance, impulse invariance, etc.)
    
    if Mode == 'difference':
        #  Converts a quadratic in s to an equivalent quadratic in z^-1, using
        #  a difference equation approximation.
        #  F(s) = A0 + A1 * s + A2 * s^2
        #  s = (1 - z^-1) / Ts  #  Ts is the sampling interval.
        #  F(z) = a0 + a1 * z^-1 + a2 * z^-2
        a0 = A0 + A1 / Ts + A2 / (Ts * Ts)
        a1 = -A1 / Ts - 2.0 * A2 / (Ts * Ts)
        a2 = A2 / (Ts * Ts)
        return a0, a1, a2


def PhaseVsFrequency(ctbq, dtbq, w0, Q, Ts, mode, frequencyrange, df):

    f = np.arange(df, frequencyrange, df)
    y1 = []
    y2 = []
    [y1.append(np.angle(ctbq.FrequencyResponse(item), deg=True)) for item in f]
    [y2.append(np.angle(dtbq.FrequencyResponse(item), deg=True)) for item in f]
    fig, ax = plt.subplots()

    # X-axis
    ax.set_xlabel('Frequency (Hz)')
    ax.xaxis.grid(which='major')
    ax.tick_params(axis='x', which='major', bottom='on')
    # ax.set_xlim()

    # Y-axis
    ax.set_ylabel('Phase (degrees)')
    ax.yaxis.grid(which='both')
    ax.tick_params(axis='y', which='minor', left='on')
    ax.set_ylim(-180.0, +180.0)

    # Curve(s)
    ax.plot(f, y1, 'r:', label='')
    ax.plot(f, y2, 'r', label='')
    # ax.semilogy(x, y, 'r', label='')

    # Caption
    line0 = mode
    line1 = ''
    if line1 == '':
        fig.subplots_adjust(top=0.92)  # Single line title.
    else:
        line0 += '\n' + line1  # Two-line title.
        fig.subplots_adjust(top=0.86)
    fig.suptitle(line0, fontsize=10, fontweight='bold')

    # Text box
    text = 'fo = ' + '{:3.1f}'.format(w0 / (2.0 * math.pi)) + ' Hz'
    text += '\nQ = ' + '{:2.1f}'.format(Q)
    ax.text(0.80, 0.95, text, fontsize=9, horizontalalignment='left',
            verticalalignment='top', transform=ax.transAxes,
            bbox=dict(facecolor='white', alpha=1.0))

    plt.show()
    print('Phase at fo:', np.angle(dtbq.FrequencyResponse(w0 / (2.0 * math.pi)), deg=True), 'degrees')


def RampFunction(start, end, Ts, rampduration):
    t = []
    y = []
    for i in range(0, math.ceil((end - start) / Ts)):
        t.append(start + i * Ts)
        if t[i] <= 0:
            y.append(0.0)
        elif t[i] < rampduration:
            y.append(t[i] / rampduration)
        else:
            y.append(1.0)
    return np.array(t), np.array(y)


def SinglePoleCoefficients(w0, Ts, mode="lowpass"):
        α1 = 1.0 / (1.0 + Ts)
        α2 = 0
        K = Ts / (1.0 + Ts)
        return α1, α2, K


def StepFunction(start, end, Ts):
    t = []
    y = []
    for i in range(0, math.ceil((end - start) / Ts)):
        t.append(start + i * Ts)
        if t[i] <= 0:
            y.append(0.0)
        else:
            y.append(1.0)
    return np.array(t), np.array(y)


def StepResponse(ctbq, dtbq, w0, Q, Ts, mode, timebase):
    # t, y1 = StepFunction(-0.03 * timebase, timebase, Ts)
    t, y1 = StepFunction(-0.00 * timebase, timebase, Ts)
    # t, y1 = KroneckerImpulse(-0.00 * timebase, timebase, Ts)
    y2 = []
    [y2.append(dtbq.Next(item)) for item in y1]
    y3 = []
    [y3.append(ctbq.Next(item)) for item in t]

    fig, ax = plt.subplots()

    # X-axis
    ax.set_xlabel('Time (s)')
    ax.xaxis.grid(which='major')
    ax.tick_params(axis='x', which='major', bottom='on')
    # ax.set_xlim()

    # Y-axis
    ax.set_ylabel('Step response')
    ax.yaxis.grid(which='both')
    ax.tick_params(axis='y', which='minor', left='on')
    # ax.set_ylim()

    # Curve(s)
    ax.plot(t, y2, 'r', label='')
    if mode == 'Lowpass':
        ax.plot(t, y3, 'rs', label='')
    # ax.semilogy(x, y, 'r', label='')

    # Caption
    line0 = mode
    line1 = ''
    if line1 == '':
        fig.subplots_adjust(top=0.92)  # Single line title.
    else:
        line0 += '\n' + line1  # Two-line title.
        fig.subplots_adjust(top=0.86)
    fig.suptitle(line0, fontsize=10, fontweight='bold')

    # Text box
    text = 'fo = ' + '{:3.1f}'.format(w0 / (2.0 * math.pi)) + ' Hz'
    text += '\nQ = ' + '{:2.1f}'.format(Q)
    # ax.text(0.80, 0.95, text, fontsize=9, horizontalalignment='left',
    #         verticalalignment='top', transform=ax.transAxes,
    #         bbox=dict(facecolor='white', alpha=1.0))

    plt.show()


def __Main__():

    # Parameters of continuous-time prototype.
    f0 = 3.0  # Characteristic frequency.
    w0 = 2.0 * math.pi * f0
    w0 = 1.0  # Assign w0, rather than f0.
    f0 = w0 / (2.0 * math.pi)
    Q = 2.0

    # Parameters of discrete-time model.
    n = 10  # Samples per period (of characteristic frequency).
    # n = 10  # Samples per period (of characteristic frequency).
    Ts = 1.0 / (n * f0)

    # Display parameters.
    timebase = 3.0
    timebase = 5 / w0  # i.e. five time constants.
    frequencyrange = n * f0  # I.e., to Nyquist frequency.
    # timebase =
    # frequencyrange =
    df = 1e-3 * frequencyrange
    # mode = 'Allpass'
    # mode = 'Bandpass'
    # mode = 'Highpass'
    mode = 'Lowpass'
    # mode = 'Notch'

    # dtbq = DiscreteTimeBiquad(w0, Q, Ts, mode)
    n0, n1, n2, d0, d1, d2 = DiscreteTimeBiquad_1(w0, Q, Ts, mode)
    dtbq = DiscreteTimeBiquad_2(n0, n1, n2, d0, d1, d2, Ts)
    ctbq = ContinuousTimeBiquad(w0, Q, Ts, mode)
    StepResponse(ctbq, dtbq, w0, Q, Ts, mode, timebase)
    AmplitudeVsFrequency(ctbq, dtbq, w0, Q, Ts, mode, frequencyrange, df)
    # PhaseVsFrequency(ctbq, dtbq, w0, Q, Ts, mode, frequencyrange, df)


__Main__()
