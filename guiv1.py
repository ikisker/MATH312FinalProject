import tkinter as tk
from tkinter import filedialog, messagebox
import wavio
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import simpleaudio as sa

k = None
t = None
audio_data = None
fig = None
ax = None
canvas = None

def load_audio():
    global k, t, audio_data, fig, ax, canvas
    audio_file = filedialog.askopenfilename(title="Select Audio File", filetypes=[("WAV Files", "*.wav")])
    if audio_file:
        k = wavio.read(audio_file)
        audio_data = k.data[:, 0]
        t = np.linspace(0, len(audio_data) / k.rate, len(audio_data))
        plot_audio()

def plot_audio():
    global fig, ax, canvas
    fig = Figure(figsize=(6, 4), dpi=100)
    ax = fig.add_subplot(111)
    ax.clear()
    ax.plot(t, audio_data)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Amplitude")
    ax.set_title("Audio Waveform")

    if canvas is None:
        canvas = FigureCanvasTkAgg(fig, master=root)
        canvas.draw()
        canvas.get_tk_widget().pack()
    else:
        canvas.draw()

def play_audio():
    if audio_data is not None:
        audio = sa.play_buffer(audio_data, 1, 2, k.rate)
        audio.wait_done()

def encrypt_audio():
    if k is None:
        messagebox.showerror("Error", "Please load an audio file first.")
        return

    wavio.write("ch1.wav", k.data[:12*44100,0], 44100, sampwidth=2)

    def Lorentz(t, a):
        u, v, w = a
        dotu = 16*(v-u)
        dotv = u*(45.6-20*w)-v
        dotw = 5*u*v - 4*w
        return np.array([dotu, dotv, dotw])

    ic = np.array([-2.3, -1.9, 2.2])
    tend = 20000
    t = np.linspace(0, tend, 12*44100)
    sol = solve_ivp(Lorentz, [0, tend], ic, t_eval=t, atol=1e-9, rtol=1e-8)
    u = sol.y[0, :]

    noisesize = np.sqrt((u**2).sum())
    musicsize = np.sqrt((k.data[:12*44100,0]**2).sum())
    signal = u + k.data[:12*44100,0]*noisesize/musicsize/1000
    signal1 = signal-signal.mean()
    signal1 = signal1 / np.abs(signal1).max()
    wavio.write("encrypted_audio.wav", signal1, 44100, sampwidth=2)

    messagebox.showinfo("Encryption Successful", "Audio file has been encrypted and saved as 'encrypted_audio.wav'")

def decrypt_audio():
    if k is None:
        messagebox.showerror("Error", "Please load an audio file first.")
        return

    signal = k.data[:, 0]

    def Lorentz2(tt, a):
        U, V, W = a
        uu = CubicSpline(t, signal)(tt)
        dotU = 16*(V-U)
        dotV = uu*(45.6-20*W)-V
        dotW = 5*uu*V - 4*W
        return np.array([dotU, dotV, dotW])

    ic = np.array([-2.3, -1.9, 2.2])
    tend = 20000
    t = np.linspace(0, tend, 12*44100)
    sol = solve_ivp(Lorentz2, [0, tend], ic, t_eval=t, atol=1e-9, rtol=1e-8)
    U = sol.y[0, :]

    recoveredmusic = U - signal
    recoveredmusic = recoveredmusic - recoveredmusic.mean()
    recoveredmusic = recoveredmusic / np.abs(recoveredmusic).max()
    wavio.write("decrypted_audio.wav", recoveredmusic, 44100, sampwidth=2)

    messagebox.showinfo("Decryption Successful", "Audio file has been decrypted and saved as 'decrypted_audio.wav'")

root = tk.Tk()
root.title("Audio Encryption/Decryption")

# Create a placeholder empty graph
fig = Figure(figsize=(6, 4), dpi=100)
ax = fig.add_subplot(111)
ax.set_xlabel("Time (s)")
ax.set_ylabel("Amplitude")
ax.set_title("Audio Waveform")
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.draw()
canvas.get_tk_widget().pack()

load_button = tk.Button(root, text="Load Audio", command=load_audio)
load_button.pack(pady=10)

play_button = tk.Button(root, text="Play Audio", command=play_audio)
play_button.pack(pady=10)

encrypt_button = tk.Button(root, text="Encrypt Audio", command=encrypt_audio)
encrypt_button.pack(pady=10)

decrypt_button = tk.Button(root, text="Decrypt Audio", command=decrypt_audio)
decrypt_button.pack(pady=10)

root.mainloop()