import tkinter as tk
from tkinter import filedialog
import wavio
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt

sample_rate = None
audio_data = None
u = None

# Existing Lorentz function from the original code
def Lorentz(t, a):
    u = a[0]  # unpack x
    v = a[1]  # unpack y
    w = a[2]  # unpack z
    dotu = 16*(v-u)
    dotv = u*(45.6-20*w)-v
    dotw = 5*u*v - 4*w
    return np.array([dotu, dotv, dotw])  # return vector of derivatives

def load_audio_file():
    global audio_data, sample_rate
    file_path = filedialog.askopenfilename(title="Load Audio File")
    if file_path:
        try:
            wav_obj = wavio.read(file_path)
            audio_data = wav_obj.data
            sample_rate = wav_obj.rate
            # Call functions to process audio data and update waveform display
            plot_waveform(input_fig, audio_data, "Input Signal")
        except Exception as e:
            print(f"Error loading audio file: {e}")

def generate_noise():
    global u
    ic = np.array([-2.3, -1.9, 2.2])  # Initial conditions
    tend = len(audio_data) / sample_rate  # Calculate the duration of the audio file
    t = np.linspace(0, tend, len(audio_data))  # Times for simulation
    sol = solve_ivp(Lorentz, [0, tend], ic, t_eval=t, atol=1e-9, rtol=1e-8)
    u = sol.y[0, :]  # Extract u(t)
    plot_waveform(noise_fig, u, "Noise Signal")

def combine_signals():
    global u
    noise_size = np.sqrt((u**2).sum())
    music_size = np.sqrt((audio_data**2).sum())
    signal = u + audio_data*noise_size/music_size/1000
    plot_waveform(combined_fig, signal, "Combined Signal")


def decode_signal():
    global u
    ic = np.array([-2.3, -1.9, 2.2])  # Initial conditions
    tend = len(audio_data) / sample_rate  # Calculate the duration of the audio file
    t = np.linspace(0, tend, len(audio_data))  # Times for simulation

    ufun = CubicSpline(t, signal)

    def Lorentz2(tt, a):
        assert np.isscalar(tt)
        U = a[0]
        V = a[1]
        W = a[2]
        uu = ufun(np.array([tt]))[0]
        dotU = 16*(V-U)
        dotV = uu*(45.6-20*W)-V
        dotW = 5*uu*V - 4*W
        return np.array([dotU, dotV, dotW])

    sol = solve_ivp(Lorentz2, [0, tend], ic, t_eval=t, atol=1e-9, rtol=1e-8)
    U = sol.y[0, :]
    recovered_music = U - u
    plot_waveform(decoded_fig, recovered_music, "Decoded Signal")

def play_audio(data):
    # Function to play the audio data
    pass

def plot_waveform(fig, data, title):
    fig.clear()
    ax = fig.add_subplot(111)
    ax.plot(data)
    ax.set_title(title)
    fig.tight_layout()
    canvas = FigureCanvasTkAgg(fig, master=graph_frame)
    canvas.draw()
    canvas.get_tk_widget().pack()

root = tk.Tk()
root.title("Audio Processing GUI")

# Create frames for different sections of the GUI
control_frame = tk.Frame(root)
graph_frame = tk.Frame(root)

# Create buttons in the control frame
load_button = tk.Button(control_frame, text="Load Audio File", command=load_audio_file)
load_button.pack()

generate_noise_button = tk.Button(control_frame, text="Generate Noise", command=generate_noise)
generate_noise_button.pack()

combine_button = tk.Button(control_frame, text="Combine Signals", command=combine_signals)
combine_button.pack()

decode_button = tk.Button(control_frame, text="Decode Signal", command=decode_signal)
decode_button.pack()

# Create graph canvases in the graph frame
input_fig = plt.figure()
noise_fig = plt.figure()
combined_fig = plt.figure()
decoded_fig = plt.figure()

# Pack the frames
control_frame.pack(side=tk.LEFT)
graph_frame.pack(side=tk.RIGHT)

root.mainloop()
