import tkinter as tk
from tkinter import filedialog
import wavio
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp
import matplotlib
import simpleaudio as sa
import librosa

file = '496127__timbre__loopable-airy-remix-of-erokias-freesound-496096.wav'

sample_wavio_obj = wavio.read(file)
sample_wav_obj = sa.WaveObject.from_wave_file(file)
audio_data = sample_wavio_obj.data
sample_rate = sample_wavio_obj.rate

SUPPORTED_SAMPLE_RATES = [8000, 11025, 22050, 44100, 48000]

# Existing Lorentz function from the original code
def Lorentz(t, a):
    u = a[0]  # unpack x
    v = a[1]  # unpack y
    w = a[2]  # unpack z
    dotu = 16*(v-u)
    dotv = u*(45.6-20*w)-v
    dotw = 5*u*v - 4*w
    return np.array([dotu, dotv, dotw])  # return vector of derivatives

# Load audio file
def load_audio_file():
    global audio_data, sample_rate
    file_path = filedialog.askopenfilename(title="Load Audio File")
    if file_path:
        try:
            sample_wavio_obj = wavio.read("496127__timbre__loopable-airy-remix-of-erokias-freesound-496096.wav")
            audio_data = sample_wavio_obj.data
            sample_rate = sample_wavio_obj.rate
            # Call functions to process audio data and update waveform display
            plot_waveform(signal_canvas, audio_data, "Input Signal")
        except Exception as e:
            print(f"Error loading audio file: {e}")


# Generate noise
def generate_noise():
    global u
    global noise1_wav_obj
    ic = np.array([-2.3, -1.9, 2.2])  # Initial conditions
    tend = len(audio_data) / sample_rate  # Calculate the duration of the audio file
    t = np.linspace(0, tend, len(audio_data))  # Times for simulation
    sol = solve_ivp(Lorentz, [0, tend], ic, t_eval=t, atol=1e-9, rtol=1e-8)
    u = sol.y[0, :]  # Extract u(t)
    plot_waveform(noise_canvas, u, "Noise Signal")
    #save the noisewavfrorm to a wav file
    wavio.write("noise.wav", u, sample_rate, sampwidth=3)

    #create a simpleaudio object for the noise
    noise1_wav_obj = sa.WaveObject.from_wave_file("noise.wav")

# Generate noise
def generate_noise2():
    global u
    global noise1_wav_obj2
    ic = np.array([-2.3, -1.9, 2.2])  # Initial conditions
    tend = len(audio_data) / sample_rate  # Calculate the duration of the audio file
    t = np.linspace(0, tend, len(audio_data))  # Times for simulation
    sol = solve_ivp(Lorentz, [0, tend], ic, t_eval=t, atol=1e-9, rtol=1e-8)
    u = sol.y[0, :]  # Extract u(t)
    plot_waveform(noise_canvas, u, "Noise Signal")
    #save the noisewavfrorm to a wav file
    wavio.write("noise2.wav", u, sample_rate, sampwidth=3)

    #create a simpleaudio object for the noise
    noise1_wav_obj2 = sa.WaveObject.from_wave_file("noise2.wav")


# Combine signal and noise
def combine_data():
    global audio_data, u, combined_data, combined_wav_obj
    
    # Check if audio data and noise data are available
    if audio_data is None or u is None:
        print("Audio data or noise data is missing.")
        return
    
    # Check if the lengths of audio data and noise data are the same
    if len(audio_data) != len(u):
        print("Audio data and noise data have different lengths.")
        return
    
    # Sum the elements along the second axis of audio_data
    audio_data_sum = np.sum(audio_data, axis=1)
    
    # Combine the audio data and noise data
    combined_data = audio_data_sum + u
    
    # Plot the combined waveform
    plot_waveform(combined_canvas, combined_data, "Combined Signal")

    # Save the combined waveform to a WAV file
    wavio.write("combined.wav", combined_data, sample_rate, sampwidth=3)

    # Create a simpleaudio object for the combined signal
    combined_wav_obj = sa.WaveObject.from_wave_file("combined.wav")

# Decode audio
def decode_audio():
    global combined_data, decoded_data, decoded_wav_obj
    
    # Check if combined data is available
    if combined_data is None:
        print("Combined data is missing.")
        return
    
    # Decode the combined data
    decoded_data = combined_data - u
    
    # Plot the decoded waveform
    plot_waveform(decoded_canvas, decoded_data, "Decoded Signal")
    
    # Save the decoded waveform to a WAV file
    wavio.write("decoded.wav", decoded_data, sample_rate, sampwidth=3)
    
    # Create a simpleaudio object for the decoded signal
    decoded_wav_obj = sa.WaveObject.from_wave_file("decoded.wav")


# Play audio data
def play_audio(audio_data, sample_rate):
    # Convert audio data to simpleaudio compatible format
    audio_data = audio_data.astype(np.float32)
    
    # Check if the sample rate is supported
    if sample_rate not in SUPPORTED_SAMPLE_RATES:
        print(f"Sample rate {sample_rate} Hz is not supported. Resampling to 44100 Hz.")
        sample_rate = 44100  # Resample to a supported rate (e.g., 44100 Hz)
        audio_data = librosa.core.resample(audio_data, orig_sr=sample_rate, target_sr=44100)
    
    # Create an instance of the simpleaudio player
    audio_player = sa.play_buffer(audio_data, 1, 2, sample_rate)
    
    # Wait for the audio to finish playing
    audio_player.wait_done()

def play_signal(wav_obj):
    play_obj = wav_obj.play()

    # filename = 'myfile.wav'
    # wave_obj = sa.WaveObject.from_wave_file(filename)
    # play_obj = wave_obj.play()
    # play_obj.wait_done()  # Wait until sound has finished playing

# Plot waveform data
def plot_waveform(canvas, audio_data, title):
    # Flatten the audio_data array if it has more than one dimension
    audio_data = audio_data.flatten()

    # Clear the canvas
    canvas.delete("all")

    # Add title to the canvas
    canvas.create_text(canvas.winfo_width()//2, 10, text=title)

    # Plot the waveform data on the canvas
    canvas_width = canvas.winfo_width()
    canvas_height = canvas.winfo_height()
    samples = len(audio_data)

    if samples > 0:
        x_scale = canvas_width / samples
        y_scale = canvas_height / 2 / np.max(np.abs(audio_data))

        canvas.create_line(0, canvas_height//2, canvas_width, canvas_height//2, fill="gray")

        x_start = 0
        y_start = int(canvas_height//2 - audio_data[0] * y_scale)

        for i in range(1, samples):
            x_end = int(i * x_scale)
            y_end = int(canvas_height//2 - audio_data[i] * y_scale)
            canvas.create_line(x_start, y_start, x_end, y_end, fill="blue")
            x_start = x_end
            y_start = y_end


root = tk.Tk()
root.title("Audio Processing GUI")

# Create frames
button_frame = tk.Frame(root)
quadrant_frame = tk.Frame(root)

# Create buttons
load_button = tk.Button(button_frame, text="Load Audio File", command=load_audio_file)
noise_button = tk.Button(button_frame, text="Generate Noise", command=generate_noise)
combine_button = tk.Button(button_frame, text="Combine Data", command=combine_data)
decode_button = tk.Button(button_frame, text="Decode Audio", command=decode_audio)

load_button.pack(side=tk.LEFT)
noise_button.pack(side=tk.LEFT)
combine_button.pack(side=tk.LEFT)
decode_button.pack(side=tk.LEFT)

# Create quadrants
signal_canvas = tk.Canvas(quadrant_frame, width=400, height=200, bg="white")
noise_canvas = tk.Canvas(quadrant_frame, width=400, height=200, bg="white")
combined_canvas = tk.Canvas(quadrant_frame, width=400, height=200, bg="white")
decoded_canvas = tk.Canvas(quadrant_frame, width=400, height=200, bg="white")

signal_canvas.grid(row=0, column=0, padx=10, pady=10)
noise_canvas.grid(row=0, column=1, padx=10, pady=10)
combined_canvas.grid(row=1, column=0, padx=10, pady=10)
decoded_canvas.grid(row=1, column=1, padx=10, pady=10)

# Create play buttons
signal_play_button = tk.Button(quadrant_frame, text="Play Signal", command=lambda: play_signal(sample_wav_obj))
noise_play_button = tk.Button(quadrant_frame, text="Play Noise", command=lambda: play_signal(noise1_wav_obj))
combined_play_button = tk.Button(quadrant_frame, text="Play Combined", command=lambda: play_signal(combined_wav_obj))
decoded_play_button = tk.Button(quadrant_frame, text="Play Decoded", command=lambda: play_signal(decoded_wav_obj))

signal_play_button.grid(row=2, column=0, padx=10, pady=10)
noise_play_button.grid(row=2, column=1, padx=10, pady=10)
combined_play_button.grid(row=3, column=0, padx=10, pady=10)
decoded_play_button.grid(row=3, column=1, padx=10, pady=10)

# Pack frames
button_frame.pack(side=tk.TOP, fill=tk.X)
quadrant_frame.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

root.mainloop()