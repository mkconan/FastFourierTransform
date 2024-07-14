import numpy as np
import matplotlib.pyplot as plt
from pprint import pprint


# input_data = [0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0]
# input_data = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# input_data = [1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1]

input_data_list = []
for i in range(16):
    input_data = [0 for _ in range(16 - i)] + [1 for _ in range(i + 1)]
    input_data_list.append(input_data)
pprint(input_data_list)

input_data_list = [[0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1]]
input_data_list = [[0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1]]

for input_data in input_data_list:
    input_data = np.array(input_data, dtype=np.float32)

    plt.plot(input_data)
    plt.savefig("input_data.png")
    plt.close()

    T = 0.1
    N = input_data.size

    fft_out = np.fft.fft(input_data)
    freq = np.fft.fftfreq(N, d=T)
    print(2.0 / N * np.abs(fft_out[2]))
    plt.plot(freq[: N // 2], 2.0 / N * np.abs(fft_out[: N // 2]))
    plt.xticks(freq[: N // 2])
    # print(freq[: N // 2])
    plt.savefig("fft.png")
    # print(fft_out)
    # print(np.abs(fft_out))
    # print(np.sum(np.abs(fft_out)[1:]))
    plt.close()
