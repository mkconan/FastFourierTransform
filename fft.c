#include <stdio.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265358979323846264338327950288
#define DATA_SIZE 16
#define CALC_REPEAT_NUM 100

// 複素数を表す構造体
typedef struct {
    double real;
    double imag;
} Complex;

// 複素数の掛け算
Complex complex_mul(Complex a, Complex b) {
    Complex result;
    result.real = a.real * b.real - a.imag * b.imag;
    result.imag = a.real * b.imag + a.imag * b.real;
    return result;
}

// 複素数の足し算
Complex complex_add(Complex a, Complex b) {
    Complex result;
    result.real = a.real + b.real;
    result.imag = a.imag + b.imag;
    return result;
}

// 複素数の引き算
Complex complex_sub(Complex a, Complex b) {
    Complex result;
    result.real = a.real - b.real;
    result.imag = a.imag - b.imag;
    return result;
}

// ビット反転順のインデックス計算
unsigned int reverse_bits(unsigned int n, unsigned int num_bits) {
    unsigned int reversed = 0;
    for (unsigned int i = 0; i < num_bits; ++i) {
        reversed |= ((n >> i) & 1) << (num_bits - 1 - i);
    }
    return reversed;
}

// 離散フーリエ変換（DFT）の実装
void dft(const Complex* input, Complex* output, int N) {
    for (int k = 0; k < N; ++k) {
        output[k].real = 0;
        output[k].imag = 0;
        for (int n = 0; n < N; ++n) {
            double angle = 2 * PI * k * n / N;
            double real_part = input[n].real * cos(angle) + input[n].imag * sin(angle);
            double imag_part = -input[n].real * sin(angle) + input[n].imag * cos(angle);
            output[k].real += real_part;
            output[k].imag += imag_part;
        }
    }
}

// 高速フーリエ変換（FFT）関数
void fft(Complex *x, int n) {
    if (n <= 1) return;

    // 偶数と奇数の部分配列を生成
    Complex even[n/2];
    Complex odd[n/2];
    for (int i = 0; i < n / 2; i++) {
        even[i] = x[i * 2];
        odd[i] = x[i * 2 + 1];
    }

    // 再帰的にFFTを計算
    fft(even, n / 2);
    fft(odd, n / 2);

    // FFTの結果を合成
    for (int k = 0; k < n / 2; k++) {
        double t = -2 * M_PI * k / n;
        Complex wk = {cos(t), sin(t)};
        Complex t1 = complex_mul(wk, odd[k]);
        x[k] = complex_add(even[k], t1);
        x[k + n / 2] = complex_sub(even[k], t1);
    }
}

void fft_butterfly(Complex* data, unsigned int n) {
    unsigned int num_bits = (unsigned int)log2(n);

    // データのビット反転順に並べ替え
    for (unsigned int i = 0; i < n; ++i) {
        unsigned int j = reverse_bits(i, num_bits);
        if (i < j) {
            Complex temp = data[i];
            data[i] = data[j];
            data[j] = temp;
        }
    }

    // FFTの計算
    for (unsigned int s = 1; s <= num_bits; ++s) {
        unsigned int m = 1 << s;
        double angle = -2.0 * PI / m;
        Complex wm = {cos(angle), sin(angle)};
        for (unsigned int k = 0; k < n; k += m) {
            Complex w = {1.0, 0.0};
            for (unsigned int j = 0; j < m / 2; ++j) {
                Complex t = complex_mul(w, data[k + j + m / 2]);
                Complex u = data[k + j];
                data[k + j] = complex_add(u, t);
                data[k + j + m / 2] = complex_sub(u, t);
                w = complex_mul(w, wm);
            }
        }
    }
    return;
}

// 2周期成分の確認および直流/交流の判定
void analyze_fft_result(Complex *x, int n) {
    double dc_component = sqrt(x[0].real * x[0].real + x[0].imag * x[0].imag);
    double ac_component_sum = 0.0;

    // 2周期成分の振幅を計算
    double two_cycle_component = sqrt(x[2].real * x[2].real + x[2].imag * x[2].imag) +
                                 sqrt(x[4].real * x[4].real + x[4].imag * x[4].imag) +
                                 sqrt(x[n-2].real * x[n-2].real + x[n-2].imag * x[n-2].imag) + 
                                 sqrt(x[n-4].real * x[n-4].real + x[n-4].imag * x[n-4].imag);

    for (int i = 1; i < n; i++) {
        ac_component_sum += sqrt(x[i].real * x[i].real + x[i].imag * x[i].imag);
    }

    printf("DC Component: %f\n", dc_component);
    printf("AC Component Sum: %f\n", ac_component_sum);
    printf("Two-cycle Component: %f\n", two_cycle_component);

    if (two_cycle_component > dc_component) {
        printf("The signal contains a significant 2-cycle component and is primarily AC.\n");
    } else {
        printf("The signal does not contain a significant 2-cycle component and may be primarily DC.\n");
    }
}

int main() {
    // 入力データ（16個の要素）
    double input[] = {0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0};
    // double input[] = {1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1};
    Complex output[DATA_SIZE] = {0};
    Complex x1[DATA_SIZE];
    for (int i = 0; i < DATA_SIZE; i++) {
        x1[i].real = input[i];
        x1[i].imag = 0.0;
    }
    Complex x2[DATA_SIZE];
    for (int i = 0; i < DATA_SIZE; i++) {
        x2[i].real = input[i];
        x2[i].imag = 0.0;
    }
    Complex x3[DATA_SIZE];
    for (int i = 0; i < DATA_SIZE; i++) {
        x3[i].real = input[i];
        x3[i].imag = 0.0;
    }

    clock_t start_time, end_time;

    // バタフライFFTの実行
    long int butterfly_fft_time = 0;
    for (int i = 0; i < CALC_REPEAT_NUM; i++){
        start_time = clock();
        fft_butterfly(x2, DATA_SIZE);
        end_time = clock();
        butterfly_fft_time += (end_time - start_time);
    }
    printf("FFT-butterfly time: %ld\n", butterfly_fft_time);

    // FFTの実行
    long int fft_time = 0;
    for (int i = 0; i < CALC_REPEAT_NUM; i++){
        start_time = clock();
        fft(x1, DATA_SIZE);
        end_time = clock();
        fft_time += (end_time - start_time);
    }
    printf("FFT time: %ld\n", fft_time);

    // DFTの実行
    long int dft_time = 0;
    for (int i = 0; i < CALC_REPEAT_NUM; i++){
        start_time = clock();
        dft(x3, output, DATA_SIZE);
        end_time = clock();
        dft_time += (end_time - start_time);
    }
    printf("DFT time: %ld\n", dft_time);

    // 直流成分と交流成分の判定
    // analyze_fft_result(x1, DATA_SIZE);

    /* printf("FT Result:\n");
    for (int i = 0; i < DATA_SIZE; i++) {
        printf("(%2.2f, %2.2f), (%2.2f, %2.2f), (%2.2f, %2.2f)\n", x1[i].real, x1[i].imag, x2[i].real, x2[i].imag, output[i].real, output[i].imag);
    } */

    return 0;
}