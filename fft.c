#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979323846264338327950288

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
    Complex x1[16];
    for (int i = 0; i < 16; i++) {
        x1[i].real = input[i];
        x1[i].imag = 0.0;
    }
    Complex x2[16];
    for (int i = 0; i < 16; i++) {
        x2[i].real = input[i];
        x2[i].imag = 0.0;
    }

    // FFTの実行
    fft(x1, 16);

    // 結果の表示
    printf("FFT Result:\n");

    // 直流成分と交流成分の判定
    // analyze_fft_result(x1, 16);

    fft_butterfly(x2, 16);

    printf("FFT-butterfly Result:\n");
    for (int i = 0; i < 16; i++) {
        printf("(%2.2f, %2.2f), (%2.2f, %2.2f)\n", x1[i].real, x1[i].imag, x2[i].real, x2[i].imag);
    }

    return 0;
}