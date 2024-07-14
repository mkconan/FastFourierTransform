#include <stdio.h>
#include <math.h>

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
    Complex x[16];
    double input[] = {0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0};
    // double input[] = {1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1};
    for (int i = 0; i < 16; i++) {
        x[i].real = input[i];
        x[i].imag = 0.0;
    }

    // FFTの実行
    fft(x, 16);

    // 結果の表示
    printf("FFT Result:\n");
    for (int i = 0; i < 16; i++) {
        printf("(%2.2f, %2.2f)\n", x[i].real, x[i].imag);
    }

    // 直流成分と交流成分の判定
    analyze_fft_result(x, 16);

    return 0;
}