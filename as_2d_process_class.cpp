#include "as_2d_process_class.h"
#include <cmath>
#include <QDebug>
#include <fftw3.h>
#include <complex.h>
#include <utils.h>

AS_2D_process_class::AS_2D_process_class(int masksize, double sf, double sc):
    masksize(masksize), sf(sf), sc(sc) {
    // define the meshgrids
    // the real size is 2*masksize+1
    QVector<double> oneline_x;
    QVector<double> oneline_y;
    for (int i = -masksize; i <= masksize; ++i) {
        oneline_x.clear();
        oneline_y.clear();
        for (int j = -masksize; j <= masksize; ++j) {
            oneline_x.push_back(static_cast<double>(i));
            oneline_y.push_back(static_cast<double>(j));
        }
        this->x_mesh.push_back(oneline_x);
        this->y_mesh.push_back(oneline_y);
    }
}

AS_2D_process_class::~AS_2D_process_class(){

}

void AS_2D_process_class::fx_2DHilbertKernel(double s,
                                             QVector<QVector<double>>& kernel_1st,
                                             QVector<QVector<double>>& kernel_2nd) {
    double ss = pow(s, 2);
    double kk;
    double d;
    QVector<double> kernel_1st_oneline;
    QVector<double> kernel_2nd_oneline;
    for (int i = -this->masksize; i <= this->masksize; ++i) {
        kernel_1st_oneline.clear();
        kernel_2nd_oneline.clear();
        for (int j = -this->masksize; j <= this->masksize; ++j) {
            kk = pow(i, 2) + pow(j, 2);
            // 1st kernel
            kernel_1st_oneline.push_back(1 / (2*M_PI*pow(ss+kk, 1.5)));
            // 2nd kernel
            d = pow(kk, 2) * pow(ss+kk, 1.5) * 2 * M_PI;
            kernel_2nd_oneline.push_back(d==0? 0:-(s * (2*ss + 3*kk) - 2*pow(ss+kk, 1.5)) / d);
        }
        kernel_1st.push_back(kernel_1st_oneline);
        kernel_2nd.push_back(kernel_2nd_oneline);
    }
}

QVector<QVector<std::complex<double>>> AS_2D_process_class::
applyFFT2D(QVector<QVector<double>>& data)
{
    QVector<QVector<std::complex<double>>> output(data.size(), QVector<std::complex<double>>(data[0].size()));
    fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * data.size() * data[0].size());
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * data.size() * data[0].size());
    fftw_plan plan = fftw_plan_dft_2d(data.size(), data[0].size(), in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    // copy input data to fftw_complex array
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[0].size(); j++) {
            in[i*data[0].size()+j][0] = data[i][j];
            in[i*data[0].size()+j][1] = 0.0;
        }
    }
    fftw_execute(plan);
    // copy output data to QVector<QVector<std::complex<double>>> format
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[0].size(); j++) {
            output[i][j] = std::complex<double>(out[i*data[0].size()+j][0], out[i*data[0].size()+j][1]);
        }
    }
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
    return output;
}

QVector<QVector<double>> AS_2D_process_class::
applyIFFT2D(QVector<QVector<std::complex<double>>>& data)
{
    QVector<QVector<double>> output(data.size(), QVector<double>(data[0].size()));
    fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * data.size() * data[0].size());
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * data.size() * data[0].size());
    fftw_plan plan = fftw_plan_dft_2d(data.size(), data[0].size(), in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    // copy input data to fftw_complex array
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[0].size(); j++) {
            in[i*data[0].size()+j][0] = data[i][j].real();
            in[i*data[0].size()+j][1] = data[i][j].imag();
        }
    }
    fftw_execute(plan);
    // copy output data to QVector<QVector<double>> format
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[0].size(); j++) {
            output[i][j] = out[i*data[0].size()+j][0] / (data.size()*data[0].size());
        }
    }
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
    return output;
}

// Compute the 2D convolution of an image and a kernel in the frequency domain
QVector<QVector<double>> AS_2D_process_class::
convolve2DFreq(const QVector<QVector<double>>& image, const QVector<QVector<double>>& kernel)
{
    int numRows = image.size();
    int numCols = image[0].size();
    int kernelNumRows = kernel.size();
    int kernelNumCols = kernel[0].size();
    // Pad the input image and the kernel to avoid circular convolution
    int paddedNumRows = numRows + kernelNumRows - 1;
    int paddedNumCols = numCols + kernelNumCols - 1;
    QVector<QVector<double>> paddedImage(paddedNumRows, QVector<double>(paddedNumCols, 0.0));
    QVector<QVector<double>> paddedKernel(paddedNumRows, QVector<double>(paddedNumCols, 0.0));
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            paddedImage[i][j] = image[i][j];
        }
    }
    for (int i = 0; i < kernelNumRows; i++) {
        for (int j = 0; j < kernelNumCols; j++) {
            paddedKernel[i][j] = kernel[i][j];
        }
    }
    // Compute the IDFT of the product
    QVector<QVector<std::complex<double>>> Image_fft = this->applyFFT2D(paddedImage);
    QVector<QVector<std::complex<double>>> Kernel_fft = this->applyFFT2D(paddedKernel);
    // Multiply the frequency-domain representations of the image and the kernel
    QVector<QVector<std::complex<double>>> fftResult(paddedNumRows,
                                                     QVector<std::complex<double>>(paddedNumCols, 0.0));
    for (int i = 0; i < paddedNumRows; i++) {
        for (int j = 0; j < paddedNumCols; j++) {
            fftResult[i][j] = Image_fft[i][j] * Kernel_fft[i][j];
        }
    }
    QVector<QVector<double>> Result = this->applyIFFT2D(fftResult);
    // Crop the result to remove the padding
    QVector<QVector<double>> croppedResult(numRows, QVector<double>(numCols, 0.0));
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            croppedResult[i][j] = Result[i + kernelNumRows / 2][j + kernelNumCols / 2];
        }
    }
    return croppedResult;
}

void AS_2D_process_class::processData(const QVector<QVector<double>>& C_scan_AS) {
    // Ensure the data is not empty
    Q_ASSERT(!C_scan_AS.isEmpty());

    int x_size = C_scan_AS.size();
    int y_size = C_scan_AS[0].size();

    this->img_Cscan = C_scan_AS;

    // 2d analytic-signals
    // calculate 2D analytic-signal
    QVector<QVector<double>> img_ori(x_size, QVector<double>(y_size, 0));
    QVector<QVector<double>> img_pha(x_size, QVector<double>(y_size, 0));
    QVector<QVector<double>> img_amp(x_size, QVector<double>(y_size, 0));
    QVector<QVector<double>> img_ape(x_size, QVector<double>(y_size, 0));
    // *************** freq. domain convolution
    QVector<QVector<double>> kernel_1st_f;
    QVector<QVector<double>> kernel_2nd_f;
    QVector<QVector<double>> kernel_1st_c;
    QVector<QVector<double>> kernel_2nd_c;
    //
    this->fx_2DHilbertKernel(this->sf, kernel_1st_f, kernel_2nd_f);
    this->fx_2DHilbertKernel(this->sc, kernel_1st_c, kernel_2nd_c);
    //
    QVector<QVector<double>> fp_kernel = kernel_1st_f * (this->sf) - kernel_1st_c * (this->sc);
    QVector<QVector<double>> q1_x_kernel = (kernel_1st_f - kernel_1st_c) * this->x_mesh;
    QVector<QVector<double>> q1_y_kernel = (kernel_1st_f - kernel_1st_c) * this->y_mesh;
    QVector<QVector<double>> q2_xx_kernel = (kernel_2nd_f - kernel_2nd_c) * this->x_mesh * this->x_mesh;
    QVector<QVector<double>> q2_xy_kernel = (kernel_2nd_f - kernel_2nd_c) * this->x_mesh * this->y_mesh;
    QVector<QVector<double>> q2_yy_kernel = (kernel_2nd_f - kernel_2nd_c) * this->y_mesh * this->y_mesh;
    //
    QVector<QVector<double>> fp = this->convolve2DFreq(this->img_Cscan, fp_kernel);
    QVector<QVector<double>> fx = this->convolve2DFreq(this->img_Cscan, q1_x_kernel);
    QVector<QVector<double>> fy = this->convolve2DFreq(this->img_Cscan, q1_y_kernel);
    QVector<QVector<double>> fxx = this->convolve2DFreq(this->img_Cscan, q2_xx_kernel);
    QVector<QVector<double>> fxy = this->convolve2DFreq(this->img_Cscan, q2_xy_kernel);
    QVector<QVector<double>> fyy = this->convolve2DFreq(this->img_Cscan, q2_yy_kernel);
    //
    QVector<QVector<double>> f_pm = (fxx-fyy) * 0.5;
    QVector<QVector<double>> f_s = fp * 0.5;
    for (int i = 0; i<x_size; ++i){
        for (int j = 0; j<y_size; ++j){
            double e = sqrt(pow(f_pm[i][j], 2) + pow(fxy[i][j], 2)) / abs(f_s[i][j]);
            double l = pow(fx[i][j], 2) + pow(fy[i][j], 2);
            double q = l*2/(1+e);
            img_pha[i][j] = atan2(sqrt(q),fp[i][j]);
            img_ori[i][j] = 0.5 * atan2(fxy[i][j], f_pm[i][j]);
            //    img_ori[i][j] = img_pha[i][j]==0?
            //         0.5 * atan2(fxy[i][j], f_pm[i][j])+M_PI/2:atan2(fy[i][j], fx[i][j]);
            img_amp[i][j] = 0.5 * sqrt(pow(fp[i][j], 2) + q);
            img_ape[i][j] = atan2(sqrt(pow(f_s[i][j],2)-pow(fxy[i][j],2)-pow(f_pm[i][j],2)),
                                  sqrt(pow(fxy[i][j],2)+pow(f_pm[i][j],2)));
            //            if (std::real(sqrt(pow(f_s[i][j],2)-pow(fxy[i][j],2)-pow(f_pm[i][j],2)))<=0)
            //                qDebug() << "negative nominator" << pow(f_s[i][j],2)-pow(fxy[i][j],2)-pow(f_pm[i][j],2);
        }
    }
    qDebug() << "Data processed. 2D AS image are generated";
}
