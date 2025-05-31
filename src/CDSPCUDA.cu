/*
 * CDSPCUDA.cpp
 *
 *  Created on: 01 янв. 2020 г.
 *      Author: user
 */

#include "CDSPCUDA.hpp"
#include <omp.h>

CDSP_CUDA::CDSP_CUDA(CLog *log) /*: CDSP(log)*/ {
    Log = log;
    //	cudaDeviceProp deviceProp;
    //	devID = gpuGetMaxGflopsDeviceId();
#ifdef QUASAR_CUDA
    devID = gpuGetMaxGflopsDeviceId();
    checkCudaErrors(cuDeviceGet(&cuDevice, devID));

    cudaGetDeviceProperties(&prop, 0);

#endif
}

CDSP_CUDA::CDSP_CUDA() {
    // TODO Auto-generated constructor stub
}

CDSP_CUDA::~CDSP_CUDA() {
    // TODO Auto-generated destructor stub
}

#ifdef QUASAR_CUDA

/*!
   \brief Вывод информации о версии библиотеки, частоте процессора, количестве
   ядер и возможностях процессора

        \return		 Функция не возвращает значения
*/

int CDSP_CUDA::printDSPInfo(void) {

    Log->LogPrintf((char *)"============\n");

    int n;
    cudaGetDeviceCount(&n);
    Log->LogPrintf("Количество вычислителей  CUDA: %d\n", n);

    char name[100];
    cuDeviceGetName(name, 100, cuDevice);
    Log->LogPrintf("Вычислитель \t\t CUDA [%d]: %s\n", devID, name);

    Log->LogPrintf("Частота вычислителя\t %u кГц\n", prop.memoryClockRate);
    Log->LogPrintf("Количество потоков\t %u\n", prop.maxThreadsPerBlock);
    Log->LogPrintf("maxThreadsDim\t %ux%ux%u\n", prop.maxThreadsDim[0],
                   prop.maxThreadsDim[1], prop.maxThreadsDim[2]);
    Log->LogPrintf("maxGridSize\t %ux%ux%u\n", prop.maxGridSize[0],
                   prop.maxGridSize[1], prop.maxGridSize[2]);

    return 0;
}

//============================
// fft()
//============================
int CDSP_CUDA::fft(Ipp32f *pImageData, int N, int M) {

    /*
        int order = log2(N);
        assert(N == pow(2,order)); // N должно быть степенью двойки

        int 		  batch = M;
        float 		 *devInputData;
        cufftComplex *devOutputData;
        checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&devInputData),
       N*batch*sizeof(float))); checkCudaErrors(cudaMalloc(reinterpret_cast<void
       **>(&devOutputData), (N/2+1)*batch*sizeof(cufftComplex)));

        checkCudaErrors(cudaMemcpy(devInputData, pImageData,
       N*batch*sizeof(float), cudaMemcpyHostToDevice));

        cufftHandle handle;
        int rank = 1;                        // --- 1D FFTs
        int n[] = { N };                	 // --- Size of the Fourier
       transform int istride = 1, ostride = 1;        // --- Distance between
       two successive input/output elements int idist = N, odist = (N / 2 + 1);
       // --- Distance between batches int inembed[] = { 0 };               //
       --- Input size with pitch (ignored for 1D transforms) int onembed[] = { 0
       };               // --- Output size with pitch (ignored for 1D
       transforms)

        checkCudaErrors(cufftPlanMany(&handle, rank, n, inembed, istride, idist,
       onembed, ostride, odist, CUFFT_R2C, batch));
        checkCudaErrors(cufftExecR2C(handle, devInputData, devOutputData));

        #pragma omp parallel for
        for(int j=0; j<M; j++){
                float *p = &pImageData[j*N];
                cufftComplex *dp = &devOutputData[j*(N/2+1)];
                checkCudaErrors(cudaMemcpy(p, dp, (N/2)*sizeof(cufftComplex),
       cudaMemcpyDeviceToHost));
        }

        checkCudaErrors(cufftDestroy(handle));
        checkCudaErrors(cudaFree(devInputData));
        checkCudaErrors(cudaFree(devOutputData));
        */

    int order = log2(N);
    assert(N == pow(2, order)); // N должно быть степенью двойки

    fftwf_plan FFTPlan;

    if (!FFTPlan) {
        FFTPlan = fftwf_plan_dft_r2c_1d(N, NULL, NULL, FFTW_ESTIMATE);
    }

#pragma omp parallel for
    for (int j = 0; j < M; j++) {
        float *p = &pImageData[j * N];
        p[0] = 0.0f;
        p[1] = 0.0f;
        fftwf_execute_dft_r2c(FFTPlan, p, (fftwf_complex *)p);
        p[0] = 0;
        p[1] = 0; // Убираем постоянную составляющую
    }

    return -1;
}

int CDSP_CUDA::fft(Ipp32fc *pImageData, int N, int M) {
    int order = log2(N);
    assert(N == pow(2, order)); // N должно быть степенью двойки

    int nbthreads = omp_get_max_threads();
    fftwf_init_threads();
    fftwf_plan_with_nthreads(nbthreads);
    omp_set_num_threads(nbthreads);

    fftwf_plan p;
    int rank = 1;
    int n[] = {(int)N};

#pragma omp critical
    p = fftwf_plan_many_dft(rank, n, M, (fftwf_complex *)pImageData, NULL, 1, N,
                            (fftwf_complex *)pImageData, NULL, 1, N,
                            FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(p);

#pragma omp critical
    fftwf_destroy_plan(p);

    /*	//#pragma omp parallel for
            for(int j=0; j<M; j++){
                    fftwf_plan p1;
                    fftwf_complex *p = (fftwf_complex*)(&pImageData[j*N]);

                    //#pragma omp critical (make_plan)
                    p1 = fftwf_plan_dft_1d(N,  p, (fftwf_complex*)p,
       FFTW_FORWARD, FFTW_ESTIMATE );

                    fftwf_execute(p1);
                    fftwf_destroy_plan(p1);

                    fftwf_cleanup_threads();
            }
    */
    return -1;
}

int CDSP_CUDA::ifft(Ipp32fc *pImageData, int N, int M) {
    int order = log2(N);
    assert(N == pow(2, order)); // N должно быть степенью двойки

    int nbthreads = omp_get_max_threads();
    fftwf_init_threads();
    fftwf_plan_with_nthreads(nbthreads);
    omp_set_num_threads(nbthreads);

    fftwf_plan p;
    int rank = 1;
    int n[] = {(int)N};

#pragma omp critical
    p = fftwf_plan_many_dft(rank, n, M, (fftwf_complex *)pImageData, NULL, 1, N,
                            (fftwf_complex *)pImageData, NULL, 1, N,
                            FFTW_BACKWARD, FFTW_ESTIMATE);
    fftwf_execute(p);

#pragma omp critical
    fftwf_destroy_plan(p);

    // #pragma omp parallel for
    //	for(int j=0; j<M; j++){
    //		fftwf_plan p1;
    //		fftwf_complex *p = (fftwf_complex*)(&pImageData[j*N]);

    // #pragma omp critical (make_plan)
    //		p1 = fftwf_plan_dft_1d(N,  p, (fftwf_complex*)p, FFTW_BACKWARD,
    // FFTW_ESTIMATE );

    //		fftwf_execute(p1);
    //		fftwf_destroy_plan(p1);
    //	}

    return -1;
}

__device__ inline void complexMul_(Ipp32fc *Dst, const Ipp32fc Src) {
    float re = 0;
    float im = 0;
    re = Dst[0].re * Src.re - Dst[0].im * Src.im;
    im = Dst[0].re * Src.im + Dst[0].im * Src.re;
    Dst[0].re = re;
    Dst[0].im = im;
}

__device__ inline void complexMul_(Ipp32fc *SrcDst1, const float Src2) {

    SrcDst1[0].re = SrcDst1[0].re * Src2;
    SrcDst1[0].im = SrcDst1[0].im * Src2;
}

__device__ inline void complexSum_(Ipp32fc *Dst, const Ipp32fc Src) {
    Dst[0].re += Src.re;
    Dst[0].im += Src.im;
}

__device__ void bpPixel(Ipp32fc *indata, float *outdata, float *V, float x,
                        float y, float *h, bool bCorr, float dt, float c1,
                        float c2, float c5, float dfr, float dr, int Np, int Ns,
                        float *devWindow) {

    float t = 0;
    float r = 0;
    int ind = 0;
    float phase = 0;

    float fp = 0;
    float fr = 0;
    float phi = 0;
    float c4 = 0;

    Ipp32fc Sum;
    Sum.re = 0;
    Sum.im = 0;

    Ipp32fc signal;
    Ipp32fc sop;
    Ipp32fc sk;

    for (int k = 0; k < Np; k++) {
        c4 = h[k] * h[k] + x * x;
        t = k * dt;
        r = sqrtf(c4 + pow((y + V[k] * t), 2));
        ind = round(r / dr);

        signal = indata[k * Ns + ind];
        phase = c1 * r;

        sop.re = cosf(phase);
        sop.im = sinf(phase);

        if (bCorr) {
            fp = (ind - 1) * dfr;

            fr = r * c5;
            phi = c2 * (fr - fp);
            sk.re = cosf(phi);
            sk.im = sinf(phi);

            complexMul_(&signal, sk);
        }

        complexMul_(&signal, devWindow[k]); // Оконное взвешивание
        complexMul_(&signal, sop);
        complexSum_(&Sum, signal);
    }

    *outdata = sqrtf(Sum.im * Sum.im + Sum.re * Sum.re);
}

__global__ void bp_(Ipp32fc *indata, float *outdata, float *V, float x0,
                    float y0, float dx, float dy, float *h, int nx, int ny,
                    bool bCorr, float dt, float c1, float c2, float c5,
                    float dfr, float dr, int Np, int Ns, float *devWindow) {
    int ind_i = blockIdx.y * blockDim.y + threadIdx.y;
    int stride_y = blockDim.y * gridDim.y;

    int ind_j = blockIdx.x * blockDim.x + threadIdx.x;
    int stride_x = blockDim.x * gridDim.x;

    float out = 0;

    for (int i = ind_i; i < ny; i += stride_y) {

        for (int j = ind_j; j < nx; j += stride_x) {
            float x = x0 + j * dx;
            float y = y0 + i * dy;
            bpPixel(indata, &out, V, x, y, h, bCorr, dt, c1, c2, c5, dfr, dr,
                    Np, Ns, devWindow);
            outdata[i * nx + j] = out;
        }
    }
}

#define BLOCK_DIVISION 16

int CDSP_CUDA::bp(CImage *image, CRadar *radar, float *outdata, float *window,
                  float *V, float x0, float y0, float dx, float dy, float *h,
                  int nx, int ny, bool bCorr, Ipp32fc *pc) {

    // float *window = (float*)malloc(image->ImageH0*sizeof(float));
    // HammingWindow(window, image->ImageH0);

    int Np = image->ImageH0;
    int Ns = image->ImageWidth;
    float Tp = radar->GetTp();
    float Mu = radar->GetMu();
    float F0 = radar->F0;
    float Fs = radar->Fs;

    float dt = Tp;
    float c1 = -4 * M_PI / SPEED_OF_LIGHT * F0;
    float c2 = -M_PI * Tp;
    float c5 = 2 * Mu / SPEED_OF_LIGHT;
    float dfr = Fs / (2 * (Ns - 1));
    float dr = Fs * SPEED_OF_LIGHT / (4 * Mu * (Ns - 1));

    Ipp32fc *indata = (Ipp32fc *)image->pImageData;

    Ipp32fc *devInData;
    float *devOutData;
    float *devWindow;
    float *devV;
    float *devH;

    cudaError_t cudaStatus;

    // In
    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&devInData),
                               Ns * Np * sizeof(Ipp32fc)));
    checkCudaErrors(cudaMemcpy(devInData, indata, Ns * Np * sizeof(Ipp32fc),
                               cudaMemcpyHostToDevice));
    // Out
    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&devOutData),
                               nx * ny * sizeof(float)));
    checkCudaErrors(cudaMemcpy(devOutData, outdata, nx * ny * sizeof(float),
                               cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemset(devOutData, 0, nx * ny * sizeof(float)));
    // Window
    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&devWindow),
                               image->ImageH0 * sizeof(float)));
    checkCudaErrors(cudaMemcpy(devWindow, window,
                               image->ImageH0 * sizeof(float),
                               cudaMemcpyHostToDevice));
    // V H
    checkCudaErrors(
        cudaMalloc(reinterpret_cast<void **>(&devV), Np * sizeof(float)));
    checkCudaErrors(
        cudaMemcpy(devV, V, Np * sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(
        cudaMalloc(reinterpret_cast<void **>(&devH), Np * sizeof(float)));
    checkCudaErrors(
        cudaMemcpy(devH, h, Np * sizeof(float), cudaMemcpyHostToDevice));

    dim3 blocks((nx + BLOCK_DIVISION - 1) / BLOCK_DIVISION,
                (ny + BLOCK_DIVISION - 1) / BLOCK_DIVISION);
    dim3 threads(BLOCK_DIVISION, BLOCK_DIVISION);

    bp_<<<blocks, threads>>>(devInData, devOutData, devV, x0, y0, dx, dy, devH,
                             nx, ny, bCorr, dt, c1, c2, c5, dfr, dr, Np, Ns,
                             devWindow);

    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "Kernel launch failed: %s\n",
                cudaGetErrorString(cudaStatus));
    }

    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr,
                "cudaDeviceSynchronize returned error code %d after launching "
                "bpKernel!\n",
                cudaStatus);
    }

    checkCudaErrors(cudaMemcpy(outdata, devOutData, nx * ny * sizeof(float),
                               cudaMemcpyDeviceToHost));

    checkCudaErrors(cudaFree(devInData));
    checkCudaErrors(cudaFree(devOutData));
    checkCudaErrors(cudaFree(devWindow));
    checkCudaErrors(cudaFree(devV));
    checkCudaErrors(cudaFree(devH));

    // free(window);

    return 0;
}

__global__ void bp_strip_(Ipp32fc *indata, float *outdata, float *V, float x0,
                          float y0, float dx, float *h, int nx, bool bCorr,
                          float dt, float c1, float c2, float c5, float dfr,
                          float dr, int Np, int Ns, float *devWindow,
                          float angle) {

    int ind_j = blockIdx.x * blockDim.x + threadIdx.x;
    int stride_x = blockDim.x * gridDim.x;

    float out = 0;

    float cosAngle = cos(angle);
    float tanAngle = tan(angle);

    for (int j = ind_j; j < nx; j += stride_x) {
        float x = (x0 + j * dx) * cosAngle;
        float y = x * tanAngle;
        bpPixel(indata, &out, V, x, y, h, bCorr, dt, c1, c2, c5, dfr, dr, Np,
                Ns, devWindow);
        outdata[j] = out;
    }
}

int CDSP_CUDA::bp_strip(CImage *image, CImage *out, CRadar *radar,
                        float *window, float *V, float x0, float y0, float dx,
                        float *h, int nx, int line, bool bCorr, Ipp32fc *pc,
                        float angle) {

    int Np = image->ImageH0;
    int Ns = image->ImageWidth;
    float Tp = radar->GetTp();
    float Mu = radar->GetMu();
    float F0 = radar->F0;
    float Fs = radar->Fs;

    float dt = Tp;
    float c1 = -4 * M_PI / SPEED_OF_LIGHT * F0;
    float c2 = -M_PI * Tp;
    float c5 = 2 * Mu / SPEED_OF_LIGHT;
    float dfr = Fs / (2 * (Ns - 1));
    float dr = Fs * SPEED_OF_LIGHT / (4 * Mu * (Ns - 1));

    Ipp32fc *indata = (Ipp32fc *)image->pImageData;
    Ipp32f *outdata = (Ipp32f *)out->pImageData + nx * line;

    Ipp32fc *devInData;
    float *devOutData;
    float *devWindow;
    float *devV;
    float *devH;

    cudaError_t cudaStatus;

    // In
    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&devInData),
                               Ns * Np * sizeof(Ipp32fc)));
    checkCudaErrors(cudaMemcpy(devInData, indata, Ns * Np * sizeof(Ipp32fc),
                               cudaMemcpyHostToDevice));
    // Out
    checkCudaErrors(
        cudaMalloc(reinterpret_cast<void **>(&devOutData), nx * sizeof(float)));
    // checkCudaErrors(cudaMemcpy(devOutData, outdata, nx*ny*sizeof(float),
    // cudaMemcpyHostToDevice)); checkCudaErrors(cudaMemset(devOutData, 0,
    // nx*ny*sizeof(float)));
    // Window
    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&devWindow),
                               image->ImageH0 * sizeof(float)));
    checkCudaErrors(cudaMemcpy(devWindow, window,
                               image->ImageH0 * sizeof(float),
                               cudaMemcpyHostToDevice));
    // V H
    checkCudaErrors(
        cudaMalloc(reinterpret_cast<void **>(&devV), Np * sizeof(float)));
    checkCudaErrors(
        cudaMemcpy(devV, V, Np * sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(
        cudaMalloc(reinterpret_cast<void **>(&devH), Np * sizeof(float)));
    checkCudaErrors(
        cudaMemcpy(devH, h, Np * sizeof(float), cudaMemcpyHostToDevice));

    dim3 blocks((nx + BLOCK_DIVISION - 1) / BLOCK_DIVISION);
    dim3 threads(BLOCK_DIVISION, BLOCK_DIVISION);
    bp_strip_<<<blocks, threads>>>(devInData, devOutData, devV, x0, y0, dx,
                                   devH, nx, bCorr, dt, c1, c2, c5, dfr, dr, Np,
                                   Ns, devWindow, angle);

    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "Kernel launch failed: %s\n",
                cudaGetErrorString(cudaStatus));
    }

    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr,
                "cudaDeviceSynchronize returned error code %d after launching "
                "bpKernel!\n",
                cudaStatus);
    }

    checkCudaErrors(cudaMemcpy(outdata, devOutData, nx * sizeof(float),
                               cudaMemcpyDeviceToHost));

    checkCudaErrors(cudaFree(devInData));
    checkCudaErrors(cudaFree(devOutData));
    checkCudaErrors(cudaFree(devWindow));
    checkCudaErrors(cudaFree(devV));
    checkCudaErrors(cudaFree(devH));

    return 0;
}

__device__ void edgeRP(float x, float phi, int *ret) {

    float r = x / cos(phi);
    float y = r * sin(phi);

    *ret = (int)(y + 0.5);
}

__global__ void bpm_(Ipp32fc *indata, float *outdata, float *V, float x0,
                     float y0, float dx, float dy, float *h, int nx, int ny,
                     bool bCorr, float dt, float c1, float c2, float c5,
                     float dfr, float dr, int Np, int Ns, float *devWindow,
                     float angle, float div) {

    int ind_i = blockIdx.y * blockDim.y + threadIdx.y;
    int stride_y = blockDim.y * gridDim.y;

    int ind_j = blockIdx.x * blockDim.x + threadIdx.x;
    int stride_x = blockDim.x * gridDim.x;

    float out = 0;

    float div2 = div / 2;

    float cosAngle = cos(angle);
    float sinAngle = sin(angle);

    for (int j = ind_j; j < nx; j += stride_x) {

        float r = x0 + j * dx;

        int b;
        int e;
        edgeRP(r, -div2, &b);
        edgeRP(r, div2, &e);

        b = ((float)b - y0) / dy;
        e = ((float)e - y0) / dy;

        if (b < 0) {
            b = 0;
        }
        if (e > ny) {
            e = ny;
        }

        if (ind_i > e || ind_i < b) {
            continue;
        }

        for (int i = ind_i; i < ny; i += stride_y) {
            float x = (x0 + j * dx) * cosAngle - (y0 + i * dy) * sinAngle;
            float y = (x0 + j * dx) * sinAngle + (y0 + i * dy) * cosAngle;
            bpPixel(indata, &out, V, x, y, h, bCorr, dt, c1, c2, c5, dfr, dr,
                    Np, Ns, devWindow);
            outdata[i * nx + j] = out;
        }
    }
}

int CDSP_CUDA::bpm(CImage *image, CRadar *radar, float *outdata, float *wnd,
                   float *V, float x0, float y0, float dx, float dy, float *h,
                   int nx, int ny, bool bCorr, float angle, float div,
                   Ipp32fc *pc) {

    float *window = (float *)malloc(image->ImageH0 * sizeof(float));
    HammingWindow(window, image->ImageH0);

    int Np = image->ImageH0;
    int Ns = image->ImageWidth;
    float Tp = radar->GetTp();
    float Mu = radar->GetMu();
    float F0 = radar->F0;
    float Fs = radar->Fs;

    float dt = Tp;
    float c1 = -4 * M_PI / SPEED_OF_LIGHT * F0;
    float c2 = -M_PI * Tp;
    float c5 = 2 * Mu / SPEED_OF_LIGHT;
    float dfr = Fs / (2 * (Ns - 1));
    float dr = Fs * SPEED_OF_LIGHT / (4 * Mu * (Ns - 1));

    Ipp32fc *indata = (Ipp32fc *)image->pImageData;

    Ipp32fc *devInData;
    float *devOutData;
    float *devWindow;
    float *devV;
    float *devH;

    cudaError_t cudaStatus;

    // In
    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&devInData),
                               Ns * Np * sizeof(Ipp32fc)));
    checkCudaErrors(cudaMemcpy(devInData, indata, Ns * Np * sizeof(Ipp32fc),
                               cudaMemcpyHostToDevice));
    // Out
    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&devOutData),
                               nx * ny * sizeof(float)));
    checkCudaErrors(cudaMemcpy(devOutData, outdata, nx * ny * sizeof(float),
                               cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemset(devOutData, 0, nx * ny * sizeof(float)));
    // Window
    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&devWindow),
                               image->ImageH0 * sizeof(float)));
    checkCudaErrors(cudaMemcpy(devWindow, window,
                               image->ImageH0 * sizeof(float),
                               cudaMemcpyHostToDevice));
    // V H
    checkCudaErrors(
        cudaMalloc(reinterpret_cast<void **>(&devV), Np * sizeof(float)));
    checkCudaErrors(
        cudaMemcpy(devV, V, Np * sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(
        cudaMalloc(reinterpret_cast<void **>(&devH), Np * sizeof(float)));
    checkCudaErrors(
        cudaMemcpy(devH, h, Np * sizeof(float), cudaMemcpyHostToDevice));

    dim3 blocks((nx + BLOCK_DIVISION - 1) / BLOCK_DIVISION,
                (ny + BLOCK_DIVISION - 1) / BLOCK_DIVISION);
    dim3 threads(BLOCK_DIVISION, BLOCK_DIVISION);
    bpm_<<<blocks, threads>>>(devInData, devOutData, devV, x0, y0, dx, dy, devH,
                              nx, ny, bCorr, dt, c1, c2, c5, dfr, dr, Np, Ns,
                              devWindow, angle, div);

    /*
    int blockSize;
int minGridSize;
int gridSize;



    cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
bpm_, 0, nx*ny); gridSize = (nx*ny + blockSize - 1) / blockSize;

    bpm_<<<gridSize, blockSize>>>(devInData, devOutData, devV, x0, y0, dx, dy,
devH, nx, ny, bCorr, dt, c1, c2, c5, dfr, dr, Np, Ns, devWindow, angle, div);
    */

    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "Kernel launch failed: %s\n",
                cudaGetErrorString(cudaStatus));
    }

    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr,
                "cudaDeviceSynchronize returned error code %d after launching "
                "bpKernel!\n",
                cudaStatus);
    }

    checkCudaErrors(cudaMemcpy(outdata, devOutData, nx * ny * sizeof(float),
                               cudaMemcpyDeviceToHost));

    checkCudaErrors(cudaFree(devInData));
    checkCudaErrors(cudaFree(devOutData));
    checkCudaErrors(cudaFree(devWindow));
    checkCudaErrors(cudaFree(devV));
    checkCudaErrors(cudaFree(devH));

    free(window);

    return 0;
}

#endif
