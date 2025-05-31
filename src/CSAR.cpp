/*
 * CSAR.cpp
 *
 * Created on: 3 янв. 2017 г.
 * Author: User
 */

#include "CSAR.hpp"
#include <tgmath.h>

CSAR::CSAR() {
    Log = NULL;
    Radar = NULL;

    dspType = DSP_IPP;
}

CSAR::CSAR(CLog *log, CRadar *radar, DSPType dsp_type)
    : Log(log), Radar(radar) {
    dspType = dsp_type;
    DSP = NULL;

    switch (dsp_type) {
    case DSP_IPP:
#ifdef QUASAR_IPP
        DSP = new CDSP_IPP(Log);
#endif
        break;
    case DSP_FFTW:
        DSP = new CDSP_FFTW(Log);
        break;
    case DSP_NEON:
#ifdef QUASAR_NEON
        DSP = new CDSP_NEON(Log);
#endif
        break;
    case DSP_CUDA:
#ifdef QUASAR_CUDA
        DSP = new CDSP_CUDA(Log);
#endif
        break;
    case DSP_AF:
#ifdef SUPPORT_AF
        DSP = new CDSP_AF(Log);
#endif
        break;

    case DSP_OCL_CPU:
    case DSP_OCL_GPU:
#ifdef SUPPORT_OCL
        DSP = new CDSP_OCL(Log, dsp_type);
#endif
        break;
    }

    assert(DSP != NULL);
}

CSAR::~CSAR() { delete (DSP); }

/*!
   \brief Вывод информации о параметрах системы ЦОС

    \return      Функция не возвращает значения
*/
void CSAR::printDSPInfo(void) { DSP->printDSPInfo(); }

/*
uint64_t CSAR::GetDSP(){
    return (uint64_t)DSP;
}
*/
/*!
   \brief           Удаление постоянной составляющей из сигнала

   \param[in]  image        Указатель на экземпляр CImage

    \return      0 - успешное выполнение;
                -1 - ошибка
*/
int CSAR::removeDC(CImage *image) {
    assert(image->ImageType == CImage::IT_FCOMPLEX);
    // assert(image->ImageHeight == 1);

    Log->LogPrintf("\t*Компенсация DC...\t\t");
    fflush(stdout);
    time_start();

    Ipp32fc *data = (Ipp32fc *)image->pImageData;
    int W0 = image->ImageW0;
    int W = image->ImageWidth;
    int H = image->ImageH0;

    DSP->removeDC(data, W0, W, H);
    /*
        #pragma omp parallel for
        for(int h=0; h<H; h++){
            float mean = DSP->mean(&data[h*image->ImageWidth], W);
            for(int w=0; w<W; w++){
                data[h*image->ImageWidth + w] -= mean;
            }
        }
    */
    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());
    return 0;
}

/*!
   \brief       Сжатие по дальности

   \param[in]  image        Указатель на экземпляр CImage

    \return      0 - успешное выполнение;
                -1 - ошибка
*/
int CSAR::rangeCompress(CImage *image) {
    Log->LogPrintf("\t*Сжатие по дальности...\t\t");
    fflush(stdout);
    time_start();

    int status = 0;
    Ipp32f *data;
    int W = image->ImageWidth;

    switch (dspType) {
    case DSP_IPP:
    case DSP_OCL_CPU:
    case DSP_OCL_GPU:
    case DSP_FFTW:
    case DSP_NEON:
    case DSP_CUDA:
        data = (Ipp32f *)image->pImageData;

        // fftwf_init_threads();
        // fftwf_plan_with_nthreads(omp_get_max_threads());
        status = DSP->fft(data, W,
                          image->ImageH0); // На выходе - односторонний спектр

        image->SetImageWidth(W / 2);
        image->SetImageType(CImage::IT_FCOMPLEX);
        break;
    case DSP_AF:
#ifdef SUPPORT_AF
        DSP->af_fft_r2c(image->afImageData); // На выходе - односторонний спектр
                                             // af::deviceGC();
#endif
        break;
    default:
        printError(Log, 0, (char *)"Неизвестный тип вычислителя",
                   (char *)__FILE__, __LINE__);
        break;
    }

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());

    return status;
}

/*!
   \brief       Транспонирование изображения

   \param[in]  image        Указатель на экземпляр CImage

    \return      0 - успешное выполнение;
                -1 - ошибка
*/
int CSAR::transpose(CImage *image) {
    Log->LogPrintf("\t*Транспонирование...\t\t");
    fflush(stdout);
    time_start();

    int status = 0;
    if ((dspType == DSP_IPP) || (dspType == DSP_OCL_CPU) ||
        (dspType == DSP_OCL_GPU) || (dspType == DSP_FFTW) ||
        (dspType == DSP_NEON)) {
        Ipp32fc *data = (Ipp32fc *)image->pImageData;
        int H = image->ImageHeight;
        int W = image->ImageWidth;

        if (image->ImageType == CImage::IT_FCOMPLEX) {
            status = DSP->transpose(data, W, H);
        } else if (image->ImageType == CImage::IT_FLOAT) {
            status = DSP->transpose((Ipp32f *)data, W, H);
        }

        image->SetImageWidth(H);
        image->SetImageHeight(W);
    } else {
#ifdef SUPPORT_AF
        status = DSP->af_transpose(image->afImageData);
        // af::deviceGC();
#endif
    }

    int tmp = image->ImageW0;
    image->ImageW0 = image->ImageH0;
    image->ImageH0 = tmp;

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());

    return status;
}

/*!
   \brief       Коррекция смещения отметок по дальности вследствие дробного
   отношения количества отсчетов АЦП за период модуляции

   \param[in]   image       Указатель на экземпляр CImage

    \return      0 - успешное выполнение;
                -1 - ошибка
*/
int CSAR::compensateShift(CImage *image) {
    Log->LogPrintf("\t*Компенсация смещения...\t");
    fflush(stdout);
    time_start();

    if ((dspType == DSP_IPP) || (dspType == DSP_OCL_CPU) ||
        (dspType == DSP_OCL_GPU)) {
        Ipp32fc *data = (Ipp32fc *)image->pImageData;
        int H = image->ImageHeight;
        int W = image->ImageWidth;

        DSP->compensateShift(data, W, H, Radar->Mp,
                             image->GetPixelSize(image->ImageType));
    } else {
#ifdef SUPPORT_AF
        DSP->af_compensateShift(image->afImageData, Radar->Mp);
        af::deviceGC();
#endif
    }

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());

    return true;
}

/*!
   \brief       Сжатие по азимуту

   \param[in]  image        Указатель на экземпляр CImage

    \return      0 - успешное выполнение;
                -1 - ошибка
*/
int CSAR::azimuthCompress(CImage *image) {
    // Сжатие по азимуту
    Log->LogPrintf("\t*Сжатие по азимуту...\t\t");
    fflush(stdout);
    time_start();

    // assert(image->ImageType == CImage::IT_FCOMPLEX);

    if ((dspType == DSP_IPP) || (dspType == DSP_OCL_CPU) ||
        (dspType == DSP_OCL_GPU) || (dspType == DSP_FFTW) ||
        (dspType == DSP_NEON)) {
        Ipp32fc *data = (Ipp32fc *)image->pImageData;
        int H = image->ImageHeight;
        int W = image->ImageWidth;

        if (image->ImageType == CImage::IT_FCOMPLEX) {
            DSP->fft(data, W, H);
            DSP->fftshift(data, W, H);
        } else if (image->ImageType == CImage::IT_FLOAT) {
            DSP->fft((Ipp32f *)data, W, H);
            image->SetImageWidth(W / 2);
            image->SetImageType(CImage::IT_FCOMPLEX);
        }
    } else {
#ifdef SUPPORT_AF
        DSP->af_fft(image->afImageData);
        DSP->af_fftshift(image->afImageData);
        // af::deviceGC();
#endif
    }

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());

    return ippStsNoErr;
}

/*!
   \brief       Коррекция миграции по дальности

   \param[in]  image        Указатель на экземпляр CImage
   \param[in]  V        Скорость носителя, м/с

    \return      0 - успешное выполнение;
                -1 - ошибка
*/
int CSAR::rmc(CImage *image, float V) {
    Log->LogPrintf("\t*Коррекция миграции...\t\t");
    fflush(stdout);
    time_start();

    int Hb = image->ImageNearRangeIndex;
    float lambda = Radar->GetLambda();
    float Tp = (float)Radar->Ns / (float)Radar->Fs;
    int status = 0;

    if ((dspType == DSP_IPP) || (dspType == DSP_OCL_CPU) ||
        (dspType == DSP_OCL_GPU) || (dspType == DSP_FFTW) ||
        (dspType == DSP_NEON)) {
        Ipp32fc *data = (Ipp32fc *)image->pImageData;
        int H = image->ImageHeight;
        int W = image->ImageWidth;
        status = DSP->rmc(data, W, H, Hb, lambda, Tp, V);
    } else {
#ifdef SUPPORT_AF
        DSP->af_rmc(image->afImageData, Hb, lambda, Tp, V);
        // af::deviceGC();
#endif
    }

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());
    return status;
}

/*!
   \brief       Коррекция миграции по дальности методом частотного
   масштабирования

   \param[in]  image        Указатель на экземпляр CImage
   \param[in]  V        Скорость носителя, м/с

    \return      0 - успешное выполнение;
                -1 - ошибка
*/
int CSAR::fsrmc(CImage *image, float V) {
    Log->LogPrintf("\t*Коррекция миграции...\t\t");
    fflush(stdout);
    time_start();

    assert(image->ImageType == CImage::IT_FCOMPLEX);

    //  int Hb = image->ImageNearRangeIndex;
    //  float lambda = Radar->GetLambda();
    //  float Tp = (float)Radar->Ns/(float)Radar->Fs;
    int status = 0;
    /*  int Ns = image->ImageWidth;
        int Np = image->ImageHeight;
        Ipp32fc* data = (Ipp32fc*)image->pImageData;
        float a = Radar->GetMu()/8.0f;
        float dt = Radar->GetTp()/image->ImageW0;

        Ipp32fc *s1 = (Ipp32fc*)malloc(Ns*2*sizeof(float));
        Ipp32fc *s2 = (Ipp32fc*)malloc(Ns*2*sizeof(float));
        Ipp32fc *s3  = (Ipp32fc*)malloc(Ns*2*sizeof(float));
        Ipp32fc *s4  = (Ipp32fc*)malloc(Ns*2*sizeof(float));

        float betta = 1.0;

        //printf("s %d %d\n", order, Ns);
        #pragma omp parallel for
        for (int p=0; p<Np; p++){
            for (int i=0; i<Ns; i++){
                float t = dt*i - Radar->GetTp()/2;
                float t2 = t*t;
                float c = betta*a;
                float b = c + a;
                float phase1 = a*t2;
                float phase2 = -b*t2;
                float phase3 = c*t2;
                float phase4 = -b*c/a*t2;
                s1[i].re = cos(phase1);
                s1[i].im = sin(phase1);
                s2[i].re = cos(phase2);
                s2[i].im = sin(phase2);
                s3[i].re = cos(phase3);
                s3[i].im = sin(phase3);
                s4[i].re = cos(phase4);
                s4[i].im = sin(phase4);
            }

            //DSP->fastconv(&data[p*Ns], s1, Ns);
            ippsMul_32fc_I(s2, &data[p*Ns], Ns);
            DSP->fastconv(&data[p*Ns], s3, Ns);
            ippsMul_32fc_I(s4, &data[p*Ns], Ns);
        }

        free(s4);
        free(s3);
        free(s2);
        free(s1);
    */
    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());
    return status;
}

/*!
   \brief       Фокусировка изображения

   \param[in]  image    Указатель на экземпляр CImage
   \param[in]  V        Скорость носителя, м/с

    \return      0 - успешное выполнение;
                -1 - ошибка
*/
int CSAR::phocus(CImage *image, float V) {
    // Фокусировка
    Log->LogPrintf("\t*Фокусировка...\t\t\t");
    fflush(stdout);
    time_start();

    float Rmax = fr2range(Radar->Fs / 2);
    float lambda = Radar->GetLambda();
    float Tp = Radar->GetTp();
    float Rb = index2range(image->ImageNearRangeIndex);
    if ((dspType == DSP_IPP) || (dspType == DSP_OCL_CPU) ||
        (dspType == DSP_OCL_GPU) || (dspType == DSP_FFTW) ||
        (dspType == DSP_NEON)) {
        assert(image->ImageType == CImage::IT_FCOMPLEX);
        Ipp32fc *data = (Ipp32fc *)image->pImageData;
        int H = image->ImageHeight;
        int W = image->ImageWidth;

        DSP->fftshift(data, W, H);
        DSP->ifft(data, W, H);

        DSP->phocus(data, W, H, Radar->Ns2 / 2, Rb, Tp, lambda, V, Rmax,
                    image->ImageW0);

        DSP->fft(data, W, H);
        DSP->fftshift(data, W, H);
    } else {
#ifdef SUPPORT_AF
        DSP->af_fftshift(image->afImageData);
        DSP->af_ifft(image->afImageData);

        DSP->af_phocus(image->afImageData, Radar->Ns2 / 2, Rb, Tp, lambda, V,
                       Rmax, image->ImageW0);

        DSP->af_fft(image->afImageData);
        DSP->af_fftshift(image->afImageData);
        // af::deviceGC();
#endif
    }

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());

    return ippStsNoErr;
}

/*!
   \brief       Вычисление энтропии изображения

   \param[in]  image    Указатель на экземпляр CImage

    \return      Значение энтропии изображения
*/
float CSAR::enthropy(CImage *image) {
    Log->LogPrintf("\t*Вычисление энтропии...\t\t");
    fflush(stdout);
    time_start();

    float E = 0;

    if (dspType != DSP_AF) {
        Ipp32fc *data = (Ipp32fc *)image->pImageData;
        int H = image->ImageHeight;
        int W = image->ImageWidth;

        if (image->ImageType == CImage::IT_FCOMPLEX) {
            E = DSP->entropy(data, W, H);
        } else if (image->ImageType == CImage::IT_FLOAT) {
            E = DSP->entropyf((Ipp32f *)data, W, H);
        } else
            E = 0;

    } else {
#ifdef SUPPORT_AF
        E = DSP->af_entropy(image->afImageData);
        // af::deviceGC();
#endif
    }

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());

    return E;
}

/*!
   \brief       Оконное взвешивание (Хэмминга)

   \param[in]  image    Указатель на экземпляр CImage

    Осуществляет построчное оконное взвешивание отсчетов изображения.

    \return      0 - успешное выполнение;
                -1 - ошибка
*/
int CSAR::hamming(CImage *image) {
    Log->LogPrintf("\t*Оконное взвешивание...\t\t");
    fflush(stdout);
    time_start();

    switch (dspType) {
    case DSP_IPP:
    case DSP_OCL_CPU:
    case DSP_OCL_GPU:
    case DSP_FFTW:
    case DSP_NEON:
    case DSP_CUDA:
        if (image->ImageType == CImage::IT_FLOAT)
            DSP->hamming((Ipp32f *)image->pImageData, image->ImageWidth,
                         image->ImageHeight, image->ImageW0,
                         image->ImageHeight);
        if (image->ImageType == CImage::IT_FCOMPLEX)
            DSP->hamming((Ipp32fc *)image->pImageData, image->ImageWidth,
                         image->ImageHeight, image->ImageW0,
                         image->ImageHeight);
        break;
    case DSP_AF:
#ifdef SUPPORT_AF
        DSP->af_hamming(image->afImageData, image->ImageW0);
        // af::deviceGC();
#endif
        break;
    default:
        printError(Log, 0, (char *)"Неизвестный тип вычислителя",
                   (char *)__FILE__, __LINE__);
        break;
    }

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());
    //  Log->LogPrintf("ОK. (%2.3f c)\n", timer::stop(start1));

    return 0;
}

/*!
   \brief       Оконное взвешивание (Хэмминга)

   \param[in]  image    Указатель на экземпляр CImage

    Осуществляет построчное оконное взвешивание отсчетов изображения.

    \return      0 - успешное выполнение;
                -1 - ошибка
*/
int CSAR::windowWeighting(CImage *image, uint64_t window_addr) {
    Log->LogPrintf("\t*Оконное взвешивание...\t\t");
    fflush(stdout);
    time_start();

    float *window = (float *)window_addr;

    // printf("pHammingWindow %p\n", window);

    switch (dspType) {
    case DSP_IPP:
    case DSP_OCL_CPU:
    case DSP_OCL_GPU:
    case DSP_FFTW:
    case DSP_NEON:
    case DSP_CUDA:
        if (image->ImageType == CImage::IT_FLOAT) {
            DSP->windowWeighting((Ipp32f *)image->pImageData, window,
                                 image->ImageWidth, image->ImageW0,
                                 image->ImageH0);
        } else if (image->ImageType == CImage::IT_FCOMPLEX) {
            DSP->windowWeighting((Ipp32fc *)image->pImageData, window,
                                 image->ImageWidth, image->ImageW0,
                                 image->ImageHeight);
        }
        break;
    case DSP_AF:
#ifdef SUPPORT_AF
        DSP->af_hamming(image->afImageData, image->ImageW0);
        // af::deviceGC();
#endif
        break;
    default:
        printError(Log, 0, (char *)"Неизвестный тип вычислителя",
                   (char *)__FILE__, __LINE__);
        break;
    }

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());
    //  Log->LogPrintf("ОK. (%2.3f c)\n", timer::stop(start1));

    return 0;
}

/*!
   \brief       Предварительное оконное взвешивание (Хэмминга)

   \param[in]  image    Указатель на экземпляр CImage

    Осуществляет построчное оконное взвешивание отсчетов изображения.

    \return      0 - успешное выполнение;
                -1 - ошибка
*/
uint64_t CSAR::createWindow(int size, WindowType w_type) {

    float *Window = NULL;
    switch (w_type) {
    case CSAR::HAMMING:
        Window = (float *)malloc(size * sizeof(float));
        DSP->HammingWindow(Window, size);
        break;
    case CSAR::BLACKMAN:
        Window = (float *)malloc(size * sizeof(float));
        DSP->BlackmanWindow(Window, size);
        break;
    case CSAR::NUTTALL:
        Window = (float *)malloc(size * sizeof(float));
        DSP->NuttallWindow(Window, size);
        break;
    case CSAR::NONE:
        Window = (float *)malloc(size * sizeof(float));
        for (int i = 0; i < size; i++) {
            Window[i] = 1;
        }
        break;
    }

    if (Window) {
        return ((uint64_t)Window);
    }

    return -1;
}

void CSAR::freeWindow(uint64_t window_addr) {
    float *window = (float *)window_addr;
    free(window);
}

/*!
   \brief       Частотная интерполяция

   \param[in]  image    Указатель на экземпляр CImage
   \param[in]  kr       Коэффициент частотной интерполяции (степень двойки)

    Осуществляется дополнение каждой строки нулями

    \return      0 - успешное выполнение;
                -1 - ошибка
*/
int CSAR::interpolate(CImage *image, int kr) {
    Log->LogPrintf("\t*Частотная интерполяция...\t");
    fflush(stdout);
    time_start();

    kr = pow(2, kr);

    if ((dspType == DSP_IPP) || (dspType == DSP_OCL_CPU) ||
        (dspType == DSP_OCL_GPU) || (dspType == DSP_FFTW) ||
        (dspType == DSP_NEON)) {
        int H = image->ImageHeight;
        int W = image->ImageWidth;
        int PS = image->GetPixelSize(image->ImageType);
        void *data = image->pImageData;

        void *newdata = aligned_alloc(64, H * W * kr * PS);
        memset(newdata, 0, H * W * kr * PS);

        // #pragma omp parallel for
        for (int i = 0; i < H; i++) {
            memcpy(&((char *)newdata)[i * W * PS * kr],
                   &((char *)data)[i * W * PS], W * PS);
        }

        free(image->pImageData);
        image->pImageData = newdata;
        image->SetImageWidth(W * kr);
    } else {
#ifdef SUPPORT_AF
        //  printf("%d %d", image->afImageData.dims(0),
        //  image->afImageData.dims(1));
        int W = image->afImageData.dims(0);
        int H = image->afImageData.dims(1);
        array Y = array(W * kr, H, c32);
        Y(seq(0, W - 1), span) = image->afImageData;
        image->afImageData = Y;
        // af::deviceGC();
#endif
    }
    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());
    return 0;
}

/*!
    \brief      Преобразование Гильберта

    \param[in]  image   Указатель на экземпляр CImage

    \return      0 - успешное выполнение;
                -1 - ошибка
*/
int CSAR::hilbert(CImage *image) {
    Log->LogPrintf("\t*Преобразование Гильберта...\t");
    fflush(stdout);
    time_start();
    assert(image->ImageType == CImage::IT_FLOAT);
    int H = image->ImageHeight;
    int W = image->ImageWidth;
    Ipp32f *pindata = (Ipp32f *)image->pImageData;

    Ipp32fc *poutdata = (Ipp32fc *)malloc(W * H * sizeof(Ipp32fc));
    DSP->hilbert(pindata, poutdata, W, H);

    free(image->pImageData);
    image->pImageData = poutdata;
    image->SetImageType(CImage::IT_FCOMPLEX);

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());
    return 0;
}

#include <signal.h> //  our new library

volatile sig_atomic_t flag = 0;
void my_function(int sig) { // can be called asynchronously
    flag = 1;               // set flag
}

int CSAR::dspProgress() {

    signal(SIGINT, my_function);

    while (1) {
        if (flag) { // my action when signal set it 1
            printf("\n Signal caught!\n");
            printf("\n default action it not termination!\n");
            flag = 0;
        }
    }

    return DSP->dspProgress();
}

int CSAR::acc_dcm(CImage *image, float Vmax, int Np) {
    Log->LogPrintf("\t*Формирование ДКМ...\t\t");
    fflush(stdout);
    time_start();

    assert(image->ImageType == CImage::IT_FCOMPLEX);

    Ipp32f *data = (Ipp32f *)image->pImageData;
    int W = image->ImageWidth;

    float mu = Radar->GetMu();
    float mue = 4 * mu * Vmax / SPEED_OF_LIGHT;
    float dt = 1 / Radar->Fs;

    int ns = image->ImageW0;
    Ipp32fc *Sop =
        (Ipp32fc *)malloc(ns * sizeof(Ipp32fc)); // Буфер под опорную функцию
    Ipp32fc *im_buf =
        (Ipp32fc *)malloc(W * Np * sizeof(Ipp32fc)); // Буфер под изображение
    memset(im_buf, 0, W * Np * sizeof(Ipp32fc));

    float dmu = 2.0f * mue / (float)(Np - 1);

#pragma omp parallel for
    for (int p = 0; p < Np; p++) {
        float mur = -mue + p * dmu;

        // Формируем опорную функцию
        for (int i = 0; i < ns; i++) {
            float t = dt * i;
            float arg = -M_PI * mur / 2.0f * pow(t, 2);
            Sop[i].re = cos(arg);
            Sop[i].im = sin(arg);
        }
        memcpy((void *)&im_buf[p * W], data, W * sizeof(Ipp32fc));
        DSP->complexMul(&im_buf[p * W], Sop, ns);
        //      im_buf[p*W].re = 0;
        //      im_buf[p*W].im = 0;
        //      im_buf[p*W+1].re = 0;
        //      im_buf[p*W+1].im = 0;
        //      im_buf[p*W+2].re = 0;
        //      im_buf[p*W+2].im = 0;
    }
    free(Sop);

    DSP->fft(im_buf, W, Np);
    // DSP->fftshift(im_buf, W, Np);

    free(image->pImageData);
    image->pImageData = im_buf;
    image->ImageH0 = Np;
    image->SetImageHeight(Np);

    image->CutW(0, W / 2); // Односторонний спектр

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());
    return 0;
}

/*!
   \brief       Пересчет значения разностной частоты к наклонной дальности

   \param[in]  fr       Значение дальностной частоты, Гц

    \return      Дальнось, м
*/
float CSAR::fr2range(float fr) {
    return fr * SPEED_OF_LIGHT / (2 * Radar->GetMu());
}

/*!
   \brief       Пересчет наклонной дальности к значению разностной частоты

   \param[in]  r Дальнось, м

    \return      Значение дальностной частоты, Гц
*/
float CSAR::range2fr(float r) {
    return 2 * Radar->GetMu() * r / SPEED_OF_LIGHT;
}

/*!
   \brief       // Пересчет значения индекса изображения к наклонной дальности

   \param[in]  ind      Значение индекса на изображении

    \return      Дальнось, м
*/
float CSAR::index2range(int ind) {
    float Rmax = fr2range(Radar->Fs / 2);
    return Rmax / (float)(Radar->Ns2 / 2) * (float)ind;
}

/*!
   \brief       // Пересчет наклонной дальности к индексу изображения

   \param[in]  ind      Значение наклонной дальности, м

    \return      Индекс
*/
// int CSAR::range2index(float range){
//   float Rmax = beat2range(Radar->Fs/2);
//   return  round(range*(float)(Radar->Ns2/2)/Rmax);
// }

/*!
   \brief       // Вычисление наклонной дальности по координатам

   \param[in]  x    Наклонная дальность на траверзе, м
   \param[in]  y    Путевая дальность, м
   \param[in]  h    Высота полета носителя, м
   \param[in]  h    Скорость полета носителя, м
   \param[in]  t    Момент времени, с

    \return      Наклонная дальность
*/
// float CSAR::get_range(float x, float y, float h, float V, float t){
//   float c1 = y + V*t;
//   return sqrt(x * x + c1 * c1 + h * h);
// }

/*!
   \brief       // Билинейное масштабирование РЛИ

   \param[in]   image   Указатель на экземпляр CImage
   \param[in]   newW    Желаемая ширина изображения, пикс.
   \param[in]   newH    Желаемая высота изображения, пикс.

    \return     0 - успешное выполнение;
                -1 - ошибка
*/
int CSAR::resize(CImage *image, int newW, int newH) {

    Log->LogPrintf("\t*Масштабирование РЛИ...\t\t");
    fflush(stdout);
    time_start();

    int oldW = image->ImageWidth;
    int oldH = image->ImageHeight;

    Ipp32f *in_matrix = NULL;
    Ipp32f *out_matrix = NULL;

    if (image->ImageType == CImage::IT_FCOMPLEX) {
        in_matrix = (Ipp32f *)malloc(oldW * oldH * sizeof(Ipp32f));
        memset(in_matrix, 0, oldW * oldH * sizeof(Ipp32f));

        image->magnitude((Ipp32fc *)image->pImageData, in_matrix,
                         image->ImageWidth, image->ImageHeight);
        free((Ipp32fc *)image->pImageData);

        image->ImageType = CImage::IT_FLOAT;
    } else {
        in_matrix = (Ipp32f *)image->pImageData;
    }

    out_matrix = (Ipp32f *)malloc(newH * newW * sizeof(Ipp32f));

#pragma omp parallel for
    for (int i = 0; i < newH; i++) {

        int l;
        int c;
        float t;
        float u;
        float tmp;
        float d1, d2, d3, d4;
        Ipp32f p1, p2, p3, p4;

        for (int j = 0; j < newW; j++) {

            tmp = (float)(i) / (float)(newH - 1) * (oldH - 1);
            l = (int)floor(tmp);
            if (l < 0) {
                l = 0;
            } else {
                if (l >= oldH - 1) {
                    l = oldH - 2;
                }
            }

            u = tmp - l;
            tmp = (float)(j) / (float)(newW - 1) * (oldW - 1);
            c = (int)floor(tmp);
            if (c < 0) {
                c = 0;
            } else {
                if (c >= oldW - 1) {
                    c = oldW - 2;
                }
            }
            t = tmp - c;

            // коэффициенты
            d1 = (1 - t) * (1 - u);
            d2 = t * (1 - u);
            d3 = t * u;
            d4 = (1 - t) * u;

            // pixels
            p1 = *((Ipp32f *)in_matrix + (l * oldW) + c);
            p2 = *((Ipp32f *)in_matrix + (l * oldW) + c + 1);
            p3 = *((Ipp32f *)in_matrix + ((l + 1) * oldW) + c + 1);
            p4 = *((Ipp32f *)in_matrix + ((l + 1) * oldW) + c);

            out_matrix[(i * newW) + j] = p1 * d1 + p2 * d2 + p3 * d3 + p4 * d4;
        }
    }

    free(in_matrix);
    image->pImageData = out_matrix;
    image->ImageWidth = newW;
    image->ImageHeight = newH;

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());

    return 0;
}

int get_extra_width(float *data, float r, float step, int width) {

    int i = 0;
    int size = 0;

    while (i < width) {

        while (r > data[i]) {
            i++;
        }
        r += step;
        size++;
    }

    return size;
}

#define getK(rmin, rmaj, r) ((rmaj - rmin) / (r - rmin))
int findK(float *data, int *index_list, float *k_list, float r, float step,
          int width) {

    int i = 0;
    int size = 0;

    while (1) {

        while (r > data[i]) {
            i++;

            if (i >= width) {
                return size;
            }
        }

        k_list[size] = getK(data[i - 1], data[i], r);
        index_list[size] = i;

        r += step;
        size++;
    }
}

/*!
   \brief       // Преобразование наклонной дальности в горизонтальную (RDA)

   \param[in]   image   Указатель на экземпляр CImage
   \param[in]   decim   Коэффициент децимации
   \param[in]   x0      Начальная координата по оси x
   \param[in]   lx      Дальность изображения
   \param[in]   H       Высота носителя

    \return     0 - успешное выполнение;
                -1 - ошибка
*/
#define getA(amin, amaj, k) ((amaj - amin) / k) + amin
int CSAR::downrange(CImage *image, float decim, float x0, float lx, float H) {

    Log->LogPrintf("\t*Преобразование\n\t наклонной дальности...\t\t");
    fflush(stdout);
    time_start();

    int width = image->ImageWidth;
    int height = image->ImageHeight;

    Ipp32f *in_matrix = NULL;
    Ipp32f *out_matrix = NULL;

    if (image->ImageType == CImage::IT_FCOMPLEX) {
        in_matrix = (Ipp32f *)malloc(width * height * sizeof(Ipp32f));
        memset(in_matrix, 0, width * height * sizeof(Ipp32f));

        image->magnitude((Ipp32fc *)image->pImageData, in_matrix,
                         image->ImageWidth, image->ImageHeight);
        free((Ipp32fc *)image->pImageData);

        image->ImageType = CImage::IT_FLOAT;
    } else {
        in_matrix = (Ipp32f *)image->pImageData;
    }

    float *rangeData = (float *)malloc(width * sizeof(float));
    float step = lx / width;

#pragma omp parallel for
    for (int i = 0; i < width; i++) {
        float R = (x0 + (step * i));
        rangeData[i] = sqrt(R * R - H * H);
    }

    x0 = sqrt(x0 * x0 - H * H);

    int extra_width = get_extra_width(rangeData, x0, step * decim, width);

    float *k_list = (float *)malloc(extra_width * sizeof(float));
    int *index_list = (int *)malloc(extra_width * sizeof(int));
    int list_size =
        findK(rangeData, index_list, k_list, x0, step * decim, width);
    free(rangeData);

    // k_list = (float*)realloc(k_list, list_size * sizeof(float));
    // index_list = (int*)realloc(index_list, list_size * sizeof(int));

    out_matrix = (float *)malloc(list_size * height * sizeof(float));

#pragma omp parallel for
    for (int i = 0; i < height; i++) {

        float *row = &in_matrix[i * width];
        row[0] = 0;
        for (int j = 1; j < list_size; j++) {
            out_matrix[i * list_size + j] =
                getA(row[index_list[j - 1]], row[index_list[j]], k_list[j]);
        }
    }

    free(in_matrix);
    free(k_list);
    free(index_list);

    image->ImageWidth = list_size;
    image->pImageData = out_matrix;

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());

    return 0;
}

/*!
   \brief       // Угловая коррекция по допплеровской частоте (RDA)

   \param[in]   image   Указатель на экземпляр CImage
   \param[in]   x0      Начальная координата по оси x
   \param[in]   lx      Дальность изображения

    \return     0 - успешное выполнение;
                -1 - ошибка
*/

inline int correctY(int x, int y, float r, float k, float x0, float dmax,
                    int height) {

    float fd = k * y * r / x0;
    return (fd * height / (2 * dmax)) + 0.5;
}

int CSAR::angleCorrection(CImage *image, float x0, float lx) {

    Log->LogPrintf("\t*Угловая коррекция...\t\t");
    fflush(stdout);
    time_start();

    int width = image->ImageWidth;
    int height = image->ImageHeight;

    Ipp32f *in_matrix = NULL;
    Ipp32f *out_matrix = NULL;

    if (image->ImageType == CImage::IT_FCOMPLEX) {
        in_matrix = (Ipp32f *)malloc(width * height * sizeof(Ipp32f));
        memset(in_matrix, 0, width * height * sizeof(Ipp32f));

        image->magnitude((Ipp32fc *)image->pImageData, in_matrix,
                         image->ImageWidth, image->ImageHeight);
        free((Ipp32fc *)image->pImageData);

        image->ImageType = CImage::IT_FLOAT;
    } else {
        in_matrix = (Ipp32f *)image->pImageData;
    }

    float Tm = (float)Radar->Ns / (float)Radar->Fs;
    float Tpp = (Tm / (float)Radar->Ns) * (Radar->Nppb + Radar->Nppe);
    float dmax = 1.0f / 2 * (Tm + Tpp);
    float k = ((2 * dmax) / (height - 1));
    float step = lx / width;
    float r = x0;
    int Hmax =
        correctY(width, height - 1, r + (width * step), k, x0, dmax, height);

    out_matrix = (float *)malloc(width * Hmax * sizeof(float));
    memset(out_matrix, 0, width * Hmax * sizeof(float));

#pragma omp parallel for
    for (int j = 0; j < width; j++) {
        int y;
        float r = x0 + (j * step);
        int offset =
            Hmax / 2 - (correctY(j, height - 1, r, k, x0, dmax, height) / 2);
        int lasty = offset;

        for (int i = 0; i < height; i++) {

            y = correctY(j, i, r, k, x0, dmax, height);
            y += offset;

            out_matrix[y * width + j] = in_matrix[i * width + j];

            int dy = y - lasty;
            if (dy > 1) {
                float dA =
                    out_matrix[y * width + j] - out_matrix[lasty * width + j];
                float stp = dA / dy;

                dA = out_matrix[lasty * width + j];
                lasty++;
                while (y > lasty) {

                    dA += stp;
                    out_matrix[lasty * width + j] = dA;

                    lasty++;
                }
            }
        }
    }

    free(in_matrix);

    image->ImageHeight = Hmax;
    image->pImageData = out_matrix;

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());
}

void smooth(float *addr_in_t, int m, int n) {

    float *addr_end = addr_in_t + m - n;
    float *tmp = addr_in_t;

    float data = 0; // регистр для вычислений
    float data0 = 0;

    for (unsigned char i = 0; i < n - 1; i++) {
        data += *(addr_in_t + i);
    }
    while (addr_in_t <= addr_end) {
        data0 = *addr_in_t;
        data += *(addr_in_t + (n - 1));
        *addr_in_t = (data) / n;

        data -= data0;
        addr_in_t++;
    }

    addr_in_t = tmp;
}

int argmax(float *arr, int n) {
    float max = 0;
    int arg = 0;
    for (int i = 0; i < n; i++) {

        if (arr[i] > max) {
            max = arr[i];
            arg = i;
        }
    }

    return arg;
}

float median_of(float newVal, float *buffer, int w_size) {
    static int count = 0;
    buffer[count] = newVal;
    if ((count < w_size - 1) and (buffer[count] > buffer[count + 1])) {
        for (int i = count; i < w_size - 1; i++) {
            if (buffer[i] > buffer[i + 1]) {
                float buff = buffer[i];
                buffer[i] = buffer[i + 1];
                buffer[i + 1] = buff;
            }
        }
    } else {
        if ((count > 0) and (buffer[count - 1] > buffer[count])) {
            for (int i = count; i > 0; i--) {
                if (buffer[i] < buffer[i - 1]) {
                    float buff = buffer[i];
                    buffer[i] = buffer[i - 1];
                    buffer[i - 1] = buff;
                }
            }
        }
    }
    if (++count >= w_size)
        count = 0;
    return buffer[(int)w_size / 2];
}

void median_filter(float *arr, int size) {
    const int buf_size = 10;
    float *buffer = (float *)malloc(buf_size * sizeof(float));
    memset(buffer, 0, buf_size);

    for (int i = 0; i < size; i++) {
        arr[i] = median_of(arr[i], buffer, buf_size);
    }
    for (int i = 0; i < buf_size; i++) {
        arr[i] = arr[buf_size + 1];
    }

    free(buffer);
}

float CSAR::driftAngle(CImage *image, float V, char *FnPlot) {

    Log->LogPrintf("\t*Вычисление угла сноса...\t");
    fflush(stdout);
    time_start();
    int width = image->ImageWidth;
    int height = image->ImageHeight;
    int x0 = width / 8;
    int x1 = width / 2;

    int w = x1 - x0;
    float Tp = (float)Radar->Ns / (float)Radar->Fs;
    float lambda = Radar->GetLambda();

    Ipp32fc *data = (Ipp32fc *)image->pImageData;

    int pow_height = 1;
    while (pow_height < height) {
        pow_height = pow_height << 1;
    }

    Ipp32fc *transposed = (Ipp32fc *)malloc(pow_height * w * sizeof(Ipp32fc));
    memset(transposed, 0, pow_height * w * sizeof(Ipp32fc));

    for (int i = 0; i < height; i++) {
        for (int j = x0; j < x1; j++) {
            transposed[(j - x0) * pow_height + i] = data[i * width + j];
        }
    }

    DSP->fft(transposed, pow_height, w);
    DSP->fftshift(transposed, pow_height, w);

    float *slices = (float *)malloc(pow_height * w * sizeof(float));

    image->magnitude(transposed, slices, pow_height, w);
    free(transposed);

    for (int j = pow_height; j < pow_height * w; j += pow_height) {
        for (int i = 0; i < pow_height; i++) {
            slices[i] += slices[i + j];
        }
    }

    median_filter(slices, pow_height);

    float fdc =
        (argmax(slices, pow_height) / (pow_height * Tp)) - (1 / (2 * Tp));

    float angle = asin(fdc * lambda / 2 / V);

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());
    Log->LogPrintf("\t\tУгол сноса = %2.3f град.\n", angle * 180.0f / M_PI);
    // return wrap_float1d(slices, pow_height);

    if (FnPlot) {
        Log->LogPrintf(
            "\t*Сохранение графика средней доплеровской частоты...\n");
        fflush(stdout);
        image->Plot.PlotChart(slices, pow_height);
        image->SaveToJpg8_(FnPlot, image->Plot.Width, image->Plot.Height,
                           image->Plot.pFrameBuffer, 80);
    }

    free(slices);

    return angle;
}
