#include "CBackProjection.hpp"
#include <limits.h>

CBackProjection::CBackProjection() {}

CBackProjection::CBackProjection(CSAR *sar, CImage *image)
    : SAR(sar), Image(image) {
    V = NULL;
    H = NULL;
    Window = (float *)SAR->createWindow(Image->ImageW0, W_type);
}

CBackProjection::~CBackProjection() {
    if (V) {
        free(V);
    }
    if (H) {
        free(H);
    }
    if (Window) {
        free(Window);
    }

    if (AccImage) {
        AccImage->Free();
        delete AccImage;
    }

    if (StripImage) {
        StripImage->Free();
        // Тут delete не нужен, т.к. StripImage передается в python и должен
        // быть уничтожен оттуда
        // delete StripImage;
    }
}

void CBackProjection::setPixelSize(float pixSize) {
    dx = pixSize;
    dy = pixSize;
}

void CBackProjection::setPixelSize(float _dx, float _dy) {
    dx = _dx;
    dy = _dy;
}

int CBackProjection::setVH(float v, float h, int H0) {

    if (H0 == -1) {
        H0 = Image->ImageH0;
    }

    if (V) {
        free(V);
    }
    if (H) {
        free(H);
    }

    int size_bytes = H0 * sizeof(float);
    V = (float *)malloc(size_bytes);
    H = (float *)malloc(size_bytes);

    for (int i = 0; i < H0; i++) {
        V[i] = v;
        H[i] = h;
    }

    return 0;
}

int CBackProjection::setVH(pybind11::list &v, pybind11::list &h) {
    return setVH(v, h, -1);
}

int CBackProjection::setVH(pybind11::list &v, pybind11::list &h, int H0) {

    if (H0 == -1) {
        H0 = Image->ImageH0;
    }

    if (V) {
        free(V);
    }
    if (H) {
        free(H);
    }

    int size_bytes = H0 * sizeof(float);

    V = (float *)malloc(size_bytes);
    H = (float *)malloc(size_bytes);
    if (v.size() == H0) {
        for (int i = 0; i < v.size(); i++) {
            V[i] = v[i].cast<float>();
            H[i] = h[i].cast<float>();
        }
    } else {
        for (int i = 0; i < H0; i++) {
            V[i] = v[0].cast<float>();
            H[i] = h[0].cast<float>();
        }
    }
    // printf("Size = %d",  v.size(), V[0]);

    return 0;
}

int CBackProjection::setWindowType(CSAR::WindowType w_type) {
    W_type = w_type;

    if (Window) {
        free(Window);
    }

    Window = (float *)SAR->createWindow(Image->ImageW0, W_type);
    return 0;
}

int CBackProjection::setProgressive(bool p) {
    progressive = p;
    return 0;
}

int CBackProjection::telescopic(float x0, float lx, float y0, float ly) {

    assert(V);
    assert(H);

    SAR->Log->LogPrintf("\t*Формирование РЛИ...\t\t");
    fflush(stdout);
    time_start();

    int nx = (int)(lx / dx);
    int ny = (int)(ly / dy);

    float *outdata = (float *)malloc(nx * ny * sizeof(float));
    memset(outdata, 0, nx * ny * sizeof(float));

    if (progressive) {
        SAR->DSP->bp_progressive(Image, SAR->Radar, outdata, Window, V, x0, y0,
                                 dx, dy, H, nx, ny, bCorr, pc);
    } else {
        SAR->DSP->bp(Image, SAR->Radar, outdata, Window, V, x0, y0, dx, dy, H,
                     nx, ny, bCorr, pc);
    }

    free(Image->pImageData);
    Image->pImageData = outdata;
    Image->ImageWidth = nx;
    Image->ImageHeight = ny;
    Image->SetImageType(CImage::IT_FLOAT);

    SAR->Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());
    return 0;
}

int CBackProjection::telescopicAngle(float x0, float lx, float angle,
                                     float div) {

    assert(V);
    assert(H);

    SAR->Log->LogPrintf("\t*Формирование РЛИ (2)...\t");
    fflush(stdout);
    time_start();

    float ly = (lx + x0) * tan(div);
    float y0 = -ly / 2;
    int nx = (int)(lx / dx);
    int ny = (int)(ly / dy);

    float *outdata = (float *)malloc(nx * ny * sizeof(float));
    memset(outdata, 0, nx * ny * sizeof(float));

    if (progressive) {
        SAR->DSP->bpm_progressive(Image, SAR->Radar, outdata, Window, V, x0, y0,
                                  dx, dy, H, nx, ny, bCorr, angle, div, pc);
    } else {
        SAR->DSP->bpm(Image, SAR->Radar, outdata, Window, V, x0, y0, dx, dy, H,
                      nx, ny, bCorr, angle, div, pc);
    }

    free(Image->pImageData);
    Image->pImageData = outdata;
    Image->ImageWidth = nx;
    Image->ImageHeight = ny;
    Image->SetImageType(CImage::IT_FLOAT);

    SAR->Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());

    return 0;
}

void unwrap(float *out, float *in, int len) {
    out[0] = in[0];
    for (int i = 1; i < len; i++) {
        float d = in[i] - in[i - 1];
        d = d > M_PI ? d - 2 * M_PI : (d < -M_PI ? d + 2 * M_PI : d);
        out[i] = out[i - 1] + d;
    }
}

void complexAngle(float *out, Ipp32fc *in, int len) {
    for (int i = 0; i < len; i++) {
        out[i] = atan2(in[i].im, in[i].re);
    }
}

int CBackProjection::phaseCorrection(float x, float y) {
    Ipp32fc *data = (Ipp32fc *)Image->pImageData;
    int Np = Image->ImageH0;
    int Ns = Image->ImageWidth;

    float Tp = SAR->Radar->GetTp();
    float Mu = SAR->Radar->GetMu();
    float F0 = SAR->Radar->F0;
    float Fs = SAR->Radar->Fs;

    float dt = Tp;
    float c1 = -4 * M_PI / SPEED_OF_LIGHT * F0;
    float c2 = -M_PI * Tp;
    float c5 = 2 * Mu / SPEED_OF_LIGHT;
    float dfr = Fs / (2 * (Ns - 1));
    float dr = Fs * SPEED_OF_LIGHT / (4 * Mu * (Ns - 1));

    Ipp32fc *signal = (Ipp32fc *)malloc(Np * 2 * sizeof(float));
    Ipp32fc *sop = (Ipp32fc *)malloc(Np * 2 * sizeof(float));
    Ipp32fc *sk = (Ipp32fc *)malloc(Np * 2 * sizeof(float));

    // #pragma omp parallel for
    // #pragma block_loop
    for (int k = 0; k < Np; k++) {
        float t = k * dt;
        float c4 = powf(x, 2) + powf(H[k], 2);
        float r = sqrtf(c4 + powf(y + V[k] * t, 2));

        int ind = roundf(r / dr);

        signal[k] = data[k * Ns + ind];
        float phase = c1 * r;

        sop[k].re = cosf(phase);
        sop[k].im = sinf(phase);

        // Коррекция скачков фазы
        if (bCorr) {
            float fp = (ind - 1) * dfr;
            // float fr = range2fr(r);
            float fr = r * c5;
            float phi = c2 * (fr - fp);

            sk[k].re = cosf(phi);
            sk[k].im = sinf(phi);
        }
    }

    if (bCorr) {
        SAR->DSP->complexMul(signal, sk, Np);
    }
    free(sk);

    SAR->DSP->complexMul(signal, Window, Np); // Оконное взвешивание
    SAR->DSP->complexMul(signal, sop, Np); // демодуляция
    free(sop);

    float *wrapped = (Ipp32f *)malloc(Np * sizeof(float));
    float *unwrapped = (Ipp32f *)malloc(Np * sizeof(float));
    // SAR->DSP->complexAngle(pc, signal, Np);

    complexAngle(wrapped, signal, Np);
    free(signal);

    unwrap(unwrapped, wrapped, Np);
    free(wrapped);

    /*
    FILE *f = fopen("angle.txt", "w");
    for(int k = 0; k < Np; k++){
        fprintf(f, "%f\n", pc[k]);
    }
    fclose(f);
    */

    pc = (Ipp32fc *)malloc(Np * 2 * sizeof(float));
    for (int k = 0; k < Np; k++) {
        pc[k].re = cosf(unwrapped[k]);
        pc[k].im = sin(unwrapped[k]);
    }

    return 0;
}

int CBackProjection::stripInit(float lx, float Tstrip, float Ts,
                               pybind11::list &v, pybind11::list &h) {

    snx = lx / dx;
    int nb = Image->ImageH0;
    int nstrip = SAR->Radar->Time2Np((Tstrip) / nb);
    int AccW = Image->GetWidth() / 2; // Complex
    int AccH = SAR->Radar->Time2Np(Ts);

    AccImage = new CImage(SAR->Log);
    StripImage = new CImage(SAR->Log);

    StripImage->Create(snx, nstrip, CImage::IT_FLOAT);
    AccImage->Create(AccW, AccH, CImage::IT_FCOMPLEX);

    if (Window) {
        free(Window);
    }

    Window = (float *)SAR->createWindow(AccImage->ImageH0, W_type);

    setVH(v, h, AccImage->ImageH0);

    return 0;
}

int CBackProjection::strip(float x0, float y0, int n, float angle) {

    // AccImage->Roll(Image->ImageH0);
    AccImage->CopyImageTop(Image);
    SAR->DSP->bp_strip(AccImage, StripImage, SAR->Radar, Window, V, x0, y0, dx,
                       H, snx, n, bCorr, NULL, angle);

    return 0;
}

CImage *CBackProjection::getStripImage() { return StripImage; }

float CBackProjection::autofocus(float v0, float v1, float point_x,
                                 float point_y, float Ts, float ang, float ls,
                                 int n, char *report) {

    assert(V);
    assert(H);

    SAR->Log->LogPrintf("\t*Автофокусировка...\t\t");
    fflush(stdout);
    time_start();

    float v_step = (v1 - v0) / n;

    if (v_step < 0.0) {
        SAR->Log->LogPrintf("\n\t\tНе правильно заданы параметры v0, v1!\n");
        fflush(stdout);
        return 0.0;
    }

    float v;
    float y_offs;
    float x_offs;
    float _V = V[0];
    float _H = H[0];

    float x = point_x * cos(ang) - point_y * sin(ang);
    float y = point_x * sin(ang) + point_y * cos(ang); // * -1 ?

    int nx = (int)(ls / dx);
    int ny = (int)(ls / dy);

    // CLog log;
    // log.SetVerbose(false);
    // CImage out(&log);
    SAR->Log->SetVerbose(false);
    CImage out(SAR->Log);
    out.Create(nx, ny, CImage::IT_FLOAT);
    out.ImageWidth = nx;
    out.ImageHeight = ny;

    char buf[32];

    float x0 = x - ls / 2;
    float y0 = y - ls / 2;

    float e = std::numeric_limits<float>::max();
    float _v = 0.0;
    float e_arr[n];

    for (int i = 0; i < n; i++) {
        v = v0 + v_step * i;
        y_offs = (_V * (Ts * _V + y) / v - Ts * _V);
        x_offs = sqrtf(x * x + y * y - y_offs * y_offs);
        y_offs -= y;
        x_offs -= x;

        setVH(v, _H, -1);

        memset(out.pImageData, 0, nx * ny * sizeof(float));
        SAR->DSP->bp(Image, SAR->Radar, (float *)out.pImageData, Window, V,
                     x0 + x_offs, y0 + y_offs, dx, dy, H, nx, ny, bCorr, pc);

        float tmp = SAR->DSP->entropyf((float *)out.pImageData, nx, ny);
        if (tmp < e) {
            e = tmp;
            _v = v;
        }

        if (report) {
            e_arr[i] = tmp;
            sprintf(buf, "%s/%0.2f km_h.jpg", report, v * 3.6);
            out.SaveToJpg8(buf, 0, 80);
        }
    }

    if (report) {
        CPlot plt = CPlot();
        plt.Initialize(800, 600);

        plt.PlotChart(v0 * 3.6, v1 * 3.6, e_arr, n);

        sprintf(buf, "%s/e.jpg", report);
        CImage::SaveToJpg8(buf, plt.Width, plt.Height, plt.pFrameBuffer, 80);
    }

    SAR->Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());

    return _v;
}
