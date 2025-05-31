/*
 * CRadarSTT.cpp
 *
 *  Created on: 3 июл. 2017 г.
 *      Author: user
 */

#include "CRadarSTT.hpp"

CRadar_STT::CRadar_STT(CLog *log, const char *address, const char *port)
    : CRadar(log) {}

CRadar_STT::CRadar_STT(CLog *log) : CRadar(log) {}

CRadar_STT::CRadar_STT() {
    // TODO Auto-generated constructor stub
}

CRadar_STT::~CRadar_STT() {
    // TODO Auto-generated destructor stub
}

/*!
   \brief Считывание данных из файла согласно переданных функцией SetParams()
   параметров.

    Порядок считывания данных:
    - открытие файла для чтения;
    - установка указателя в соответствии с параметром tshift;
    - выделение объема памяти под считываемые данные, кратного степени двойки;
    - считывание данных, нормировка к единичной шкале и преобразование отсчетов
   к типу float;

    Функция осуществляет выделение памяти под считываемые данные. Освобождение
   памяти осуществляется программистом. Количество периодов и считываемых за
   период отсчетов дополняется до степени двойки.

    \param[out] *data   Указатель на область памяти, в которую будет произведено
   считывание данных \param[in]  np      Количество периодов, подлежащих
   считыванию \param[in]  tshift  Сдвиг по времени от начала голограммы, с
    \param[in]  syncpos Смещение до синхропаузы
    \param[in]  kr      Коэффициент частотной интерполяции (степень двойки). При
   kr = -1 дополнение нулями до степени двойки не осуществляется \param[in]
   comp2   Признак дополнения отсчетов до степени двойки в строках (по
   количеству зондирований) \param[in]  rmdist  Признак удаления искажений в
   начале и в конце периода

    \return      0 - успешное выполнение;
                -1 - ошибка
*/
int CRadar_STT::ReadDataFromFile(char *filename, CImage *image, int np,
                                 float tshift, int syncpos, int kr, int comp2,
                                 int rmdist) {
    Log->LogPrintf("\t*Считывание голограммы...\t");
    fflush(stdout);
    time_start();

    int Np2 = 0;

    if (comp2 == 0) {
        Np2 = np;
    } else {
        Np2 = nearest_power_of_two(np);
    }

    if (kr < 0) {
        //      Np2 = np;
        Ns2 = Ns;
    } else {
        //      Np2 = nearest_power_of_two(np) * pow(2, kr);
        Ns2 = nearest_power_of_two(Ns) * pow(2, kr);
    }

    float bytespersample =
        (float)(Sample_Size) / 8.0f; // Количество байт на один отсчет

    // Открываем файл для чтения
    FILE *fid = fopen(filename, "r");
    if (fid == NULL) {
        Log->LogPrintf("\e[31m CRadar_STT::ReadDataFromFile(). Ошибка чтения "
                       "файла голограммы \"%s\" \e[0m \n",
                       filename, "");
        return -1;
    }

    // Сдвигаем указатель начала файла на нужную позицию
    int Shift = tshift * Fs * bytespersample;
    if (bytespersample > SAMPLE_SIZE_8BIT)
        Shift = (Shift / 2) * 2; // Начинаем с четного байта
    fseek(fid, Shift + syncpos * bytespersample, SEEK_SET);

    // Создаем изображение
    image->Create(Ns2, Np2, CImage::IT_FLOAT);
    float *data = (float *)image->pImageData;

    // Считываем данные попериодно (для ускорения процесса считывания) и
    // преобразовываем их к типу float
    size_t rb = 0;
    float m = 0;
    void *samplebuf = aligned_alloc(
        64, Ns * bytespersample); // Буфер для попериодного считывания
    for (int i = 0; i < np; i++) {
        m = 0;
        switch (Sample_Size) {
        case SAMPLE_SIZE_UNDEFINED:
            break;
        case SAMPLE_SIZE_1BIT:
        case SAMPLE_SIZE_2BIT:
        case SAMPLE_SIZE_4BIT:
            rb = fread(samplebuf, sizeof(uint8_t), Ns * bytespersample, fid);
            if (rb != (size_t)(Ns * bytespersample)) {
                Log->LogPrintf("ReadDataFromFile(). Достигнут конец файла.");
            }

            // убираем переходные процессы в начале и конце периода
            if (rmdist) {
                memset(samplebuf, 8, (int)(Nppb * bytespersample));
                memset(samplebuf + Ns - Nppe, 8, (int)(Nppe * bytespersample));
            }

            // Преобразование к типу float, компенсация постоянной составляющей
            // и нормировка
            for (int j = 0; j < Ns; j += 2) {
                data[i * Ns2 + j] =
                    ((float)(((uint8_t *)samplebuf)[j]) - m) / 16.0f - 0.5f;
                data[i * Ns2 + j + 1] =
                    ((float)(((uint8_t *)samplebuf)[j + 1]) - m) / 16.0f - 0.5f;
            }
            break;

        case SAMPLE_SIZE_8BIT:
            rb = fread(samplebuf, sizeof(uint8_t), Ns, fid);
            if (rb != (size_t)Ns) {
                Log->LogPrintf("ReadDataFromFile(). Достигнут конец файла.");
            }

            // убираем переходные процессы в начале и конце периода
            if (rmdist) {
                /*
                for(int i=0; i<Nppb; i++)
                    ((uint8_t*)samplebuf)[i] = 127;

                for(int i=0; i<Nppe; i++)
                    ((uint8_t*)samplebuf)[Ns-Nppe+i] = 127;
                */

                memset(samplebuf, 127, Nppb);
                memset(samplebuf + Ns - Nppe, 127, Nppe);
            }
            /*
            for(int j=0; j<Ns; j++){
                m = m + (((uint8_t*)samplebuf)[j]); // Преобразование к типу
            float и нормировка
            }
            m = (m/(float)Ns);
            //m = 0;
            */
            // Преобразование к типу float, компенсация постоянной составляющей
            // и нормировка
            for (int j = 0; j < Ns; j++)
                data[i * Ns2 + j] =
                    ((float)(((uint8_t *)samplebuf)[j]) - m) / 255.0f - 0.5f;
            break;
        case SAMPLE_SIZE_16BIT:
            rb = fread(samplebuf, sizeof(int16_t), Ns, fid);
            if (rb != (size_t)Ns) {
                Log->LogPrintf("ReadDataFromFile(). Достигнут конец файла.");
            }
            for (int j = 0; j < Ns; j++) {
                data[i * Ns2 + j] =
                    (float)(((int16_t *)samplebuf)[j]) /
                    65536.0f; // Преобразование к типу float и нормировка
            }
            break;
        }
    }
    free(samplebuf);
    fclose(fid);

    // Установка параметров изображения
    image->ImageW0 = Ns;
    image->ImageH0 = np;
    image->SetImageWidth(Ns2);
    image->SetImageHeight(Np2);
    image->SetImageType(CImage::IT_FLOAT);

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());

    return 0;
}

/*!
   \brief Поблочное считывание данных.

    \param[in]  filename    Имя файла
    \param[in] *image   Указатель  на CImage
    \param[in]  Np      Количество периодов, подлежащих считыванию
    \param[in]  kr      Коэффициент частотной интерполяции (степень двойки). При
   kr = -1 дополнение нулями до степени двойки не осуществляется

    \return      0 - успешное выполнение;
                -1 - ошибка
*/
int CRadar_STT::ReadBlockFromFile(char *filename, CImage *image, int np,
                                  int kr) {
    Log->LogPrintf("\t*Считывание блока...\t\t");
    fflush(stdout);
    time_start();

    int Np2 = np;

    if (kr < 0) {
        Ns2 = Ns;
    } else {
        Ns2 = nearest_power_of_two(Ns) * pow(2, kr);
    }

    float bytespersample =
        (float)(Sample_Size) / 8.0f; // Количество байт на один отсчет

    // Открываем файл для чтения
    FILE *fid = fopen(filename, "r");
    if (fid == NULL) {
        Log->LogPrintf("\e[31m CRadar_STT::ReadDataFromFile(). Ошибка чтения "
                       "файла голограммы \"%s\" \e[0m \n",
                       filename, "");
        return -1;
    }

    fseek(fid, NSeek * bytespersample, SEEK_SET);

    // Создаем изображение
    image->Create(Ns2, Np2, CImage::IT_FLOAT);
    float *data = (float *)image->pImageData;
    //  #pragma omp parallel for
    //  for (int i=0; i<Ns2*Np2; i++){
    //      data[i] = 0.5f;
    //  }

    // Считываем данные попериодно (для ускорения процесса считывания) и
    // преобразовываем их к типу float
    size_t rb = 0;
    float m = 0;
    void *samplebuf = aligned_alloc(
        64, Ns * bytespersample); // Буфер для попериодного считывания
    for (int i = np - 1; i >= 0; i--) {
        m = 0;
        switch (Sample_Size) {
        case SAMPLE_SIZE_UNDEFINED:
            break;
        case SAMPLE_SIZE_1BIT:
        case SAMPLE_SIZE_2BIT:
        case SAMPLE_SIZE_4BIT:
        case SAMPLE_SIZE_8BIT:
        case SAMPLE_SIZE_16BIT:
            rb = fread(samplebuf, sizeof(uint8_t), Ns, fid);
            if (rb != (size_t)Ns) {
                Log->LogPrintf("ReadDataFromFile(). Достигнут конец файла.");
            }

            // убираем переходные процессы в начале и конце периода
            memset(samplebuf, 127, Nppb);
            memset(samplebuf + Ns - Nppe, 127, Nppe);

            // Преобразование к типу float, компенсация постоянной составляющей
            // и нормировка
            for (int j = 0; j < Ns; j++)
                data[i * Ns2 + j] =
                    ((float)(((uint8_t *)samplebuf)[j]) - m) / 255.0f - 0.5f;

            break;
        }
    }
    free(samplebuf);
    fclose(fid);

    NSeek += Ns * np;

    // Установка параметров изображения
    image->ImageW0 = Ns;
    image->ImageH0 = np;
    image->SetImageWidth(Ns2);
    image->SetImageHeight(Np2);
    image->SetImageType(CImage::IT_FLOAT);

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());

    return 0;
}

/*!
   \brief Удаление переходных процессов.

    \param[in] *image   Указатель  на CImage
    \return      0 - успешное выполнение;
                -1 - ошибка
*/
int CRadar_STT::RemoveTransients(CImage *image) {

    int W = image->GetWidth();
    int H = image->GetHeight();
    int transient_begin = (int)(Fs * GetTm());
    int transient_size =
        (image->GetWidth() - (int)(Fs * GetTm())) * sizeof(float);

    float *data = (float *)image->pImageData;

    for (int i = 0; i < H; i++) {

        memset(data, 0, Nsync * sizeof(float));
        memset(&data[transient_begin], 0, transient_size);

        data = &data[W];
    }

    return 0;
}

int CRadar_STT::ReadBlockSetPos(int syncpos, float tshift) {
    float bytespersample = (float)(Sample_Size) / 8.0f;
    int Shift = tshift * Fs * bytespersample;
    NSeek = Shift + syncpos * bytespersample;

    return 0;
};

int CRadar_STT::ReadBlockResetPos(void) {
    NSeek = 0;
    return 0;
};

/*!
   \brief Обертка функции FindSyncPos() для Python

    \param[in]  image       Указатель на экземпляр CImage

    \return      Начало и длительность синхропаузы
*/
pybind11::tuple CRadar_STT::FindSyncPosPy(CImage *image) {
    int sync_begin;
    int sync_width;
    FindSyncPos(image, &sync_begin, &sync_width);
    return pybind11::make_tuple(sync_begin, sync_width);
}

/*!
   \brief Поиск положения синхропаузы

    \param[in]  image       Указатель на экземпляр CImage
    \param[out] sync_begin  Начало синхропаузы
    \param[out] sync_width  Длительность синхропаузы в отсчетах

    \return      0 - успешное выполнение;
                -1 - ошибка
*/
int CRadar_STT::FindSyncPos(CImage *image, int *sync_begin, int *sync_width) {
    Log->LogPrintf("\t*Поиск начала зондирования...\t");
    fflush(stdout);
    time_start();

    // Вычисляем порог
    float tresh = -0.4f;
    float *data = (float *)image->pImageData;
    int len = image->ImageWidth * image->ImageHeight;
    // for(int i = 0; i < len; i++){
    //   tresh += fabs(data[i]);
    // }
    // tresh /= 3*len;

    // Поиск синхропаузы
    int sync_is_found = 0; // Признак того, что найдена синхропауза
    *sync_begin = -1; // Отсчет начала синхропаузы
    *sync_width = 0;  // Длительность синхропаузы

    for (int i = 0; i < (len - Nsync); i++) {
        // Если синхропауза была найдена, то вычисление ее длительности
        if ((fabs(data[i]) > tresh) && (sync_is_found == 1)) {
            sync_is_found = 0;
            *sync_width = i - *sync_begin;
            break;
        }

        // Если синхропауза не найдена, то производится поиск ее начала
        // для этого двигается окно длительностьтю Nsync, если в пределах
        // этого окна сигнал хотя бы раз превышает порог, то считается, что
        // синхропаузы нет
        int j = 0;
        if (sync_is_found == 0) {
            for (j = i; j < (i + Nsync); j++) {
                if (data[j] > tresh) {
                    break;
                }
            }
        }

        // Если в пределах окна небыло превышения порога, то считается, что
        // синхропауза найдена
        if (j == (i + Nsync)) {
            sync_is_found = 1;
            *sync_begin = i;
        }
    }

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());
    return 0;
}

/*!
   \brief Возвращает длительность файла радиоголограммы в секундах

    \param[in]  filename    Имя файла РГ

    \return      длительность файла РГ в секундах
*/
float CRadar_STT::GetFileDuration(char *filename) {
    // Открываем файл для чтения
    FILE *fid = fopen(filename, "r");
    if (fid == NULL) {
        Log->LogPrintf("\e[31m CSAR::LoadDataFromFile(). Ошибка чтения файла "
                       "голограммы. \"%s\" \e[0m \n",
                       filename, "");
        return 0;
    }

    // Вычисляем длительность
    fseek(fid, 0L, SEEK_END);
    int size_of_file = ftell(fid);
    fseek(fid, 0L, SEEK_SET);

    float bitspersample = (float)Sample_Size / 8;
    return (float)size_of_file / (float)Fs /
           bitspersample; // Длительность голограммы
}

int CRadar_STT::OpenHologram(char *fn, CImage *image, int np, int syncpos,
                             int kr, int comp2) {

    Log->LogPrintf("\t*Открытие голограммы...\t");
    fflush(stdout);
    time_start();

    int Np2 = 0;

    if (comp2 == 0) {
        Np2 = np;
    } else {
        Np2 = nearest_power_of_two(np);
    }

    if (kr < 0) {
        Ns2 = Ns;
    } else {
        Ns2 = nearest_power_of_two(Ns) * pow(2, kr);
    }

    // Создание изображения
    image->Create(Ns2, Np2, CImage::IT_FLOAT);
    image->ImageW0 = Ns;
    image->ImageH0 = np;
    image->SetImageWidth(Ns2);
    image->SetImageHeight(Np2);
    image->SetImageType(CImage::IT_FLOAT);

    // Создание буффера голограммы
    FILE *fid = fopen(fn, "rb");
    if (fid == NULL) {
        Log->LogPrintf("\e[31m CRadar_STT::OpenHologram(). Ошибка чтения файла "
                       "голограммы \"%s\" \e[0m \n",
                       fn, "");
        return -1;
    }

    fseek(fid, 0, SEEK_END);
    long fsize = ftell(fid) - syncpos;
    fseek(fid, syncpos, SEEK_SET);

    HologramData.BlockSize = Ns * np;
    HologramData.Np2 = Np2;
    HologramData.Np = np;
    HologramData.pImageDataSize = (Ns2 * Np2 + 2) * 4;
    HologramData.pHologram = (uint8_t *)malloc(fsize);
    fread(HologramData.pHologram, 1, fsize, fid);
    fclose(fid);

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());

    return 0;
}

int CRadar_STT::ReadBlock(CImage *image, int block) {

    // image->Free();
    // image->Create(Ns2, HologramData.Np2, CImage::IT_FLOAT);

    memset(image->pImageData, 0, HologramData.pImageDataSize);
    image->SetImageType(CImage::IT_FLOAT);
    image->ImageW0 = Ns;
    image->ImageH0 = HologramData.Np;
    image->SetImageWidth(Ns2);
    image->SetImageHeight(HologramData.Np2);

    int offset = HologramData.BlockSize * block;
    float *data = (float *)image->pImageData;
    uint8_t *tmpHologram = HologramData.pHologram + offset;

    tmpHologram += Ns * (image->ImageH0 - 1);

    for (int i = 0; i < image->ImageH0; i++) {

        // убираем переходные процессы в начале и конце периода
        memset(tmpHologram, 127, Nppb);
        memset(tmpHologram + Ns - Nppe, 127, Nppe);

        for (int j = 0; j < Ns; j++) {
            data[j] = ((float)tmpHologram[j]) / 255.0f;
        }
        data += Ns2;
        tmpHologram -= Ns;
    }

    return 0;
}

void CRadar_STT::FreeHologram() {
    Log->LogPrintf("\t*Освобождение памяти голограммы...\n");
    fflush(stdout);
    free(HologramData.pHologram);
}

int CRadar_STT::OpenHologramSocket(char *ip_addr, CImage *image, int np,
                                   int syncpos, int kr, int comp2) {

    Log->LogPrintf("\t*Открытие сокета...\t");
    fflush(stdout);
    time_start();

    int Np2 = 0;

    if (comp2 == 0) {
        Np2 = np;
    } else {
        Np2 = nearest_power_of_two(np);
    }

    if (kr < 0) {
        Ns2 = Ns;
    } else {
        Ns2 = nearest_power_of_two(Ns) * pow(2, kr);
    }

    // Создание изображения
    image->Create(Ns2, Np2, CImage::IT_FLOAT);
    image->ImageW0 = Ns;
    image->ImageH0 = np;
    image->SetImageWidth(Ns2);
    image->SetImageHeight(Np2);
    image->SetImageType(CImage::IT_FLOAT);

    HologramData.BlockSize = Ns * np;
    HologramData.Np2 = Np2;
    HologramData.Np = np;
    HologramData.pImageDataSize = (Ns2 * Np2 + 2) * 4;
    HologramData.pHologram = (uint8_t *)malloc(HologramData.BlockSize);

    tcp_srv = new CTCPServer();
    int error = tcp_srv->bind(ip_addr);

    if (error) {
        Log->LogPrintf("Ошибка. (%2.3f c)\n", time_stop());
    } else {
        Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());
    }

    return error;
}

int CRadar_STT::ReadBlockFromSocket(CImage *image) {

    memset(image->pImageData, 0, HologramData.pImageDataSize);
    image->SetImageType(CImage::IT_FLOAT);
    image->ImageW0 = Ns;
    image->ImageH0 = HologramData.Np;
    image->SetImageWidth(Ns2);
    image->SetImageHeight(HologramData.Np2);

    float *data = (float *)image->pImageData;
    uint8_t *tmpHologram = HologramData.pHologram;

    tcp_srv->recv(tmpHologram, HologramData.BlockSize);

    tmpHologram += Ns * (image->ImageH0 - 1);

    for (int i = 0; i < image->ImageH0; i++) {

        // убираем переходные процессы в начале и конце периода
        memset(tmpHologram, 127, Nppb);
        memset(tmpHologram + Ns - Nppe, 127, Nppe);

        for (int j = 0; j < Ns; j++) {
            data[j] = ((float)tmpHologram[j]) / 255.0f;
        }
        data += Ns2;
        tmpHologram -= Ns;
    }

    return 0;
}

int CRadar_STT::InitDataFile(char *_fn, CImage *image, int np, int kr,
                             int comp2) {

    Log->LogPrintf("\t*Открытие голограммы...\t");
    fflush(stdout);
    time_start();

    int Np2 = 0;

    if (comp2 == 0) {
        Np2 = np;
    } else {
        Np2 = nearest_power_of_two(np);
    }

    if (kr < 0) {
        Ns2 = Ns;
    } else {
        Ns2 = nearest_power_of_two(Ns) * pow(2, kr);
    }

    // Создание изображения
    image->Create(Ns2, Np2, CImage::IT_FLOAT);
    image->ImageW0 = Ns;
    image->ImageH0 = np;
    image->SetImageWidth(Ns2);
    image->SetImageHeight(Np2);
    image->SetImageType(CImage::IT_FLOAT);

    char *fn = (char *)malloc(strlen(_fn));
    strcpy(fn, _fn);

    HologramData.BlockSize = Ns * np;
    HologramData.Np2 = Np2;
    HologramData.Np = np;
    HologramData.pImageDataSize = (Ns2 * Np2 + 2) * 4;
    HologramData.pHologram = (uint8_t *)fn;

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());

    return 0;
}

int CRadar_STT::ReadBlockFromInitializedFile(CImage *image, int syncpos,
                                             int block) {

    memset(image->pImageData, 0, HologramData.pImageDataSize);
    image->SetImageType(CImage::IT_FLOAT);
    image->ImageW0 = Ns;
    image->ImageH0 = HologramData.Np;
    image->SetImageWidth(Ns2);
    image->SetImageHeight(HologramData.Np2);

    FILE *f = fopen((char *)HologramData.pHologram, "r");

    // Ожидание следующего пакета
    while (1) {
        fseek(f, 0L, SEEK_END);
        size_t file_size = ftell(f);
        size_t required =
            HologramData.BlockSize * block + syncpos + HologramData.BlockSize;
        if (file_size >= required) {
            break;
        }
        usleep(100000);
    }

    uint8_t *tmpHologram = (uint8_t *)malloc(HologramData.BlockSize);
    uint8_t *tmp = tmpHologram;
    fseek(f, HologramData.BlockSize * block + syncpos, SEEK_SET);
    fread(tmpHologram, 1, HologramData.BlockSize, f);
    fclose(f);

    // int offset = HologramData.BlockSize * block;
    float *data = (float *)image->pImageData;
    // uint8_t* tmpHologram = HologramData.pHologram + offset;

    tmpHologram += Ns * (image->ImageH0 - 1);

    for (int i = 0; i < image->ImageH0; i++) {

        // убираем переходные процессы в начале и конце периода
        memset(tmpHologram, 127, Nppb);
        memset(tmpHologram + Ns - Nppe, 127, Nppe);

        for (int j = 0; j < Ns; j++) {
            data[j] = ((float)tmpHologram[j]) / 255.0f;
        }
        data += Ns2;
        tmpHologram -= Ns;
    }

    free(tmp);
    return 0;
}
