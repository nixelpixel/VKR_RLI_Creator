/*
 * CRadar.cpp
 *
 *  Created on: 31 мая 2017 г.
 *      Author: Рязанцев Леонид
 */

#include "CRadar.hpp"
#include "helpers.hpp"

CRadar::CRadar() : Log(NULL) {
    DF = 0;
    F0 = 0;
    Fs = 0;
    Ns = 0;
    Nsync = 0;
    Sample_Size = SAMPLE_SIZE_UNDEFINED;
    Ramp_Type = RAMP_UNDEFINED;
    NSeek = 0;
}

/*!
   \brief Конструктор

    \param[in]  log     Указатель на экземпляр класса, реализующего ведение
   журнала работы программы

    \return      Функция не возвращает значения
*/
CRadar::CRadar(CLog *log) : Log(log) {
    DF = 0;
    F0 = 0;
    Fs = 0;
    Ns = 0;
    Nsync = 0;
    Nppb = 0;
    Nppe = 0;
    Sample_Size = SAMPLE_SIZE_UNDEFINED;
    Ramp_Type = RAMP_UNDEFINED;
    Tm = 0;
    NSeek = 0;
}

CRadar::CRadar(CLog *log, const char *address, const char *port) {
    DF = 0;
    F0 = 0;
    Fs = 0;
    Ns = 0;
    Nsync = 0;
    Nppb = 0;
    Nppe = 0;
    Sample_Size = SAMPLE_SIZE_UNDEFINED;
    Ramp_Type = RAMP_UNDEFINED;
    Tm = 0;
    NSeek = 0;
}

CRadar::~CRadar() {
    // TODO Auto-generated destructor stub
}

/*!
   \brief Установка параметров сигналов РЛС и выдаваемых данных

    \param[in]  f0              Начальная частота ЛЧМ зондирующего сигнала АЦП,
   Гц \param[in]  dF              Девиация ЛЧМ зондирующего сигнала, Гц
    \param[in]  fs              Частота дискретизации АЦП, Гц
    \param[in]  ns              Количество отсчетов за период зондирования
    \param[in]  nsync           Количество отсчетов в синхропаузе
    \param[in]  nppb            Количество отсчетов переходных процессов в
   начале периода \param[in]  nppe            количество отсчетов переходных
   процессов в конце периода \param[in]  mp              Коэффициент коррекции
   смещения отметок по дальности вследствие дробного отношения количества
   отсчетов АЦП за период модуляции \param[in]  sample_size     Количество бит,
   отводимых на один отсчет \param[in]  ramp_type       Тип модуляции сигнала
   (симметричный, несимметричный)

    \return      Функция не возвращает значения
*/
void CRadar::SetParams(float f0, float df, int fs, int ns, int nsync, int nppb,
                       int nppe, float mp, float tm, SampleSize sample_size,
                       RampType ramp_type) {
    F0 = f0;
    DF = df;
    Fs = fs;
    Ns = ns;
    Nsync = nsync;
    Nppb = nppb;
    Nppe = nppe;
    Mp = mp;
    Sample_Size = sample_size;
    Ramp_Type = ramp_type;
    Tm = tm;
}

int CRadar::SetParamsFromFile(char *path) {
    header_t header;
    uint32_t key;
    uint32_t size;
    FILE *f = fopen(path, "r");
    fread(&key, sizeof(uint32_t), 1, f);

    if (key == HEADER_KEY) {
        fread(&size, sizeof(uint32_t), 1, f);

        if (size != sizeof(header_t)) {
            Log->LogPrintf(
                "ВНИМАНИЕ! Несоответствие структуры заголовков файла.\n");
            Log->LogPrintf("Header size = %d.\n", size);
            Log->LogPrintf("Struct size = %d.\n", sizeof(header_t));
        }
        rewind(f);
        fread(&header, sizeof(header_t), 1, f);
        fclose(f);

        /*
        printf("Ns %d\n", header.Ns);
        printf("Fs %d\n", header.Fs);
        printf("Nsync %d\n", header.Nsync);
        printf("Nppb %d\n", header.Nppb);
        printf("Nppe %d\n", header.Nppe);

        printf("Tm %f\n", header.Tm);
        printf("F0 %f\n", header.F0);
        printf("DF %f\n", header.DF);
        printf("Mp %f\n", header.Mp);

        printf("RampType %d\n", header.RampType);
        printf("SampleSize %d\n", header.SampleSize);
        */

        F0 = header.F0;
        DF = header.DF;
        Fs = header.Fs;
        Ns = header.Ns;
        Nsync = header.Nsync;
        Nppb = header.Nppb;
        Nppe = header.Nppe;
        Mp = header.Mp;
        Sample_Size = (SampleSize)header.SampleSize;
        Ramp_Type = (RampType)header.RampType;
        Tm = header.Tm;

        return size;
    }

    Log->LogPrintf("ВНИМАНИЕ! Файл %s не имеет заголовка!\n", path);
    return -1;
}

/*!
   \brief Определение количества периодов за время синтезирования

    \param[in]  Ts              Время синтезирования, с

    \return     Количество периодов
*/
int CRadar::Time2Np(float Ts) {
    return round(Fs * Ts / (float)Ns); // Количество периодов в обработке
}

/*!
   \brief Возвращает длинну волны зондирующего сигнала

    \return      Длинна волны зондирующего сигнала.
*/
float CRadar::GetLambda(void) { return SPEED_OF_LIGHT / (float)F0; }

/*!
   \brief Возвращает период повторения зондирующего сигнала

    \return      Период повторения зондирующего сигнала.
*/
float CRadar::GetTp(void) { return (float)Ns / (float)Fs; }

/*!
   \brief Возвращает период модуляции зондирующего сигнала

    \return      Период модуляции зондирующего сигнала.
*/
float CRadar::GetTm(void) {
    return Tm; //(float)Ns/(float)Fs;
}

/*!
   \brief Возвращает крутизну модуляции ЛЧМ зондирующего сигнала

    \return      Крутизна модуляции ЛЧМ зондирующего сигнала.
*/
float CRadar::GetMu(void) { return (float)DF / GetTm(); }
