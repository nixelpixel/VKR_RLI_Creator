/*
 * CImage.cpp
 *
 *  Created on: 19 июн. 2017 г.
 *      Author: user
 */

#include "CImage.hpp"

CImage::CImage() {
    Log = NULL;
    //  DSP = new CDSP_IPP(Log);
    Plot.Initialize(DEF_PLOT_WIDTH, DEF_PLOT_HEIGHT);
    pImageData = NULL;
    ImageWidth = 0;
    ImageHeight = 0;
    ImageType = IT_UNDEFINED;
    ImageW0 = 0;
    ImageH0 = 0;
    ImageNearRangeIndex = 0;
    actualIndex = -1;
#ifdef SUPPORT_AF
    afImageData = 0;
#endif
}

CImage::CImage(CLog *log) : Log(log) {
    //  DSP = new CDSP_IPP(Log);
    Plot.Initialize(DEF_PLOT_WIDTH, DEF_PLOT_HEIGHT);
    pImageData = NULL;
    ImageWidth = 0;
    ImageHeight = 0;
    ImageType = IT_UNDEFINED;
    ImageW0 = 0;
    ImageH0 = 0;
    ImageNearRangeIndex = 0;
    actualIndex = -1;
#ifdef SUPPORT_AF
    afImageData = 0;
#endif
}

CImage::~CImage() {
    if (pImageData) {
        free(pImageData);
    }
    //  delete(DSP);

    //  if(afImageData){
    //      delete(afImageData);
    //      afImageData = NULL;
    //  }
    // TODO Auto-generated destructor stub
}

/*!
    \brief  Преобразует данные для дальнейшего использования библиотеки
   ArrayFire

    \return      Функция не возвращает значения.
*/

void CImage::PrepareForAF(void) {
    Log->LogPrintf("\t*Подготовка данных для AF...\t");
    fflush(stdout);
    time_start();

    //  Ipp32f* data = (Ipp32f*)pImageData;
    //  for(int i=0; i<16; i++) data[i] = i;
    //  ImageWidth = 16;
    //  ImageHeight = 1;

#ifdef SUPPORT_AF
    afImageData =
        af::array(ImageWidth, ImageHeight, (float *)pImageData, afHost);
#endif
    Free();

    // printf("%d", afImageData.isdouble());

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());
}

/*!
    \brief  Клонирование экземпляра класса

    \param[in]  image   Указатель на экземпляр CImage

    \return      Функция не возвращает значения.
*/
void CImage::CloneFrom(CImage *image) {
    ImageWidth = image->ImageWidth;
    ImageHeight = image->ImageHeight;
    ImageType = image->ImageType;
    ImageW0 = image->ImageW0;
    ImageH0 = image->ImageH0;
    ImageNearRangeIndex = image->ImageNearRangeIndex;
    pImageData = malloc(ImageWidth * ImageHeight * GetPixelSize(ImageType));
    memcpy(pImageData, image->pImageData,
           ImageWidth * ImageHeight * GetPixelSize(ImageType));
}

/*!
   \brief Установка ширины изображения (количества отсчетов в строке).

    \param[in]  width   Количество отсчетов в строке.

    \return      Функция не возвращает значения.
*/
void CImage::SetImageWidth(int width) { ImageWidth = width; }

/*!
   \brief Установка высоты изображения (количества строк).

    \param[in]  height  Количество строк.

    \return      Функция не возвращает значения.
*/
void CImage::SetImageHeight(int height) { ImageHeight = height; }

/*!
   \brief Установка типа изображения.

    \param[in]  type    Тип изображения.

    \return      Функция не возвращает значения.
*/
void CImage::SetImageType(eImageType type) { ImageType = type; }

/*!
   \brief Сохранение изображения в бинарный файл.

    \param[in]  filename    Имя файла

    \return      0 - успешное выполнение;
                -1 - ошибка
*/
int CImage::SaveToBin(char *fileName) {
    Log->LogPrintf("\t*Сохранение данных...\t\t");
    fflush(stdout);
    time_start();
    Log->LogPrintf("%dx%d ", ImageWidth, ImageHeight);

    // Открываем файл для записи
    FILE *fid = fopen(fileName, "w");
    if (fid == NULL) {
        Log->LogPrintf(
            "CSAR::LoadDataFromFile(). Ошибка создания файла. \"%s\"\n",
            fileName, "");
        return -1;
    }

    fwrite(pImageData, GetPixelSize(ImageType), ImageWidth * ImageHeight, fid);

    fclose(fid);

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());

    return 0;
}

int CImage::LoadBin(char *fileName, int width, int height, eImageType type) {

    Log->LogPrintf("\t*Загрузка данных...\t\t");
    fflush(stdout);
    time_start();
    Log->LogPrintf("%dx%d ", width, height);

    FILE *fid = fopen(fileName, "r");
    if (fid == NULL) {
        Log->LogPrintf(
            "CSAR::LoadDataFromFile(). Ошибка создания файла. \"%s\"\n",
            fileName, "");
        return -1;
    }

    Create(width, height, type);

    fread(pImageData, GetPixelSize(ImageType), ImageWidth * ImageHeight, fid);

    fclose(fid);

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());
}

//==================================
// magnitude()
//==================================
int CImage::magnitude(Ipp32fc *pInData, float *pOutData, int N, int M) {
#pragma omp parallel for
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < M; i++) {
            pOutData[j * M + i] =
                sqrt(pInData[j * M + i].im * pInData[j * M + i].im +
                     pInData[j * M + i].re * pInData[j * M + i].re);
            // ippsMagnitude_32fc(&pInData[j*M], &pOutData[j*M], M);
        }
    }

    return 0;
}

/*!
   \brief Сохранение амплитудного РЛИ в *.png файл c предварительной
   нормировкой.

    \param[in]  fileName    Имя файла
    \param[in]  brightness  Коэффициент увеличения яркости изображения
    \param[in]  colormap    Цветовая схема

    \return      0 - успешное выполнение;
                -1 - ошибка
*/
int CImage::SaveToPng8(char *fileName, float brightness, ColorMap colormap) {
    Log->LogPrintf("\t*Сохранение изображения...\t");
    fflush(stdout);
    time_start();

    // Взятие абсолютного значения
    float *magbuff;
    float *pbuff;
    int len = ImageWidth * ImageHeight;

    if (ImageType == IT_FCOMPLEX) {
        magbuff = (float *)malloc(len * sizeof(float));
        memset(magbuff, 0, len * sizeof(float));
        magnitude((Ipp32fc *)pImageData, magbuff, ImageWidth, ImageHeight);
        pbuff = magbuff;
    } else {
        pbuff = (float *)pImageData;
    }

    float maxval = 0;
    float minval = pbuff[0];

    // Нахождение максимального значения
    // #pragma omp parallel for
    for (int i = 0; i < len; i++) {
        if (pbuff[i] > maxval)
            maxval = pbuff[i];
        if (pbuff[i] < minval)
            minval = pbuff[i];
    }
    //  ippsMax_32f(pbuff, len, &maxval);
    //  ippsMin_32f(pbuff, len, &minval);

    // Нормировка
    uint8_t *savebuff;
    savebuff = (uint8_t *)malloc(len * sizeof(uint8_t));

    for (int i = 0; i < len; i++) {
        float val =
            (pbuff[i] - minval) / (maxval - minval) * 255.0f * brightness;
        savebuff[i] = val < 255.0f ? (uint8_t)val : 255;
    }

    // Сохранение в файл
    SaveToPng8_(fileName, ImageWidth, ImageHeight, savebuff, colormap);

    free(savebuff);

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());

    return true;
}

/*!
   \brief Сохранение РЛИ в *.png файл.

    \param[in]  fileName    Имя файла
    \param[in]  Width       Ширина изображения
    \param[in]  Height      Высота изображения
    \param[in]  pData       Указатель намассив данных
    \param[in]  colormap    Цветовая схема

    \return      0 - успешное выполнение;
                -1 - ошибка
*/
int CImage::SaveToPng8_(char *fileName, int Width, int Height, uint8_t *pData,
                        ColorMap colormap) {
    FILE *fp = NULL;
    png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
    const int depth = 8;
    png_colorp palette;
    int i;
    png_byte **rows;

    fp = fopen(fileName, "wb");
    if (!fp) {
        printf("Can`t open file for writing: %s \n", fileName);
        return -1;
    }

    rows = (png_byte **)aligned_alloc(64, Height * sizeof(png_byte *));

    for (i = 0; i < Height; i++) {
        rows[i] = (uint8_t *)pData + i * Width;
    }

    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (png_ptr == NULL) {
        printf("CSARProc::SaveToPng()->png_create_write_struct() FAILED \n");
        return -1;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == NULL) {
        printf("CSARProc::SaveToPng()->png_create_info_struct() FAILED \n");
        png_destroy_write_struct(&png_ptr, &info_ptr);
        return -1;
    }

    palette = (png_colorp)png_malloc(png_ptr, PNG_MAX_PALETTE_LENGTH *
                                                  (sizeof(png_color)));
    png_color color;

    for (unsigned int i = 0; i < PNG_MAX_PALETTE_LENGTH; i++) {
        switch (colormap) {
        case COLORMAP_GRAY:
            color.red = i;
            color.green = i;
            color.blue = i;
            break;
        case COLORMAP_JET:
            color = GetJetColor((float)i / PNG_MAX_PALETTE_LENGTH);
            break;
        }
        palette[i] = color;
    }

    png_set_PLTE(png_ptr, info_ptr, palette, PNG_MAX_PALETTE_LENGTH);

    png_init_io(png_ptr, fp);
    // Запись заголовка
    png_set_IHDR(png_ptr, info_ptr, Width, Height, depth,
                 PNG_COLOR_TYPE_PALETTE, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

    png_write_info(png_ptr, info_ptr);
    png_write_image(png_ptr, rows);
    png_write_end(png_ptr, info_ptr);

    png_destroy_write_struct(&png_ptr, &info_ptr);

    free(rows);
    png_free(png_ptr, palette);

    fclose(fp);

    return true;
}

/*!
   \brief Сохранение амплитудного РЛИ в *.jpg файл c предварительной
   нормировкой.

    \param[in]  fileName    Имя файла
    \param[in]  brightness  Коэффициент увеличения яркости изображения
    \param[in]  colormap    Цветовая схема

    \return      0 - успешное выполнение;
                -1 - ошибка
*/
int CImage::SaveToJpg8(char *fileName, float brightness, int quality) {
    return SaveToJpg8(fileName, brightness, quality, NULL);
}

/*!
   \brief Сохранение амплитудного РЛИ в *.jpg файл c предварительной
   нормировкой.

    \param[in]  fileName    Имя файла
    \param[in]  brightness  Коэффициент увеличения яркости изображения
    \param[in]  colormap    Цветовая схема
    \param[in]  label       Пометки

    \return      0 - успешное выполнение;
                -1 - ошибка
*/

int CImage::SaveToJpg8(char *fileName, float brightness, int quality,
                       char *label) {

    uint8_t *savebuff =
        (uint8_t *)malloc(ImageWidth * ImageHeight * sizeof(uint8_t));

    Prepocessing(savebuff, brightness, quality, label);

    SaveToJpg8_(fileName, ImageWidth, ImageHeight, savebuff, quality);

    free(savebuff);

    return 0;
}
/*

int CImage::CompressJpg8(float brightness, int quality, char *label){


    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;

    uint8_t* pData = Prepocessing(brightness, quality, label);

    JSAMPROW rowpointer[1];
    int row_stride = Width;

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);

    //jpeg_stdio_dest(&cinfo, fp);

    cinfo.image_width = Width;
    cinfo.image_height = Height;
    cinfo.input_components = 1;
    cinfo.in_color_space = JCS_GRAYSCALE;

    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, quality, TRUE);
    jpeg_start_compress(&cinfo, TRUE);

    while(cinfo.next_scanline < cinfo.image_height){
        rowpointer[0] = &pData[cinfo.next_scanline * row_stride];
        (void) jpeg_write_scanlines(&cinfo, rowpointer, 1);
    }

    jpeg_finish_compress(&cinfo);


    jpeg_destroy_compress(&cinfo);


}

*/
uint8_t *CImage::Prepocessing(uint8_t *savebuff, float brightness, int quality,
                              char *label) {
    // Log->LogPrintf("\t*Сохранение изображения...\t"); fflush(stdout);
    // time_start();
    Log->LogPrintf("\t*Сохранение изображения...\n");
    // Взятие абсолютного значения
    float *magbuff;
    float *pbuff;
    int len = ImageWidth * ImageHeight;

    if (ImageType == IT_FCOMPLEX) {
        magbuff = (float *)malloc(len * sizeof(float));
        memset(magbuff, 0, len * sizeof(float));
        magnitude((Ipp32fc *)pImageData, magbuff, ImageWidth, ImageHeight);
        pbuff = magbuff;
    } else {
        pbuff = (float *)pImageData;
    }

    float maxval = 0;
    float minval = pbuff[0];

    // Нахождение максимального значения
    // #pragma omp parallel for
    for (int i = 0; i < len; i++) {
        if (pbuff[i] > maxval)
            maxval = pbuff[i];
        if (pbuff[i] < minval)
            minval = pbuff[i];
    }
    //  ippsMax_32f(pbuff, len, &maxval);
    //  ippsMin_32f(pbuff, len, &minval);

    // Нормировка

    // savebuff = (uint8_t*)malloc(len * sizeof(uint8_t));

    for (int i = 0; i < len; i++) {
        // float val = (pbuff[i] - minval)/(maxval - minval)*255.0f*brightness;
        float val = (pbuff[i] - minval) / (maxval - minval) * 255.0f;
        pbuff[i] = val < 255.0f ? val : 255.0f;
    }

    if (brightness == 0.0f) {
        brightness = AutoBrigtness(pbuff);
    }

    for (int i = 0; i < len; i++) {
        float val = pbuff[i] * brightness;
        savebuff[i] = val < 255.0f ? (uint8_t)val : 255;
    }

    if (label) {
        AddLabel(savebuff, label);
    }

    // Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());

    return savebuff;
}

// #define powerX(r)  ( 1.0 / ( ( (4 * 3.14) * (4 * 3.14) ) * (r*r) ) )
// #define powerX(r)  ( 1.0 / ( (r*r) ) )
#define powerX(r) pow(r, 1.0 / 2)

int CImage::RangeBrightness(float x0, float dx) {

    Log->LogPrintf("\t*Яркостная коррекция по дальности...\n");

    float *pbuff = (float *)pImageData;
    int len = ImageWidth * ImageHeight;

    float p0 = powerX(x0);
    float px;

    for (int x = 0; x < ImageWidth; x++) {
        px = powerX(((x + 1) * dx));
        for (int y = 0; y < ImageHeight; y++) {

            pbuff[y * ImageWidth + x] *= (px / p0);
        }
    }

    return 0;
}

void CImage::JpegCompress(unsigned char *image, int width, int height,
                          int quality, unsigned long *jpegSize,
                          unsigned char **jpegBuf, CompressType compress_type) {
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    JSAMPROW row_pointer[1];
    int row_stride;

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);

    cinfo.image_width = width;
    cinfo.image_height = height;

    if (compress_type == STANDART) {
        cinfo.progressive_mode = 0;
    }

    cinfo.input_components = 1;
    cinfo.in_color_space = JCS_GRAYSCALE;
    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, quality, TRUE);

    jpeg_mem_dest(&cinfo, jpegBuf, jpegSize);
    jpeg_start_compress(&cinfo, TRUE);

    row_stride = width;

    while (cinfo.next_scanline < cinfo.image_height) {
        row_pointer[0] = &image[cinfo.next_scanline * row_stride];
        jpeg_write_scanlines(&cinfo, row_pointer, 1);
    }

    jpeg_finish_compress(&cinfo);
    jpeg_destroy_compress(&cinfo);
}

/*!
   \brief Сохранение РЛИ в *.jpg файл.

    \param[in]  fileName    Имя файла
    \param[in]  Width       Ширина изображения
    \param[in]  Height      Высота изображения
    \param[in]  pData       Указатель намассив данных
    \param[in]  colormap    Цветовая схема

    \return      0 - успешное выполнение;
                -1 - ошибка
*/

#define JPEG_HEADER_SIZE 20
int CImage::SaveToJpg8_(char *fileName, int Width, int Height, uint8_t *pData,
                        int quality) {
    FILE *fp = NULL;
    unsigned char *mem = NULL;
    unsigned long mem_size = 0;

    JpegCompress(pData, Width, Height, quality, &mem_size, &mem);

    fp = fopen(fileName, "wb");
    if (!fp) {
        printf("Can`t open file for writing: %s \n", fileName);
        return -1;
    }

    uint16_t meta_size = sizeof(metadata);
    uint16_t meta_size_swp =
        __builtin_bswap16(meta_size + 2); // metadata + (uint16_t)size
    uint16_t meta_marker = __builtin_bswap16(0xFFE1);

    metadata.crc =
        crc16_bitwise((uint8_t *)&metadata, meta_size - sizeof(metadata.crc));

    fwrite(mem, JPEG_HEADER_SIZE, 1, fp);
    fwrite(&meta_marker, 2, 1, fp);
    fwrite(&meta_size_swp, 2, 1, fp);
    fwrite(&metadata, meta_size, 1, fp);
    fwrite(mem + JPEG_HEADER_SIZE, mem_size - JPEG_HEADER_SIZE, 1, fp);

    fclose(fp);
    free(mem);

    return 0;
}

int CImage::SaveToJpg8(char *fileName, int Width, int Height, uint8_t *pData,
                       int quality) {
    FILE *fp = NULL;
    unsigned char *mem = NULL;
    unsigned long mem_size = 0;

    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;

    JSAMPROW rowpointer[1];
    int row_stride = Width;

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    // jpeg_stdio_dest(&cinfo, fp);
    jpeg_mem_dest(&cinfo, &mem, &mem_size);

    cinfo.image_width = Width;
    cinfo.image_height = Height;
    cinfo.input_components = 1;
    cinfo.in_color_space = JCS_GRAYSCALE;

    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, quality, TRUE);

    jpeg_simple_progression(&cinfo);
    jpeg_start_compress(&cinfo, TRUE);

    while (cinfo.next_scanline < cinfo.image_height) {
        rowpointer[0] = &pData[cinfo.next_scanline * row_stride];
        (void)jpeg_write_scanlines(&cinfo, rowpointer, 1);
    }

    jpeg_finish_compress(&cinfo);
    jpeg_destroy_compress(&cinfo);

    fp = fopen(fileName, "wb");
    if (!fp) {
        printf("Can`t open file for writing: %s \n", fileName);
        return -1;
    }

    fwrite(mem, mem_size, 1, fp);

    fclose(fp);

    return 0;
}
/*!
   \brief Преобразование оттенков серого в цветовую схему JET.

    \param[in]  v   Яркось пикселя


    \return     Значение пикселя rgb
*/
png_color CImage::GetJetColor(float v) {
    png_color c = {0, 0, 0};

    if (v < 0.125f) {
        c.red = 0;
        c.green = 0;
        c.blue = (4.0f * v + 0.5f) * 255;
    } else if (v < 0.375f) {
        c.red = 0;
        c.green = (4.0f * v - 0.5f) * 255;
        c.blue = 255;
    } else if (v < 0.625f) {
        c.red = (4.0f * v - 1.5f) * 255;
        c.green = 255;
        c.blue = (-4.0f * v + 2.5f) * 255;
    } else if (v < 0.875f) {
        c.red = 255;
        c.green = (-4.0f * v + 3.5f) * 255;
        c.blue = 0;
    } else {
        c.red = (-4.0f * v + 4.5f) * 255;
        c.green = 0;
        c.blue = 0;
    }
    return (c);
}

/*!
   \brief Инициализация экземпляра и выделение памяти

    \param[in]  width       Ширина изображения
    \param[in]  height      Высота изображения
    \param[in]  type        Тип данных изображения

    \return      0 - успешное выполнение;
                -1 - ошибка
*/
int CImage::Create(int width, int height, eImageType type) {
    // pImageData = malloc((width * height + 2) * GetPixelSize(type)); //2
    // комплексных отсчета добавляется после комплексного преобразования фурье
    pImageData =
        calloc((width * height + 2),
               GetPixelSize(type)); // 2 комплексных отсчета добавляется после
                                    // комплексного преобразования фурье
    
    memset(pImageData, 0, (width * height + 2) * GetPixelSize(type));

    ImageW0 = width;
    ImageH0 = height;
    SetImageWidth(width);
    SetImageHeight(height);
    SetImageType(type);

    return 0;
}

/*!
   \brief Прокрутка РЛИ

    \param[in]  n       Количество строк, на которое осуществляется сдвиг РЛИ

    \return      0 - успешное выполнение;
                -1 - ошибка
*/
int CImage::Roll(int n) {
    Log->LogPrintf("\t*Сдвиг изображения...\t\t");
    fflush(stdout);
    time_start();
    //  assert(ImageType == IT_FLOAT);

    int8_t *srcdata = (int8_t *)pImageData;
    int8_t *dstdata = (int8_t *)pImageData;

    for (int y = ImageHeight - 1; y >= n; y--) {
        int len = ImageWidth * GetPixelSize(ImageType);
        memcpy((void *)&dstdata[y * len], (void *)&srcdata[(y - n) * len], len);
    }

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());

    return 0;
}

/*
int CImage::CopyImageTop(CImage *srcimg) {
    // assert(srcimg->ImageType == IT_FCOMPLEX);

    Log->LogPrintf("\t*Копирование изображения...\t");
    fflush(stdout);
    time_start();

    void *srcdata = (void *)srcimg->pImageData;
    void *dstdata = (void *)pImageData;

    int len = srcimg->ImageHeight * srcimg->ImageWidth *
              GetPixelSize(srcimg->ImageType);
    memcpy(dstdata, srcdata, len);

    Log->LogPrintf("ОK. (%2.3f c)\n", time_stop());
    return 0;
}*/

int CImage::CopyImageTop(CImage *srcimg) {

    long int size_buffer = ImageHeight * ImageWidth * GetPixelSize(ImageType);
    
    if (actualIndex == -1){
        actualIndex = size_buffer; 
    }
    void *srcdata = (void *)srcimg->pImageData;
    void *dstdata = (void *)pImageData;         

    int len = srcimg->ImageHeight * srcimg->ImageWidth *
              GetPixelSize(srcimg->ImageType);
    
    int left_part_size = 0;
    int right_part_size = 0;
    
    

    if (len <= actualIndex ){
        memcpy(dstdata + actualIndex - len, srcdata, len);
        actualIndex = actualIndex - len;
    }
    else{

        right_part_size = actualIndex;
        left_part_size = len - right_part_size;
        
        memcpy(dstdata, srcdata + left_part_size, right_part_size);
        memcpy(dstdata + size_buffer - left_part_size, srcdata, left_part_size);
        actualIndex = size_buffer - left_part_size;
    }
    
    return 0;
}

#define THRESHOLD 12

float CImage::AutoBrigtnessFull(float *rawimg) {

    int len = ImageHeight * ImageWidth;
    int j = 0;
    float mean_rawimg = 0;
    float mean_buf = 0;

    for (int i = 0; i < len; i++) {

        if (rawimg[i] > THRESHOLD) {
            mean_buf += rawimg[i];
            j++;
        }

        mean_rawimg += rawimg[i];
    }

    float b = (float)(((float)mean_buf / (float)j) /
                      ((float)mean_rawimg / (float)len));

    Log->LogPrintf("\t*Яркость = %f\n", b);
    return b;
}

int radiation_pattern_edge(float x, float phi) {

    float r = x / cos(phi);
    float y = r * sin(phi);

    return (int)(y + 0.5);
}

float CImage::AutoBrigtness(float *rawimg) {
    if (metadata.div == 0.0f) {
        Log->LogPrintf("\t*Угол расхожения не задан (metadata)!\n");
        return AutoBrigtnessFull(rawimg);
    }
    if (metadata.dx == 0.0f) {
        Log->LogPrintf("\t*Шаг по дальности не задан (metadata)!\n");
        return AutoBrigtnessFull(rawimg);
    }
    if (metadata.dy == 0.0f) {
        Log->LogPrintf("\t*Шаг по путевой дальности не задан (metadata)!\n");
        return AutoBrigtnessFull(rawimg);
    }
    if (metadata.x0 == 0.0f) {
        Log->LogPrintf("\t*Начальная дальность не задана (metadata)!\n");
        return AutoBrigtnessFull(rawimg);
    }

    float div2 = (metadata.div / 2) * (M_PI / 180.0f);
    int j = 0;
    float mean_rawimg = 0;
    float mean_buf = 0;

    int len = 0;
    for (int x = 0; x < ImageWidth; x++) {
        float r = metadata.x0 + x * metadata.dx;
        float dy = metadata.dy;
        float y0 = ImageHeight * dy / 2; // metadata.y0;

        int b = (radiation_pattern_edge(r, -div2) + y0) / dy;
        int e = (radiation_pattern_edge(r, div2) + y0) / dy;
        if (b < 0) {
            b = 0;
        }
        if (e > ImageHeight) {
            e = ImageHeight;
        }

        for (int y = b; y < e; y++) {
            int index = y * ImageWidth + x;
            if (rawimg[index] > THRESHOLD) {
                mean_buf += rawimg[index];
                j++;
            }
            mean_rawimg += rawimg[index];
            len++;
        }
    }

    float b = (float)(((float)mean_buf / (float)j) /
                      ((float)mean_rawimg / (float)len));

    Log->LogPrintf("\t*Яркость = %f\n", b);
    return b;
}

/*!
   \brief Освобождение памяти, выделенной под изображение.

    \return     Функция не возвращает значения
*/
void CImage::Free(void) {
    if (pImageData) {
        free(pImageData);
        pImageData = NULL;
    }

    ImageWidth = 0;
    ImageHeight = 0;
    ImageType = IT_UNDEFINED;
    ImageNearRangeIndex = 0;
}

int CImage::GetWidth(void) { return ImageWidth; }

int CImage::GetHeight(void) { return ImageHeight; }

/*!
   \brief Функция вырезает часть изображения в соответствии с номером начальной
   и конечной строки

    \param[in]  nb  Номер начальной строки
    \param[in]  ne  Номер конечной строки

    \return     Функция не возвращает значения
*/
void CImage::CutH(int hb, int he) {
    assert(!(he > ImageHeight));
    assert(!(he < hb));

    // Выделяем память под новое изображение
    int len = (he - hb) * ImageWidth * GetPixelSize(ImageType);
    void *buf = malloc(len);

    // Копируем часть изображения в выделенный буфер
    memcpy(buf,
           &((int8_t *)pImageData)[hb * ImageWidth * GetPixelSize(ImageType)],
           len);

    // Переприсваиваем указатель на буфер
    free(pImageData);
    pImageData = buf;
    ImageHeight = he - hb;
    ImageNearRangeIndex = hb;
}

/*!
   \brief Функция вырезает часть изображения в соответствии с номером начального
   и конечного столбца

    \param[in]  nb  Номер начального стролбца
    \param[in]  ne  Номер конечного столбца

    \return     Функция не возвращает значения
*/
/*
void CImage::CutW(int wb, int we){
    if(pImageData){
        assert(!(we > ImageWidth));
        assert(!(we < wb));

        // Выделяем память под новое изображение
        int len  = (we - wb) * ImageHeight * GetPixelSize(ImageType);
        void* buf = malloc(len);

        // Копируем часть изображения в выделенный буфер
        for(int i=0; i<ImageHeight; i++){
            char* pdst = &((char*)buf)[i*(we - wb) * GetPixelSize(ImageType)];
            char* psrc = &((char*)pImageData)[(ImageWidth * i + wb) *
GetPixelSize(ImageType)]; memcpy(pdst, psrc, (we - wb) *
GetPixelSize(ImageType));
        }
        // Переприсваиваем указатель на буфер
        free(pImageData);
        pImageData = buf;
        ImageWidth = we - wb;
    }else{
#ifdef SUPPORT_AF
        afImageData = afImageData(af::seq(wb,we), af::span);
#endif
    }

    ImageNearRangeIndex = wb;
}
*/

void CImage::CutW(int wb, int we) {
    if (pImageData) {
        assert(!(we > ImageWidth));
        assert(!(we < wb));

        // Выделяем память под новое изображение
        int len = (we - wb) * ImageHeight * GetPixelSize(ImageType);
        int llen = len / ImageHeight;

        void *buf = malloc(len);

        // Копируем часть изображения в выделенный буфер
        for (int i = 0; i < ImageHeight; i++) {
            char *pdst = &((char *)buf)[i * llen];
            char *psrc = &((char *)pImageData)[(ImageWidth * i + wb) *
                                               GetPixelSize(ImageType)];
            memcpy(pdst, psrc, llen);
        }
        // Переприсваиваем указатель на буфер
        free(pImageData);
        pImageData = buf;
        ImageWidth = we - wb;
    } else {
#ifdef SUPPORT_AF
        afImageData = afImageData(af::seq(wb, we), af::span);
#endif
    }

    ImageNearRangeIndex = wb;
}

/*!
   \brief Установка пикселя на изображении

    \param[in]  w   Номер пикселя в строке
    \param[in]  h   Номер столбца пикселя

    \return     Функция не возвращает значения
*/
int CImage::SetPixel(int w, int h, float value) {
    int H = ImageHeight;
    int W = ImageWidth;
    Ipp32f *data = (Ipp32f *)pImageData;

    assert(ImageType == IT_FLOAT);
    assert(!(h > H));
    assert(!(w > W));

    data[h * W + w] = value;

    return 0;
}

/*!
   \brief Возвращает количество байт, занимаемое одним отсчетом РЛИ в памяти

    \param[in]  eImageType  тип РЛИ

    \return     Количество байт, занимаемое одним отсчетом РЛИ в памяти
*/
int CImage::GetPixelSize(eImageType type) {
    switch (type) {
    case IT_INT16:
        return sizeof(int16_t);
    case IT_FLOAT:
        return sizeof(Ipp32f);
    case IT_FCOMPLEX:
        return sizeof(Ipp32fc);
    case IT_UNDEFINED:
    default:
        // Log->LogPrintf("\e[31m CImage::GetPixelSize(). Неверный вормат РЛИ
        // \"%d\" \e[0m \n", type, "");
        return 0;
    }
}

pybind11::object CImage::GetData(int begin, int cnt) {
    float *data = (float *)pImageData;
    pyValues.clear();

    int len = ImageWidth * ImageHeight;
    if (cnt > len) {
        Log->LogPrintf(
            "\e[31m CImage::GetData(). cnt>len, cnt=%d, len=%d \e[0m \n", cnt,
            len, "");
        cnt = len;
    } else if (!len) {
        cnt = len;
    }

    for (int i = 0; i < cnt; i++) {
        pyValues.push_back((double)(data[begin + i]));
    }

    return wrap_double1d(pyValues.data(), pyValues.size());
}

pybind11::object CImage::GetImage(void) {
    float *data = (float *)pImageData;
    pyValues.clear();

    int len = ImageWidth * ImageHeight;

    for (int i = 0; i < len; i++) {
        pyValues.push_back((double)(data[i]));
    }

    return wrap_double2d(pyValues.data(), ImageWidth, ImageHeight);
}

void CImage::AddLabel(uint8_t *img, char *label) {

    uint8_t *data = img;
    uint8_t *c;
    int x = 0;
    int len = strlen(label);

    c = (uint8_t *)malloc(8 * 8);

    for (int i = 0; i < len; i++) {

        if (label[i] == '\n') {
            data += ImageWidth * 8;
            x = 0;
            continue;
        }

        Plot.RenderChar(c, (int)label[i]);
        for (int y = 0; y < 8; y++) {
            memcpy(data + (ImageWidth * y) + x, c + (8 * y), 8);
        }
        x += 8;
    }

    free(c);
}

void CImage::VerticalMirror() {

    int LineSize = GetPixelSize(ImageType) * ImageWidth;
    uint8_t *data = (uint8_t *)pImageData;
    uint8_t *buf = (uint8_t *)malloc(ImageHeight * LineSize);

    buf += ImageHeight * LineSize;

    for (int i = 0; i < ImageHeight; i++) {
        buf -= LineSize;
        memcpy(buf, data, LineSize);
        data += LineSize;
    }

    free(pImageData);
    pImageData = buf;
}

int CImage::Norm() {

    float *img = (float *)pImageData;

    float mean = 0;
    for (int i = 0; i < ImageHeight * ImageWidth; i++) {
        mean += img[i];
    }

    mean /= ImageHeight * ImageWidth;

    for (int i = 0; i < ImageHeight * ImageWidth; i++) {
        // img[i] /= log(img[i]);
        if (img[i] > mean) {
            img[i] = 0;
        }
    }

    return 0;
}

void CImage::Rotate(float angle) {

    float *newImg = (float *)malloc(ImageWidth * ImageHeight * sizeof(float));
    memset(newImg, 0, ImageWidth * ImageHeight * sizeof(float));

    float *data = (float *)pImageData;

    float b;

    // значения полуширины и полувысоты и изображений (исходного и итогового)
    // расшифровка имен
    // первая буква: h - half; вторая: w - width, h - height;
    // третья буква: d - dest (целевая картинка), s - source (исходная)
    int hwd = ImageWidth / 2;
    int hhd = ImageHeight / 2;
    int hws = ImageWidth / 2;
    int hhs = ImageHeight / 2;

    float r = sqrt(hws * hws + hhs * hhs);

    b = atan2(1.0 * hhs, hws);

    for (int i = 0 - hwd; i < ImageWidth - hwd; i++) {
        for (int j = 0 - hhd; j < ImageHeight - hhd; j++) {
            // I и J - координаты точек исходной картинки
            // правая часть формул (до r) нужна для поворота
            // картинки вокруг начальной точки (левый верхний угол)
            // слагаемые r*... нужны для того, чтобы вращение происходило
            // вокруг центра изображения

            int I = i * cos(angle) - j * sin(angle) + r * cos(b);
            int J = i * sin(angle) + j * cos(angle) + r * sin(b);
            // проверяем, не выходят ли точки за границы исходной картинки
            if (I < 2 * hws && I >= 0 && J < 2 * hhs && J >= 0) {
                newImg[(i + hwd) + (j + hhd) * ImageWidth] =
                    data[I + J * ImageWidth];
            }
        }
    }

    pImageData = newImg;
    free(data);
}

void CImage::PlotSignal(char *fileName) {

    Log->LogPrintf("\t*Сохранение 'осциллограммы' сигнала...\n");

    float *f = (float *)pImageData;

    Plot.PlotChart(f, ImageW0);
    SaveToJpg8_(fileName, Plot.Width, Plot.Height, Plot.pFrameBuffer, 80);
}

int CImage::ReadMetadata(char *fileName) {

    FILE *fp = fopen(fileName, "rb");
    if (!fp) {
        printf("Can`t open file for read metadata: %s \n", fileName);
        // fclose(fp);
        return -1;
    }

    // metadata.crc = crc16_bitwise( (uint8_t *)&metadata,
    // meta_size-sizeof(metadata.crc) );

    uint8_t jpeg_head[JPEG_HEADER_SIZE];
    uint16_t meta_marker;
    uint16_t meta_size = sizeof(metadata);
    uint16_t meta_size_rd = 0;

    fread(jpeg_head, JPEG_HEADER_SIZE, 1, fp);
    fread(&meta_marker, 2, 1, fp);

    if (meta_marker != __builtin_bswap16(0xFFE1)) {
        printf("Meta marker not found: %s \n", fileName);
        fclose(fp);
        return -2;
    }

    fread(&meta_size_rd, 2, 1, fp);
    meta_size_rd = __builtin_bswap16(meta_size_rd) - sizeof(uint16_t);

    fread(&metadata, meta_size, 1, fp);

    fclose(fp);

    return 0;
}
