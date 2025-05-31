LIBS := -lpng16 -lgomp -lfftw3f_threads -ljpeg -lfftw3f -lconfig

INCLUDES := \
-I/home/dmitry/vkr/quasar-dsp_refactor/venv/lib/python3.13/site-packages/pybind11/include \
$(shell pkg-config --cflags-only-I python3)\
$(shell pkg-config --cflags-only-I libpng16)

CPP_SRCS += \
../src/CLog.cpp \
../src/helpers.cpp \
../src/CPlot.cpp \
../src/CImage.cpp \
../src/CRadar.cpp \
../src/CRadarSTT.cpp \
../src/CConfig.cpp \
../src/CSAR.cpp \
../src/CBackProjection.cpp \
../src/CDSP.cpp \
../src/CDSPFFTW.cpp \
../src/CStripStream.cpp \
../src/CNav.cpp \
../src/CLinkedImage.cpp \
../src/CTCPServer.cpp \
../src/pywrapper.cpp 

OBJS += \
./src/CLog.o \
./src/helpers.o \
./src/CPlot.o \
./src/CImage.o \
./src/CRadar.o \
./src/CRadarSTT.o \
./src/CConfig.o \
./src/CSAR.o \
./src/CBackProjection.o \
./src/CDSP.o \
./src/CDSPFFTW.o \
./src/CStripStream.o \
./src/CNav.o \
./src/CLinkedImage.o \
./src/CTCPServer.o \
./src/pywrapper.o 

CPP_DEPS += \
./src/CLog.d \
./src/helpers.d \
./src/CPlot.d \
./src/CImage.d \
./src/CRadar.d \
./src/CRadarSTT.d \
./src/CConfig.d \
./src/CSAR.d \
./src/CBackProjection.d \
./src/CDSP.d \
./src/CDSPFFTW.d \
./src/CStripStream.d \
./src/CNav.d \
./src/CLinkedImage.d \
./src/CTCPServer.d \
./src/pywrapper.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++14 -DQUASAR_FFTW $(INCLUDES) -O3 -g3 -Wall -c -fmessage-length=0 -fPIC -fopenmp -ffast-math -fabi-version=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


