# Путь к CUDA 12.6
CUDA_HOME := /usr/local/cuda-12.6
NVCC := $(CUDA_HOME)/bin/nvcc

# Библиотеки и пути к include
LIBS := -lpng16 -lcufft -lcuda -lgomp -lfftw3f_threads -ljpeg -lfftw3f -lconfig

INCLUDES := \
-I../src/pybind11_legacy \
$(shell pkg-config --cflags-only-I python3) \
-I$(shell python3 -c "import numpy; print(numpy.get_include())") \
$(shell pkg-config --cflags-only-I libpng16)

GENCODE := -gencode arch=compute_89,code=sm_89

# Исходные файлы
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
../src/CDSPCUDA.cu \
../src/CDSPFFTW.cpp \
../src/CStripStream.cpp \
../src/CNav.cpp \
../src/CLinkedImage.cpp \
../src/CTCPServer.cpp \
../src/pywrapper.cpp

# Object файлы
OBJS := $(CPP_SRCS:../src/%.cpp=./src/%.o)
OBJS := $(OBJS:../src/%.cu=./src/%.o)

# Dependency файлы
CPP_DEPS := $(OBJS:.o=.d)

# Компиляция .cpp
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	$(NVCC) -allow-unsupported-compiler -std=c++14 -DQUASAR_CUDA $(INCLUDES) -O3 --use_fast_math \
		-Xcompiler -fPIC -Xcompiler -fopenmp -Xcompiler -ffast-math \
		$(GENCODE) -M -o "$(@:%.o=%.d)" "$<"
	$(NVCC) -allow-unsupported-compiler -std=c++14 -DQUASAR_CUDA $(INCLUDES) -O3 --use_fast_math \
		-Xcompiler -fPIC -Xcompiler -fopenmp -Xcompiler -ffast-math \
		--compile -x c++ -o "$@" "$<"
	@echo 'Finished building: $<'

# Компиляция .cu
src/%.o: ../src/%.cu
	@echo 'Building file: $<'
	$(NVCC) -allow-unsupported-compiler -DQUASAR_CUDA $(INCLUDES) -O3 --use_fast_math \
		-Xcompiler -fPIC -Xcompiler -fopenmp -Xcompiler -ffast-math \
		$(GENCODE) -M -o "$(@:%.o=%.d)" "$<"
	$(NVCC) -allow-unsupported-compiler -DQUASAR_CUDA $(INCLUDES) -O3 --use_fast_math \
		-Xcompiler -fPIC -Xcompiler -fopenmp -Xcompiler -ffast-math \
		--compile --relocatable-device-code=false $(GENCODE) -x cu -o "$@" "$<"
	@echo 'Finished building: $<'
