#ifndef VECTOR_TO_NPARRAY
#define VECTOR_TO_NPARRAY

// Отключить устаревшие API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

// Подключаем Python.h до NumPy!
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include <numpy/ndarrayobject.h>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <vector>

// wrap c++ vector as numpy array
static pybind11::array_t<double> wrap_double1d(double* data, npy_intp size) {
    return pybind11::array_t<double>(
        {size}, // shape
        {sizeof(double)}, // C-style contiguous strides for double
        data); // numpy array references this parent
}

static pybind11::array_t<float> wrap_float1d(float* data, npy_intp size) {
    return pybind11::array_t<float>(
        {size}, // shape
        {sizeof(float)}, // C-style contiguous strides for float
        data); // numpy array references this parent
}

static pybind11::array_t<double> wrap_double2d(double* data, npy_intp sizex, npy_intp sizey) {
    return pybind11::array_t<double>(
        {sizey, sizex}, // shape
        {sizeof(double), sizeof(double)}, // C-style contiguous strides
        data); // numpy array references this parent
}

#endif // VECTOR_TO_NPARRAY
