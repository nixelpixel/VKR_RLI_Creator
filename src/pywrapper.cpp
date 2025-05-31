/*
 * pywrapper.cpp
 *
 *  Created on: 21 июн. 2017 г.
 *      Author: user
 */

#include "CRadar.hpp"
// #include "CRadarACC.hpp"
#include "CConfig.hpp"
#include "CImage.hpp"
#include "CLog.hpp"
#include "CRadarSTT.hpp"
#include "CSAR.hpp"
// #include "CSim.hpp"
#include "CBackProjection.hpp"
#include "CLinkedImage.hpp"
#include "CNav.hpp"
#include "CStripStream.hpp"
#include "base.hpp"

#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(libQuaSAR, m) {

    m.def("initQuaSAR", &initQuaSAR);

    py::enum_<CImage::ColorMap>(m, "ColorMap")
        .value("COLORMAP_GRAY", CImage::COLORMAP_GRAY)
        .value("COLORMAP_JET", CImage::COLORMAP_JET);

    py::enum_<CImage::eImageType>(m, "ImageType")
        .value("IT_INT16", CImage::IT_INT16)
        .value("IT_FLOAT", CImage::IT_FLOAT)
        .value("IT_FCOMPLEX", CImage::IT_FCOMPLEX);

    py::enum_<CRadar::RampType>(m, "RampType")
        .value("RAMP_SAW", CRadar::RAMP_SAW)
        .value("RAMP_TRI", CRadar::RAMP_TRI);

    py::enum_<CRadar::SampleSize>(m, "SampleSize")
        .value("SAMPLE_SIZE_1BIT", CRadar::SAMPLE_SIZE_1BIT)
        .value("SAMPLE_SIZE_2BIT", CRadar::SAMPLE_SIZE_2BIT)
        .value("SAMPLE_SIZE_4BIT", CRadar::SAMPLE_SIZE_4BIT)
        .value("SAMPLE_SIZE_8BIT", CRadar::SAMPLE_SIZE_8BIT)
        .value("SAMPLE_SIZE_16BIT", CRadar::SAMPLE_SIZE_16BIT);

    py::enum_<DSPType>(m, "DSPType")
        .value("DSP_IPP", DSP_IPP)
        .value("DSP_AF", DSP_AF)
        .value("DSP_OCL_CPU", DSP_OCL_CPU)
        .value("DSP_OCL_GPU", DSP_OCL_GPU)
        .value("DSP_FFTW", DSP_FFTW)
        .value("DSP_NEON", DSP_NEON)
        .value("DSP_CUDA", DSP_CUDA);

    py::enum_<CSAR::WindowType>(m, "WindowType")
        .value("HAMMING", CSAR::HAMMING)
        .value("BLACKMAN", CSAR::BLACKMAN)
        .value("NUTTALL", CSAR::NUTTALL)
        .value("NONE", CSAR::NONE);

    py::class_<CLog>(m, "CLog")
        .def(py::init<const char *>())
        .def("SetVerbose", &CLog::SetVerbose, py::arg("verbose"));

    py::class_<CImage>(m, "CImage")
        .def(py::init<CLog *>())

        .def("Create", &CImage::Create)

        .def("GetData", &CImage::GetData)
        .def("GetImage", &CImage::GetImage)

        .def("SaveToBin", &CImage::SaveToBin)
        .def("SaveToPng8", &CImage::SaveToPng8)
        .def("SaveToJpg8",
             static_cast<int (CImage::*)(char *, float, int, char *)>(
                 &CImage::SaveToJpg8))
        .def("SaveToJpg8", static_cast<int (CImage::*)(char *, float, int)>(
                               &CImage::SaveToJpg8))
        .def("RangeBrightness", &CImage::RangeBrightness)

        .def("LoadBin", &CImage::LoadBin)

        .def("Free", &CImage::Free)
        .def("CloneFrom", &CImage::CloneFrom)
        .def("CutH", &CImage::CutH)
        .def("CutW", &CImage::CutW)

        .def("SetPixel", &CImage::SetPixel)
        .def("PrepareForAF", &CImage::PrepareForAF)

        .def("Roll", &CImage::Roll)
        .def("GetWidth", &CImage::GetWidth)
        .def("GetHeight", &CImage::GetHeight)
        .def("CopyImageTop", &CImage::CopyImageTop)

        .def("VerticalMirror", &CImage::VerticalMirror)

        .def("Norm", &CImage::Norm)
        .def("Rotate", &CImage::Rotate)
        .def("PlotSignal", &CImage::PlotSignal)

        .def("ReadMetadata", &CImage::ReadMetadata)

        .def_readonly("ImageH0", &CImage::ImageH0)
        .def_readonly("ImageW0", &CImage::ImageW0)

        //.def_property("metadata", &CImage::metadata)
        .def_readwrite("metadata", &CImage::metadata)

        ;

    py::class_<CImage::metadata_t>(m, "metadata_t")
        .def_readwrite("lat", &CImage::metadata_t::lat)
        .def_readwrite("lon", &CImage::metadata_t::lon)
        .def_readwrite("dx", &CImage::metadata_t::dx)
        .def_readwrite("dy", &CImage::metadata_t::dy)
        .def_readwrite("x0", &CImage::metadata_t::x0)
        .def_readwrite("y0", &CImage::metadata_t::y0)
        .def_readwrite("ang", &CImage::metadata_t::ang)
        .def_readwrite("driftAngle", &CImage::metadata_t::driftAngle)
        .def_readwrite("lx", &CImage::metadata_t::lx)
        .def_readwrite("ly", &CImage::metadata_t::ly)
        .def_readwrite("div", &CImage::metadata_t::div)
        .def_readwrite("v", &CImage::metadata_t::v)
        .def_readwrite("h", &CImage::metadata_t::h)
        .def_readwrite("kr", &CImage::metadata_t::kr)
        .def_readwrite("t", &CImage::metadata_t::t)
        .def_readwrite("ts", &CImage::metadata_t::ts)
        .def_readwrite("mode", &CImage::metadata_t::mode)
        .def_readwrite("type", &CImage::metadata_t::type);

    py::class_<CRadar>(m, "CRadar")
        .def(py::init<CLog *>())
        .def("SetParams", &CRadar::SetParams)
        .def("SetParamsFromFile", &CRadar::SetParamsFromFile);

    py::class_<CRadar_STT, CRadar>(m, "CRadar_STT")
        .def(py::init<CLog *>())
        .def("SetParams", &CRadar_STT::SetParams)
        .def("ReadDataFromFile", &CRadar_STT::ReadDataFromFile)
        .def("ReadBlockFromFile", &CRadar_STT::ReadBlockFromFile)
        .def("ReadBlockSetPos", &CRadar_STT::ReadBlockSetPos)
        .def("ReadBlockResetPos", &CRadar_STT::ReadBlockResetPos)
        .def("FindSyncPos", &CRadar_STT::FindSyncPosPy)
        .def("GetFileDuration", &CRadar_STT::GetFileDuration)
        .def("Time2Np", &CRadar_STT::Time2Np)

        .def("OpenHologram", &CRadar_STT::OpenHologram)
        .def("FreeHologram", &CRadar_STT::FreeHologram)
        .def("ReadBlock", static_cast<int (CRadar_STT::*)(CImage *, int)>(
                              &CRadar_STT::ReadBlock))

        .def("RemoveTransients", &CRadar_STT::RemoveTransients)

        .def("OpenHologramSocket", &CRadar_STT::OpenHologramSocket)
        .def("ReadBlockFromSocket", &CRadar_STT::ReadBlockFromSocket)

        .def("InitDataFile", &CRadar_STT::InitDataFile)
        .def("ReadBlockFromInitializedFile",
             &CRadar_STT::ReadBlockFromInitializedFile);

    py::class_<CConfig>(m, "CConfig")
        .def(py::init<CLog *>())
        .def("GetValueInt", &CConfig::GetValueInt)
        .def("GetValueDouble", &CConfig::GetValueDouble);

    py::class_<CSAR>(m, "CSAR")
        .def(py::init<CLog *, CRadar *, DSPType>())
        .def("printDSPInfo", &CSAR::printDSPInfo)
        .def("removeDC", &CSAR::removeDC)
        .def("rangeCompress", &CSAR::rangeCompress)
        .def("transpose", &CSAR::transpose)
        .def("compensateShift", &CSAR::compensateShift)
        .def("azimuthCompress", &CSAR::azimuthCompress)
        .def("rmc", &CSAR::rmc)
        .def("fsrmc", &CSAR::fsrmc)
        .def("phocus", &CSAR::phocus)
        .def("enthropy", &CSAR::enthropy)
        .def("index2range", &CSAR::index2range)
        .def("range2index", &CSAR::range2index)
        .def("hamming", &CSAR::hamming)
        .def("windowWeighting", &CSAR::windowWeighting)

        .def("fr2range", &CSAR::fr2range)
        .def("interpolate", &CSAR::interpolate)
        .def("hilbert", &CSAR::hilbert)

        .def("dspProgress", &CSAR::dspProgress)

        .def("acc_dcm", &CSAR::acc_dcm)

        .def("createWindow",
             static_cast<uint64_t (CSAR::*)(int, CSAR::WindowType)>(
                 &CSAR::createWindow))
        .def("freeWindow", &CSAR::freeWindow)

        .def("resize", &CSAR::resize)
        .def("downrange", &CSAR::downrange)
        .def("angleCorrection", &CSAR::angleCorrection)
        //.def( "driftAngle", &CSAR::driftAngle, CSAR_overloads(args("image",
        //"V", "FnPlot"), "Высчисление угла сноса") )
        .def("driftAngle", &CSAR::driftAngle, "Высчисление угла сноса",
             py::arg("Image"), py::arg("V"), py::arg("FnPlot") = (nullptr))

        ;

    py::class_<CBackProjection>(m, "CBackProjection")
        .def(py::init<CSAR *, CImage *>())

        .def("setPixelSize", static_cast<void (CBackProjection::*)(float)>(
                                 &CBackProjection::setPixelSize))
        .def("setPixelSize",
             static_cast<void (CBackProjection::*)(float, float)>(
                 &CBackProjection::setPixelSize))
        .def("setVH",
             static_cast<int (CBackProjection::*)(py::list &, py::list &, int)>(
                 &CBackProjection::setVH))
        .def("setVH",
             static_cast<int (CBackProjection::*)(py::list &, py::list &)>(
                 &CBackProjection::setVH))

        .def("setProgressive", &CBackProjection::setProgressive)
        .def("telescopic", &CBackProjection::telescopic)
        .def("telescopicAngle", &CBackProjection::telescopicAngle)
        .def("phaseCorrection", &CBackProjection::phaseCorrection)

        .def("stripInit", &CBackProjection::stripInit)
        .def("strip", &CBackProjection::strip)
        .def("getStripImage", &CBackProjection::getStripImage)
        .def("autofocus", &CBackProjection::autofocus,
             "Автофокусировка на основе нахождения минимума энтропии",
             py::arg("v0"), py::arg("v1"), py::arg("point_x"),
             py::arg("point_y"), py::arg("Ts"), py::arg("ang") = (0.0),
             py::arg("ls") = (40.0), py::arg("n") = (20),
             py::arg("report") = (nullptr))

        ;

    py::enum_<CStripStream::DataFormat>(m, "DataFormat")
        .value("UINT8", CStripStream::UINT8)
        .value("UINT16", CStripStream::UINT16)
        .value("COMPRESSED", CStripStream::COMPRESSED);

    py::class_<CStripStream>(m, "CStripStream")
        .def(py::init<CSAR *, CImage *, float, float, float>())
        .def("strip", &CStripStream::strip)
        .def("set_format", &CStripStream::set_format,
             "Установить формат передаваемой посылки", py::arg("format"),
             py::arg("chunk_height") = 8)
        .def("set_compress_quality", &CStripStream::set_compress_quality)
        .def("pack_len", &CStripStream::pack_len)
        .def("set_nav_data", &CStripStream::set_nav_data,
             "Установка навигационных данных носителя", py::arg("lat"),
             py::arg("lon"), py::arg("velocity"), py::arg("ele"),
             py::arg("course"), py::arg("drift_ang"), py::arg("pitch") = (0.0),
             py::arg("roll") = (0.0))

        //.def("DataFormat", &CStripStream::DataFormat)
        ;

    py::class_<CNav>(m, "CNav")
        .def(py::init<CLog *, const char *>())
        .def("read_msec", &CNav::read_msec)
        .def("read_msec_from", &CNav::read_msec_from)
        .def("read_sec", &CNav::read_sec)
        .def("read_sec_from", &CNav::read_sec_from)
        .def("velocity", &CNav::velocity)
        .def("elevation", &CNav::elevation)
        .def("size", &CNav::size)
        .def("gyro", &CNav::gyro)
        .def("accel", &CNav::accel)
        .def("gyro", &CNav::compass)

        .def("lat", &CNav::lat, "Широта", py::arg("n") = (0))
        .def("lon", &CNav::lon, "Долгота", py::arg("n") = (0))
        .def("course", &CNav::course, "Курс", py::arg("n") = (0));

    py::class_<ParameterBase>(m, "ParameterBase").def(py::init<CNav *>())
        //.def( "mean", &ParameterBase::mean)
        ;

    py::class_<Velocity, ParameterBase>(m, "Velocity")
        .def(py::init<CNav *>())
        //.def( "mean_kmh", &Velocity::mean_kmh)
        .def("mean", &Velocity::mean)
        .def("as_kmh", &Velocity::as_kmh)
        .def("get", &Velocity::get);

    py::class_<Elevation, ParameterBase>(m, "Elevation")
        .def(py::init<CNav *>())
        .def("mean_above", &Elevation::mean_above)
        .def("mean", &Elevation::mean);

    py::class_<CLinkedImage>(m, "CLinkedImage")
        .def(py::init())
        .def("Create", &CLinkedImage::Create)
        .def("SaveToBin", &CLinkedImage::SaveToBin)
        .def_readonly("ImageH0", &CLinkedImage::ImageH0)
        .def_readonly("ImageW0", &CLinkedImage::ImageW0);

    py::class_<Dim3>(m, "Dim3")
        .def(py::init())
        .def_readonly("x", &Dim3::x)
        .def_readonly("y", &Dim3::y)
        .def_readonly("z", &Dim3::z);
}

/*
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(CSAR_overloads, CSAR::driftAngle, 2, 3)

BOOST_PYTHON_MODULE( libQuaSAR ){
//PYBIND11_MODULE( libQuaSAR ){




    class_<CRadar_ACC, bases<CRadar> >( "CRadar_ACC" )
            .def( init<CLog*>( args( "log" ) ) )
            .def( init<CLog*, const char*, const char*>( args( "log", "address",
"port" ) ) ) .def( "SetParams", &CRadar_ACC::SetParams, args( "f0", "df", "fs",
"ns", "nsync", "mp", "sample_size", "ramp_type" ) ) .def( "ReadDataFromSocket",
&CRadar_ACC::ReadDataFromSocket, args( "image1", "imqge2", "sync_width", "kr") )
            .def( "Time2Np", &CRadar_ACC::Time2Np, args( "ts" ) )
    ;





    class_<CSim>( "CSim" )
            .def( init<CLog*, CRadar*, DSPType>( args( "log", "radar",
"dsp_type" ) ) ) .def( "create", &CSim::create, args( "image", "Ts", "ImageType"
) ) .def( "writeToFile", &CSim::writeToFile, args( "image", "filename" ) ) .def(
"addPointTargetSignal", &CSim::addPointTargetSignal, args( "image", "V", "H",
"X", "Y", "A" ) ) .def( "addSync", &CSim::addSync, args( "image", "Nsync" ) )
    ;

}
*/
