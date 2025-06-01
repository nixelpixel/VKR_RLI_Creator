#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import time
import numpy
from termcolor import colored
import quasar
import pyquasar
import QuaSAR_misc


class StripNavigationParameters:
    def __init__(self, sar: quasar.SAR, args: pyquasar.args.CommandLineArgs, index_begin: int):
        nav_gpx = quasar.CNavGPX(args.paths[args.nt])
        gpx_time = nav_gpx.time()

        # Получение отсчетов времени для каждого периода зондирования
        interpolation_timestamps = numpy.linspace(
            0,
            gpx_time[-1],
            sar.radar.period_count(gpx_time[-1])
        )
        interpolation_timestamps_list = interpolation_timestamps.tolist()

        # Интерполяция навигационных координат
        elevation = nav_gpx.GetEleAtTime(interpolation_timestamps_list)
        speed = nav_gpx.GetSpeedAtTime(interpolation_timestamps_list)
        course = nav_gpx.GetDirectionAtTime(interpolation_timestamps_list)
        x = nav_gpx.GetXAtTime(interpolation_timestamps_list)
        y = nav_gpx.GetYAtTime(interpolation_timestamps_list)
        self.elevation = [x - args.altitude_shift for x in elevation]
        self.latitude = nav_gpx.GetLatAtTime(interpolation_timestamps_list)[index_begin]
        self.longitude = nav_gpx.GetLonAtTime(interpolation_timestamps_list)[index_begin]

        # Оценка траекторных нестабильностей
        self.course_initial = course[int(len(course) / 2)]
        self.course = [x - self.course_initial for x in course]
        azimuth = nav_gpx.GetAzimuthAtTime(interpolation_timestamps_list)
        azimuth_initial = azimuth[int(len(azimuth) / 2)]
        self.azimuth = [x - azimuth_initial for x in azimuth]
        self.x_difference = []
        y_lin = []  # Уравнение прямой относительно которой рассчитываются отклонения по х
        x_raw = nav_gpx.x()
        y_raw = nav_gpx.y()
        cx = x_raw[int(len(x_raw) / 2)]
        cy = y_raw[int(len(y_raw) / 2)]
        phi = cy / cx
        for i in range(len(x)):
            y_lin.append(phi * x[i])
            self.x_difference.append((cy * x[i] - cx * y[i]) / numpy.sqrt(cy * cy + cx * cx))
        self.speed = [s / 3.6 * numpy.cos(c) for s, c in zip(speed, self.course)]

        # Оценка степени отклонения от прямолинейной траектории
        #print_once = False
        #for xi in range(len(self.x_difference)):
            #if self.x_difference[xi] > args.x0:
                #self.x_difference[xi] = args.x0
                #if print_once is False:
                    #print(colored(" Предельное отклонение траектории. Увеличь ближнюю границу кадра!\t ".expandtabs(4), "yellow", attrs=["bold"]))
                    #print_once = True

        # Вывод графиков
        QuaSAR_misc.plot_nav(
            args.filename,
            interpolation_timestamps,
            self.elevation,
            speed,
            self.course,
            x,
            y,
            y_lin,
            self.x_difference,
            self.azimuth
        )


def load_and_estimate_drift_angle(sar: quasar.SAR, args: pyquasar.args.CommandLineArgs, pos: int) -> float:
    temp_image = sar.radar.read_image(sar.radar.period_count(args.synthesis_time), args.time_shift, pos, args.fic, False, True)
    window = quasar.Window(args.dsp, pyquasar.window_cast(args.window_function), temp_image.width)
    sar.dsp.window_weighting(temp_image, window)
    sar.dsp.range_compress(temp_image)
    temp_image.move_to(quasar.DSP.FFTW)         # ВРЕМЕННО (quasar.DSP.FFTW)
    sar.radar_mut.type = quasar.DSP.FFTW        # ВРЕМЕННО (quasar.DSP.FFTW)

    if numpy.isnan(args.drift_angle):
        return sar.driftAngle(temp_image, args.velocity)
    return numpy.deg2rad(args.drift_angle)


def main():
    # 1. Исходные данные
    trash_offset = 500
    vh_correction = True
    ang_cor = 3  # Коррекция ошибки курса, град

    # 2 Инициализация
    args = pyquasar.args.CommandLineArgs()
    logger = pyquasar.logger.init_default_logger(args.log_level)

    if args.strip_time < 0:
        args.strip_time = QuaSAR_misc.CalcTstrip(args.paths['data'], args.synthesis_time, args.time_shift)

    #nav_delay = QuaSAR_misc.GetNavDelay(args.paths['raw'])  # для старых записей используйте nav_delay = 2.5
    nav_delay = 2.5
    print(f"Задержка между навигацией и данными: {nav_delay:0.2f} c")

    sar = quasar.SAR(
        quasar.RadarSTT(
            args.paths['config'],
            args.paths['data'],
            args.dsp
        )
    )
    nav = quasar.CNav(args.paths['nav'])
    nav.read_sec_from(args.time_shift, args.synthesis_time)

    if args.velocity < 0:
        args.velocity = nav.velocity().mean()
    else:
        args.velocity /= 3.6
    if args.altitude < 0:
        args.altitude = nav.elevation().mean()
    args.altitude -= args.altitude_shift

    sar.dsp.print_description()

    # 3. Загрузка кусочка файла голограммы и поиск начала зондирования
    pos, width = sar.radar.sync_pos(args.time_shift, trash_offset)

    # 4. Определение угла сноса
    args.drift_angle = load_and_estimate_drift_angle(sar, args, pos)

    # 5. Подготовка к поблочному считыванию и формированию РЛИ
    period_count_per_block = sar.radar.period_count(args.dy / args.velocity)  # Количество периодов в блоке
    block_image = sar.radar.read_image(period_count_per_block, args.time_shift, pos, args.fic, False, True)
    bp = quasar.CBackProjection(sar, block_image)
    bp.stripInit(args.dx, args.lx, args.strip_time, args.synthesis_time)
    logger.level = quasar.LogLevel.Warn

    block_nav_index_begin = sar.radar.period_count(args.time_shift + nav_delay)                 # Индекс отсчета навигации, соответствующий началу блока
    strip_image_max_index = sar.radar.period_count(args.strip_time) + block_nav_index_begin     # Индекс последнего отсчета навигации, соответствующий концу полосового РЛИ
    nav_params = StripNavigationParameters(sar, args, index_begin=block_nav_index_begin)

    t_read = 0.0
    t_ww = 0.0
    t_rc = 0.0
    t_strip = 0.0

    sar.radar_mut.type = args.dsp  # ВРЕМЕННО
    window = quasar.Window(args.dsp, pyquasar.window_cast(args.window_function), block_image.width)
    block_count = int(sar.radar.period_count(args.strip_time) / period_count_per_block)
    pb = quasar.ProgressBar('Формирование РЛИ', block_count)
    for k in range(block_count):
        np = sar.radar.period_count(args.synthesis_time)
        if block_nav_index_begin + np >= len(nav_params.speed):
            print(colored(" Достигнут конец файла РГ\t ".expandtabs(4), "yellow", attrs=["bold"]))
            break
        if vh_correction:
            period_count_per_block = sar.radar.period_count(args.dy / nav_params.speed[block_nav_index_begin])
        t1 = time.time()
        block_chunk = sar.radar_mut.block_read(period_count_per_block, args.fic)
        # -- temporary --
        #block_chunk.save_signal_as_csv(f"{args.paths['image']}_block_{k}.csv")
        #QuaSAR_misc.PlotSignalFromCSV(f"{args.paths['image']}_block_{k}.csv", f"Block_{k}", "Signal", "Block no#", xlims=[], [-1, 1])
        # -- temporary --
        t2 = time.time()
        bp.set_image(block_chunk)
        t3 = time.time()
        sar.dsp.window_weighting(block_chunk, window)
        t4 = time.time()
        sar.dsp.range_compress(block_chunk)
        t5 = time.time()
        if vh_correction:
            speed_slice = nav_params.speed[block_nav_index_begin:block_nav_index_begin + np]
            ele_slice = nav_params.elevation[block_nav_index_begin:block_nav_index_begin + np]
            diffx_slice = nav_params.x_difference[block_nav_index_begin:block_nav_index_begin + np]
            block_nav_index_begin += period_count_per_block
            if block_nav_index_begin > strip_image_max_index:
                break
            t6 = time.time()
            bp.strip(args.x0, args.y0, k, args.drift_angle, speed_slice, ele_slice, diffx_slice)
            t7 = time.time()
        else:
            t6 = time.time()
            bp.strip(args.x0, args.y0, k, args.drift_angle, [args.velocity], [args.altitude], [0])
            t7 = time.time()

        t_read  += t2 - t1
        t_ww    += t4 - t3
        t_rc    += t5 - t4
        t_strip += t7 - t6

        pb.update(k)
    pb.finish()
    logger.level = args.log_level

    t_read  /= block_count
    t_ww    /= block_count
    t_rc    /= block_count
    t_strip /= block_count

    print("ZERO_COPY_READ")
    print('t_read :' + str(t_read))
    print('t_ww   :' + str(t_ww))
    print('t_rc   :' + str(t_rc))
    print('t_strip:' + str(t_strip))

    image = bp.retrieve_strip_image()
    image.VerticalMirror()
    if args.brightness_correction:
        image.RangeBrightness(args.x0, args.dx)
    args.populate_metadata(image)  # В полосовом режиме driftAngle в наземке не учитывается
    image.metadata_mut().latitude = nav_params.latitude
    image.metadata_mut().longitude = nav_params.longitude
    image.metadata_mut().azimuth = nav_params.course_initial + args.drift_angle + numpy.deg2rad(ang_cor)
    image.metadata_mut().image_type = quasar.ImageKind.Strip
    image.save(args.paths['image'], quasar.ImageFormat.JPG, args.brightness)

    args.paths['image'] = args.paths['image'] + '.jpg'

    pyquasar.send.send_file(args.paths['image'], args.remote)
    #pyquasar.send.send_file(args.paths['image'][:-4] + '-meta.json', args.remote)
    if args.detection:
        pyquasar.detection.neural_detect(args.paths['image'], args.remote)


if __name__ == "__main__":
    main()
