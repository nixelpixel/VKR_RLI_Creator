#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import numpy
from termcolor import colored
import quasar
import pyquasar.logger
import pyquasar.args
import QuaSAR_misc


def main():
    args = pyquasar.args.CommandLineArgs()
    logger = pyquasar.logger.init_default_logger(args.log_level)
    args.altitude -= args.altitude_shift
    
    nav_delay = QuaSAR_misc.GetNavDelay(args.paths['raw'])  # для старых записей используйте nav_delay = 2.5
    #nav_delay = 2.5  # Задержка между стартом записи навигации и записью голограммы
    trash_offset = 500

    # Инициализация
    sar = quasar.SAR(
        quasar.RadarSTT(
            args.paths['config'],
            args.paths['data'],
            args.dsp
        )
    )

    # Считывание навигационных данных
    nav = quasar.CNav(args.paths['nav'])
    nav.read_sec_from(args.time_shift + nav_delay, args.synthesis_time)

    if args.velocity < 0:
        args.velocity = nav.velocity().mean()
    else:
        args.velocity = args.velocity / 3.6

    if args.altitude < 0:
        args.altitude = nav.elevation().mean()

    print(f"Velocity:  {args.velocity:.1f} m/s")
    print(f"Elevation: {args.altitude:.1f} m\n")

    sar.dsp.print_description()

    # 0. Загрузка кусочка файла голограммы и поиск начала зондирования
    pos, width = sar.radar.sync_pos(args.time_shift, trash_offset)

    # 1. Загрузка файла голограммы в соответствии с заданным args.synthesis_time и найденной позицией начала зондирования
    image = sar.radar.read_image(sar.radar.period_count(args.synthesis_time), args.time_shift, pos, args.fic, False, True)

    # 1.5 Сохранение осциллограммы сигнала
    filename = args.paths['image'] + '-s'
    image.save_signal_as_csv(filename + '.csv')
    QuaSAR_misc.PlotSignalFromCSV(filename, "Сигнал", "Номер отсчета", "Амплитуда", [], [-1, 1])

    # 2. Сжатие по дальности
    window = quasar.Window(args.dsp, pyquasar.window_cast(args.window_function), image.width)
    sar.dsp.window_weighting(image, window)
    sar.dsp.range_compress(image)
    image.move_to(quasar.DSP.FFTW)  # ВРЕМЕННО!

    filename = args.paths['image'] + '-d'
    drift_angle = sar.driftAngle(image, args.velocity, filename + '.csv')
    QuaSAR_misc.PlotSignalFromCSV(filename, f"Угол сноса = {numpy.rad2deg(drift_angle):.2f} град.", "Угол, град.", "Амплитуда", [], [])

    if not drift_angle:
        drift_angle = 0

    # 3. Выполнение bp
    image.move_to(args.dsp)  # ВРЕМЕННО!
    bp = quasar.CBackProjection(sar, image)
    bp.setPixelSize(args.dx, args.dy)
    bp.telescopic(args.x0, args.lx, drift_angle, numpy.deg2rad(args.div), [args.velocity], [args.altitude], [0])

    mid = int(nav.size() / 2)
    lat = nav.lat(mid)
    lon = nav.lon(mid)
    course = nav.course(mid)

    if args.brightness_correction:
        image.RangeBrightness(args.x0, args.dx)

    args.populate_metadata(image)
    image.metadata_mut().latitude = numpy.rad2deg(lat)
    image.metadata_mut().longitude = numpy.rad2deg(lon)
    image.metadata_mut().azimuth = course
    image.metadata_mut().image_type = quasar.ImageKind.Telescopic
    image.save(args.paths['image'], quasar.ImageFormat.JPG, args.brightness)

    args.paths['image'] = args.paths['image'] + '.jpg'
    pyquasar.send.send_file(args.paths['image'], args.remote)
    #pyquasar.send.send_file(args.paths['image'][:-4] + '-meta.json', args.remote)
    if args.detection:
        pyquasar.detection.neural_detect(args.paths['image'], args.remote)


if __name__ == '__main__':
    main()
