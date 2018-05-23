from distutils.core import setup, Extension

c_ext = Extension("ads1256", ["wrapper.c", "RadiometerRun.c"], libraries = ['bcm2835'])

setup(
    ext_modules=[c_ext],
)
