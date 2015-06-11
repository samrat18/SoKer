from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy
 
extensions=[
    Extension("tensorproduct",
            ["tensorproduct.pyx"],            
            include_dirs=[numpy.get_include()],
            extra_compile_args=["-w"]
            ),
    Extension("PLegendre",
            ["PLegendre.pyx"],
            include_dirs=[numpy.get_include()],
            extra_compile_args=["-w"]
            ),
    Extension("mi_cython",
            ["mi_cython.pyx"],
            include_dirs=[numpy.get_include()]
            ),
    Extension("power_spectrum",
			["power_spectrum.pyx"],
			include_dirs=[numpy.get_include()],
			extra_compile_args=["-w"]
			)
]

setup(
    ext_modules=cythonize(extensions),
)
