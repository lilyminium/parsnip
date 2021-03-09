"""
parsnip
Build molecular topologies.
"""
import sys
import os
import platform
from distutils.ccompiler import new_compiler
from distutils.sysconfig import customize_compiler
from setuptools import setup, find_packages, Extension

try:
    from Cython.Build import cythonize
except ImportError:
    cython_found = False
    cython_linetrace = False
else:
    cython_found = True
    cython_linetrace = bool(os.environ.get('CYTHON_TRACE_NOGIL', False))
import versioneer
from subprocess import getoutput

is_release = False

short_description = __doc__.split("\n")

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except:
    long_description = "\n".join(short_description[2:])

def using_clang():
    """Will we be using a clang compiler?"""
    compiler = new_compiler()
    customize_compiler(compiler)
    compiler_ver = getoutput("{0} -v".format(compiler.compiler[0]))
    return 'clang' in compiler_ver

def get_rdkit_include():
    try:
        import rdkit
    except ImportError:
        print('*** package "rdkit" not found ***')
        sys.exit(-1)
    prefix = os.environ["CONDA_PREFIX"]
    return os.path.join(prefix, "include/rdkit")

def get_lib():
    prefix = os.environ["CONDA_PREFIX"]
    return os.path.join(prefix, "lib")

def get_numpy_include():
    # Obtain the numpy include directory. This logic works across numpy
    # versions.
    # setuptools forgets to unset numpy's setup flag and we get a crippled
    # version of it unless we do it ourselves.
    try:
        # Python 3 renamed the ``__builin__`` module into ``builtins``.
        # Here we import the python 2 or the python 3 version of the module
        # with the python 3 name. This could be done with ``six`` but that
        # module may not be installed at that point.
        import builtins
    except ImportError:
        import __builtin__ as builtins
    builtins.__NUMPY_SETUP__ = False
    try:
        import numpy as np
    except ImportError:
        print('*** package "numpy" not found ***')
        sys.exit(-1)
    return np.get_include()

# shamelessly copied from MDA.
def extensions(use_cython=True, debug_cflags=False):
    # dev installs must build their own cythonized files.

    include_dirs = [get_numpy_include(), get_rdkit_include()]
    cpp_extra_compile_args = ["-std=c++11"]
    cpp_extra_link_args = ["-L"+get_lib(), "-lRDKitGraphMol",
                            "-lRDKitMolAlign",
                            "-lRDKitSmilesParse",
                            "-lRDKitFileParsers",
                            "-lRDKitRDGeneral",
                            "-lRDKitForceFieldHelpers",
                            "-lRDKitMolStandardize"]

    define_macros = []
    if debug_cflags:
        cpp_extra_compile_args.extend(['-Wall', '-pedantic'])
        define_macros.extend([('DEBUG', '1')])

    # needed to specify c++ runtime library on OSX
    if platform.system() == 'Darwin' and using_clang():
        cpp_extra_compile_args.append('-stdlib=libc++')
        cpp_extra_compile_args.append('-mmacosx-version-min=10.9')
        cpp_extra_link_args.append('-stdlib=libc++')
        cpp_extra_link_args.append('-mmacosx-version-min=10.9')

    # Needed for large-file seeking under 32bit systems (for xtc/trr indexing
    # and access).
    largefile_macros = [
        ('_LARGEFILE_SOURCE', None),
        ('_LARGEFILE64_SOURCE', None),
        ('_FILE_OFFSET_BITS', '64')
    ]

    if use_cython:
        print('Will attempt to use Cython.')
        if not cython_found:
            print("Couldn't find a Cython installation. "
                  "Not recompiling cython extensions.")
            use_cython = False
    else:
        print('Will not attempt to use Cython.')

    cpp_source_suffix = '.pyx' if use_cython else '.cpp'

    # Windows automatically handles math library linking
    # and will not build MDAnalysis if we try to specify one
    if os.name == 'nt':
        mathlib = []
    else:
        mathlib = ['m']

    if cython_linetrace:
        cpp_extra_compile_args.append("-DCYTHON_TRACE_NOGIL")

    mon = Extension('parsnip.lib.pymonomer',
                    ['parsnip/lib/pymonomer' + cpp_source_suffix,
                     "parsnip/lib/src/utils.cpp",
                     "parsnip/lib/src/atom.cpp",
                     "parsnip/lib/src/tag.cpp",],
                    language="c++",
                    include_dirs=include_dirs + ['parsnip/lib/include'],
                    define_macros=define_macros,
                    extra_link_args= cpp_extra_link_args,
                    extra_compile_args=cpp_extra_compile_args)
    
    pol = Extension('parsnip.lib.pypolymer',
                    ['parsnip/lib/pypolymer' + cpp_source_suffix],
                    language="c++",
                    include_dirs=include_dirs + ['parsnip/lib/include'],
                    define_macros=define_macros,
                    extra_link_args= cpp_extra_link_args,
                    extra_compile_args=cpp_extra_compile_args)

    pre_exts = [mon]


    cython_generated = []
    if use_cython:
        extensions = cythonize(
            pre_exts,
            compiler_directives={'linetrace' : cython_linetrace,
                                 'embedsignature' : False},
        )
        if cython_linetrace:
            print("Cython coverage will be enabled")
        for pre_ext, post_ext in zip(pre_exts, extensions):
            for source in post_ext.sources:
                if source not in pre_ext.sources:
                    cython_generated.append(source)
    else:
        #Let's check early for missing .c files
        extensions = pre_exts
        for ext in extensions:
            for source in ext.sources:
                if not (os.path.isfile(source) and
                        os.access(source, os.R_OK)):
                    raise IOError("Source file '{}' not found. This might be "
                                "caused by a missing Cython install, or a "
                                "failed/disabled Cython build.".format(source))
    return extensions, cython_generated




if __name__ == "__main__":
    exts, cythonfiles = extensions(use_cython=not is_release)


    setup(
        # Self-descriptive entries which should always be present
        name='parsnip',
        author='Lily Wang',
        author_email='lily.wang@anu.edu.au',
        description=short_description[0],
        long_description=long_description,
        long_description_content_type="text/markdown",
        version=versioneer.get_version(),
        cmdclass=versioneer.get_cmdclass(),
        license='LGPLv3',
        ext_modules = exts,

        # Which Python importable modules should be included when your package is installed
        # Handled automatically by setuptools. Use 'exclude' to prevent some specific
        # subpackage(s) from being added, if needed
        packages=find_packages(),

        # Optional include package data to ship with your package
        # Customize MANIFEST.in if the general case does not suit your needs
        # Comment out this line to prevent the files from being packaged with your software
        include_package_data=True,

        # Allows `setup.py test` to work correctly with pytest
        setup_requires=[] + pytest_runner,

        # Additional entries you may want simply uncomment the lines you want and fill in the data
        # url='http://www.my_package.com',  # Website
        # install_requires=[],              # Required packages, pulls from pip if needed; do not use for Conda deployment
        # platforms=['Linux',
        #            'Mac OS-X',
        #            'Unix',
        #            'Windows'],            # Valid platforms your code works on, adjust to your flavor
        # python_requires=">=3.5",          # Python version restrictions

        # Manual control if final package is compressible or not, set False to prevent the .egg from being made
        # zip_safe=False,

    )

    # if not config.get('keep_cythonized', default=is_release) and not cython_linetrace:
    #     for cythonized in cythonfiles:
    #         try:
    #             os.unlink(cythonized)
    #         except OSError as err:
    #             print("Warning: failed to delete cythonized file {0}: {1}. "
    #                 "Moving on.".format(cythonized, err.strerror))
