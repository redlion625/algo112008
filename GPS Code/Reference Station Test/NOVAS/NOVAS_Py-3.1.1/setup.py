# -*- coding: utf-8 -*-

import os
import platform
from distutils.core import setup, Extension
from distutils.core import Command
from distutils.command.build_clib import build_clib
from distutils.sysconfig import get_python_lib
from distutils import log

class build_dynamic_clib(build_clib):
    def finalize_options (self):
        self.set_undefined_options('build',
                                   ('build_lib', 'build_clib'),
                                   ('build_temp', 'build_temp'),
                                   ('compiler', 'compiler'),
                                   ('debug', 'debug'),
                                   ('force', 'force'))

        self.libraries = self.distribution.libraries
        if self.libraries:
            self.check_library_list(self.libraries)

        if self.include_dirs is None:
            self.include_dirs = self.distribution.include_dirs or []
        if isinstance(self.include_dirs, str):
            self.include_dirs = self.include_dirs.split(os.pathsep)

    def build_libraries (self, libraries):
        for (lib_name, build_info) in libraries:
            sources = build_info.get('sources')
            if sources is None or not isinstance(sources, (list, tuple)):
                raise DistutilsSetupError(
                       "in 'libraries' option (library '%s'), "
                       "'sources' must be present and must be "
                       "a list of source filenames" % lib_name)
            sources = list(sources)

            log.info("building '%s' library", lib_name)

            macros = build_info.get('macros')
            include_dirs = build_info.get('include_dirs')
            objects = self.compiler.compile(sources,
                output_dir=self.build_temp,
                macros=macros,
                include_dirs=include_dirs,
                extra_postargs=build_info.get('extra_compile_args', []),
                debug=self.debug)

            package = build_info.get('package', '')
            self.compiler.link_shared_lib(
                objects, lib_name,
                output_dir=os.path.join(self.build_clib, package),
                extra_postargs=build_info.get('extra_link_args', []),
                debug=self.debug,)

    def run(self):
        log.info('running build_dynamic_clib')
        build_clib.run(self)

c_sources = [
    'Cdist/solsys1.c',
    'Cdist/readeph0.c',
    'Cdist/eph_manager.c',
    'Cdist/nutation.c',
    'Cdist/novascon.c',
    'Cdist/novas.c'
]

calling_dir = None

def main():
    global calling_dir
    calling_dir = os.getcwd()

    cwd = os.path.dirname(__file__)
    if cwd:
        os.chdir(cwd)

    system = platform.system().lower()

    if 'darwin' in system:
        novaslib = [(
            'novas', {
                'package': 'novas',
                'sources': c_sources,
                'include_dirs': ['Cdist'],
                'extra_compile_args': ['-arch', 'i386',
                                       '-arch', 'x86_64',
                                       '-O2', '-Wall', '-fPIC']
            }
        )]
    elif 'windows' in system:
        raise OSError("Operating system not supported at this time")
    else:
        novaslib = [(
            'novas', {
                'package': 'novas',
                'sources': c_sources,
                'include_dirs': ['Cdist'],
                'extra_compile_args': ['-O2', '-Wall', '-fPIC']
            }
        )]

    options = {
        'name': 'NOVAS_Py',
        'version': '3.1.1',
        'description': "Python wrappers for the US Naval Observatory's \
                        NOVAS-C package.",
        'author': 'Eric G. Barron',
        'author_email': "%(firstdotlast)s@%(place)s" %
                        {'firstdotlast': 'eric.barron',
                         'place': 'usno.navy.mil'},
        'maintainer': 'Eric G. Barron',
        'maintainer_email': "%(firstdotlast)s@%(place)s" %
                            {'firstdotlast': 'eric.barron',
                             'place': 'usno.navy.mil'},
        'url': 'http://www.usno.navy.mil/USNO/astronomical-applications/software-products/novas',
        'download_url' : 'http://www.usno.navy.mil/USNO/astronomical-applications/software-products/novas',
        'platforms': ['macosx', 'linux'],
        'classifiers': [
            'Development Status :: 5 - Production/Stable',
            'Intended Audience :: Science/Research',
            'Natural Language :: English',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: POSIX :: Linux',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3.3',
            'Programming Language :: Python :: 3.4',
            'Topic :: Scientific/Engineering :: Astronomy',
        ],

        'packages': ['novas', 'novas.compat'],
        'package_dir': {
            'novas': 'novas_py',
            'novas.compat': 'compat'
        },
        'libraries': novaslib,
        'cmdclass': {
            'build_clib': build_dynamic_clib,
        }
    }

    setup(**options)

if __name__ == '__main__':
    main()
