# -*- coding: utf-8 -*-
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import amas

with open('README.md') as readme_file:
    readme = readme_file.read()

setup(
    name='amas',
    version=amas.__version__,
    description="Calculate various summary statistics on a multiple sequence alignment",
    long_description=readme,
    author=amas.__author__,
    author_email=amas.__email__,
    url='https://github.com/marekborowiec/AMAS',
    packages=[
        'amas',
    ],
    package_dir={'amas':
                 'amas'},
    include_package_data=True,
    license="GNU PL v3",
    zip_safe=False,
    keywords='amas',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.3',
    ],
    test_suite='tests',
)
