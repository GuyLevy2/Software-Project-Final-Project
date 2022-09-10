from setuptools import setup, find_packages, Extension

setup(
    name='spkmeansmodule',
    version='0.1.0',
    author="Liad and Guy",
    author_email="sample@example.com",
    description="C-API Extansion",
    install_requires=['invoke'],
    packages=find_packages(),

    license='GPL-2',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: Implementation :: CPython',
    ],
    ext_modules=[
        Extension(
            'spkmeansmodule',
            ['spkmeansmodule.c', 'spkmeans.c'],
        ),
    ]
)