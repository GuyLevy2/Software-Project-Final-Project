from setuptools import setup, find_packages, Extension

setup(
    name='capi_demo1',
    version='0.1.0',
    author="Example Author",
    author_email="liad@example.com",
    description="A sample C-API",
    install_reqires=['invoke'],
    packages=find_packages(),


    license='GPL-2',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: Implementation :: CPython',
    ],
    ext_modules=[Extension(
            'capi_demo1',
            ['capi1.c'],
        ),
    ]
)