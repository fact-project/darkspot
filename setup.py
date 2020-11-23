from setuptools import setup


setup(
    name='darkspot',
    version='0.3.0',
    description='Find a dark spot in the sky e. g. for a Ratescan',
    url='http://github.com/fact-project/darkspot',
    author='Maximilian Noethe',
    author_email='maximilian.noethe@tu-dortmund.de',
    license='MIT',
    packages=['darkspot'],
    package_data={'darkspot': ['hipparcos_vmag_less_10.tsv']},
    entry_points={
        'console_scripts': ['find_darkspot = darkspot.__main__:main']
    },
    python_requires='>= 3.6',
    install_requires=[
        'numpy',
        'ephem',
        'tqdm',
        'click',
        'astropy>=4',
    ],
    extras_require={
        'plot': [
            'matplotlib',
            'cartopy',
        ]
    },
    zip_safe=False,
)
