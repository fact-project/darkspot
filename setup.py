from setuptools import setup


setup(
    name='darkspot',
    version='0.2',
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
    install_requires=[
        'numpy',
        'ephem',
        'pandas',
        'tqdm',
        'docopt',
        'blessings',
        'numexpr',
    ],
    zip_safe=False,
)
