from setuptools import setup, find_packages
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

with (this_directory / "antares" / "version.py").open() as f:
    version = f.read().split(" = '")[1].split("'\n")[0]


setup(
    name='antares_hep',
    version=version,
    license='GNU General Public License v3.0',
    description='Automated Numerical To Analytical REconstruction Software',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Giuseppe De Laurentis, Daniel Maitre',
    author_email='g.dl@hotmail.it',
    url='https://github.com/GDeLaurentis/antares',
    download_url=f'https://github.com/GDeLaurentis/antares/archive/v{version}.tar.gz',
    project_urls={
        'Documentation': 'https://gdelaurentis.github.io/antares/',
        'Issues': 'https://github.com/GDeLaurentis/antares/issues',
    },
    keywords=['Analytic Reconstruction', 'Spinor Helicity', 'Scattering Amplitudes', 'QFT'],
    packages=find_packages(),
    include_package_data=True,
    data_files=[],
    install_requires=[
        'lips>=0.5.0',
        'pyadic>=0.2.4',
        'syngular>=0.5.0',
        # 'linac',
        'pyyaml',
        'pandas',
        'multiset',
        # 'ortools',
    ],
    entry_points={
        'console_scripts': [
            'SpinorLatexCompiler=antares.scripts.SpinorLatexCompiler:main',  # Define the entry point
        ],
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
    ],
)
