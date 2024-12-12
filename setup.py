from setuptools import setup, find_packages


setup(
    name='antares',
    version='0.1.0',
    author='Giuseppe De Laurentis',
    author_email='g.dl@hotmail.it',
    description='Automated Numerical To Analytical Rational coefficients Extraction for Spinor helicity scattering amplitudes',
    packages=find_packages(),
    include_package_data=True,
    data_files=[],
    install_requires=[
        'lips',
        'pyadic',
        'syngular',
        # 'linac',
        'pyyaml',
        'pandas',
        'multiset',
        'ortools',
    ],
    entry_points={
        'console_scripts': [
            'SpinorLatexCompiler=antares.scripts.SpinorLatexCompiler:main',  # Define the entry point
        ],
    },
)
