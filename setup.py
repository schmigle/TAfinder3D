from setuptools import setup, find_packages

setup(
    name='tafinder3d',
    version='1.0.0',
    author='Moshe Steyn',
    author_email='moshe.steyn@pnnl.gov',
    description='Toxin-antitoxin system identification with 3D structure-based searching',
    long_description=open('README.md').read() if os.path.exists('README.md') else '',
    long_description_content_type='text/markdown',
    url='https://github.com/schmigle/TAfinder3D',
    py_modules=[
        'config',
        'utils',
        'input_and_neighborhoods',
        'combined_search',
        'result_processing',
        'output_generation',
        'tafinder'
    ],
    entry_points={
        'console_scripts': [
            'tafinder=tafinder:main',
        ],
    },
    python_requires='>=3.9',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],
)
