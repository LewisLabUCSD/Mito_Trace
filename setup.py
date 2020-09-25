from setuptools import find_packages, setup

setup(
    name='src',
    packages=find_packages(),
    version='0.1.0',
    description='Lineage Tracing with mitochondrion reads',
    author='Isaac Shamie',
    license='MIT',
    include_package_data=True,
    install_requires=[
            'Click',
        ],
    entry_points='''
            [console_scripts]
            bam_barcodes=src.bam_barcodes_function:main
        '''
)
