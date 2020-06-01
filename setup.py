# coding: utf-8

from setuptools import setup

setup(
    name='PyQTST',
    version='2.0.4',
    author='linqiaosong',
    author_email='linqiaosong@outlook.com',
    url='https://github.com/Linqiaosong/PyQTST',
    description=u'Software package for calculating chemical reaction rate constant by using transition state theory',
    packages=['PyQTST'],
    install_requires=[
        'numpy>=1.15.0',
        'matplotlib>=3.0.0'       
        ],
    entry_points={
        'console_scripts': [
            
        ]
    }
)