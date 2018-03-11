from setuptools import setup

setup(
    name='open_metal_detector',
    version='2018.12.3',
    packages=['omsdetector'],
    install_requires=["pandas>=0.22.0",
                      "matplotlib>=2.1.2",
                      "numpy>=1.14.1",
                      "pymatgen>=2018.2.13"],
    url='',
    license='MIT',
    author='Emmanuel Haldoupis',
    author_email='emmhald@gmail.com',
    description='A tool to detect open metal sites in collections of MOFs.',
    python_requires='>=3'
)
