import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

REQUIRED_PACKAGES=['numba>=0.50.1', 'numpy>=1.19.1', 'pandas>=1.1.1']    

setuptools.setup(
    name="pyteiser",
    version="0.0.1",
    install_requires=REQUIRED_PACKAGES,
    author="Matvei Khoroshkin, Hani Goodarzi",
    author_email="khorms21@gmail.com",
    description="Python implementation of TEISER",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/goodarzilab/pyteiser",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    scripts = ['bin/pyteiser_pipeline'],
)