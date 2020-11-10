import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyteiser",
    version="0.0.1",
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