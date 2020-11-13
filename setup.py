import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="JAS",
    version="0.1.4",
    author="Juliane Oliveira and Moreno Rodrigues",
    author_email="julianlanzin@gmail.com, rodriguesmsb@gmail.com",
    description="A python package that could be used to  run compartimental models on epidemic data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
