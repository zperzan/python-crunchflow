from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="crunchflow", 
    version="1.0.2",
    author="Zach Perzan",
    author_email="zperzan@stanford.edu",
    description="Python toolbox for working with the CrunchFlow reactive transport model",
    license='GPL',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/zperzan/python-crunchflow",
    packages=find_packages(),
    install_requires=[
        'numpy',
        'matplotlib',
        'pandas',
        'shlex',
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)"
    ],
    python_requires='>=3.6',
    include_package_data=True,
    zip_safe=False,
)
