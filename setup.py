import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="py2vyu",
    version="0.1.0",
    author="Shohan Hasan",
    author_email="shohan.hasan@nyu.edu",
    description="Python2 backport of pyvyu",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    python_requires='>=2.7',
)
