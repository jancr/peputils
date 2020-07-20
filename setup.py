import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()
    setuptools.setup(
        name='peputils',
        version='0.1',
        scripts=[],
        author="Jan Christian Refsgaard",
        author_email="jcr@intomics.com",
        description="Peptide Utility Scripts",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/jancr/peputils",
        #  packages=setuptools.find_packages(),
        packages=['peputils'],
        classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
        ],
    )
