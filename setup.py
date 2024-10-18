from setuptools import find_packages, setup

version = "0.4.1"

# Dynamically read dependencies from requirements file
with open("requirements.txt") as f:
    requirements = f.readlines()

if __name__ == "__main__":
    setup(version=version, install_requires=requirements, packages=find_packages())
