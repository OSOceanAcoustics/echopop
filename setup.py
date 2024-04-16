from setuptools import find_packages, setup

# Dynamically read dependencies from requirements file
with open("requirements.txt") as f:
    requirements = f.readlines()

if __name__ == "__main__":
    setup(install_requires=requirements, packages=find_packages())
