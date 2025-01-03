from setuptools import setup

# Dynamically read dependencies from requirements file
with open("requirements.txt") as f:
    requirements = f.readlines()

if __name__ == "__main__":
    setup(install_requires=requirements)
