from setuptools import find_packages, setup

export_packages = find_packages(exclude=["tests"])

setup(
    name="backpressure_report",
    version="0.1",
    author="Batch Tasks",
    packages=export_packages,
    python_requires=">=3.8",
    install_requires=[
        "python-dateutil",
    ],
)
