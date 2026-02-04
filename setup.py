from setuptools import setup, find_packages

setup(
    name="triangles",
    version="1.0.0",
    description="Mesh generator for SWEpy",
    author="Juan A. Fuenzalida A.",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy",
        # cupy is optional, so not strictly required here, but good to mention in README
    ],
    py_modules=["main", "cases", "filesave", "grids", "overlap", "triangles", "backend"],
    entry_points={
        'console_scripts': [
            'triangles-mesh=run:main', 
        ],
    },
)
