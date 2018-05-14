from setuptools import setup,find_packages
setup(name='mk_mass',
      version='1.0',
      py_modules=['mk_mass'],
      packages=find_packages() + ['resources'],
      include_package_data=True
      )