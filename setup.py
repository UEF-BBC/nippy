from setuptools import setup


setup(name='nippy',
      version='1.0',
      description='Semi-automatic NIR preprocessing module',
      author='Jari Torniainen, Isaac Afara, Biophysics of Bone and Cartilage, Applied Physics, University of Eastern Finland',
      author_email='jari.torniainen@uef.fi, isaac.afara@uef.fi',
      packages=['nippy'],
      package_dir={'nippy': 'nippy'},
      include_package_data=False,
      url='https://github.com/uef-bbc/nippy',
      install_requires=['numpy>=1.13.3', 'scipy>=0.19.1', 'sklearn>=0.0']
      )
