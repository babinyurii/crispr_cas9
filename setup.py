# this import is a quick fix of warning 'install_requires"
# https://stackoverflow.com/questions/8295644/pypi-userwarning-unknown-distribution-option-install-requires

from distutils.core import setup
import setuptools

setup(
  name = 'crispr_cas9',
  packages = ['crispr_cas9'],   
  version = '0.3',      
  license='MIT',        
  description = 'script for indels counting in sequencing data',   
  author = 'Yuriy Babin',                  
  author_email = 'babin.yurii@gmail.com',      
  url = 'https://github.com/babinyurii/crispr_cas9', 
  download_url = 'https://github.com/babinyurii/crispr_cas9/archive/v_0.3.tar.gz',
  keywords = ['bioinformatics', 'indels', 'sequencing', 'NGS'],   
  install_requires=[            
          'pandas==0.23.0',
          'biopython==1.73',
          'matplotlib==2.2.2',
          'numpy==1.11.3',
          'seaborn==0.8.1'
      ],
  classifiers=[
    'Development Status :: 4 - Beta',      
    'Intended Audience :: Developers',      
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',   
    'Programming Language :: Python :: 3',      
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
  ],
)
