from distutils.core import setup

setup(
  name = 'crispr_cas9',
  packages = ['crispr_cas9'],   
  version = '0.1',      
  license='MIT',        
  description = 'script for indels counting in deep sequencing data',   
  author = 'Yuriy Babin',                  
  author_email = 'babin.yurii@gmail.com',      
  url = 'https://github.com/babinyurii/crispr_cas9', 
  download_url = 'https://github.com/babinyurii/crispr_cas9/archive/v_0.1.tar.gz',
  keywords = ['bioinformatics', 'SNP'],   
  install_requires=[            
          'pandas',
          'biopython',
          'matplotlib',
          'numpy',
          'seaborn'
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
