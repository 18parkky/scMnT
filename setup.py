from setuptools import setup

setup(
    name='scMnT',
    version='0.1.1',    
    description='A collection of Python scripts for analyzing microsatellite instability in single-cell resolution',
    url='https://github.com/18parkky/scMnT',
    author='Gyumin Park',
    author_email='18parkky@gm.gist.ac.kr',
    license='BSD 2-clause',
    packages=['scMnT'],
    entry_points={
        'console_scripts' : [ 'getAlleleTable = scMnT.getAlleleTable:main', 'scMnT-score = scMnT.scMnT_score:main', 'scMnT-find = scMnT.scMnT_find:main', 'scMnT = scMnT.scMnT:main',                            
] 
    },
    install_requires=['pysam>=0.20.0',
                      'Levenshtein>=0.21.1',     
                      'numpy==1.26.4',
                      'pandas>=2.0.0',
                      'scipy>=1.7.1',
                      'matplotlib>=3.7.1',
                      'seaborn>=0.13.0',
                      'scanpy>=1.11.4',           
                      ],

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',  
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 3.12.5',
    ],
)
