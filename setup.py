import setuptools

setuptools.setup(author='Chris Rosenthal',
                 author_email='crosenth@gmail.com',
                 classifiers=[
                     'License :: OSI Approved :: '
                     'GNU General Public License v3 (GPLv3)',
                     'Environment :: Console',
                     'Operating System :: OS Independent',
                     'Intended Audience :: End Users/Desktop',
                     'License :: OSI Approved :: '
                     'GNU General Public License v3 (GPLv3)',
                     'Programming Language :: Python :: 3',
                     ],
                 description='Alignment based taxonomic classifier',
                 entry_points={
                     'console_scripts': {'classify=classifier.classify:main'}},
                 install_requires=['pandas>=2.0.2'],
                 keywords=[
                     'ncbi', 'blast', 'classifier', 'genetics', 'genomics',
                     'dna', 'rna', 'bioinformatics'],
                 license='GPLv3',
                 name='moose_classifier',
                 packages=setuptools.find_packages(exclude=['tests']),
                 python_requires='>=3.7',
                 url='https://github.com/crosenth/moose',
                 version=0.9
                 )
