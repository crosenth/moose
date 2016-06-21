import classifier
import setuptools

setuptools.setup(author='Chris Rosenthal',
                 author_email='crosenth@uw.edu',
                 description='Alignment based taxonomic classifier',
                 name='classifier',
                 packages=setuptools.find_packages(exclude=['tests']),
                 entry_points={'console_scripts':
                               {'classify = classifier:main'}},
                 version=classifier.version.version(),
                 url='https://bitbucket.org/uwlabmed/classifier',
                 package_data={'classifier': ['data/*']},
                 install_requires=['pandas>=0.17.1'],
                 license='GPLv3',
                 classifiers=[
                     'Development Status :: 4 - Beta'
                     'Environment :: Console',
                     'Operating System :: OS Independent',
                     ('License :: OSI Approved :: '
                      'GNU General Public License v3 (GPLv3)'),
                     'Programming Language :: Python',
                     'Programming Language :: Python :: 2',
                     'Programming Language :: Python :: 3',
                     'Programming Language :: Python :: 2.7',
                     'Programming Language :: Python :: 3.4',
                     'Programming Language :: Python :: 3.5'])
