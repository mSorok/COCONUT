from setuptools import setup

setup(
    name='pyclassyfire',

    version='0.1',

    description='A python client for the ClassyFire API',

    url='https://github.com/jamesjeffryes/pyclassyfire',

    author='James Jeffryes',
    author_email='jamesgjeffryes@gmail.com',

    license='MIT',

    classifiers=[
        'Development Status :: 3 - Alpha',

        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Chemistry',

        'License :: OSI Approved :: MIT License',

        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],

    keywords='Cheminformatics ClassyFire ontology',

    packages=['pyclassyfire'],

    install_requires=['requests'],

    extras_require={
        'sdf_query': ['rdkit'],
    },
)
