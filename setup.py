from setuptools import setup, find_packages

setup(
    name='RulesTools',
    version='0.1.0',
    packages=find_packages(include=['RulesTools', 'RulesTools.database_processing', 'RulesTools.filters',
                                    'RulesTools.transformations']),
    python_requires='>=3.9',
    install_requires=['CGRtools', 'tqdm', 'StructureFingerprint', 'numpy', 'py-mini-racer>=0.4.0'],
    license='MIT'
)
