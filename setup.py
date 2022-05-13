from setuptools import setup, find_packages

setup(
    name='RulesTools',
    version='0.1.0',
    packages=find_packages(include=['RulesTools', 'RulesTools.database_processing', 'RulesTools.filters',
                                    'RulesTools.transformations']),
    install_requires=['CGRtools', 'tqdm', 'StructureFingerprint', 'numpy', 'python>=3.8'],
    license='MIT'
)
