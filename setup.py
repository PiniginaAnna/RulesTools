from setuptools import setup, find_packages

setup(
    name='RulesTools',
    version='0.1.0',
    packages=find_packages(include=['RulesTools']),
    install_requires=['CGRtools', 'tqdm'],
    license='MIT',
)



# from setuptools import find_packages, setup
#
# setup(
#     # name='mypythonlib',
#     packages=find_packages(include=['mypythonlib']),
#     # version='0.1.0',
#     # description='My first Python library',
#     # author='Me',
#     # license='MIT',
# )