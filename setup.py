from setuptools import setup, find_packages


def readme():
    with open('README.md') as f:
        return f.read()


setup(
    name='BioCatalyzer',
    version='0.0.1',
    python_requires='>=3.7',
    package_dir={'': 'src'},
    packages=find_packages('src'),
    include_package_data=True,
    zip_safe=False,
    url='https://github.com/jcorreia11/BioCatalyzer',
    license='MIT',
    author='jcorreia',
    author_email='jfscorreia95@gmail.com',
    description='...',
    keywords='...',
    long_description=readme(),
    long_description_content_type='text/markdown',
    entry_points={
        'console_scripts': [
            'biocatalyzer=biocatalyzer.cli:main',
        ]
    },
)


