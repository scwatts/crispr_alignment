#!/usr/bin/env python3
import setuptools


# Set package details
package_name = 'crispr_align'
package_description = 'Create strong ordering of CRISPR spacers using graphical topology sorting'
package_version = '0.0.1'
authors = 'Stephen Watts, Alex Tokolyi, Kat Holt'
license = 'GPLv3'


# Call setup
setuptools.setup(
        name=package_name,
        version=package_version,
        license=license,
        author=authors,
        test_suite='tests',
        packages=setuptools.find_packages(),
        scripts=['crispr_spacer_alignment.py']
        )
