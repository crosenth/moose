#!/usr/bin/env python
"""
execute the classifier program without setuptools installation
"""

import classifier
import sys

if __name__ == '__main__':
    sys.exit(classifier.main(sys.argv[1:]))
