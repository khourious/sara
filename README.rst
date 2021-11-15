====
sara
====

**s**\emi\  **a**\utomatic\  **R**\NA-Seq\  **a**\nalysis

This repository contains a workflow to performs a RNA-Seq analysis created by Laise de Moraes and Joyce Silva (FIOCRUZ-IGM).

***********************
Setting up the pipeline
***********************

Download and install the pipeline from the GitHub repo:

.. code:: bash

    git clone --recursive https://github.com/khourious/sara.git; cd sara
    chmod 700 -R DEPENDENCIES
    bash DEPENDENCIES

****************************
How to use the sara pipeline
****************************

It is necessary to create the sample sheet (.csv). You need create in ``sheets`` directory.

The csv file name **corresponds to the library name** and contains in the first line: xxx,xxx,xxx and the next lines: xxxxxxxxxxxxxxxxxxxxxxxxxxxxx -- ATTENTION: **NO HEADER!!**

.. code-block:: text

    xxx,xxx,xxx

.. code-block:: text

    semi automatic RNA-Seq analysis

    Usage: sara -c <config file name>

    -c  Name of CSV file that contains a list of the softwares and the experimental groups (default path: /sara/sheets).
    -t  Max number of threads (default: all cores).
