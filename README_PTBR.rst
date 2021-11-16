====
sara
====

**s**\emi\  **a**\utomatic\  **R**\NA-Seq\  **a**\nalysis

Este repositório contém um fluxo de trabalho para realizar análise de RNA-Seq criado por Laise de Moraes e Joyce Silva (FIOCRUZ-IGM).

********************************
Configurando o fluxo de trabalho
********************************

Baixe e instale o fluxo de trabalho a partir do repositório do GitHub:

.. code:: bash

    git clone --recursive https://github.com/khourious/sara.git; cd sara
    chmod 700 -R DEPENDENCIES
    bash DEPENDENCIES

********************
Como utilizar o sara
********************

Será necessário criar uma planilha em formato CSV. Esta planilha precisa estar armazenada no diretório ``SHEETS``.

O nome do arquivo CSV irá **corresponder ao nome da bibliotera de RNA-Seq** e conterá na primeira linha: xxx,xxx,xxx e nas próximas linhas: xxxxxxxxxxxxxxxxxxxxxxxxxxxxx -- ATENÇÃO: **SEM CABEÇALHO!!**

.. code-block:: text

    xxx,xxx,xxx

.. code-block:: text

    semi automatic RNA-Seq analysis <<< análise semi automática de RNA-Seq >>>

    Uso: sara -c <nome do arquivo de configuração>

    -c  Nome do arquivo CSV que contém a lista de programas e os grupos experimentais.
    -t  Número máximo de threads (padrão: todos os núcleos).
