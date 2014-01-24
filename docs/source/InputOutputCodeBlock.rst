Input Output Code Block
=======================

This document is a proof of concept example of using Sphinx/reST syntax to create a code block with a different colored background.

Input Code Block
^^^^^^^^^^^^^^^^^

    An input code block:
   
    .. code-block:: r
         :linenos:
        
         print("just a test")
         print(8/2)

Output Code Block
^^^^^^^^^^^^^^^^^

.. cssclass:: output
..

   An output code block:
   
   .. code-block:: r
        :linenos:
       
        just a test
        4
