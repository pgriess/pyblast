=======
PyBLAST
=======

Run `NCBI BLAST`_ with an easy-to-use Pythonic API.

Running NCBI BLAST manually is of course not rocket science, but this
module provides several benefits over doing so:

* Automatically runs a BLAST process for each CPU on the system. This
  achieves far better throughput than the ``-num_threads`` option.
* Provides an iterator API that emits native Python objects for each
  BLAST result, obviating the need for manually parsing results

.. _`NCBI BLAST`: http://blast.ncbi.nlm.nih.gov/

API
---

``blastn(input_file, *args, **kwargs)``

Iterator to process the contents of the FASTA ``input_file`` using
``blastn``.

The ``*args`` and ``**kwargs`` arguments control how `blastn` is
invoked. The former are passed as options without values, while the
latter are passed as options with values. For example,
``blastn(some_file, 'ungapped', db='foo/bar')`` will run ``blastn``
with the ``-ungapped -db foo/bar`` options.

In addition, the following keyword arguments are handled specially and
are not passed on to BLAST:

* ``pb_num_processes`` -- number of BLAST processes to spawn; default is ``sysconf(SC_NPROCESSORS_ONLN)``

``blastp(input_file, *args, **kwargs)``

See documentation for ``blastn``.

``blastx(input_file, *args, **kwargs)``

See documentation for ``blastn``.
