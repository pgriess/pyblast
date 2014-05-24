=======
PyBLAST
=======

Run `NCBI BLAST`_ with an easy-to-use Pythonic API.

Running NCBI BLAST manually is of course not rocket science, but this
module provides several benefits over doing so:

* Automatically runs a BLAST process for each CPU on the system;
  achieves far better throughput than the ``-num_threads`` option
* Provides an iterator API that emits native Python objects for each
  BLAST result `as they're produced`, rather than at the end
* ``Result`` and ``Hit`` objects obviate the need for manually parsing
  results; all values represented by their native Python types (e.g.
  ``Hit.evalue`` is a ``float``, etc)

.. _`NCBI BLAST`: http://blast.ncbi.nlm.nih.gov/

Example
-------

Here's a simple example with comments hilighting some relevant
features

::

    import pyblast
    
    with open('data.fasta') as f:
        # Use the pyblast.blastx() iterator function
        for r in pyblast.blastx(f, db='/path/to/swissprot'):
            msg = 'query {} has {} hits'.format(r.query_id, len(r.hits))
            if r.hits:
                # Use Hit.evalue as a float for comparison
                min_evalue = sorted([h.evalue for h in r.hits])[0]
                msg += '; minimum evalue {:f}'.format(min_evalue)

            print msg

This will produce output like the following

::

    query M00181:167:000000000-A4VBV:1:1101:11880:1874 1:N:0:6 has 6 hits; minimum evalue 0.310000
    query M00181:167:000000000-A4VBV:1:1101:17067:1875 1:N:0:6 has 14 hits; minimum evalue 0.200000
    query M00181:167:000000000-A4VBV:1:1101:15039:1878 1:N:0:6 has 4 hits; minimum evalue 4.400000
    query M00181:167:000000000-A4VBV:1:1101:17090:1895 1:N:0:6 has 6 hits; minimum evalue 1.700000
    query M00181:167:000000000-A4VBV:1:1101:15843:1907 1:N:0:6 has 2 hits; minimum evalue 1.800000


API
---

.. _`blastn`:

``blastn(input_file, *args, **kwargs)``

Iterator to process the contents of the FASTA ``input_file`` using
``blastn``; yields `Result`_ objects.

The ``*args`` and ``**kwargs`` arguments control how `blastn` is
invoked. The former are passed as options without values, while the
latter are passed as options with values. For example,
``blastn(some_file, 'ungapped', db='foo/bar')`` will run ``blastn``
with the ``-ungapped -db foo/bar`` options.

In addition, the following keyword arguments are handled specially and
are not passed on to BLAST:

- ``pb_num_processes``: number of BLAST processes to spawn; default is ``sysconf(SC_NPROCESSORS_ONLN)``
- ``pb_fields``: iterable of field names to retrieve for each hit; default is `DEFAULT_HIT_FIELDS`_. The list of valid field names (and their meanings) can be found in the ``*** Formatting options`` section of ``blastn -help``.

``blastp(input_file, *args, **kwargs)``

See documentation for `blastn`_.

``blastx(input_file, *args, **kwargs)``

See documentation for `blastn`_.

.. _`Result`:

``Result``

The result of BLAST processing a single query sequence. The set of
attributes on this object are:

- ``id``: identifier for the query sequence; can be ``None``
- ``description``: textual description of the query sequence; can be ``None``
- ``hits``: array of `Hit`_ objects

.. _`Hit`:

``Hit``

A single sequence hit in a `Result`_ object.

The attributes of this object are the names of the fields requested of
BLAST. For example, if ``blastn`` was run with ``pb_fields=['qseqid',
...]`` then one could access the ``qseqid`` value of the ``Hit``
object ``h`` like so: ``h.qseqid``. Fields referenced that were not
requested of BLAST have a ``None`` value.

In addition, BLAST fields are converted to their native Python types.
For example, ``evalue`` fields are automatically converted to floating
point values.

.. _`DEFAULT_HIT_FIELDS`:

``DEFAULT_HIT_FIELDS``

The default fields returned for each ``Hit`` object.
