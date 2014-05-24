'''
Run `NCBI BLAST`_ with an easy-to-use Pythonic API.

See https://github.com/pgriess/pyblast for further details.
'''

import errno
import fcntl
import os
import select
import StringIO
import subprocess


VERSION = '0.1'
'''
Version of pyblast.
'''


DEFAULT_HIT_FIELDS = [
    'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
    'sstart', 'send', 'evalue', 'bitscore']
'''
The default set of fields available in Hit objects.
'''


class Result(object):
    '''
    The result of a BLAST processing a single query sequence.

    Each Result object contains some number of Hit objects
    representing hits in the searched database.
    '''

    id = None
    '''
    Identifier for the query corresponding to this result. Can be
    None if no identifier was provided.
    '''

    description = None
    '''
    Description of the query corresponding to this result. Can be
    None if no description was provided.
    '''

    hits = None
    '''
    Iterable of of Hit objects.
    '''

    def __init__(self):
        self.hits = []


class Hit(object):
    '''
    A single hit in the BLAST database.

    Fields are referenced as attributes on the instance (e.g.
    h.qseqid to get the 'qseqid' field of the 'h' Hit object). If the
    given field is not present in the hit, None is returned.
    '''

    __TYPES = [int, float, str]

    def __init__(self, fields):
        self.__fields = {}

        for k, v in fields.iteritems():
            for t in Hit.__TYPES:
                try:
                    self.__fields[k] = t(v)
                    break
                except:
                    pass

    def __getattr__(self, fn):
        return self.__fields.get(fn, None)

    def __dir__(self):
        return dir(type(self)) + self.__fields.keys()


def __read_single_fasta_query_lines(f):
    '''
    Read and return sequence of lines (including newlines) that
    represent a single FASTA query record. The provided file is
    expected to be blocking.

    Returns None if there are no more query sequences in the file.
    '''

    def readline():
        l = f.readline()
        if l == '':
            raise EOFError()
        return l

    rec = None
    try:
        l = readline()
        assert l.startswith('>')

        rec = [l]
        while True:
            pos = f.tell()
            l = readline()
            if l.startswith('>'):
                f.seek(pos, 0)
                break
            rec += [l]
    except EOFError:
        pass

    return rec


def __read_single_query_result(rs, field_names):
    '''
    Read the result of a single query from the given string,
    returning a tuple of (record, remaining-string). If no complete
    record could be read, the first element of the tuple is None and
    the second element is the original imput string.
    '''

    rf = StringIO.StringIO(rs)

    def readline():
        l = rf.readline()
        if not l.endswith('\n'):
            raise EOFError()
        return l.strip()

    result = Result()

    try:
        l = readline()
        assert l.startswith('# BLAST')

        l = readline()
        assert l.startswith('# Query: ')
        query_str = l[len('# Query: '):].strip()
        if query_str:
            if ' ' in query_str:
                result.id, result.description = [
                    s.strip() for s in query_str.split(' ', 1)]
            else:
                result.id = query_str

        l = readline()
        assert l.startswith('# Database: ')

        l = readline()
        if l.startswith('# Fields: '):
            fns = l[len('# Fields: '):].split(', ')
            assert len(field_names) == len(fns)
            l = readline()

        assert l.endswith(' hits found')
        nhits = int(l[len('# '):-1 * len(' hits found')])

        while nhits > 0:
            l = readline()
            field_vals = l.split('\t')
            assert len(field_vals) == len(field_names)

            fields = dict(zip(field_names, field_vals))
            result.hits.append(Hit(fields))
            nhits -= 1

        return result, rf.read()
    except EOFError:
        return None, rs


def __run_blast_select_loop(input_file, popens, fields):
    '''
    Run the select(2) loop to handle blast I/O to the given set of Popen
    objects.

    Yields records back that have been read from blast processes.
    '''

    def make_nonblocking(f):
        fl = fcntl.fcntl(f.fileno(), fcntl.F_GETFL)
        fl |= os.O_NONBLOCK
        fcntl.fcntl(f.fileno(), fcntl.F_SETFL, fl)

    rfds = set()
    wfds = set()
    fd_map = {}

    for p in popens:
        make_nonblocking(p.stdout)
        rfds.add(p.stdout.fileno())
        fd_map[p.stdout.fileno()] = {
            'popen': p,
            'query_buffer': '',
            'result_buffer': ''}

        make_nonblocking(p.stdin)
        wfds.add(p.stdin.fileno())
        fd_map[p.stdin.fileno()] = fd_map[p.stdout.fileno()]

    while len(rfds) + len(wfds) > 0:
        # XXX: Should we be tracking excepted file descriptors as well?
        rl, wl, _ = select.select(rfds, wfds, [])

        # For each of our readable blast processes, read response
        # records and emit them
        for fd in rl:
            rs = fd_map[fd]['result_buffer']
            rbuf = os.read(fd, select.PIPE_BUF)

            # The blast process has finished emitting records. Stop
            # attempting to read from or write to it. If we have
            # excess data in our result_buffer, c'est la vie.
            if rbuf == '':
                p = fd_map[fd]['popen']

                rfds.remove(p.stdout.fileno())
                p.stdout.close()

                if not p.stdin.closed:
                    wfds.remove(p.stdin.fileno())
                    p.stdin.close()
                
                continue

            rs += rbuf
            while True:
                rec, rs = __read_single_query_result(rs, fields)
                if rec is None:
                    break
                yield rec

            fd_map[fd]['result_buffer'] = rs

        # For each of our writable blast processes, grab a new query
        # sequence and send it off to them.
        for fd in wl:
            qs = fd_map[fd]['query_buffer']
            if not qs:
                ql = __read_single_fasta_query_lines(input_file)

                # No more input records available. Close the pipe to
                # signal this to the blast process.
                if ql is None:
                    p = fd_map[fd]['popen']

                    wfds.remove(p.stdin.fileno())
                    p.stdin.close()

                    continue

                qs = ''.join(ql)

            # XXX: For some reason, despite select(2) indicating that
            #      this file descriptor is writable, writes can fail
            #      with EWOULDBLOCK. Handle this gracefully.
            try:
                written = os.write(fd, qs)
                qs = qs[written:]
            except OSError, e:
                assert e.errno == errno.EWOULDBLOCK

            fd_map[fd]['query_buffer'] = qs


def __run_blast(blast_command, input_file, *args, **kwargs):
    '''
    Run a blast variant on the given input file.
    '''

    # XXX: Eventually, translate results on the fly as requested? Or
    #      just always use our parsed object?
    if 'outfmt' in kwargs:
        raise Exception('Use of the -outfmt option is not supported')

    num_processes = kwargs.get(
        'pb_num_processes', os.sysconf('SC_NPROCESSORS_ONLN'))
    fields = kwargs.get('pb_fields', DEFAULT_HIT_FIELDS)

    blast_args = [blast_command]
    blast_args += ['-outfmt', '7 {}'.format(' '.join(fields))]
    for a in args:
        blast_args += ['-' + a]
    for k, v in kwargs.iteritems():
        if not k.startswith('pb_'):
            blast_args += ['-' + k, str(v)]

    popens = []
    for _ in range(num_processes):
        popens.append(
            subprocess.Popen(
                args=blast_args, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                stderr=None, close_fds=True))

    try:
        for r in __run_blast_select_loop(input_file, popens, fields):
            yield r
    finally:
        for p in popens:
            if p.poll() is None:
                p.terminate()
            p.wait()


def blastn(input_file, *args, **kwargs):
    '''
    Run blastn on the given input file.

    All positional arguments are passed directly to blastn as commandline
    options without values. For example, passing 'ungapped' will invoke blastn
    with the '-ungapped' option. Likewise, all keyword arguments are passed to
    blastn as options with values. For example, db='foo/bar' will result in
    blastn being invoked with '-db foo/bar'.

    Some extra keyword arguments are supported which are not passed to blast:

      pb_num_processes      number of blast processes to spawn; defaults to
                            sysconf(SC_NPROCESSORS_ONLN)
      pb_fields             iterable of field names to retrieve for each hit
    '''

    return __run_blast('blastn', input_file, *args, **kwargs)

def blastp(input_file, *args, **kwargs):
    return __run_blast('blastp', input_file, *args, **kwargs)

def blastx(input_file, *args, **kwargs):
    return __run_blast('blastx', input_file, *args, **kwargs)

# Crib documentation strings from blastn
blastp.__doc__ = blastn.__doc__.replace('blastn', 'blastp')
blastx.__doc__ = blastn.__doc__.replace('blastn', 'blastx')
