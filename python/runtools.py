#!/usr/bin/env python
import os
import sys
import time
import numpy as N

atmos_dir = '/mydata'
seq_file = 'RH_SEQUENCE'
outsuff = 'output/output_ray_CaII_PRD_s%03i.ncdf'
mpicmd = 'mpiexec'  # 'aprun'
bin = './rh15d_ray_pool'
defkey = 'keyword.save'
log = 'rh_running.log'
tm = 40   # timeout in minutes to stop a single run
rerun = True
rerun_opt = [{'NG_DELAY': 40, 'NG_PERIOD': 20, '15D_DEPTH_REFINE': 'FALSE',
              '15D_ZCUT': 'TRUE', 'N_MAX_ITER': 250, 'PRD_SWITCH': 0.002},
             {'NG_DELAY': 120, 'NG_PERIOD': 100, '15D_DEPTH_REFINE': 'FALSE',
              '15D_ZCUT': 'FALSE', 'N_MAX_ITER': 350, 'PRD_SWITCH': 0.001,
              'ITER_LIMIT': '1.0E-2'}]
rerun_max = len(rerun_opt)


def run_sequence(seq_file, atmos_dir='.'):
    """
    Runs run_single for a sequence written to file seq_file, which
    contains one atmos file per line, optionally followed by a subset of
    its snapshot numbers.
    """
    retries = 0
    while True:
        # Try I/O operations three times to ensure Pleiades hiccups don't
        # cause problems with the run
        try:
            job = getjob(seq_file, fdir=atmos_dir)
            retries = 0
        except IOError, ee:
            print(ee)
            if retries == 3:
                break
            sleep_time = 60 * 2 ** (retries + 1)
            print('Sleeping for %i minutes and retrying...' % (sleep_time/60))
            time.sleep(sleep_time)
            retries += 1
            continue
        if job is None:
            break
        atmos_file, snaps, snap_nums = job
        nt = len(snaps)
        print('--- ' + os.path.split(atmos_file)[1])
        for i, n, s in zip(range(nt), snaps, snap_nums):
            print('     running snapshot %s (%i/%i)' % (s, i + 1, nt))
            run_single(atmos_file, n, snapnum=s)
    print('*** Finished gracefully.')
    return


def run_single(atmos_file, snap, snapnum=None):
    if snapnum is None:
        snapnum = snap
    outfile = outsuff % (snapnum)
    # Option for hexagon
    if mpicmd == 'aprun':     # Do NOT put time before commands! (Wrong kill)
        runcmd = '%s -B %s' % (mpicmd, bin)
    else:
        runcmd = '%s %s' % (mpicmd, bin)
    nv = {'SNAPSHOT': snap, 'ATMOS_FILE': atmos_file, '15D_RERUN': 'FALSE'}
    keywordin_update(defkey, 'keyword.input', nv)

    if os.path.isfile(outfile):
        print('File %s exists, skipping' % outfile)
        return
    # Retry a few times to avoid Pleiades timeouts
    if run_timeout(runcmd, timeout=tm, log=log):
        print('ERROR: unclean exit status, sleeping 5 minutes and retrying...')
        time.sleep(60 * 5)
        nv['15D_RERUN'] = 'TRUE'
        keywordin_update(defkey, 'keyword.input', nv)
        if run_timeout(runcmd, timeout=tm, log=log):
            print('ERROR: unclean exit status, sleeping 15 minutes and'
                  ' retrying...')
            time.sleep(60 * 15)
            if run_timeout(runcmd, timeout=tm, log=log):
                print('ERROR: unclean exit status, skipping.')
                return
    # re-run to obtain convergence
    if rerun:
        import netCDF4
        nv['15D_RERUN'] = 'TRUE'
        for i in range(rerun_max):
            # find out how many non-converged items
            try:
                f = netCDF4.Dataset('output/output_indata.ncdf', 'r')
                conv = f.groups['mpi'].variables['convergence'][:]
                f.close()
            except:
                print('WARNING: output_indata.ncdf not found or lacking'
                      'convergence info. Aborting rerun.')
                break
            nremain = N.sum(conv < 1)
            if N.ma.isMaskedArray(conv):
                masked_remain = N.sum(conv.mask)
                if masked_remain == conv.size:
                    print('WARNING: no columns were fully run in previous try')
                    nremain = 0
                nremain += masked_remain
            if nremain == 0:
                break
            print('Remaining columns: %i' % nremain)
            # change some stuff on keyword.input (including rerun = TRUE)
            nn = nv.copy()
            nn.update(rerun_opt[i])
            keywordin_update(defkey, 'keyword.input', nn)
            # Run with retries
            if run_timeout(runcmd, timeout=tm, log=log):
                print('ERROR: unclean exit status, sleeping 5 minutes and'
                      ' retrying...')
                time.sleep(60 * 5)
                if run_timeout(runcmd, timeout=tm, log=log):
                    print('ERROR: unclean exit status, sleeping 15 minutes and'
                          ' retrying...')
                    time.sleep(60 * 15)
                    if run_timeout(runcmd, timeout=tm, log=log):
                        print('ERROR: unclean exit status, skipping'
                              ' but copying file.')
    print('mv output/output_ray.ncdf ' + outfile)
    if os.system('mv output/output_ray.ncdf ' + outfile):
        print('ERROR: unclean exit status, sleeping 5 minutes and'
              ' retrying...')
        time.sleep(60 * 5)
        if os.system('mv output/output_ray.ncdf ' + outfile):
            print('ERROR: unclean exit status, sleeping 15 minutes and'
                 ' retrying...')
            time.sleep(60 * 15)
            if os.system('mv output/output_ray.ncdf ' + outfile):
                print('ERROR: file copy failed.')
    return


def run_timeout(cmd, timeout=0, log='rh_running.log'):
    '''
    Runs cmd and stops execution if nothing is written to stdout
    within timeout minutes (if zero, there is no timeout).
    '''
    import sys
    import subprocess
    import signal
    import shlex
    import time
    print(cmd)
    logfile = open(log, 'w', 1)
    args = shlex.split(cmd)   # tokenise args list
    p = subprocess.Popen(args, shell=False, bufsize=0, stdin=None,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    def _handler(signum, frame):
        print('Timeout of %s min reached, stopping execution' % timeout)
        p.terminate()
        kill_mpirun()
        time.sleep(30)  # to ensure no ghost process is left running
        raise RuntimeError('Timeout')

    signal.signal(signal.SIGALRM, _handler)
    try:
        while True:
            signal.alarm(int(timeout * 60))
            inline = p.stdout.readline()
            if not inline:
                break
            logfile.write(inline)
            signal.alarm(0)
    except RuntimeError:
        logfile.close()
        return 0

    p.communicate()   # wait for process to finish, get return code
    logfile.close()
    return p.returncode


def kill_mpirun():
    """
    Tries to kill the actual mpirun process (not the wrapper mpiexec) in
    systems like pleiades
    """
    import os
    import subprocess
    import signal
    p = subprocess.Popen(['ps', 'x'], stdout=subprocess.PIPE)
    out, err = p.communicate()
    for line in out.splitlines():
        if ('mpiexec.params' in line) or ('mpirun' in line):
            pid = int(line.split(None, 1)[0])
            os.kill(pid, signal.SIGKILL)
    return


def getjob(filename, fdir='.'):
    """
    Gets RH job from filename. Job contains input atmos file and which
    snapshots to run. After read, the first line is deleted from the
    file. When there are no more entries in the file, it returns None.

    Parameters:
    -----------
    filename - string with name of input file
    fdir - string with directory name

    Returns:
    -------
    result - tuple
        First element contains the input file name, the second a list
        with the snapshot indexes, and the third a list with the snapshot
        numbers.
    """
    import os
    import netCDF4
    # read first line, write rest to the file
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()
    try:
        curr = lines[0]
        j = 1
        while len(curr.split()) == 0:   # skip blank lines
            curr = lines[j]
            j += 1
    except IndexError:
        return None
    f = open(filename, 'w')
    for l in lines[j:]:
        f.write(l)
    f.close()
    result = curr.split()
    result = (os.path.join(fdir, result[0]), result[1:])
    for retries in range(3):
        if retries == 2:
            raise IOError('File %s not found, skipping.' % result[0])
        if os.path.isfile(result[0]):
            break
        print('File %s not found, sleeping %i minutes and '
              'retrying...' % (result[0], retries + 2))
        time.sleep((retries + 1) * 60)
    snaps = get_snaps(result[0])
    if len(result[1]) == 0:
        selected = range(len(snaps))
        selected_snaps = snaps
    else:
        selected = []
        selected_snaps = []
        for s in result[1]:
            try:
                selected.append(snaps.index(int(s)))
                selected_snaps.append(int(s))
            except ValueError:
                print('WARNING: snapshot %s not in %s.' % (s, result[0]))
    result = (result[0], selected, selected_snaps)
    return result


def get_snaps(ncdf_file):
    ''' Gets snapshot numbers from netCDF file'''
    import netCDF4
    f = netCDF4.Dataset(ncdf_file, 'r')
    snaps = f.variables['snapshot_number'][:]
    f.close()
    return list(snaps)


def keywordin_update(infile, outfile, new_values):
    ''' Updates a given number of fields with values on a keywords.input file.
        These are given in a dictionary: fvalues = {field: value}.
        Reads from infile and writes into outfile.'''
    out = open(outfile, 'w')
    for line in file(infile):
        if line[0] == '#':
            out.write(line)
        elif line.find('=') < 0:
            out.write(line)
        else:
            ss = line.split('=')[0]
            ssv = ss.strip().upper()
            if ssv in new_values.keys():
                out.write('%s= %s\n' % (ss, str(new_values[ssv])))
            else:
                out.write(line)
    return

if __name__ == '__main__':
    run_sequence(seq_file, atmos_dir=atmos_dir)
