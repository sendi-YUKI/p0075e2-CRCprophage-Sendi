#!/usr/bin/env python3
"""Restrict native child processes to data writes beneath explicit task directories.

Landlock stacks restrictions on existing system policy without privileges or namespaces.
This is data-write confinement, not a container, network sandbox or metadata-chmod guard.
"""
import argparse
import ctypes
import os
from pathlib import Path
import platform
import stat
import sys

CREATE, ADD, RESTRICT = 444, 445, 446
WRITE_FILE=1<<1
REMOVE_DIR=1<<4
REMOVE_FILE=1<<5
MAKE_CHAR=1<<6
MAKE_DIR=1<<7
MAKE_REG=1<<8
MAKE_SOCK=1<<9
MAKE_FIFO=1<<10
MAKE_BLOCK=1<<11
MAKE_SYM=1<<12
REFER=1<<13
TRUNCATE=1<<14
WRITE_MASK=WRITE_FILE|REMOVE_DIR|REMOVE_FILE|MAKE_CHAR|MAKE_DIR|MAKE_REG|MAKE_SOCK|MAKE_FIFO|MAKE_BLOCK|MAKE_SYM|REFER|TRUNCATE
ORIGINAL_HOME=Path.home().resolve()
IPC_PATH=Path('/dev/shm')
IPC_MASK=WRITE_FILE|TRUNCATE|MAKE_REG|REMOVE_FILE
IPC_POLICY={
    'version':'posix_ipc_shm_v1', 'path':'/dev/shm',
    'allowed_access':['WRITE_FILE','TRUNCATE','MAKE_REG','REMOVE_FILE'],
    'denied_additional_access':['MAKE_DIR','REMOVE_DIR','MAKE_SYM','MAKE_CHAR','MAKE_BLOCK','MAKE_SOCK','MAKE_FIFO','REFER'],
    'purpose':'standard POSIX shared memory and named semaphores; scientific inputs and databases forbidden here',
    'dac_uid_and_sticky_bit_preserved':True,
    'same_uid_ipc_isolation':False, 'ipc_namespace_isolation':False,
    'limitation':'Landlock cannot restrict this grant by generated IPC filename prefix; existing DAC-permitted same-user shm is reachable',
}

def reject_ipc_storage(paths):
    for value in paths:
        if value is None or value == '': continue
        path=Path(value).expanduser().resolve()
        if path==IPC_PATH or path.is_relative_to(IPC_PATH):
            raise ValueError('Scientific inputs/databases/output/work cannot be stored in POSIX IPC /dev/shm: '+str(path))

def verify_ipc_directory():
    path=IPC_PATH.resolve(strict=True);info=path.stat()
    if path!=IPC_PATH or not stat.S_ISDIR(info.st_mode) or info.st_uid!=0 or stat.S_IMODE(info.st_mode)!=0o1777:
        raise RuntimeError('Standard /dev/shm must remain root-owned mode 01777; no permissions are changed')
    mounts=[line.split() for line in Path('/proc/self/mountinfo').read_text().splitlines()]
    entries=[row for row in mounts if row[4]=='/dev/shm']
    if len(entries)!=1 or entries[0][entries[0].index('-')+1]!='tmpfs' or not {'rw','nosuid','nodev'}.issubset(entries[0][5].split(',')):
        raise RuntimeError('Require existing separate rw,nosuid,nodev tmpfs at /dev/shm')
    return path

class Ruleset(ctypes.Structure):
    _fields_=[('handled_access_fs',ctypes.c_uint64)]

class Beneath(ctypes.Structure):
    _pack_=1
    _fields_=[('allowed_access',ctypes.c_uint64),('parent_fd',ctypes.c_int32)]

def apply_write_guard(writable):
    if platform.machine() not in {'aarch64','x86_64'}:
        raise RuntimeError('Unsupported Landlock syscall architecture')
    libc=ctypes.CDLL(None,use_errno=True)
    libc.syscall.restype=ctypes.c_long
    abi=libc.syscall(CREATE,ctypes.c_void_p(),0,1)
    if abi<3:
        raise RuntimeError(f'Landlock ABI >=3 required, observed {abi}, errno={ctypes.get_errno()}')
    paths=[]
    for item in writable:
        p=Path(item).resolve(strict=True)
        reject_ipc_storage([p])
        if not p.is_dir() or p in {Path('/'),Path('/home'),ORIGINAL_HOME}:
            raise ValueError(f'Refusing overbroad or non-directory write grant: {p}')
        paths.append(p)
    if not paths: raise ValueError('Explicit task write directory required')
    attr=Ruleset(WRITE_MASK)
    fd=libc.syscall(CREATE,ctypes.byref(attr),ctypes.sizeof(attr),0)
    if fd<0: raise OSError(ctypes.get_errno(),'landlock_create_ruleset')
    try:
        grants=[(p,WRITE_MASK) for p in paths]
        # This project rule retains ordinary POSIX IPC; it does not grant system
        # permissions or relax protection of inputs, databases or environments.
        grants.append((verify_ipc_directory(),IPC_MASK))
        grants.extend((Path(p),WRITE_FILE) for p in ['/dev/null','/dev/tty'] if Path(p).exists())
        for path,access in grants:
            parent=os.open(path,os.O_PATH|os.O_CLOEXEC)
            try:
                rule=Beneath(access,parent)
                if libc.syscall(ADD,fd,1,ctypes.byref(rule),0)!=0:
                    raise OSError(ctypes.get_errno(),f'landlock_add_rule: {path}')
            finally: os.close(parent)
        if libc.prctl(38,1,0,0,0)!=0:
            raise OSError(ctypes.get_errno(),'PR_SET_NO_NEW_PRIVS')
        if libc.syscall(RESTRICT,fd,0)!=0:
            raise OSError(ctypes.get_errno(),'landlock_restrict_self')
    finally: os.close(fd)
    return abi

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--writable',action='append',required=True)
    parser.add_argument('command',nargs=argparse.REMAINDER)
    args=parser.parse_args()
    command=args.command[1:] if args.command[:1]==['--'] else args.command
    if not command: parser.error('A child command is required')
    apply_write_guard(args.writable)
    os.execvpe(command[0],command,os.environ)

if __name__=='__main__':
    try: main()
    except (OSError,ValueError,RuntimeError) as error:
        print(f'spark_fs_guard: {error}',file=sys.stderr)
        raise SystemExit(78)
