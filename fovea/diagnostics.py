"""
Structured logging, noSQL database, and diagnostic administration
"""
from __future__ import division, absolute_import

import structlog
from structlog import BoundLoggerBase, PrintLogger, wrap_logger
from structlog.processors import JSONRenderer, KeyValueRenderer
from pprint import pprint
from tinydb import TinyDB, where
from tinydb.storages import MemoryStorage, JSONStorage
from tinydb.middlewares import CachingMiddleware

from collections import OrderedDict
from PyDSTool import info  # use for viewing log dicts
import json

import datetime
import time
import hashlib
import sys, os

__all__ = ['diagnostic_manager', 'get_unique_name', 'info', 'load_log',
           'where', 'filter_log']


global name_registry
name_registry = {}

def get_unique_name(name, start=None):
    """
    Use optional start integer to begin numbering of the name with that value
    """
    try:
        count = name_registry[name]
    except KeyError:
        if start is not None:
            start = int(start)
            start_str = str(start)
            name_registry[name] = start
            out_name = name + '_' + start_str
        else:
            name_registry[name] = 0
            out_name = name
    else:
        out_name = name + '_%i' % (count + 1)
        name_registry[name] = count + 1
    return out_name


class diagnostic_manager(object):
    def __init__(self, name, make_log=False, dbfilepath=None, dirpath='.'):
        """
        By default, does not create a log for diagnostics unless make_log=True.
        Log can be set with make_log() method after object created.
        """
        self.name = name
        self.global_count = 0
        self._dirpath = dirpath
        if make_log:
            self.log = wrap_logger(FoveaPrintLogger(dbfilepath=\
                                             os.path.join(dirpath,dbfilepath)),
                                   wrapper_class=SemanticLogger,
                                   processors=[ev_store, KeyValueRenderer()],
                       ) #JSONRenderer(indent=1, sort_keys=True)])
            # reference to the database
            self.db = self.log._logger.db
        else:
            self.log = None
            self.db = None
        # store any metadata associated with log table entries
        self.log_items_digest = {}
        self.name_to_digest = {}

    def use_dir(self, dirpath):
        """Uses directory from given path (relative or absolute) for
        log files, etc. Directory is created if it doesn't already exist.
        """
        import os
        if os.path.isfile(dirpath):
            # exists as a file already
            raise ValueError("Path exists but is not a directory")
        elif not os.path.isdir(dirpath):
            os.mkdir(dirpath)
        self._dirpath = dirpath

    def make_log(self, dbfilepath=None):
        """
        Resets and recreates log if one not already created, using
        optional dbfilepath for file output.
        """
        self.log = wrap_logger(FoveaPrintLogger(dbfilepath=\
                                       os.path.join(self._dirpath,dbfilepath)),
                               wrapper_class=SemanticLogger,
                               processors=[ev_store, KeyValueRenderer()],
                               ) #JSONRenderer(indent=1, sort_keys=True)])
        # reference to the database
        self.db = self.log._logger.db

    def get_events(self):
        """
        Convenience method to fetch "invisible" event_list obscured
        inside the _logger attribute (due to mixin-like inheritance magic
        in structlog)
        """
        try:
            return self.log._logger.event_list
        except AttributeError:
            raise AttributeError("Log not yet created")

    def get_unique_name(self, name):
        """
        Convenience method for access to function
        """
        return get_unique_name(name)

    def attach_obj(self, obj, name):
        # sha the name and obj repr
        digest = unique_sha(repr(obj)+name)
        # add obj to log_items
        try:
            self.log_items_digest[digest] = obj
        except AttributeError:
            raise AttributeError("Log not yet created")
        if name in self.name_to_digest:
            raise ValueError()
        else:
            self.name_to_digest[name] = digest
        self.log._attach_obj(obj_name=name, obj_digest=digest)


class counter_util(object):
    def __init__(self):
        # will start at 1 when used
        self.n = 0

    def get_count(self):
        self.n += 1
        return self.n


def ev_store(logger, log_method, event_dict):
    logger.event = event_dict
    return event_dict


class SemanticLogger(BoundLoggerBase):
    def __init__(self, logger, processors, context):
        self._logger = logger
        self._processors = processors
        self._context = context

    def get_DB(self):
        """Return the connection and cursor for the database
        """
        return self._logger.db

    def msg(self, event, **kw):
        if not 'status' in kw:
            res = self._proxy_to_logger('msg', event, status='ok', **kw)
        else:
            res = self._proxy_to_logger('msg', event, **kw)
        return res

    def _attach_obj(self, obj_name, obj_digest):
        self.msg('add_obj', name=obj_name, obj_digest=obj_digest)

    def user_error(self, event, **kw):
        self.msg(event, status='user_error', **kw)

    def user_action(self, event, **kw):
        self.msg(event, status='user_action', **kw)

    def dump_events(self):
        return self._logger.event_list


class FoveaPrintLogger(object):
    """
    (Non-thread safe version of structlog.PrintLogger)
    Prints events into a file AND store event structure internally
    as `event_list` attribute and an SQL database (accessible via get_DB method).

    :param filestream file: File (stream) to print to. (default: stdout)
    :param dbfilepath: tinydb file path or None (default) for in-memory only
    """
    def __init__(self, filestream=None, dbfilepath=None):
        self._file = filestream or sys.stdout
        self._write = self._file.write
        self._flush = self._file.flush
        # permanent
        self.event_list = []
        self.counter = counter_util()
        # ephemeral
        self.event = {}

        # DB
        self._dbfilepath = dbfilepath
        self.db = make_DB(dbfilepath)

    def set_DB_filepath(self, dbfilepath):
        """
        Reset file output and *RESET DATABASE* (clears any previous values).
        Typically, use before logging starts but after diagnostic manager
        object created.
        """
        self._dbfilepath = dbfilepath
        self.db = make_DB(dbfilepath)

    def __repr__(self):
        return '<FoveaPrintLogger(file={0!r},dbfilepath={1!r})>'.format(self._file,self._dbfilepath)

    def msg(self, message):
        """
        Print *message* and commit insertion to TinyDB nosql database.
        """
        self._write(message + '\n')
        self._flush()
        self.event_list.append(self.event.copy())
        # this id will always be the same as tinydb's internal eid
        id = self.counter.get_count()
        db_dict = {'id': id}
        db_dict.update(self.event)
        self.db.insert(db_dict)

    err = debug = info = warning = error = critical = log = msg


def filter_log(logdict, event, invert=False):
    """
    Filter a dictionary of log entries (loaded by load_log)
    by the named event, inverting (logical NOT) if the optional
    invert argument is True.
    """
    d = {}
    for id, val in logdict.items():
        if val['event'] == event:
            if invert:
                continue
            else:
                d[id] = val
        else:
            if invert:
                d[id] = val
            else:
                continue
    return OrderedDict(d)


def load_log(logpath):
    """
    Load and format stored JSON log to use integer keys in an
    OrderedDict.

    Note: You can view the returned log dictionary with info(<logdict>)
    """
    with open(logpath) as f:
        j=json.load(f)['_default']

    # convert keys to integers
    d = {}
    for istr, val in j.items():
        d[int(istr)] = val

    return OrderedDict(d)

def unique_sha(string=''):
    """
    Uses optional string ID combined with current time.
    Returns 40-char hex string digest
    """
    return hashlib.sha1(string+repr(time.time()).replace('.','')).hexdigest()


def make_DB(filepath=None):
    """
    *filepath* to json document (e.g. 'path/to/db.json') else in-memory if
    none is provided.
    """
    if filepath:
        return TinyDB(filepath, storage=CachingMiddleware(JSONStorage))
    else:
        return TinyDB(storage=MemoryStorage)
