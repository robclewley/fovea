"""
Structured logging, noSQL database, and diagnostic administration
"""
import structlog
from structlog import BoundLoggerBase, PrintLogger, wrap_logger
from structlog.processors import JSONRenderer, KeyValueRenderer
from pprint import pprint
from tinydb import TinyDB, where
from tinydb.storages import MemoryStorage

import datetime
import time
import hashlib
import sys

__all__ = ['diagnostic_manager', 'get_unique_name']

def get_unique_name(name):
    try:
        count = name_registry[name]
    except KeyError:
        name_registry[name] = 0
        out_name = name
    else:
        out_name = name + '_%i' % (count + 1)
        name_registry[name] = count + 1
    return out_name


class diagnostic_manager(object):
    def __init__(self, name, dbfilepath=None):
        self.name = name
        self.global_count = 0
        self.log = wrap_logger(FoveaPrintLogger(dbfilepath=dbfilepath),
                       wrapper_class=SemanticLogger,
                       processors=[ev_store, KeyValueRenderer()],
                  ) #JSONRenderer(indent=1, sort_keys=True)])
        # reference to the database
        self.db = self.log._logger.db
        # store any metadata associated with log table entries
        self.log_items_digest = {}
        self.name_to_digest = {}

    def get_events(self):
        """
        Convenience method to fetch "invisible" event_list obscured
        inside the _logger attribute (due to mixin-like inheritance magic
        in structlog)
        """
        return self.log._logger.event_list

    def get_unique_name(self, name):
        """
        Convenience method for access to function
        """
        return get_unique_name(name)

    def attach_obj(self, obj, name):
        # sha the name and obj repr
        digest = unique_sha(repr(obj)+name)
        # add obj to log_items
        self.log_items_digest[digest] = obj
        if name in self.name_to_digest:
            raise ValueError()
        else:
            self.name_to_digest[name] = digest
        self.log._attach_obj(obj_name=name, obj_digest=digest)


global name_registry
name_registry = {}

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
        self._logger.db

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

    def get_DB(self):
        """Return the connection and cursor for the database
        """
        return self._logger.db


class FoveaPrintLogger(object):
    """
    (Non-thread safe version of structlog.PrintLogger)
    Prints events into a file AND store event structure internally
    as `event_list` attribute and an SQL database (accessible via get_DB method).

    :param file file: File to print to. (default: stdout)
    :param dbfilepath: tinydb file path or None (default) for in-memory only
    """
    def __init__(self, file=None, dbfilepath=None):
        self._file = file or sys.stdout
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


def unique_sha(string=''):
    """
    uses optional string ID with current time
    Returns 40-char hex string digest
    """
    return hashlib.sha1(string+repr(time.time()).replace('.','')).hexdigest()


def make_DB(filepath=None):
    """
    *filepath* to json document (e.g. 'path/to/db.json') else in-memory if none provided.
    """
    if filepath:
        return TinyDB(filepath)
    else:
        return TinyDB(storage=MemoryStorage)
