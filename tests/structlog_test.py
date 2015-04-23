import structlog
from structlog import BoundLoggerBase, PrintLogger, wrap_logger
from structlog.processors import JSONRenderer, KeyValueRenderer
from pprint import pprint
import sqlite3

global current_ID
# will start at 1 when used
current_ID = 0

def get_id():
    global current_ID
    current_ID += 1
    return current_ID

def ev_store(logger, log_method, event_dict):
    logger.event = event_dict
    return event_dict

class SemanticLogger(BoundLoggerBase):
    def __init__(self, logger, processors, context):
        self._logger = logger
        self._processors = processors
        self._context = context

    def SQLdump(self):
        self._logger.db

    def msg(self, event, **kw):
        #print("Enter SemanticLogger.msg: event list is %s" % str(self._logger.event))
        if not 'status' in kw:
            res = self._proxy_to_logger('msg', event, status='ok', **kw)
        else:
            res = self._proxy_to_logger('msg', event, **kw)
        #print("Exit SemanticLogger.msg: event list is %s" % str(self._logger.event))
        #print("\n")
        return res

    def user_error(self, event, **kw):
        self.msg(event, status='user_error', **kw)

    def user_action(self, event, **kw):
        self.msg(event, status='user_action', **kw)

    def dump_events(self):
        return self._logger.event_list

    def get_DB(self):
        """Return the connection and cursor for the SQL database
        """
        con = self._logger.db
        cur = con.cursor()
        return con, cur

# for debugging purposes only
import sys
global count
count = 0

class MyPrintLogger(object):
    """
    (Non-thread safe version of structlog.PrintLogger)
    Prints events into a file AND store event structure internally
    as `event_list` attribute and an SQL database (accessible via get_DB method).

    :param file file: File to print to. (default: stdout)
    """
    def __init__(self, file=None):
        self._file = file or sys.stdout
        self._write = self._file.write
        self._flush = self._file.flush
        # permanent
        self.event_list = []
        # ephemeral
        self.event = {}
        # SQL
        self.db = sqlite3.connect(':memory:')
        self.db.row_factory = sqlite3.Row
        cur = self.db.cursor()
        cur.execute("CREATE TABLE Events(Id INTEGER PRIMARY KEY, Event TEXT, User TEXT, Status TEXT, Value FLOAT);")
        # debug -- only gets called once!
        global count
        self.count = count
        count += 1

    def __repr__(self):
        return '<MyPrintLogger(file={0!r})-count{1!r}>'.format(self._file,self.count)

    def msg(self, message):
        """
        Print *message*.
        """
        self._write(message + '\n')
        self._flush()
        self.event_list.append(self.event.copy())
        id = get_id()
        eventname = self.event['event']
        user = self.event['user']
        status = self.event['status']
        val = self.event['val']
        cur = self.db.cursor()
        cur.execute("INSERT INTO Events VALUES (%i, '%s', '%s', '%s', %f);" % \
                    (id, eventname, user, status, val))
        self.db.commit()
        #print("In MyPL.msg: event_list = %s, event = %s" % (self.event_list,
        #                                                     self.event))

    err = debug = info = warning = error = critical = log = msg


log = wrap_logger(MyPrintLogger(), wrapper_class=SemanticLogger,
                  processors=[ev_store, KeyValueRenderer()],
                  #context_class=structlog.threadlocal.wrap_dict(dict) #dict
                  ) #JSONRenderer(indent=1, sort_keys=True)])

# ----------------------------

log = log.bind(user='fprefect', val=3.1416)
log.msg('begin', status='sofarsogood')
log.user_error('user.forgot_towel')
log.msg('done', status='nope')

log = log.bind(user='zaphod', val=-1.01)
log.msg('begin', status='yarp')
log.user_action('user.remembered_towel')
log.msg('done')

print("\n")
pprint(log.dump_events())

con, cur = log.get_DB()

cur.execute("SELECT * From Events")
#con.commit() # needed ?
rows = cur.fetchall()

print(rows[0]['User']) # 'user' also works

cur.execute("SELECT User, Status FROM Events WHERE Id=3")
print(cur.fetchone())