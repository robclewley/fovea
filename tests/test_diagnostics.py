import fovea.diagnostics as dgn
from PyDSTool import args
from pprint import pprint

dm = dgn.diagnostic_manager('test')

log = dm.log

log = log.bind(user='fprefect', val=3.1416)
log.msg('begin', status='sofarsogood')
log.user_error('user.forgot_towel')
log.msg('done', status='nope')

test_obj = args(name='test_obj', contents='stuff')

dm.attach_obj(test_obj, 'dummy_objname')

log = log.bind(user='zaphod', val=-1.01)
log.msg('begin', status='yarp')
log.user_action('user.remembered_towel')
log.msg('done')

print("\n")
pprint(log.dump_events())

pprint(dm.log_items_digest[dm.name_to_digest['dummy_objname']])

print("\n")
print(dm.db.all())